#include "twoPhaseModel.H"
#include "fixedValueFvsPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "surfaceInterpolate.H"
#include "MULES.H"
#include "subCycle.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcCurl.H"
#include "fvmDdt.H"
#include "fvmDiv.H"

Foam::PhaseModel::PhaseModel(const Foam::fvMesh & mesh, const Foam::dictionary & dict, const Foam::word & phaseName)
:
    name_(phaseName),
    rho0_("rho", dimDensity, dict.subDict(name_)),
    mesh_(mesh),
    alpha_
    (
        IOobject
        (
            IOobject::groupName("alpha", name_),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0)
    ),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", name_),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
    ),
    viscosityModel_(ViscosityModel::New(dict.subDict(name_), *this))
{
    const word phiName = IOobject::groupName("phi", name_);

    IOobject phiHeader
    (
        phiName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ
    );

    if (phiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        Info<< "Reading face flux field " << phiName << endl;

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U_.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U_.boundaryField(), i) {
            if(
                isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i])
             || isA<slipFvPatchVectorField>(U_.boundaryField()[i])
             || isA<partialSlipFvPatchVectorField>(U_.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvsPatchScalarField::typeName;
            }
        }

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::flux(U_),
                phiTypes
            )
        );
    }
}

void Foam::PhaseModel::correctInflowOutflow(surfaceScalarField& alphaPhi) const
{
    surfaceScalarField::Boundary& alphaPhiBf = alphaPhi.boundaryFieldRef();
    const volScalarField::Boundary& alphaBf = alpha_.boundaryField();
    const surfaceScalarField::Boundary& phiBf = phi().boundaryField();

    forAll(alphaPhiBf, patchi) {
        fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];
        if (!alphaPhip.coupled())
            alphaPhip = phiBf[patchi]*alphaBf[patchi];
    }
}

Foam::twoPhaseModel::twoPhaseModel(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    phase1_(mesh, *this, wordList(lookup("phases"))[0]),
    phase2_(mesh, *this, wordList(lookup("phases"))[1]),
    d_("d", dimLength, *this),
    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->calcPhi()
    )
{
    
}

Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseModel::calcPhi() const
{
    return
        fvc::interpolate(phase1_.alpha())*phase1_.phi()
      + fvc::interpolate(phase2_.alpha())*phase2_.phi();
}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseModel::rho() const
{
    return tmp<volScalarField>(
            new volScalarField("rho", phase1_.alpha()*phase1_.rho() + phase2_.alpha()*phase2_.rho())
        );
}

void Foam::twoPhaseModel::correct()
{
    phase1_.viscosityModel().correct();
    phase2_.viscosityModel().correct();
}

void Foam::twoPhaseModel::solve()
{
    const Time& runTime = mesh_.time();

    volScalarField& alpha1 = phase1_.alpha();
    volScalarField& alpha2 = phase2_.alpha();

    const dictionary& alphaControls = mesh_.solverDict
    (
        alpha1.name()
    );

    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    alpha1.correctBoundaryConditions();

    surfaceScalarField phic("phic", phi());
    surfaceScalarField phir("phir", phase1().phi() - phase2().phi());

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        /*surfaceScalarField phiTot(
              phic - fvc::flux(-phir, alpha2, alpharScheme)
        );
        phase1_.correctInflowOutflow(phiTot);
        fvScalarMatrix alpha1Eqn(
            fvm::ddt(alpha1)
            + fvm::div(phiTot, alpha1, alphaScheme)
        );
        alpha1Eqn.relax();
        alpha1Eqn.solve();*/

        surfaceScalarField alphaPhic1
        (
            fvc::flux
            (
                phic,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

		// Correct inflow and outflow boundary conditions for alpha1
		phase1().correctInflowOutflow(alphaPhic1);

        // Limit the flux
        if (nAlphaSubCycles > 1) {
            for(subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles); !(++alphaSubCycle).end(); ) {
                surfaceScalarField alphaPhic10(alphaPhic1);
                MULES::explicitSolve(
                    geometricOneField(), 
                    alpha1, 
                    phi_, 
                    alphaPhic10, 
                    zeroField(), 
                    zeroField(),
                    1, 0);

                if (alphaSubCycle.index() == 1)
                    phase1_.alphaPhi() = alphaPhic10;
                else
                    phase1_.alphaPhi() += alphaPhic10;
            }

            phase1_.alphaPhi() /= nAlphaSubCycles;
        } else {
            MULES::explicitSolve(
                    geometricOneField(), 
                    alpha1, 
                    phi_, 
                    alphaPhic1, 
                    zeroField(), 
                    zeroField(),
                    1, 0);
            phase1_.alphaPhi() = alphaPhic1;
        }

        phase2_.alphaPhi() = phi() - phase1().alphaPhi();
        //phase2_.correctInflowOutflow(phase2().alphaPhi());

        Info<< alpha1.name() << " volume fraction = "
            << alpha1.weightedAverage(mesh_.V()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        // Ensure the phase-fractions are bounded
        alpha1.max(0);
        alpha1.min(1);

        alpha2 = scalar(1) - alpha1;
    }
}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseModel::dragForce() const
{
    volScalarField limitedAlpha(min(1e-4, phase1_.alpha()));

    tmp<volScalarField> Kd(
        new volScalarField(
            IOobject(
                "Kd",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            9 * phase2_.viscosityModel().mu() * limitedAlpha * (pow(limitedAlpha, 0.43) + exp(2.68*limitedAlpha)) / (2.*pow(d_, 2)),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField::Boundary& Kbf = Kd.ref().boundaryFieldRef();

    // Remove drag at fixed-flux boundaries
    forAll(phase1_.phi().boundaryField(), patchi) {
        if(isA<fixedValueFvsPatchScalarField>(phase1_.phi().boundaryField()[patchi]))
            Kbf[patchi] = 0.0;
    }

    Kd.ref().correctBoundaryConditions();

    return Kd;
}