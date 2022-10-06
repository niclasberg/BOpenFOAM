surfaceScalarField alphaf1("alphaf1", fvc::interpolate(alpha1));
surfaceScalarField alphaf2("alphaf2", scalar(1) - alphaf1);

volScalarField rAU1
(
    IOobject::groupName("rAU", phase1.name()),
    1.0
   /(
        U1Eqn.A()
      + max(phase1.residualAlpha() - alpha1, scalar(0))
       *rho1/runTime.deltaT()
    )
);
volScalarField rAU2
(
    IOobject::groupName("rAU", phase2.name()),
    1.0
   /(
        U2Eqn.A()
      + max(phase2.residualAlpha() - alpha2, scalar(0))
       *rho2/runTime.deltaT()
    )
);

surfaceScalarField alpharAUf1
(
    fvc::interpolate(max(alpha1, phase1.residualAlpha())*rAU1)
);
surfaceScalarField alpharAUf2
(
    fvc::interpolate(max(alpha2, phase2.residualAlpha())*rAU2)
);

// Lift fluxes
tmp<surfaceScalarField> phiF1(fvc::flux(rAU1*fluid.F()));
tmp<surfaceScalarField> phiF2(-fvc::flux(rAU2*fluid.F()));
{
    // Phase-1 dispersion, lift and wall-lubrication flux
    phiF1 = ;

    // Phase-1 dispersion, lift and wall-lubrication flux
    phiF2 = ;
}


// --- Pressure corrector loop
while (pimple.correct())
{
    // Update continuity errors due to temperature changes
    #include "correctContErrs.H"

    volScalarField rho("rho", fluid.rho());

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    //MRF.correctBoundaryFlux(U1, phi1);
    //MRF.correctBoundaryFlux(U2, phi2);

    volVectorField HbyA1
    (
        IOobject::groupName("HbyA", phase1.name()),
        U1
    );
    HbyA1 =
        rAU1
       *(
            U1Eqn.H()
          + max(phase1.residualAlpha() - alpha1, scalar(0))
           *rho1*U1.oldTime()/runTime.deltaT()
        );

    volVectorField HbyA2
    (
        IOobject::groupName("HbyA", phase2.name()),
        U2
    );
    HbyA2 =
        rAU2
       *(
            U2Eqn.H()
         +  max(phase2.residualAlpha() - alpha2, scalar(0))
           *rho2*U2.oldTime()/runTime.deltaT()
        );

    surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    surfaceScalarField phig1
    (
        alpharAUf1
       *(
           ghSnGradRho
         - alphaf2*fvc::interpolate(rho1 - rho2)*(g & mesh.Sf())
        )
    );

    surfaceScalarField phig2
    (
        alpharAUf2
       *(
           ghSnGradRho
         - alphaf1*fvc::interpolate(rho2 - rho1)*(g & mesh.Sf())
        )
    );


    // ddtPhiCorr filter -- only apply in pure(ish) phases
    surfaceScalarField alphaf1Bar(fvc::interpolate(fvc::average(alphaf1)));
    surfaceScalarField phiCorrCoeff1(pos0(alphaf1Bar - 0.99));
    surfaceScalarField phiCorrCoeff2(pos0(0.01 - alphaf1Bar));

    {
        surfaceScalarField::Boundary& phiCorrCoeff1Bf =
            phiCorrCoeff1.boundaryFieldRef();

        surfaceScalarField::Boundary& phiCorrCoeff2Bf =
            phiCorrCoeff2.boundaryFieldRef();

        forAll(mesh.boundary(), patchi)
        {
            // Set ddtPhiCorr to 0 on non-coupled boundaries
            if
            (
               !mesh.boundary()[patchi].coupled()
             || isA<cyclicAMIFvPatch>(mesh.boundary()[patchi])
            )
            {
                phiCorrCoeff1Bf[patchi] = 0;
                phiCorrCoeff2Bf[patchi] = 0;
            }
        }
    }

    // Phase-1 predicted flux
    surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1.name()),
        fvc::flux(HbyA1)
      + phiCorrCoeff1*fvc::interpolate(alpha1.oldTime()*rho1.oldTime()*rAU1)
       *(
            phi1.oldTime() - fvc::flux(U1.oldTime())
        )/runTime.deltaT()
      - phiF1()
      - phig1
    );

    // Phase-2 predicted flux
    surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", phase2.name()),
        fvc::flux(HbyA2)
      + phiCorrCoeff2*fvc::interpolate(alpha2.oldTime()*rho2.oldTime()*rAU2)
       *(
            phi2.oldTime() - fvc::flux(U2.oldTime())
        )/runTime.deltaT()
      - phiF2()
      - phig2
    );

    // Face-drag coefficients
    surfaceScalarField rAUKd1(fvc::interpolate(rAU1*Kd));
    surfaceScalarField rAUKd2(fvc::interpolate(rAU2*Kd));

    // Construct the mean predicted flux
    // including explicit drag contributions based on absolute fluxes
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alphaf1*(phiHbyA1 + rAUKd1*phi2)
      + alphaf2*(phiHbyA2 + rAUKd2*phi1)
    );

    // Construct pressure "diffusivity"
    surfaceScalarField rAUf
    (
        "rAUf",
        mag(alphaf1*alpharAUf1 + alphaf2*alpharAUf2)
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - (
                alphaf1.boundaryField()*phi1.boundaryField()
              + alphaf2.boundaryField()*phi2.boundaryField()
            )
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    // Iterate over the pressure equation to correct for non-orthogonality
    while (pimple.correctNonOrthogonal())
    {
        // Construct the transport part of the pressure equation
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        solve
        (
            pEqnIncomp,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

        // Correct fluxes and velocities on last non-orthogonal iteration
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqnIncomp.flux();

            surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);

            // Partial-elimination phase-flux corrector
            {
                surfaceScalarField phi1s
                (
                    phiHbyA1 + alpharAUf1*mSfGradp
                );

                surfaceScalarField phi2s
                (
                    phiHbyA2 + alpharAUf2*mSfGradp
                );

                surfaceScalarField phir
                (
                    ((phi1s + rAUKd1*phi2s) - (phi2s + rAUKd2*phi1s))
                   /(1 - rAUKd1*rAUKd2)
                );

                phi1 = phi + alphaf2*phir;
                phi2 = phi - alphaf1*phir;
            }

            // Optionally relax pressure for velocity correction
            p_rgh.relax();

            mSfGradp = pEqnIncomp.flux()/rAUf;

            // Partial-elimination phase-velocity corrector
            {
                volVectorField Us1
                (
                    HbyA1
                  + fvc::reconstruct(alpharAUf1*mSfGradp - phiF1() - phig1)
                );

                volVectorField Us2
                (
                    HbyA2
                  + fvc::reconstruct(alpharAUf2*mSfGradp - phiF2() - phig2)
                );

                volScalarField D1(rAU1*Kd);
                volScalarField D2(rAU2*Kd);

                U = alpha1*(Us1 + D1*U2) + alpha2*(Us2 + D2*U1);
                volVectorField Ur(((1 - D2)*Us1 - (1 - D1)*Us2)/(1 - D1*D2));

                U1 = U + alpha2*Ur;
                U1.correctBoundaryConditions();
                //fvOptions.correct(U1);

                U2 = U - alpha1*Ur;
                U2.correctBoundaryConditions();
                //fvOptions.correct(U2);

                U = fluid.U();
            }
        }
    }

    // Update and limit the static pressure
    p = p_rgh + rho*gh;

    // Limit p_rgh
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();
}

