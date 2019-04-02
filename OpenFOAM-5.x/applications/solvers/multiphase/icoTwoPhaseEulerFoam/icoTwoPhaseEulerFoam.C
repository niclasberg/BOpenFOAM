    #include "fvCFD.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedValueFvsPatchFields.H"
#include "twoPhaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "CourantNos.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            Info << max(alpha1) << endl;
			// Solve transport equation for phase fraction
			fluid.solve();

			// Update viscosities
            fluid.correct();
            /*alphaPhi1.write();
            alphaPhi2.write();
            if(runTime.timeIndex() >= 2)
                return 0;*/
            
            // Small number
            //scalar tol = 1e-3;

            // Evaluate viscosities
            Info << "computing viscosities" << endl;
            volScalarField nu1 = phase1.viscosityModel().mu() / phase1.rho();
            volScalarField mu2 = phase2.viscosityModel().mu();
            volScalarField nu2 = mu2 / phase2.rho();

            // Evaluate drag force coefficient
            Info << "drag coefficients" << endl;
            volScalarField Kd(
                IOobject(
                    "Kd",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                9 * mu2 * alpha1 * (pow(alpha1, 0.43) + exp(2.68*alpha1)) / (2.*pow(fluid.d(), 2))/*,
                zeroGradientFvPatchScalarField::typeName*/
            );

            // Evaluate lift force coefficient
            /*volSymmTensorField gamma2(symm(fvc::grad(U2)));
            volSymmTensorField Kl(
                IOobject(
                    "Kl",
                    mesh.time().timeName(),
                    mesh
                ),
                19.38*sqrt(phase2.rho() * mu2) * pow(max(tr(gamma2), dimensionedScalar("smallShear", dimLess/dimTime, tol), -0.25) * gamma2 / (4. * 3.1415927 * pow(fluid.d(), 2)),
                zeroGradientFvPatchScalarField::typeName
            );*/

            // Remove drag and lift at fixed-flux boundaries
            volScalarField::Boundary & Kdb = Kd.boundaryFieldRef();
            //volScalarField::Boundary & Klb = Kl.boundaryFieldRef();
            forAll(phase1.phi().boundaryField(), patchi)
                if(isA<fixedValueFvsPatchScalarField>(phase1.phi().boundaryField()[patchi])) {
                    Kdb[patchi] = 0.0;
                    //Klb[patchi] = 0.0;
                }

            //Kd.correctBoundaryConditions();
            //Kl.correctBoundaryConditions();

            // Construct momentum equation matrices
            Info<< "Constructing momentum equations" << endl;

            fvVectorMatrix U1Eqn
            (
                // Material derivative
                fvm::ddt(alpha1, U1)            
                + fvm::div(alphaPhi1, U1)

                // Diffusion term
                - fvm::laplacian(alpha1*nu1, U1)
                - fvc::div
                    (
                        alpha1*(nu1*dev(T(fvc::grad(U1))) /*- ((2.0/3.0)*I)*k*/),
                        "div(Rc)"
                    )
                
                == 
                - fvm::Sp(Kd / phase1.rho(), U1)    // Implicit part of the drag term
                //- fvm::Sp(Kl / phase1.rho(), U1)    // Implicit part of the lift term
            );
            //U1Eqn.relax();
            //U1.correctBoundaryConditions();
            
            fvVectorMatrix U2Eqn
            (
                // Material derivative
                fvm::ddt(alpha2, U2) 
                + fvm::div(alphaPhi2, U2)

                // Diffusion term
                - fvm::laplacian(alpha2*nu2, U2)
                - fvc::div
                    (
                        alpha2*(nu2*dev(T(fvc::grad(U2))) /*- ((2.0/3.0)*I)*k*/),
                        "div(Rc)"
                    )
                
                ==
                - fvm::Sp(Kd / phase2.rho(), U2)
            );
            //U2Eqn.relax();
            //U2.correctBoundaryConditions();

            // Evaluate density gradient
            volScalarField rho("rho", fluid.rho());
            //surfaceScalarField ghSnGradRho(ghf*fvc::snGrad(rho)*mesh.magSf());

            surfaceScalarField alphaf1("alphaf1", fvc::interpolate(alpha1));
            surfaceScalarField alphaf2("alphaf2", scalar(1) - alphaf1);

            volScalarField rAU1(1.0 / U1Eqn.A());
            volScalarField rAU2(1.0 / U2Eqn.A());

            surfaceScalarField rAUf1(alphaf1 / fvc::interpolate(U1Eqn.A()));
            surfaceScalarField rAUf2(alphaf2 / fvc::interpolate(U2Eqn.A()));

            volVectorField HbyA1 = rAU1*U1Eqn.H();
            volVectorField HbyA2 = rAU2*U2Eqn.H();

            // Explicit source terms in the pressure equation
            surfaceScalarField phiDrag1(
                fvc::interpolate(Kd * rAU1) * phi2 / phase1.rho()
            );

            surfaceScalarField phiDrag2(
                fvc::interpolate(Kd * rAU2) * phi1 / phase2.rho()
            );

            // Remove drag fluxes at the fixed flux boundaries
            /*forAll(p_rgh.boundaryField(), patchi) {
                if(isA<zeroGradientFvPatchScalarField>(p_rgh.boundaryField()[patchi]))
                {
                    phiDrag1.boundaryFieldRef()[patchi] = 0.0;
                    phiDrag2.boundaryFieldRef()[patchi] = 0.0;
                }
            }*/

            // Lift force

            // Predict fluxes
            surfaceScalarField phiHbyA1(
                fvc::flux(HbyA1)  // replace by (fvc::interpolate(HbyA1) & mesh.Sh())?
              + rAUf1*fvc::ddtCorr(U1, phi1)
              + phiDrag1
            );

            surfaceScalarField phiHbyA2(
                fvc::flux(HbyA2) 
              + rAUf2*fvc::ddtCorr(U2, phi2)
              + phiDrag2
            );      

            HbyA1 += Kd * rAU1 * U2 / phase1.rho();
            HbyA2 += Kd * rAU2 * U1 / phase2.rho();

            //phi = alphaf1*phiHbyA1 + alphaf2*phiHbyA2;
            surfaceScalarField phiHbyA(alphaf1*phiHbyA1 + alphaf2*phiHbyA2);
            surfaceScalarField rAUf(alphaf1*rAUf1/phase1.rho() + alphaf2*rAUf2/phase2.rho() );

            // Update the fixedFluxPressure BCs to ensure flux consistency
            {
                surfaceScalarField::Boundary phib(phi.boundaryField());
                phib = 
                    alphaf1.boundaryField() * (mesh.Sf().boundaryField() & U1.boundaryField())
                  + alphaf2.boundaryField() * (mesh.Sf().boundaryField() & U2.boundaryField());

                setSnGrad<fixedFluxPressureFvPatchScalarField>
                (
                    p_rgh.boundaryFieldRef(),
                    (phiHbyA.boundaryField() - phib) / (mesh.magSf().boundaryField()*rAUf.boundaryField())
                );
            }

            while (pimple.correctNonOrthogonal()) {
                fvScalarMatrix pEqnIncomp(fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh) );
                pEqnIncomp.setReference(pRefCell, pRefValue);

                solve(pEqnIncomp, mesh.solver(p_rgh.select(pimple.finalInnerIter())));

                if (pimple.finalNonOrthogonalIter()) {
                    surfaceScalarField SfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);

                    // Correct fluxes
                    /*phi1 = phiHbyA1 + rAUf1 * SfGradp / phase1.rho();
                    phi2 = phiHbyA2 + rAUf2 * SfGradp / phase2.rho();
                    phi = alphaf1*phiHbyA1 + alphaf2*phiHbyA2;

                    //phase.phi() = phiHbyAs[phasei] + rAlphaAUfs[phasei]*mSfGradp/phase.rho();
                    //phi += alphafs[phasei]*phiHbyAs[phasei] + mag(alphafs[phasei]*rAlphaAUfs[phasei]) *mSfGradp/phase.rho();

                    // Relax pressure
                    p_rgh.relax();
                    p = p_rgh + rho*gh;
                    SfGradp = pEqnIncomp.flux()/rAUf;

                    // Correct velocities
                    U1 = HbyA1 + fvc::reconstruct(rAUf1 * SfGradp / phase1.rho());
                    U2 = HbyA2 + fvc::reconstruct(rAUf2 * SfGradp / phase2.rho());
                    U1.correctBoundaryConditions();
                    U2.correctBoundaryConditions();*/

                    // Evaluate mixture velocity
                    U = fluid.U();
                }
            }
        }

    

        #include "write.H"

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
