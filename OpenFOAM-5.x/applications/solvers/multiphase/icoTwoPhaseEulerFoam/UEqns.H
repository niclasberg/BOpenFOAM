Info<< "Constructing momentum equations" << endl;

/*MRF.correctBoundaryVelocity(U1);
MRF.correctBoundaryVelocity(U2);
MRF.correctBoundaryVelocity(U);*/

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

// Evaluate drag force coefficient
volScalarField Kd(
	9 * mu1 * alpha1 * (pow(alpha1, 0.43) + exp(2.68*alpha1)) / (2.*a*a)
);

{
    //volScalarField Vm(fluid.Vm());

    {
        U1Eqn =
        (
            fvm::ddt(alpha1, rho1, U1) 
 		  + fvm::div(alphaRhoPhi1, U1)
          + phase1.turbulence().divDevRhoReff(U1)
        );
        U1Eqn.relax();
        U1Eqn += fvm::Sp(Kd, U1);		// Drag force
        fvOptions.constrain(U1Eqn);
        U1.correctBoundaryConditions();
    }

    {
        U2Eqn =
        (
            fvm::ddt(alpha2, rho2, U2) + fvm::div(alphaRhoPhi2, U2)
          - fvm::Sp(contErr2, U2)
          ///+ MRF.DDt(alpha2*rho2 + Vm, U2)
          + phase2.turbulence().divDevRhoReff(U2)
         /*==
          - Vm
           *(
                fvm::ddt(U2)
              + fvm::div(phi2, U2)
              - fvm::Sp(fvc::div(phi2), U2)
              - DDtU1
            )
          + fvOptions(alpha2, rho2, U2)*/
        );
        U2Eqn.relax();
        U2Eqn += fvm::Sp(Kd, U2);		// Drag force
        fvOptions.constrain(U2Eqn);
        U2.correctBoundaryConditions();
        //fvOptions.correct(U2);
    }
}
