// Note: the convective part of the convection-diffusion equation for alpha has 
// been solved before this part of the code is executed.
// The purpose of the following code is to correct for the diffusive part
{
	// RBC radius
	scalar a = 4e-6;

	// Evaluate diffusivity based on current shear rate and hematocrit
	// Note: We have to prefix the sqrt function with a namespace specifier to avoid compiler confusion.
	volScalarField strainRate = Foam::sqrt(2.0) * mag(symm(fvc::grad(U)));
	volScalarField D = Dab + 0.15 * a*a * strainRate * alpha1 * pow(1 - alpha1, 0.8);

	// Form equation system
    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1)
      - fvc::ddt(alpha1)
      - fvm::laplacian
        (
            volScalarField("Dab", D + alphatab*turbulence->nut()),
            alpha1
        )
    );

    alpha1Eqn.solve();

    alpha2 = 1.0 - alpha1;
    rhoPhi += alpha1Eqn.flux()*(rho1 - rho2);
}

rho = alpha1*rho1 + alpha2*rho2;
