{
    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    alpha1.correctBoundaryConditions();

    surfaceScalarField phic("phic", phi);
    surfaceScalarField phir("phir", phi1 - phi2);

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
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
		{
		    surfaceScalarField::Boundary& alphaPhiBf = alphaPhic1.boundaryFieldRef();
			const volScalarField::Boundary& alphaBf = boundaryField();
			const surfaceScalarField::Boundary& phiBf = phi1.boundaryField();

			forAll(alphaPhiBf, patchi) {
				fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];
				if (!alphaPhip.coupled())
				    alphaPhip = phiBf[patchi]*alphaBf[patchi];
			}
		}

        if (nAlphaSubCycles > 1)
        {
            for
            (
                subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
                !(++alphaSubCycle).end();
            )
            {
                surfaceScalarField alphaPhic10(alphaPhic1);

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phi,
                    alphaPhic10,
                    1,
                    0
                );

                if (alphaSubCycle.index() == 1)
                {
                    alphaPhi1 = alphaPhic10;
                }
                else
                {
                    alphaPhi1 += alphaPhic10;
                }
            }

            alphaPhi1 /= nAlphaSubCycles;
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi_,
                alphaPhic1,
                1,
                0
            );

            alphaPhi1 = alphaPhic1;
        }

        //alphaRhoPhi1 = fvc::interpolate(phase1_.rho())*phase1_.alphaPhi();

        alphaPhi2 = phi - alphaPhi1;
        {
		    surfaceScalarField::Boundary& alphaPhiBf = alphaPhi2.boundaryFieldRef();
			const volScalarField::Boundary& alphaBf = alpha2.boundaryField();
			const surfaceScalarField::Boundary& phiBf = phi2.boundaryField();

			forAll(alphaPhiBf, patchi) {
				fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];
				if (!alphaPhip.coupled())
				    alphaPhip = phiBf[patchi]*alphaBf[patchi];
			}
		}
        //alphaRhoPhi2 = fvc::interpolate(phase2_.rho())*phase2_.alphaPhi();

        Info<< alpha1.name() << " volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        // Ensure the phase-fractions are bounded
        alpha1.max(0);
        alpha1.min(1);

        alpha2 = scalar(1) - alpha1;
    }
}
