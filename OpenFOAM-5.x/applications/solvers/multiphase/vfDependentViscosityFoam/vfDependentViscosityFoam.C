/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    vfDependentViscosityFoam

Description
    Solver for mixing 2 incompressible fluids.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "transportModel/vfDependentViscosityTwoPhaseMixture/vfDependentViscosityTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include <cassert>
#include "fieldInterpolator.H"
#include "fvMeshFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	Foam::argList::addBoolOption
	(
		argList::postProcessOptionName,
		"Execute functionObjects only"
	);

	// Run in post-process mode if the -postProcess flag has been supplied
	if (argList::postProcess(argc, argv))
	{
		Foam::timeSelector::addOptions();
		#include "addRegionOption.H"
		#include "addFunctionObjectOptions.H"

		// Set functionObject post-processing mode
		functionObject::postProcess = true;

		#include "setRootCase.H"

		if (args.optionFound("list"))
		{
		    functionObjectList::list();
		    return 0;
		}

		#include "createTime.H"
		Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
		#include "createMesh.H"
		//#include "createControl.H"
		#include "createFields.H"

		// Externally stored dictionary for functionObjectList if not constructed from runTime
		dictionary functionsDict;

		// Construct functionObjectList
		HashSet<word> selectedFields;
		autoPtr<functionObjectList> functionsPtr
		(
		    functionObjectList::New(args, runTime, functionsDict, selectedFields)
		);

		// Read number of substeps
		int substeps = 4;

		// Create interpolators
		FieldInterpolator<volVectorField> uInterp(mesh, "U");
		FieldInterpolator<volScalarField> p_rghInterp(mesh, "p_rgh");
		FieldInterpolator<volScalarField> alphaInterp(mesh, "alpha.RBC");

		for(int timei = 0; timei < timeDirs.size(); ++timei) 
		{
		    runTime.setTime(timeDirs[timei], timei);
		    Info<< "Time = " << runTime.timeName() << ", deltaT = " << runTime.deltaTValue() << endl;

		    if (mesh.readUpdate() != polyMesh::UNCHANGED)
		    {
		        // Update functionObjects if mesh changes
		        functionsPtr = functionObjectList::New
		        (
		            args,
		            runTime,
		            functionsDict,
		            selectedFields
		        );
		    }

			// Read data
			uInterp.readTime(timeDirs[timei].name());
			p_rghInterp.readTime(timeDirs[timei].name());
			alphaInterp.readTime(timeDirs[timei].name());

			FatalIOError.throwExceptions();

			if(timei > 0) {
				scalar t0 = timeDirs[timei-1].value();	// Old time
				scalar t1 = timeDirs[timei].value();	// New time

				scalar dataDt = t1 - t0;
				runTime.setTime(timeDirs[timei], timei);
				runTime.setDeltaT(dataDt);

				// Excecute subcycles
				for(subCycleTime subIter = subCycleTime(runTime, substeps); !subIter.end(); ++subIter) {
					scalar timeStepFraction = (runTime.value() - t0) / dataDt;
					Info << "Executing subcycle, current time = " << runTime.timeName() << ", deltaT = " << runTime.deltaTValue() << ", fraction of timestep  = " << timeStepFraction << endl;

					// Interpolate fields
					U = uInterp.interpolate(timeStepFraction, dataDt);
					p_rgh = p_rghInterp.interpolate(timeStepFraction, dataDt);
					alpha1 = alphaInterp.interpolate(timeStepFraction, dataDt);

					// Evaluate dependent fields
					alpha2 = 1. - alpha1;
					rho = alpha1*rho1 + alpha2*rho2;
					p = p_rgh + rho*gh;

					// Update viscosity
				    twoPhaseProperties.correct();

					try
					{
						// Execute function objects
						if(functionsPtr->status())
							forAll(functionsPtr(), objectI)
        						functionsPtr()[objectI].execute();
					}
					catch (IOerror& err)
					{
						Warning<< err << endl;
					}
				}
			} 

			// Execute the write operation for each function object
			if(functionsPtr->status())
				forAll(functionsPtr(), objectI)
					functionsPtr()[objectI].write();

			// Execute the functionObject 'end()' function for the last time
			if (timei == timeDirs.size()-1)
			{
			    functionsPtr->end();
			}
		}

		Info<< "End\n" << endl;

	} else {

		#include "setRootCase.H"
		#include "createTime.H"
		#include "createMesh.H"
		#include "createControl.H"
		#include "initContinuityErrs.H"
		#include "createFields.H"
		#include "createTimeControls.H"
		#include "CourantNo.H"
		#include "setInitialDeltaT.H"

		turbulence->validate();

		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

		Info<< "\nStarting time loop\n" << endl;

		while (runTime.run())
		{
		    #include "readTimeControls.H"
		    #include "CourantNo.H"
		    #include "alphaCourantNo.H"
		    #include "setDeltaT.H"

		    runTime++;

		    Info<< "Time = " << runTime.timeName() << nl << endl;

		    twoPhaseProperties.correct();

		    #include "alphaEqnSubCycle.H"
		    #include "alphaDiffusionEqn.H"

		    // --- Pressure-velocity PIMPLE corrector loop
		    while (pimple.loop())
		    {
		        #include "UEqn.H"

		        // --- Pressure corrector loop
		        while (pimple.correct())
		        {
		        	#include "pEqn.H"
		        }

		        if (pimple.turbCorr())
		        {
		            turbulence->correct();
		        }
		    }

		    runTime.write();

		    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
		        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
		        << nl << endl;
		}

		Info<< "End\n" << endl;
	}

	return 0;
}


// ************************************************************************* //
