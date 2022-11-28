#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "vfDependentViscosityTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include <cassert>
#include "../fieldInterpolator.H"
#include "addToRunTimeSelectionTable.H"
#include "plateletLpt.H"
#include <iostream>
#include "OSspecific.H"
#include "writeFuns.H"

using namespace Foam;

void writeFieldHeader(std::ostream & os, word fieldName, int numTuples, int numComponents, word dataType)
{
    os << fieldName << " " << numComponents << " " << numTuples << " " << dataType << std::endl;
}

void writeCloudAsVTK(const basicPlateletCloud & cloud, word timeName)
{
    // Create file name
    fileName outputFolder = "VTK" / cloud.name();
    if( ! isDir(outputFolder))
        mkDir(outputFolder);
    word vtkFileName = outputFolder + "/" + cloud.name() + timeName + ".vtk";
    word csvFileName = outputFolder + "/" + cloud.name() + timeName + ".csv";
    Info << "Writing to " << vtkFileName << endl;
    Info << "Writing CSV file to " << csvFileName << endl;

    // Create output stream
    std::ofstream os(vtkFileName.c_str(), std::iostream::binary);
    // std::ofstream osCsv(csvFileName.c_str(), std::iostream::out);

    // Write header
    bool isBinary = true;
    writeFuns::writeHeader(os, isBinary, cloud.name());
    os << "DATASET POLYDATA" << std::endl;

    // Extract data
    DynamicList<floatScalar> position; 
    DynamicList<floatScalar> U;
    DynamicList<label> active;
    DynamicList<label> origId;
    DynamicList<floatScalar> age;
    DynamicList<floatScalar> stressHistory;
    DynamicList<floatScalar> pas;
    DynamicList<floatScalar> stressRateHistory;
    DynamicList<floatScalar> tau;
    DynamicList<floatScalar> Uc;
    DynamicList<floatScalar> Re;
    DynamicList<floatScalar> d;
    DynamicList<floatScalar> rho;
    DynamicList<floatScalar> cD;

    label i = 0;
    forAllConstIter(basicPlateletCloud, cloud, iter)
    {
        const basicPlateletParcel & p = iter();
        active.append(p.active());
        writeFuns::insert(p.position(), position);
        writeFuns::insert(p.U(), U);
        origId.append(p.origId());
        writeFuns::insert(p.age(), age);
		writeFuns::insert(p.stressHistory(), stressHistory);
        writeFuns::insert(p.pas(), pas);
        floatScalar pRe = p.Re(p.U(), p.d(), p.rhoc(), p.muc());
        writeFuns::insert(p.stressRateHistory(), stressRateHistory);
        writeFuns::insert(p.tau(), tau);
        writeFuns::insert(p.Uc(), Uc);
        writeFuns::insert(pRe, Re);
        writeFuns::insert(p.d(), d);
        writeFuns::insert(p.rho(), rho);
        floatScalar pcD = pRe > 1000 ? 0.424 : 24.0 / pRe *(1.0 + 1.0/6.0*pow(pRe, 2.0/3.0));
        writeFuns::insert(pcD, cD);
        i++;
    }

    // Write points
    os << "POINTS " << cloud.size() << " float" << std::endl;
    writeFuns::write(os, isBinary, position);

    // Write point data
    if(cloud.size() > 0) {
        label numFields = 13;
        writeFuns::writePointDataHeader(os, cloud.size(), numFields);
        
        // Write fields
        writeFieldHeader(os, "U", cloud.size(), 3, "float");
        writeFuns::write(os, isBinary, U);

        writeFieldHeader(os, "active", cloud.size(), 1, "int");
        writeFuns::write(os, isBinary, active);

        writeFieldHeader(os, "id", cloud.size(), 1, "int");
        writeFuns::write(os, isBinary, origId);

        writeFieldHeader(os, "age", cloud.size(), 1, "float");
        writeFuns::write(os, isBinary, age);

        writeFieldHeader(os, "stressHistory", cloud.size(), 1, "float");
        writeFuns::write(os, isBinary, stressHistory);

        writeFieldHeader(os, "pas", cloud.size(), 1, "float");
        writeFuns::write(os, isBinary, pas);

        writeFieldHeader(os, "stressRateHistory", cloud.size(), 1, "float");
        writeFuns::write(os, isBinary, stressRateHistory);

        writeFieldHeader(os, "tau", cloud.size(), 6, "float");
        writeFuns::write(os, isBinary, tau);

        writeFieldHeader(os, "fluidVelocity", cloud.size(), 3, "float");
        writeFuns::write(os, isBinary, Uc);

        writeFieldHeader(os, "Rep", cloud.size(), 1, "float");
        writeFuns::write(os, isBinary, Re);

        writeFieldHeader(os, "d", cloud.size(), 1, "float");
        writeFuns::write(os, isBinary, d);

        writeFieldHeader(os, "rho", cloud.size(), 1, "float");
        writeFuns::write(os, isBinary, rho);

        writeFieldHeader(os, "cD", cloud.size(), 1, "float");
        writeFuns::write(os, isBinary, cD);

    }
    
    // osCsv << "x,y,z,U,active,id,age,stressHistory,pas,stressRateHistory,tau,fluidVelocity,Rep,d,rho,cD" << endl;
    // forAll(cD, i) {
    //     os  << fField[i];

    //         if (i > 0 && (i % 10) == 0)
    //         {
    //             os  << std::endl;
    //         }
    //         else
    //         {
    //             os  << ' ';
    //         }
    //     }
    //     os << std::endl;
    // }

    os.close();
}

struct timePair {
    label timeIndex;
    scalar deltaT;
    instant inputFileTime;
    instant simulationTime;
};

List<timePair> buildTimeList(const instantList & timeDirs, label repeat)
{
    label inputIterationCount = timeDirs.size();
    label totalIterationCount = repeat*(inputIterationCount - 1) + 1;
    scalar periodTime = timeDirs.last().value() - timeDirs.first().value();

    List<timePair> times(totalIterationCount);
    label timeIndex = 0;
    for(int rIt = 0; rIt < repeat; ++rIt) {
        // If we're repeating, only include the last time instant in the last cycle
        int numberOfTimesToInclude;
        if(rIt == (repeat-1)) 
            numberOfTimesToInclude = timeDirs.size();
        else
            numberOfTimesToInclude = timeDirs.size() - 1;

        for(int timeI = 0; timeI < numberOfTimesToInclude; ++timeI) {
            timePair currentTime;
            currentTime.timeIndex = timeIndex;
            currentTime.inputFileTime = timeDirs[timeI];
            currentTime.simulationTime = instant(rIt*periodTime + timeDirs[timeI].value() - timeDirs[0].value());

            if(timeIndex > 0)
                currentTime.deltaT = currentTime.simulationTime.value() - times[timeIndex-1].simulationTime.value();
            else
                currentTime.deltaT = timeDirs[timeI+1].value() - timeDirs[timeI].value();

            times[timeIndex] = currentTime;
            timeIndex += 1;
        }
    }
    return times;
}

int main(int argc, char *argv[])
{
    // Add command line parameters
    Foam::argList::addOption
    (
        "subCycle",
        "Integer",
        "Number of subcycles for each timestep"
    );

    Foam::argList::addOption
    (
        "cloudName",
        "word",
        "Name of the cloud"
    );

    Foam::argList::addOption
    (
        "repeat",
        "Integer",
        "Number of times to repeat the input data (in time)"
    );

    Foam::timeSelector::addOptions();
    //#include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    #include "createMesh.H"

    int substeps = args.optionLookupOrDefault<int>("subCycle", 1);
    int repeat = args.optionLookupOrDefault<int>("repeat", 1);
    word cloudName = args.optionLookupOrDefault<word>("cloudName", "plateletCloud");

    #include "../createFields.H"

    volScalarField mu("mu", rho*twoPhaseProperties.nu());
    volSymmTensorField tau("tau", mu * dev(twoSymm(fvc::grad(U))));
    basicPlateletCloud plateletCloud(cloudName, rho, U, mu, g, tau);

    List<timePair> timeInstants = buildTimeList(timeDirs, repeat);
    /*forAllConstIter(List<timePair>, timeInstants, iter) {
        Info << iter->timeIndex << ": " << iter->simulationTime.name() << " -> " << iter->inputFileTime.name() << endl;
    }*/

   if(substeps == 1) {
       forAllConstIter(List<timePair>, timeInstants, iter) {
            runTime.setTime(iter->simulationTime, iter->timeIndex);
            runTime.setDeltaT(iter->deltaT, false);
            Info<< "Time = " << runTime.timeName() << ", deltaT = " << runTime.deltaTValue() << endl;

            if(iter->timeIndex > 0) {
                // Read fields
                U = volVectorField(IOobject(U.name(), iter->inputFileTime.name(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
                p_rgh = volScalarField(IOobject(p_rgh.name(), iter->inputFileTime.name(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
                alpha1 = volScalarField(IOobject(alpha1.name(), iter->inputFileTime.name(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
                phi = surfaceScalarField(IOobject(phi.name(), iter->inputFileTime.name(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), fvc::flux(U));

                // Evaluate dependent fields
                alpha2 = 1. - alpha1;
                rho = alpha1*rho1 + alpha2*rho2;
                p = p_rgh + rho*gh;
            }

            twoPhaseProperties.correct();
            mu = rho * twoPhaseProperties.nu();
            tau = mu * dev(twoSymm(fvc::grad(U)));

            plateletCloud.evolve();
            writeCloudAsVTK(plateletCloud, std::to_string(iter->timeIndex));
        }
    } else {
        FieldInterpolator<volVectorField> uInterp(mesh, U.name());
        FieldInterpolator<volScalarField> p_rghInterp(mesh, p_rgh.name());
        FieldInterpolator<volScalarField> alphaInterp(mesh, alpha1.name());
        FieldInterpolator<surfaceScalarField> phiInterp(mesh, phi.name());
        timePair tLast;

        forAllConstIter(List<timePair>, timeInstants, iter) {
            runTime.setTime(iter->simulationTime, iter->timeIndex);
            Info<< "Time = " << runTime.timeName() << ", deltaT = " << runTime.deltaTValue() << endl;

            // Read data
            uInterp.readTime(iter->inputFileTime.name());
            p_rghInterp.readTime(iter->inputFileTime.name());
            alphaInterp.readTime(iter->inputFileTime.name());
            phiInterp.readTime(iter->inputFileTime.name(), fvc::flux(U));

            // Throw exceptions if IO errors occured
            FatalIOError.throwExceptions();

            if(iter->timeIndex > 0) {
                runTime.setDeltaT(iter->deltaT, false);	// The false flag here prevents the runtime from adjusting the timestep we set

                // Excecute subcycles
                for(subCycleTime subIter = subCycleTime(runTime, substeps); !(++subIter).end(); ) {
                    scalar timeStepFraction = (runTime.value() - tLast.simulationTime.value()) / iter->deltaT;
                    Info << "Executing subcycle, current time = " << runTime.timeName() << ", deltaT = " 
                        << runTime.deltaTValue() << ", fraction of timestep  = " << timeStepFraction << endl;

                    // Interpolate fields
                    U = uInterp.interpolate(timeStepFraction, iter->deltaT);
                    p_rgh = p_rghInterp.interpolate(timeStepFraction, iter->deltaT);
                    alpha1 = alphaInterp.interpolate(timeStepFraction, iter->deltaT);
                    phi = phiInterp.interpolate(timeStepFraction, iter->deltaT);

                    // Correct boundary conditions after interpolation
                    U.correctBoundaryConditions();
                    p_rgh.correctBoundaryConditions();
                    alpha1.correctBoundaryConditions();

                    // Evaluate dependent fields
                    alpha2 = 1. - alpha1;
                    rho = alpha1*rho1 + alpha2*rho2;
                    p = p_rgh + rho*gh;

                    // Update viscosity and shear stress
                    twoPhaseProperties.correct();
                    mu = rho * twoPhaseProperties.nu();
                    tau = mu * dev(twoSymm(fvc::grad(U)));

                    // Update platelet cloud
                    plateletCloud.evolve();
                }
            }
            tLast = *iter;
            writeCloudAsVTK(plateletCloud, std::to_string(iter->timeIndex));
        }
    }

    Info<< "End\n" << endl;
}