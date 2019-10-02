#include "addToRunTimeSelectionTable.H"
#include "vorticitySources.H"
#include "fvCFD.H"
#include "turbulentTransportModel.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(vorticitySources, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        vorticitySources,
        dictionary
    );
}
}

Foam::functionObjects::vorticitySources::vorticitySources
(
    const word & name,
    const Time& runTime,
    const dictionary& dict
)
: 
    fvMeshFunctionObject(name, runTime, dict),
    writeLocalObjects(obr_, log),
    uName_(dict.lookupOrDefault<word>("U", "U")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    pName_(dict.lookupOrDefault<word>("p", "p"))
{
    // Allocate arrays
    wordList typeNames;
    typeNames.append(vorticityArrayName);
    mesh_.objectRegistry::store(
        new volVectorField
        (
            IOobject
            (
                vorticityArrayName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("0", dimless/dimTime, Zero)
        )
    );

    typeNames.append(ddtVorticityArrayName);
    typeNames.append(convectiveTermName);
    typeNames.append(vortexStretchingTermName);
    typeNames.append(divTermName);
    typeNames.append(baroclinicTermName);
    typeNames.append(diffusiveFluxTermName);

    // Allocate the rest of the arrays
    // All of these arrays have the same dimension
    // Loop starts at i=1 (we skip recreating the vorticity array)
    for(int i = 1; i < typeNames.size(); ++i) {
        mesh_.objectRegistry::store(
            new volVectorField
            (
                IOobject
                (
                    typeNames[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("0", dimless/sqr(dimTime), Zero)
            )
        );
    }

    read(dict);
    resetLocalObjectNames(typeNames);
}

bool Foam::functionObjects::vorticitySources::read(const dictionary & dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    uName_ = dict.lookupOrDefault<word>("U", "U");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    pName_ = dict.lookupOrDefault<word>("p", "p");

    return true;
}

bool Foam::functionObjects::vorticitySources::execute() 
{
    // Lookup fields
    typedef incompressible::turbulenceModel icoTurb;
    const volVectorField & U = mesh_.lookupObject<volVectorField>(uName_);
    const volScalarField & p = mesh_.lookupObject<volScalarField>(pName_);
    const volScalarField & rho = mesh_.lookupObject<volScalarField>(rhoName_);
    const icoTurb & turb = mesh_.lookupObject<icoTurb>(turbulenceModel::propertiesName);

    // Lookup output arrays
    volVectorField & vorticity = mesh_.lookupObjectRef<volVectorField>(vorticityArrayName);
    volVectorField & ddtVorticity = mesh_.lookupObjectRef<volVectorField>(ddtVorticityArrayName);

    // Compute vorticity for current and last timestep
    vorticity = fvc::curl(U); 
    vorticity.oldTime() = fvc::curl(U.oldTime());
    ddtVorticity = fvc::ddt(vorticity);

    // Evaluate the terms in the vorticity equation
    mesh_.lookupObjectRef<volVectorField>(convectiveTermName) = U & fvc::grad(vorticity);
    mesh_.lookupObjectRef<volVectorField>(vortexStretchingTermName) = vorticity & fvc::grad(U);
    mesh_.lookupObjectRef<volVectorField>(divTermName) = fvc::div(U) * vorticity;
    mesh_.lookupObjectRef<volVectorField>(baroclinicTermName) = (fvc::grad(rho) ^ fvc::grad(p)) / sqr(rho);
    mesh_.lookupObjectRef<volVectorField>(diffusiveFluxTermName) = fvc::curl(fvc::div(turb.devReff()));

    return true;
}

bool Foam::functionObjects::vorticitySources::write()
{
    Info << "WRITE" << endl;
    writeLocalObjects::write();
    return true;
}