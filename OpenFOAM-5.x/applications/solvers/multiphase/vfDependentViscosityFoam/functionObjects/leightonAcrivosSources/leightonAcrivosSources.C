#include "addToRunTimeSelectionTable.H"
#include "leightonAcrivosSources.H"
#include "fvCFD.H"
#include "turbulentTransportModel.H"
#include "vfDependentViscosityTwoPhaseMixture.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(leightonAcrivosSources, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        leightonAcrivosSources,
        dictionary
    );
}
}

Foam::functionObjects::leightonAcrivosSources::leightonAcrivosSources
(
    const word & name,
    const Time& runTime,
    const dictionary& dict
)
: 
    fvMeshFunctionObject(name, runTime, dict),
    writeLocalObjects(obr_, log),
    uName_(dict.lookupOrDefault<word>("U", "U")),
    Db_("Db", dimLength*dimLength/dimTime, dict),
    a_("a", dimLength, dict),
    Kgamma_("Kgamma", dimless, dict),
    Kalpha_("Kalpha", dimless, dict),
    Kmu_("Kmu", dimless, dict)
{
    // Allocate arrays
    wordList typeNames;
    typeNames.append(alphaTermArrayName);
    typeNames.append(gammaTermArrayName);
    typeNames.append(muTermArrayName);

    for(int i = 0; i < typeNames.size(); ++i) {
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
                dimensionedVector("0", dimLength/dimTime, Zero)
            )
        );
    }

    read(dict);
    resetLocalObjectNames(typeNames);
}

bool Foam::functionObjects::leightonAcrivosSources::read(const dictionary & dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    uName_ = dict.lookupOrDefault<word>("U", "U");
    Db_ = dimensionedScalar("Db", dimLength*dimLength/dimTime, dict);
    a_ = dimensionedScalar("a", dimLength, dict);
    Kgamma_ = dimensionedScalar("Kgamma", dimless, dict);
    Kalpha_ = dimensionedScalar("Kalpha", dimless, dict);
    Kmu_ = dimensionedScalar("Kmu", dimless, dict);

    return true;
}

bool Foam::functionObjects::leightonAcrivosSources::execute() 
{
    // Lookup fields
    const volVectorField & U = mesh_.lookupObject<volVectorField>(uName_);
    const vfDependentViscosityTwoPhaseMixture & transport = 
        mesh_.lookupObject<vfDependentViscosityTwoPhaseMixture>("transportProperties");
    const viscosityModelC & viscModel = transport.muModel();
    const volScalarField & alpha = transport.alpha1();
    const volScalarField gamma = viscModel.strainRate(U);
 
    // Compute alpha-gradient dependent term
    mesh_.lookupObjectRef<volVectorField>(alphaTermArrayName) = 
        (Db_ + a_*a_*alpha*gamma*Kalpha_) * fvc::grad(alpha); 
    mesh_.lookupObjectRef<volVectorField>(gammaTermArrayName) = 
        a_*a_*alpha*alpha*Kgamma_*fvc::grad(gamma);
    mesh_.lookupObjectRef<volVectorField>(muTermArrayName) = 
        (Kmu_*a_*a_*gamma*alpha*alpha/viscModel.mu()) * (viscModel.dMuDgamma() * fvc::grad(gamma) + viscModel.dMuDalpha() * fvc::grad(alpha));

    return true;
}

bool Foam::functionObjects::leightonAcrivosSources::write()
{
    writeLocalObjects::write();
    return true;
}