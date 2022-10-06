#include "meshQualityGradient.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"




namespace Foam
{

namespace functionObjects
{
    defineTypeNameAndDebug(meshQualityGradient,0);

    addToRunTimeSelectionTable
    (
        functionObject,
        meshQualityGradient,
        dictionary
    );

}
}


Foam::functionObjects::meshQualityGradient::meshQualityGradient
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}

Foam::functionObjects::meshQualityGradient::~meshQualityGradient()
{}

bool Foam::functionObjects::meshQualityGradient::read(const dictionary& dict)
{
	fvMeshFunctionObject::read(dict);
	return true;
}

bool Foam::functionObjects::meshQualityGradient::execute()
{
    Info << type() << " write:" << endl;

    volScalarField V
    (
        IOobject
        (
            mesh_.V().name(),
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(mesh_.V().name(), mesh_.V().dimensions(), 0),
        calculatedFvPatchField<scalar>::typeName
    );

    V.ref() = mesh_.V();
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    volScalarField lengthScale = mag(U)/max(mag(fvc::grad(U)), dimensionedScalar("small", dimensionSet(0, 0, -1, 0, 0, 0, 0), VSMALL));
    //volScalarField quality = lengthScale / max(pow(V, 1.0/3.0), dimensionedScalar("small", dimLength, SMALL));
    //quality.write();
    volScalarField quality
    (
     	IOobject
	(
	 	"quality",
		mesh_.time().timeName(),
		mesh_,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	lengthScale / max(pow(V, 1.0/3.0), dimensionedScalar("small", dimLength, VSMALL))
    );
    quality.write();
    return true;

}

bool Foam::functionObjects::meshQualityGradient::write()
{
    return true;
}
