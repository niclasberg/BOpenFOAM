#include "diffusionModel.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffusionModel, 0);
    defineRunTimeSelectionTable(diffusionModel, diffusion);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusionModel::diffusionModel
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const viscosityModelC & viscosityModel
) 
:
    U_(U),
    phi_(phi),
    viscosityModel_(viscosityModel)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diffusionModel::shearRate() const
{
    return sqrt(2.0)*mag(symm(fvc::grad(U_)));
}

const Foam::viscosityModelC & Foam::diffusionModel::viscosityModel() const
{
    return viscosityModel_;
}

/*bool Foam::diffusionModel::read(const dictionary& diffusionProperties)
{
    diffusionProperties_ = diffusionProperties;

    return true;
}*/

Foam::autoPtr<Foam::diffusionModel> Foam::diffusionModel::New
(
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const viscosityModelC & viscosityModel
)
{
    const word modelType(dict.lookup("diffusionModel"));
	dictionary subDict = dict.subDict(modelType + "Coeffs");

    diffusionConstructorTable::iterator cstrIter =
        diffusionConstructorTablePtr_->find(modelType);

    if (cstrIter == diffusionConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "diffusionModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown diffusionModel type "
            << modelType << nl << nl
            << "Valid diffusionModels are : " << endl
            << diffusionConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<diffusionModel>(cstrIter()(subDict, U, phi, viscosityModel));
}


// ************************************************************************* //
