#include "noFluxFvPatchField.H"
#include "vfDependentViscosityTwoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
    makePatchTypeField(
        fvPatchScalarField,
        noFluxFvPatchField
    );
}

Foam::noFluxFvPatchField::noFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::noFluxFvPatchField::noFluxFvPatchField
(
    const noFluxFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::noFluxFvPatchField::noFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
{
    fvPatchField<scalar>::operator=(Field<scalar>("value", dict, p.size()));

    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::noFluxFvPatchField::noFluxFvPatchField
(
    const noFluxFvPatchField& ptf
)
:
    mixedFvPatchScalarField(ptf)
{}


Foam::noFluxFvPatchField::noFluxFvPatchField
(
    const noFluxFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::noFluxFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Lookup the diffusion model
    const diffusionModel & diffModel = 
        this->db().lookupObject<vfDependentViscosityTwoPhaseMixture>("transportProperties").diffModel();
    
    this->valueFraction() = diffModel.noFluxBoundaryWeights(this->patch());

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::noFluxFvPatchField::write(Ostream& os) const
{
    // Only write the values at the boundary
    fvPatchField<scalar>::write(os);
    this->writeEntry("value", os);
}