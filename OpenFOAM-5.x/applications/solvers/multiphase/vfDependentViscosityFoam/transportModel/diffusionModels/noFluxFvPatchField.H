#ifndef NOFLUXBOUNDARYCONDITION_H_
#define NOFLUXBOUNDARYCONDITION_H_
#include "mixedFvPatchFields.H"

namespace Foam {

class noFluxFvPatchField 
: public mixedFvPatchScalarField
{
public:
    TypeName("noFlux");

    noFluxFvPatchField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    noFluxFvPatchField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given inletOutletFvPatchField onto a new patch
    noFluxFvPatchField
    (
        const noFluxFvPatchField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    noFluxFvPatchField
    (
        const noFluxFvPatchField &
    );

    //- Construct and return a clone
    virtual tmp<fvPatchField<scalar>> clone() const
    {
        return tmp<fvPatchField<scalar>>
        (
            new noFluxFvPatchField(*this)
        );
    }

    //- Construct as copy setting internal field reference
    noFluxFvPatchField
    (
        const noFluxFvPatchField &,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<scalar>> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchField<scalar>>
        (
            new noFluxFvPatchField(*this, iF)
        );
    }


// Member functions

    // Attributes

        //- Return true: this patch field is altered by assignment
        virtual bool assignable() const
        {
            return true;
        }


    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
    
private:

};

}


#endif