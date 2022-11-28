#ifndef passiveScalar_H
#define passiveScalar_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "fvOptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class scalarTransport Declaration
\*---------------------------------------------------------------------------*/

class passiveScalar
:
    public fvMeshFunctionObject
{
    // Private data

        //- Name of field to process
        word fieldName_;

        //- Name of flux field (optional)
        word phiName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Diffusion coefficient (optional)
        scalar D_;

        //- Flag to indicate whether a constant, uniform D_ is specified
        bool constantD_;

        //- Laminar diffusion coefficient (optional)
        scalar alphaD_;

        //- Turbulent diffusion coefficient (optional)
        scalar alphaDt_;

        //- Number of corrector iterations (optional)
        label nCorr_;

        //- Name of field whose schemes are used (optional)
        word schemesField_;

        //- Run-time selectable finite volume options, e.g. sources, constraints
        fv::optionList fvOptions_;

        //- The scalar field
        volScalarField s_;

		//- Diffusion model
		diffusionModel DModel_;


    // Private Member Functions

        //- Return the diffusivity field
        tmp<volScalarField> D(const surfaceScalarField& phi) const;

        //- Disallow default bitwise copy construct
        passiveScalar(const passiveScalar&);

        //- Disallow default bitwise assignment
        void operator=(const passiveScalar&);


public:

    //- Runtime type information
    TypeName("passiveScalar");


    // Constructors

        //- Construct from Time and dictionary
        passiveScalar
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~passiveScalar();


    // Member Functions

        //- Read the scalarTransport data
        virtual bool read(const dictionary&);

        //- Calculate the scalarTransport
        virtual bool execute();

        //- Do nothing.
        //  The volScalarField is registered and written automatically
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
