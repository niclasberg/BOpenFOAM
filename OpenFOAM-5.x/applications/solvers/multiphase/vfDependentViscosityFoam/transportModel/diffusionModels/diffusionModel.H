#ifndef diffusionModel_H
#define diffusionModel_H

#include "dictionary.H"
#include "volFields.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "fvMatricesFwd.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class diffusionModel
{

protected:

    // Protected data

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        //volScalarField mu_;

    // Private Member Functions

        //- Disallow copy construct
        diffusionModel(const diffusionModel&);

        //- Disallow default bitwise assignment
        void operator=(const diffusionModel&);


public:

    //- Runtime type information
    TypeName("diffusionModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            diffusionModel,
            diffusion,
            (
                const dictionary& diffusionProperties,
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (diffusionProperties, U, phi)
        );


    // Selectors

        static autoPtr<diffusionModel> New
        (
            const dictionary& diffusionProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        diffusionModel
        (
            //const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~diffusionModel()
    {}


    // Member Functions

        //- Correct the diffusion model
        virtual void correct() = 0;

		//- Get the flux
		virtual tmp<volVectorField> flux(const volScalarField &) const = 0;

		//- Divergence of flux
		virtual tmp<fvScalarMatrix> divFlux(const volScalarField &) const = 0;

        //- Return the strain rate
        tmp<volScalarField> shearRate() const;
};

} // End namespace Foam

#endif
