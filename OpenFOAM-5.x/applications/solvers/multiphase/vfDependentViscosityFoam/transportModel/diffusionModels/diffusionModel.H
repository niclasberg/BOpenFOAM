#ifndef diffusionModel_H
#define diffusionModel_H

#include "dictionary.H"
#include "volFields.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "fvMatricesFwd.H"
#include "runTimeSelectionTables.H"
#include "viscosityModelC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class diffusionModel
{

protected:

    // Protected data

        const volVectorField& U_;
        const surfaceScalarField& phi_;
        const viscosityModelC & viscosityModel_;

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
                const surfaceScalarField& phi,
                const viscosityModelC & viscosityModel
            ),
            (diffusionProperties, U, phi, viscosityModel)
        );


    // Selectors

        static autoPtr<diffusionModel> New
        (
            const dictionary& diffusionProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const viscosityModelC & viscosityModel
        );


    // Constructors

        //- Construct from components
        diffusionModel
        (
            //const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const viscosityModelC & viscosityModel
        );


    //- Destructor
    virtual ~diffusionModel()
    {}


    // Member Functions

        //- Implementation functions for no flux boundary condition
        virtual scalarField noFluxBoundaryWeights(const fvPatch & patch) const = 0;

        //- Correct the diffusion model
        virtual void correct(const Foam::volScalarField & alpha) = 0;

		//- Get the flux
		virtual tmp<volVectorField> flux(const volScalarField &) const = 0;

		//- Divergence of flux
		virtual tmp<fvScalarMatrix> divFlux(const volScalarField &) const = 0;

        //- Return the strain rate
        tmp<volScalarField> shearRate() const;

        const viscosityModelC & viscosityModel() const;
};

} // End namespace Foam

#endif
