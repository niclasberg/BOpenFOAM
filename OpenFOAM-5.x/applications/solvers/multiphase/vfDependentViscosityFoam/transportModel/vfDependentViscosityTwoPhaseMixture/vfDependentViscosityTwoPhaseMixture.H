#ifndef vfDependentViscosityTwoPhaseMixture_H
#define vfDependentViscosityTwoPhaseMixture_H

#include "transportModel.H"
#include "viscosityModelC.H"
#include "diffusionModel.H"
#include "twoPhaseMixture.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class vfDependentViscosityTwoPhaseMixture
:
    public IOdictionary,
    public transportModel,
    public twoPhaseMixture
{
protected:

    // Protected data

        autoPtr<viscosityModelC> muModel_;
        autoPtr<diffusionModel> diffusionModel_;

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        volScalarField nu_;


    // Protected Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();


public:

    TypeName("vfDependentViscosityTwoPhaseMixture");


    // Constructors

        //- Construct from components
    	vfDependentViscosityTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~vfDependentViscosityTwoPhaseMixture()
    {}


    // Member Functions

        //- Return const-access to phase1 viscosityModel
        const viscosityModelC& muModel() const
        {
            return muModel_();
        }

		const diffusionModel & diffModel() const
		{
			return diffusionModel_();
		}

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const
		{
        	return muModel_->mu();
		}

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated kinematic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //- Correct the laminar viscosity
        virtual void correct()
        {
            calcNu();
            diffusionModel_->correct(alpha1_);
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
