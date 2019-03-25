#ifndef functionObjects_plateletLpt_H
#define functionObjects_plateletLpt_H

#include "fvMeshFunctionObject.H"
#include "cloud/basicPlateletCloud.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class vfDependentViscosityTwoPhaseMixture;

namespace functionObjects
{

class plateletLpt
:
    public fvMeshFunctionObject
{
    // Private member data

        //- Optional acceleration due to gravity
        uniformDimensionedVectorField g_;

        //- Reference to the laminar transport model
        const vfDependentViscosityTwoPhaseMixture & transportModel_;

        //- Density field needed to construct plateletCloud
        const volScalarField & rho_;

        //- Reference to the velocity field
        const volVectorField& U_;

        //- Dynamic viscosity field needed to construct kinematicCloud
        volScalarField mu_;

		//- Shear stress field
		volSymmTensorField tau_;

        //- Name of the cloud
        word plateletCloudName_;

        //- The kinematic cloud
        basicPlateletCloud plateletCloud_;


    // Private member functions

        //- Disallow default bitwise copy construct
        plateletLpt(const plateletLpt&);

        //- Disallow default bitwise assignment
        void operator==(const plateletLpt&);


public:

    //- Runtime type information
    TypeName("plateletLpt");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        plateletLpt
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~plateletLpt();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Track the cloud
        virtual bool execute();

        //- Write the cloud
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
