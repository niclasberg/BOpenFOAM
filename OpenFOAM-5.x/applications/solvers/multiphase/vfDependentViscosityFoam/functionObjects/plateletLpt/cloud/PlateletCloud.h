/*---------------------------------------------------------------------------*\

Author: Niclas Berg (niber@mech.kth.se)

Class
    Foam::PlateletCloud

Description
    Adds platelet activation modelling to KinematicCloud

SourceFiles
    PlateletCloudI.H
    PlateletCloud.C

\*---------------------------------------------------------------------------*/

#ifndef PlateletCloud_H
#define PlateletCloud_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class PlateletCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PlateletCloud
:
    public CloudType
{
public:

	//- Runtime type information
        //TypeName("PlateletCloud");

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef PlateletCloud<CloudType> plateletCloudType;

private:

    // Private data

        //- Cloud copy pointer
        autoPtr<PlateletCloud<CloudType> > cloudCopyPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        PlateletCloud(const PlateletCloud&);

        //- Disallow default bitwise assignment
        void operator=(const PlateletCloud&);


protected:

    // Protected data

		 // References to the carrier gas fields
		const volSymmTensorField& tau_;

        //- Parcel constant properties
        typename parcelType::constantProperties constProps_;

            //- Reset state of cloud
            void cloudReset(PlateletCloud<CloudType>& c);

public:

	// Constructors

        //- Construct given carrier gas fields
        PlateletCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
		    const volScalarField& mu,
            const dimensionedVector& g,
			const volSymmTensorField& tau,
            bool readFields = true
        );

        //- Copy constructor with new name
        PlateletCloud(PlateletCloud<CloudType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        PlateletCloud
        (
            const fvMesh& mesh,
            const word& name,
            const PlateletCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > clone(const word& name)
        {
            return autoPtr<Cloud<parcelType> >
            (
                new PlateletCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType> >
            (
                new PlateletCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~PlateletCloud() { }


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const PlateletCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;

            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();

			//- Return access to the shear stress field
			const volSymmTensorField& tau() const { return tau_; }

		// Cloud evolution

			//- Evolve the cloud
	        void evolve();

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

        // I-O

            //- Print cloud information
            void info();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PlateletCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PlateletCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
