#ifndef wallShearStressGradient_H
#define wallShearStressGradient_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class wallShearStress Declaration
\*---------------------------------------------------------------------------*/

class wallShearStressGradient
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{

protected:
    static const constexpr char * WssArrayName = "wallShearStress";
    static const constexpr char * WssDdtArrayName = "ddt(wallShearStress)";


    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the shear-stress
        void projectShearStress
        (
            const volSymmTensorField& Reff,
            volVectorField& shearStress
        );


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        wallShearStressGradient(const wallShearStressGradient&);

        //- Disallow default bitwise assignment
        void operator=(const wallShearStressGradient&);


public:

    //- Runtime type information
    TypeName("wallShearStressGradient");


    // Constructors

        //- Construct from Time and dictionary
        wallShearStressGradient
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~wallShearStressGradient();


    // Member Functions

        //- Read the wallShearStress data
        virtual bool read(const dictionary&);

        //- Calculate the wall shear-stress
        virtual bool execute();

        //- Write the wall shear-stress
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
