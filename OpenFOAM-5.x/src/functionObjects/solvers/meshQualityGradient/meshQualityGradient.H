#ifndef functionObjects_meshQualityGradient_H
#define functionObjects_meshQualityGradient_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
//#include "fieldExpression.H"

namespace Foam
{
namespace functionObjects
{
/*------------------------------------------------------------------------- *\
            Class meshQualityGradient Declaration
\*--------------------------------------------------------------------------*/

class meshQualityGradient : public fvMeshFunctionObject
{
    word outputFieldName_;


    meshQualityGradient(const meshQualityGradient&);

    void operator=(const meshQualityGradient&);


public:

    TypeName("meshQualityGradient");

    //Constructor

    meshQualityGradient
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    //Destructor
    virtual ~meshQualityGradient();

    //Member functions
    virtual bool read(const dictionary&);

    virtual bool execute();

    virtual bool write();

};



} //End namespace functionObjects

} //End namespace Foam





#endif
