#ifndef ACTIVATION_MODEL_H_
#define ACTIVATION_MODEL_H_

namespace Foam {

template<class CloudType>
class activationModel {
public:
    TypeName("activationModel");
    
    enum StressNorm { 
        Frobenius, 
        vonMises
    };

    declareRunTimeSelectionTable
    (
        autoPtr,
        activationModel,
        dictionary,
        ( const dictionary& dict ),
        (dict)
    );

    // Virtual constructor
    static autoPtr<activationModel> New(const dictionary& dict );

    activationModel(const dictionary &);
    virtual ~activationModel() {}

    virtual void calc(typename CloudType::parcelType &, scalar dt) const = 0;
    virtual autoPtr<activationModel> clone() const = 0;

    scalar stressNorm(const symmTensor &) const;
private:
    static StressNorm wordToStressNorm(const word &);

    StressNorm stressNormType_;
};

}

#include "activationModel.C"

#define makeParticleActivationModel(CloudType)                                      \
                                                                               \
    typedef Foam::CloudType::plateletCloudType plateletCloudType;            \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::activationModel<plateletCloudType>, 0);                          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            activationModel<plateletCloudType>,                                 \
            dictionary                                                         \
        );                                                                     \
    }


#define makeParticleActivationModelType(SS, CloudType)                              \
                                                                               \
    typedef Foam::CloudType::plateletCloudType plateletCloudType;            \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<plateletCloudType>, 0);      \
                                                                               \
    Foam::activationModel<plateletCloudType>::                                  \
        adddictionaryConstructorToTable<Foam::SS<plateletCloudType>>          \
            add##SS##CloudType##plateletCloudType##ConstructorToTable_;

#endif