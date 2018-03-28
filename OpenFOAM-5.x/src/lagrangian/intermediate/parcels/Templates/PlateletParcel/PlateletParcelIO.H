#include "CollidingParcel.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::PlateletParcel<ParcelType>::propertyList_ =
    Foam::PlateletParcel<ParcelType>::propertyList();


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::PlateletParcel<ParcelType>::PlateletParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
	mechanicalDose_(0.0),
	pas_(0.0),
    tau_(symmTensor::zero),
	tau_last_(symmTensor::zero)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> mechanicalDose_;
            is >> pas_;
            is >> tau_;
			is >> tau_last_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&mechanicalDose_),
              + sizeof(mechanicalDose_)
              + sizeof(pas_)
              + sizeof(tau_)
			  + sizeof(tau_last_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "PlateletParcel<ParcelType>::Collisions"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::PlateletParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);

    IOField<scalar> mechanicalDose(c.fieldIOobject("mechanicalDose", IOobject::MUST_READ));
    c.checkFieldIOobject(c, mechanicalDose);

    IOField<scalar> pas(c.fieldIOobject("pas", IOobject::MUST_READ));
    c.checkFieldIOobject(c, pas);

    IOField<symmTensor> tau(c.fieldIOobject("tau", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tau);

    IOField<symmTensor> tau_last(c.fieldIOobject("tau_last", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tau_last);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        PlateletParcel<ParcelType>& p = iter();

		p.mechanicalDose_ = mechanicalDose[i];
        p.pas_ = pas[i];
        p.tau_ = tau[i];
		p.tau_last_ = tau_last[i];

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::PlateletParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np =  c.size();

    IOField<scalar> mechanicalDose(c.fieldIOobject("mechanicalDose", IOobject::NO_READ), np);
    IOField<scalar> pas(c.fieldIOobject("pas", IOobject::NO_READ), np);
    IOField<symmTensor> tau(c.fieldIOobject("tau", IOobject::NO_READ), np);
    IOField<symmTensor> tau_last(c.fieldIOobject("tau_last", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const PlateletParcel<ParcelType>& p = iter();

		mechanicalDose[i] = p.mechanicalDose_;
        pas[i] = p.pas_;
        tau[i] = p.tau_;
        tau_last[i] = p.tau_last_;

        i++;
    }

	const bool valid = (np > 0);

    mechanicalDose.write(valid);
    pas.write(valid);
    tau.write(valid);
	tau_last.write(valid);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PlateletParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.mechanicalDose()
            << token::SPACE << p.pas()
            << token::SPACE << p.tau_
            << token::SPACE << p.tau_last_;
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.mechanicalDose_),
            sizeof(p.mechanicalDose())
          + sizeof(p.pas())
          + sizeof(p.tau_)
          + sizeof(p.tau_last_)
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const PlateletParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
