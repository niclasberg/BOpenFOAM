
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PlateletCloud<CloudType>::PlateletCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
	const volSymmTensorField& tau,
    bool readFields
)
:
    CloudType(cloudName, rho, U, mu, g, false),
	tau_(tau),
    constProps_(this->particleProperties())
{
    if (this->solution().steadyState())
    {
        FatalErrorIn
        (
            "Foam::PlateletCloud<CloudType>::PlateletCloud"
            "("
                "const word&, "
                "const volScalarField&, "
                "const volVectorField&, "
                "const volScalarField&, "
                "const dimensionedVector&, "
				"volSymmTensorField&,"
                "bool"
            ")"
        )   << "Platelet modelling not currently available for steady state "
            << "calculations" << exit(FatalError);
    }

    if (this->solution().active())
    {

        if (readFields)
        {
            parcelType::readFields(*this);
        }
    }
}


template<class CloudType>
Foam::PlateletCloud<CloudType>::PlateletCloud
(
    PlateletCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
	tau_(c.tau_)
{}


template<class CloudType>
Foam::PlateletCloud<CloudType>::PlateletCloud
(
    const fvMesh& mesh,
    const word& name,
    const PlateletCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    tau_(c.tau_)
{}

template<class CloudType>
void Foam::PlateletCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<PlateletCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}

template<class CloudType>
void Foam::PlateletCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}

template<class CloudType>
void Foam::PlateletCloud<CloudType>::cloudReset(PlateletCloud<CloudType>& c)
{
    CloudType::cloudReset(c);
}

template<class CloudType>
void Foam::PlateletCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::template
            TrackingData<PlateletCloud<CloudType> > td(*this);

        this->solve(td);
    }
}

template<class CloudType>
void Foam::PlateletCloud<CloudType>::info()
{
    CloudType::info();
}

