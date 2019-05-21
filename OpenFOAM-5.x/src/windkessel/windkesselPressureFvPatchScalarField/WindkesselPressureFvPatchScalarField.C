#include "WindkesselPressureFvPatchScalarField.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"

Foam::WindkesselPressureFvPatchScalarField::WindkesselPressureFvPatchScalarField(
    const Foam::fvPatch & patch, 
    const Foam::DimensionedField<Foam::scalar, Foam::volMesh> & iF
) 
: 
    fixedValueFvPatchScalarField(patch, iF),
    Uname_("U"),
    lastTimeIndex_(this->db().time().timeIndex()),
    qs_("Q"),
    ps_("p"),
    L_(0),
    C_(0),
    Rp_(1),
    Rc_(0)
{
    // Set intial condition
    ps_.shiftAndInsert(0);
}

Foam::WindkesselPressureFvPatchScalarField::WindkesselPressureFvPatchScalarField(
    const Foam::fvPatch & patch, 
    const Foam::DimensionedField<Foam::scalar, Foam::volMesh> & iF,
    const Foam::dictionary & dict
)
: 
    fixedValueFvPatchScalarField(patch, iF, dict),
    Uname_(dict.lookupOrDefault<word>("U", "U")),
    lastTimeIndex_(this->db().time().timeIndex()), 
    qs_("Q", dict),
    ps_("p", dict),
    L_(dict.lookupOrDefault<scalar>("L", 0)),
    C_(dict.lookupOrDefault<scalar>("C", 0)),
    Rp_(dict.lookupOrDefault<scalar>("Rp", 1)),
    Rc_(dict.lookupOrDefault<scalar>("Rc", 0))
{
    // Evaluate initial condition if not present
    if(ps_.size() == 0) {
        // Compute average pressure over the boundary
        ps_.shiftAndInsert(
            gSum(this->patch().magSf() * this->patchInternalField()) 
            / gSum(this->patch().magSf()) 
        );
    }
}

Foam::WindkesselPressureFvPatchScalarField::WindkesselPressureFvPatchScalarField
(
    const Foam::WindkesselPressureFvPatchScalarField& rhs,
    const Foam::fvPatch& patch,
    const Foam::DimensionedField<scalar, volMesh>& iF,
    const Foam::fvPatchFieldMapper& mapper
)
: 
    fixedValueFvPatchScalarField(rhs, patch, iF, mapper),
    Uname_(rhs.Uname_),
    lastTimeIndex_(rhs.lastTimeIndex_), 
    qs_(rhs.qs_),
    ps_(rhs.ps_),
    L_(rhs.L_),
    C_(rhs.C_),
    Rp_(rhs.Rp_),
    Rc_(rhs.Rc_)
{

}

// Copy constructors
Foam::WindkesselPressureFvPatchScalarField::WindkesselPressureFvPatchScalarField(
    const Foam::WindkesselPressureFvPatchScalarField & rhs
)
: 
    fixedValueFvPatchScalarField(rhs),
    Uname_(rhs.Uname_),
    lastTimeIndex_(rhs.lastTimeIndex_), 
    qs_(rhs.qs_),
    ps_(rhs.ps_),
    L_(rhs.L_),
    C_(rhs.C_),
    Rp_(rhs.Rp_),
    Rc_(rhs.Rc_)
{

}

Foam::WindkesselPressureFvPatchScalarField::WindkesselPressureFvPatchScalarField(
    const Foam::WindkesselPressureFvPatchScalarField & rhs,
    const Foam::DimensionedField<Foam::scalar, Foam::volMesh> & iF
)
:
    fixedValueFvPatchScalarField(rhs, iF),
    Uname_(rhs.Uname_),
    lastTimeIndex_(rhs.lastTimeIndex_), 
    qs_(rhs.qs_),
    ps_(rhs.ps_),
    L_(rhs.L_),
    C_(rhs.C_),
    Rp_(rhs.Rp_),
    Rc_(rhs.Rc_)
{

}

void Foam::WindkesselPressureFvPatchScalarField::updateCoeffs()
{
    // The ODE has the form
    //  a0 dp/dt + p = a1 Q + a2 dQ/dt + a3 d2Q/dt2
    // Evaluate coefficients
    scalar a0 = Rp_ * C_;
    scalar a1 = Rc_ + Rp_;
    scalar a2 = C_ * Rc_ * Rp_;
    scalar a3 = L_ * C_ * Rp_;

    // Compute old flow rate if not yet computed
    if(qs_.size() == 0) {
        qs_.shiftAndInsert(this->computePatchFlowRate(this->U().oldTime()));
    }

    // Check if we're on a new time step
    label timeIndex = this->db().time().timeIndex();
    if(timeIndex != lastTimeIndex_) {
        // Add dummy values for the flow rate and pressure (will be replaced later)
        qs_.shiftAndInsert(0);
        ps_.shiftAndInsert(0);
        lastTimeIndex_ = timeIndex;
    }

    //Compute flow rate
    qs_.currentValue() = this->computePatchFlowRate(this->U());

    // Get current timestep
    scalar dt = this->db().time().deltaTValue();

    if(qs_.size() < 3) {
        // First timestep, integrate pressure with implicit Euler
        ps_.currentValue() = (a0*ps_.oldValue() + a1*dt*qs_.currentValue() + a2*(qs_.currentValue() - qs_.oldValue())) / (a0 + dt);
    } else {
        // Get old timestep
        scalar dt0 = this->db().time().deltaT0Value();

        // Evaluate derivative weights
        // The derivative approximations are backwards finite differences, expressed such that:
        //  df/dt^(n+1) = d*f^(n+1) + d0*f^(n) + d00*f^(n-1)
        //  d2f/dt2^(n+1) = dd*f^(n+1) + dd0*f^(n) + dd00*f^(n-1)
        // First order derivatives
        scalar d00 = dt / (dt0 * (dt + dt0));
        scalar d0 = -(dt + dt0) / (dt*dt0);
        scalar d = (2*dt + dt0) / (dt*(dt + dt0));

        // Second order derivatives
        scalar dd00 = 2. / (dt0 * (dt + dt0));
        scalar dd0 = -2. / (dt*dt0);
        scalar dd = 2. / (dt * (dt + dt0));

        // Integrate the equation for the pressure
        ps_.currentValue() = (
                    -a0*d0*ps_.oldValue() - a0*d00*ps_.oldOldValue() + 
                    (a1 + a2*d + a3*dd)*qs_.currentValue() + 
                    (a2*d0+a3*dd0)*qs_.oldValue() + 
                    (a2*d00 + a3*dd00)*qs_.oldOldValue()
            ) / (1 + a0*d);
    }

    //ps_.write(Info);
    //qs_.write(Info);

    // Apply values to the faces
    operator==(ps_.currentValue());

    fixedValueFvPatchField<scalar>::updateCoeffs();
}

Foam::scalar Foam::WindkesselPressureFvPatchScalarField::computePatchFlowRate(const volVectorField & U) const
{
    label patchId = this->patch().index();
    return gSum(this->U().boundaryField()[patchId] & this->patch().Sf());
}

void Foam::WindkesselPressureFvPatchScalarField::write(Foam::Ostream & os) const
{
    fvPatchScalarField::write(os);
    if (Uname_ != "U"){
        os.writeKeyword("U") << Uname_ << token::END_STATEMENT << nl;
    }
    this->writeEntryIfDifferent(os, "Rp", scalar(1), Rp_);
    this->writeEntryIfDifferent(os, "Rc", scalar(0), Rc_);
    this->writeEntryIfDifferent(os, "L", scalar(0), L_);
    this->writeEntryIfDifferent(os, "C", scalar(0), C_);
    qs_.write(os);
    ps_.write(os);
    this->writeEntry("value", os);
}

namespace Foam {
    makePatchTypeField(
        fvPatchScalarField,
        WindkesselPressureFvPatchScalarField
    );
}