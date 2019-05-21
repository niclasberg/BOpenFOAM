#ifndef WINDKESSELPRESSUREFVPATCHSCALARFIELD_H_
#define WINDKESSELPRESSUREFVPATCHSCALARFIELD_H_
#include "fixedValueFvPatchFields.H"

namespace Foam {

class WindkesselPressureFvPatchScalarField : public fixedValueFvPatchScalarField
{
public:
    TypeName("WindkesselPressure");

    // Object to store the current value and two previous values for flow rate and pressure
    class CircularTruple {
    public:
        CircularTruple(word name) : name_(name) { }
        CircularTruple(word name, const dictionary & dict) : name_(name)
        {
            wordList names = variableNames();
            scalar tmp;
            for(label i = 0; i < names.size(); ++i) {
                if(dict.readIfPresent(name, tmp)) {
                    values_.append(tmp);
                }
            }
        }

        wordList variableNames() const 
        {
            wordList names;
            names.append(name_);
            names.append(name_ + "0");
            names.append(name_ + "00");
            return names;
        }

        void write(Ostream & os) const
        {
            wordList names = variableNames();
            for(label i = 0; i < this->size(); ++i)
                os.writeKeyword(names[i]) << values_[i] << token::END_STATEMENT << nl;
        }

        scalar & operator[](label i) { return values_[i]; }
        const scalar & operator[](label i) const { return values_[i]; }

        scalar & currentValue() { return values_[0]; }
        const scalar & currentValue() const { return values_[0]; }
        scalar & oldValue() { return values_[1]; }
        const scalar & oldValue() const { return values_[1]; }
        scalar & oldOldValue() { return values_[2]; }
        const scalar & oldOldValue() const { return values_[2]; }

        void shiftAndInsert(const scalar & v)
        {   
            //Grow the list if element count is less than 3
            if(values_.size() < 3)
                values_.append(0);
            // Shift
            for(int i = (values_.size()-1); i > 0; --i)
                values_[i] = values_[i-1];
            // Insert
            values_[0] = v;
        }

        label size() const { return values_.size(); }
    private:
        word name_;
        scalarList values_;
    };

    // Constructors
    WindkesselPressureFvPatchScalarField(
        const fvPatch &, 
        const DimensionedField<scalar, volMesh> &
    );

    WindkesselPressureFvPatchScalarField(
        const fvPatch &, 
        const DimensionedField<scalar, volMesh> &,
        const dictionary &
    );

    WindkesselPressureFvPatchScalarField
    (
        const WindkesselPressureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    // Copy constructors
    WindkesselPressureFvPatchScalarField(
        const WindkesselPressureFvPatchScalarField &
    );

    WindkesselPressureFvPatchScalarField(
        const WindkesselPressureFvPatchScalarField &,
        const DimensionedField<scalar, volMesh> & iF
    );

    // Clone operators
    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>(
            new WindkesselPressureFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone(const DimensionedField<scalar, volMesh> & iF) const
    {
        return tmp<fvPatchScalarField>(
            new WindkesselPressureFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();
    virtual void write(Ostream &) const;

    // Getters for model parameters
    scalar & Rp() { return Rp_; }
    const scalar & Rp() const { return Rp_; }
    scalar & Rc() { return Rc_; }
    const scalar & Rc() const { return Rc_; }
    scalar & L() { return L_; }
    const scalar & L() const { return L_; }
    scalar & C() { return C_; }
    const scalar & C() const { return C_; }

protected:
    // Access to the velocity field
    const volVectorField & U() const {
        return this->db().lookupObject<volVectorField>(Uname_);
    }

    scalar computePatchFlowRate(const volVectorField &) const;

private:
    word Uname_;

    // Field values
    label lastTimeIndex_;
    CircularTruple qs_;
    CircularTruple ps_;
    bool hasDataAtT0_;
    bool hasDataAtT00_;

    // Model parameters
    scalar L_;
    scalar C_;
    scalar Rp_;
    scalar Rc_;
};

}

#endif