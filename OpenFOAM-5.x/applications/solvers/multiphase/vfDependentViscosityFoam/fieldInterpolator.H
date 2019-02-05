template<class FieldType>
class FieldInterpolator {
public:
	FieldInterpolator(fvMesh & mesh, word fieldName) 
	: mesh_(mesh), fieldName_(fieldName), f0_(0), df0dt_(0), f1_(0), df1dt_(0) 
	{ 

	}

	~FieldInterpolator() 
	{
		if(f0_) delete f0_;
		if(f1_) delete f1_;
		if(df0dt_) delete df0dt_;
		if(df1dt_) delete df1dt_;
	}	

	word fieldName() const { return fieldName_;}
	fvMesh & mesh() { return mesh_; }
	const fvMesh & mesh() const { return mesh_; }
	FieldType & f0() { assert(f0_); return *f0_; }
	FieldType & f1() { assert(f1_); return *f1_; }
	FieldType & df0dt() { assert(df0dt_); return *df0dt_; }
	FieldType & df1dt() { assert(df1dt_); return *df1dt_; }

	void readTime(word timeName) 
	{
		// Read data
		FieldType * f = new FieldType
		(
			IOobject
			(
				fieldName(),
				timeName,
				mesh_,
				IOobject::MUST_READ,
				IOobject::NO_WRITE,
				false				// Don't register in the object database
			),
			mesh()
		);

		// Evaluate time derivative
		FieldType * dfdt = 0;
		/*if(f->nOldTimes() > 0) {
			df0dt_ = new FieldType(fvc::ddt(*f));
		}*/

		if( ! f0_) {
			// No data read
			f0_ = f;
			if(dfdt) df0dt_ = dfdt;
		} else if( ! f1_) {
			// Only one timestep read
			f1_ = f;
			if(dfdt) df1dt_ = dfdt;
		} else {
			// Two timesteps already read
			// Deallocate f0
			delete f0_;
			if(df0dt_) delete df0dt_;

			// Swap pointers
			f0_ = f1_;
			df0dt_ = df1dt_;

			// Assign new data
			f1_ = f;
			if(dfdt) df1dt_ = dfdt;
		}
	}

	// Interpolate the field
	tmp<FieldType> interpolate(scalar tsFrac, scalar deltaT) {
		scalar tsFracSqr = tsFrac*tsFrac;

		if(hasTimeDerivative()) {
			// Fourth order interpolation
			return (1. + tsFracSqr*(4. - 3.*tsFrac)) * f0() +
	 			   tsFrac * (1. - 4.*tsFrac + 3.*tsFracSqr)*deltaT * df0dt() +
				   tsFracSqr * (3.*tsFrac - 2.) * f1() +
				   6. * tsFracSqr * (1. - tsFrac) * deltaT * df1dt();
		} else {
			// Second order interpolation (linear)
			return (1. - tsFrac) * f0() + tsFrac * f1();
		}
	}

	bool hasTimeDerivative() const { return df0dt_ && df1dt_; }

private:
	// Disable copy construtor and assignment operator
	FieldInterpolator(const FieldInterpolator &);
	FieldInterpolator & operator=(const FieldInterpolator &);

	fvMesh & mesh_;
	word fieldName_;
	FieldType * f0_;
	FieldType * df0dt_;
	FieldType * f1_;
	FieldType * df1dt_;
};
