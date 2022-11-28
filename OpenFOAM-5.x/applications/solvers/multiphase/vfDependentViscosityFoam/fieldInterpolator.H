template<class FieldType>
class FieldInterpolator {
public:
	FieldInterpolator(fvMesh & mesh, word fieldName) 
	: mesh_(mesh), fieldName_(fieldName), f0_(0), f1_(0) 
	{ 

	}

	~FieldInterpolator() 
	{
		if(f0_) delete f0_;
		if(f1_) delete f1_;
	}	

	word fieldName() const { return fieldName_;}
	fvMesh & mesh() { return mesh_; }
	const fvMesh & mesh() const { return mesh_; }
	FieldType & f0() { assert(f0_); return *f0_; }
	FieldType & f1() { assert(f1_); return *f1_; }

	void readTime(word timeName) 
	{
		// Read data
		setField(
			new FieldType(
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
			)
		);
	}

	void setField(FieldType * f)
	{
		if( ! f0_) {
			// No data read yet
			f0_ = f;
		} else if( ! f1_) {
			// Only one timestep read
			f1_ = f;
		} else {
			// Two timesteps already read
			// Deallocate f0
			delete f0_;

			// Swap pointers
			f0_ = f1_;
			f1_ = f;
		}
	}

	void readTime(word timeName, const FieldType & defaultValue) 
	{
		// Read data
		setField(
			new FieldType(
				IOobject
				(
					fieldName(),
					timeName,
					mesh_,
					IOobject::READ_IF_PRESENT,
					IOobject::NO_WRITE,
					false				// Don't register in the object database
				),
				defaultValue
			)
		);
	}

	// Interpolate the field
	tmp<FieldType> interpolate(scalar tsFrac, scalar deltaT) {
		return (1. - tsFrac) * f0() + tsFrac * f1();
	}

private:
	// Disable copy construtor and assignment operator
	FieldInterpolator(const FieldInterpolator &);
	FieldInterpolator & operator=(const FieldInterpolator &);

	fvMesh & mesh_;
	word fieldName_;
	FieldType * f0_;
	FieldType * f1_;
};
