/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Namespace
    Foam::viscosityModels

Description
    A namespace for various incompressible viscosityModel implementations
    Coupled (C) to a suspended phase (alpha1).

Class
    Foam::viscosityModelC

Description
    An abstract base class for incompressible viscosityModels.

    The strain rate is defined by:

       sqrt(2) * mag(symm(grad(U)))


SourceFiles
    viscosityModelC.C
    viscosityModelNewC.C

\*---------------------------------------------------------------------------*/

#ifndef viscosityModelC_H
#define viscosityModelC_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class viscosityModelC Declaration
\*---------------------------------------------------------------------------*/

class viscosityModelC
{

protected:

    // Protected data

        word name_;
        dictionary viscosityProperties_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;
        const volScalarField& alpha1_;


    // Private Member Functions

        //- Disallow copy construct
        viscosityModelC(const viscosityModelC&);

        //- Disallow default bitwise assignment
        void operator=(const viscosityModelC&);


public:

    //- Runtime type information
    TypeName("viscosityModelC");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            viscosityModelC,
            dictionary,
            (
                const word& name,
                const dictionary& viscosityProperties,
                const volVectorField& U,
                const surfaceScalarField& phi,
    		const volScalarField& alpha1
            ),
            (name, viscosityProperties, U, phi, alpha1)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<viscosityModelC> New
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
    	    const volScalarField& alpha1
        );


    // Constructors

        //- Construct from components
        viscosityModelC
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
    	    const volScalarField& alpha1
        );


    //- Destructor
    virtual ~viscosityModelC()
    {}


    // Member Functions

        //- Return the phase transport properties dictionary
        const dictionary& viscosityProperties() const
        {
            return viscosityProperties_;
        }

        //- Return the strain rate
        tmp<volScalarField> strainRate() const;

        //- Return the Haematocrit (VF = Volume Fraction = alpha1)
        tmp<volScalarField> VF() const;

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const = 0;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const = 0;

        //- Correct the laminar viscosity
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
