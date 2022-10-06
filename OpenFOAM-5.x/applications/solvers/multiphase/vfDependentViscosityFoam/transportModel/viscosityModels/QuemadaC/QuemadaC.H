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

Class
    Foam::viscosityModels::QuemadaC

Description
    An incompressible Quemada non-Newtonian viscosity model Coupled (C) to the
    suspended phase (alpha1).

SourceFiles
    QuemadaC.C

\*---------------------------------------------------------------------------*/

#ifndef QuemadaC_H
#define QuemadaC_H

#include "viscosityModelC.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class QuemadaC Declaration
\*---------------------------------------------------------------------------*/

class QuemadaC
:
    public viscosityModelC
{
    // Private data

        dictionary QuemadaCCoeffs_;

        dimensionedScalar mup_;
        dimensionedScalar muMax_;

public:

    //- Runtime type information
    TypeName("QuemadaC");


    // Constructors

        //- Construct from components
        QuemadaC
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
	        const volScalarField& alpha1
        );


    //- Destructor
    ~QuemadaC()
    {}


    // Member Functions

        //- Correct the laminar viscosity
        virtual tmp<volScalarField> calcMu(const volScalarField & alpha, const volVectorField & U) const;

        tmp<volScalarField> dMuDalpha() const override;
        tmp<volScalarField> dMuDgamma() const override;

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
