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
    Foam::viscosityModels::CassonC

Description
    An incompressible Casson non-Newtonian viscosity model Coupled (C) to the
    suspended phase (alpha1).

SourceFiles
    CassonC.C

\*---------------------------------------------------------------------------*/

#ifndef CassonC_H
#define CassonC_H

#include "viscosityModelC.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class CassonC Declaration
\*---------------------------------------------------------------------------*/

class CassonC
:
    public viscosityModelC
{
    // Private data

        dictionary CassonCCoeffs_;

        dimensionedScalar mup_;
        dimensionedScalar muMax_;
        dimensionedScalar A_;
        dimensionedScalar B_;


public:

    //- Runtime type information
    TypeName("CassonC");


    // Constructors

        //- Construct from components
        CassonC
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& alpha1
        );


    //- Destructor
    ~CassonC()
    {}


    // Member Functions

        //- Correct the laminar viscosity
        virtual tmp<volScalarField> calcMu(const volScalarField & alpha, const volVectorField & U) const;

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);

        tmp<volScalarField> dMuDalpha() const override
        {
            notImplemented("CassonC::dMuDalpha()");
        }

        tmp<volScalarField> dMuDgamma() const override
        {
            notImplemented("CassonC::dMuDgamma()");
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
