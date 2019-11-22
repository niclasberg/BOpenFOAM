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
    Foam::viscosityModels::WalburnSchneckC

Description
    An incompressible WalburnSchneck non-Newtonian viscosity model Coupled (C) to the
    suspended phase (alpha1).

SourceFiles
    WalburnSchneckC.C

\*---------------------------------------------------------------------------*/

#ifndef WalburnSchneckC_H
#define WalburnSchneckC_H

#include "../viscosityModelC/viscosityModelC.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class WalburnSchneckC Declaration
\*---------------------------------------------------------------------------*/

class WalburnSchneckC
:
    public viscosityModelC
{
    // Private data

        dictionary WalburnSchneckCCoeffs_;

        dimensionedScalar C1_;
        dimensionedScalar TPMA_;
        dimensionedScalar muMin_;
        dimensionedScalar muMax_;

public:

    //- Runtime type information
    TypeName("WalburnSchneckC");


    // Constructors

        //- Construct from components
        WalburnSchneckC
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& alpha1
        );


    //- Destructor
    ~WalburnSchneckC()
    {}


    // Member Functions

        //- Correct the laminar viscosity
        virtual tmp<volScalarField> calcMu(const volScalarField & alpha, const volVectorField & U) const;

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);

        tmp<volScalarField> dMuDalpha() const override
        {
            notImplemented("WalburnSchneckC::dMuDalpha()");
        }

        tmp<volScalarField> dMuDgamma() const override
        {
            notImplemented("WalburnSchneckC::dMuDgamma()");
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
