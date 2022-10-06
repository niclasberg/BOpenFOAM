/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
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
    Foam::functionObjects::powerLawHemolysis

Group
    grpSolversFunctionObjects

Description
    Calculates the hemolysis state field using shear stress (history).

    - A constant diffusivity may be specified with the \c D entry,

    - Alternatively if a turbulence model is available a turbulent diffusivity
      may be constructed from the laminar and turbulent viscosities using the
      optional diffusivity coefficients \c alphaD and \c alphaDt (which default
      to 1):
      \verbatim
          D = alphaD*nu + alphaDt*nut
      \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    powerLawHemolysis.C
\*--------------------------------------------------------------------------*/

#ifndef functionObjects_powerLawHemolysis_H
#define functionObjects_powerLawHemolysis_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "fvOptionList.H"

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class powerLawHemolysis Declaration
\*---------------------------------------------------------------------------*/

class powerLawHemolysis : public fvMeshFunctionObject
{
    // Private stuff

    //- Model constants
    scalar C_, alpha_, beta_, criticalTau_, rhoValue_;

    //- Name of field to process
    word fieldName_;

    //- Number of corrector iterations (optional)
    label nCorr_;

    //- Runtime selectable finite volume options, e.g. sources, constraints
    fv::optionList fvOptions_;

    //- The scalar field
    volScalarField s_;

    //- Name of field whose schemes are used (optional)
    word schemesField_;

    //- Name of flux field (optional)
    word phiName_;

    //- Name of density field (optional)
    word rhoName_;

    //Private member functions

    //- Return diffusivity field

    //- Disallow default bitwise copy construct
    powerLawHemolysis(const powerLawHemolysis&);

    //- Disallow default bitwise assignment
    void operator=(const powerLawHemolysis&);

public:

    //- Runtime type information
    TypeName("powerLawHemolysis");

    //Constructors
    
    //- Construct from time and dictionary
    powerLawHemolysis
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    //- Destructor
    virtual ~powerLawHemolysis();

    // Member Functions

    //- Read the powerLawHemolysis data
    virtual bool read(const dictionary&);

    //- Solve the equations
    virtual bool execute();

    //- Do nothing. 
    //- The volScalarField is registered and written automatically
    //- WHAT DOES THIS MEAN???
    virtual bool write();
};

/*---------------------------------------------------------------------------*/

} // End namespace functionObjects
} // End namespace Foam

#endif

/*---------------------------------------------------------------------------*/
