/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "wallShearStressGradient.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "vfDependentViscosityTwoPhaseMixture.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallShearStressGradient, 0);
    addToRunTimeSelectionTable(functionObject, wallShearStressGradient, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::wallShearStressGradient::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Wall shear stress");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "wssMin");
    writeTabbed(file(), "wssMax");
    writeTabbed(file(), "ddtWssMin");
    writeTabbed(file(), "ddtWssMax");
    file() << endl;
}


void Foam::functionObjects::wallShearStressGradient::projectShearStress
(
    const volSymmTensorField& Reff,
    volVectorField& shearStress
)
{
    shearStress.dimensions().reset(Reff.dimensions());

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

        vectorField& ssp = shearStress.boundaryFieldRef()[patchi];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
        const symmTensorField& Reffp = Reff.boundaryField()[patchi];

        ssp = (-Sfp/magSfp) & Reffp;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallShearStressGradient::wallShearStressGradient
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    patchSet_()
{
    mesh_.objectRegistry::store(
        new volVectorField
        (
            IOobject
            (
                WssArrayName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "0",
                sqr(dimLength)/sqr(dimTime),
                Zero
            )
        )
    );
    mesh_.objectRegistry::store(
        new volVectorField
        (
            IOobject
            (
                WssDdtArrayName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "0",
                sqr(dimLength)/sqr(dimTime),
                Zero
            )
        )
    );

    wordList typeNames;
    typeNames.append(WssArrayName);
    typeNames.append(WssDdtArrayName);

    read(dict);
    resetName(typeName);
    resetLocalObjectNames(typeNames);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::functionObjects::wallShearStressGradient::~wallShearStressGradient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallShearStressGradient::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall shear stress on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::wallShearStressGradient::execute()
{
    // Get the transport model and viscosity field
    const vfDependentViscosityTwoPhaseMixture & transportModel(
        mesh_.lookupObject<vfDependentViscosityTwoPhaseMixture>("transportProperties")
    );

    // Lookup current and old U fields
    const volVectorField & U = mesh_.lookupObject<volVectorField>("U");
    const volVectorField & U0 = U.oldTime();
    const volScalarField & alpha = transportModel.alpha1();
    const volScalarField & alpha0 = alpha.oldTime();

    // Compute shear stress tensor at current and last timestep
    volSymmTensorField wss = transportModel.muModel().calcMu(alpha, U) * dev(twoSymm(fvc::grad(U)));
    wss.oldTime() = transportModel.muModel().calcMu(alpha0, U0) * dev(twoSymm(fvc::grad(U0)));
    const volSymmTensorField ddtWss = fvc::ddt(wss);

    projectShearStress(wss, mesh_.lookupObjectRef<volVectorField>(WssArrayName));
    projectShearStress(ddtWss, mesh_.lookupObjectRef<volVectorField>(WssDdtArrayName));

    return true;
}


bool Foam::functionObjects::wallShearStressGradient::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volVectorField& wallShearStress =
        obr_.lookupObject<volVectorField>(WssArrayName);
    const volVectorField& ddtWallShearStress =
        obr_.lookupObject<volVectorField>(WssDdtArrayName);

    const fvPatchList& patches = mesh_.boundary();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const vectorField& ssp = wallShearStress.boundaryField()[patchi];
        const vectorField& dssp = ddtWallShearStress.boundaryField()[patchi];

        vector minSsp = gMin(ssp);
        vector maxSsp = gMax(ssp);
        vector minDdtSsp = gMin(dssp);
        vector maxDdtSsp = gMax(dssp);

        if (Pstream::master())
        {
            file() << mesh_.time().value()
                << token::TAB << pp.name()
                << token::TAB << minSsp
                << token::TAB << maxSsp
                << token::TAB << minDdtSsp
                << token::TAB << maxDdtSsp
                << endl;
        }

        Log << "    min/max(" << pp.name() << ") = "
            << minSsp << ", " << maxSsp << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
