/*---------------------------------------------------------------------------*\
 ##   ####  ######     |
 ##  ##     ##         | Copyright: ICE Stroemungsfoschungs GmbH
 ##  ##     ####       |
 ##  ##     ##         | http://www.ice-sf.at
 ##   ####  ######     |
-------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Contributors/Copyright:
    2010-2014, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2013 Bruno Santos <wyldckat@gmail.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#include "manipulateField.H"

#include "FieldValueExpressionDriver.H"

namespace Foam {
    defineTypeNameAndDebug(manipulateField,0);
}

Foam::manipulateField::manipulateField
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    active_(true),
    writeManipulated_(false),
    obr_(obr),
    dict_(dict)
{
    if (!isA<fvMesh>(obr_))
    {
        active_=false;
        WarningIn("manipulateField::manipulateField")
            << "Not a fvMesh. Nothing I can do"
                << endl;
    }
    read(dict);
    bool writeAtStart=true;
    if(dict.found("manipulateAtStartup")) {
        writeAtStart=readBool(dict.lookup("manipulateAtStartup"));
    } else {
        WarningIn("Foam::manipulateField::manipulateField")
            << "'manipulateAtStartup' not set in " << dict.name() << nl
                << "Assuming: true"
                << endl;
    }
    if(writeAtStart) {
        write();
    }
}

Foam::manipulateField::~manipulateField()
{}

template<class TData,class TMask>
void Foam::manipulateField::manipulate(
    const TData &data,
    const TMask &mask,
    const word entity
)
{
    TData &original=const_cast<TData &>(obr_.lookupObject<TData>(name_));
    label cnt=0;
    forAll(original,cellI) {
        if(mask[cellI]>SMALL) {
            cnt++;
            original[cellI]=data[cellI];
        }
    }

    reduce(cnt,plusOp<label>());
    Info << "Manipulated field " << name_ << " in " << cnt
        << " " << entity << " with the expression " << expression_ << endl;
    original.correctBoundaryConditions();

    if(
        obr_.time().outputTime()
        &&
        original.writeOpt()==IOobject::AUTO_WRITE
    ) {
        if(this->writeManipulated_) {
            Info << "Rewriting manipulated field " << original.name() << endl;

            original.write();
        } else {
            Info << "Manipulated field " << original.name()
                << " not rewritten. Set 'writeManipulated'" << endl;
        }
    }
}

template<class TData,class TMask>
void Foam::manipulateField::manipulateSurface(
    const TData &data,
    const TMask &mask
)
{
    TData &original=const_cast<TData &>(obr_.lookupObject<TData>(name_));
    label cnt=0;
    forAll(original,cellI) {
        if(mask[cellI]>SMALL) {
            cnt++;
            original[cellI]=data[cellI];
        }
    }

    reduce(cnt,plusOp<label>());
    Info << "Manipulated field " << name_ << " on " << cnt
        << " faces with the expression " << expression_ << endl;

    // this does not work for surface fields
    //    original.correctBoundaryConditions();

    if(
        obr_.time().outputTime()
        &&
        original.writeOpt()==IOobject::AUTO_WRITE
    ) {
        if(this->writeManipulated_) {
            Info << "Rewriting manipulated field " << original.name() << endl;

            original.write();
        } else {
            Info << "Manipulated field " << original.name()
                << " not rewritten. Set 'writeManipulated'" << endl;
        }
    }
}

void Foam::manipulateField::timeSet()
{
    // Do nothing
}

void Foam::manipulateField::read(const dictionary& dict)
{
    if(active_) {
        name_=word(dict.lookup("fieldName"));
        expression_=exprString(
            dict.lookup("expression"),
            dict
        );
        maskExpression_=exprString(
            dict.lookup("mask"),
            dict
        );
        writeManipulated_=dict.lookupOrDefault<bool>("writeManipulated",false);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        driver_.set(
            new FieldValueExpressionDriver(
                mesh.time().timeName(),
                mesh.time(),
                mesh,
                false, // no caching. No need
                true,  // search fields in memory
                false  // don't look up files in memory
            )
        );

        driver_->readVariablesAndTables(dict_);

        // this might not work when rereading ... but what is consistent in that case?
        driver_->createWriterAndRead(name_+"_"+type());
    }
}

#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
bool
#else
void
#endif
Foam::manipulateField::write()
{
    if(active_) {
        FieldValueExpressionDriver &driver=driver_();

        driver.clearVariables();

        driver.parse(maskExpression_);

        if(!driver.isLogical()) {
            FatalErrorIn("manipulateField::execute()")
                << maskExpression_ << " does not evaluate to a logical expression"
                    << endl
                    << exit(FatalError);
        }

        if(driver.resultIsTyp<volScalarField>(true)) {
            volScalarField conditionField(driver.getResult<volScalarField>());

            driver.parse(expression_);

            if(driver.resultIsTyp<volVectorField>()) {
                manipulate(
                    driver.getResult<volVectorField>(),
                    conditionField
                );

            } else if(driver.resultIsTyp<volScalarField>()) {
                manipulate(
                    driver.getResult<volScalarField>(),
                    conditionField
                );
            } else if(driver.resultIsTyp<volTensorField>()) {
                manipulate(
                    driver.getResult<volTensorField>(),
                    conditionField
                );
            } else if(driver.resultIsTyp<volSymmTensorField>()) {
                manipulate(
                    driver.getResult<volSymmTensorField>(),
                    conditionField
                );
            } else if(driver.resultIsTyp<volSphericalTensorField>()) {
                manipulate(
                    driver.getResult<volSphericalTensorField>(),
                    conditionField
                );
            } else {
                WarningIn("Foam::manipulateField::execute()")
                    << "Expression '" << expression_
                        << "' evaluated to an unsupported type "
                        << driver.typ() << " that is incompatible with a mask defined on cells"
                        << endl;
            }
        } else if(driver.resultIsTyp<surfaceScalarField>(true)) {
            surfaceScalarField conditionField(driver.getResult<surfaceScalarField>());

            driver.parse(expression_);

            if(driver.resultIsTyp<surfaceVectorField>()) {
                manipulateSurface(
                    driver.getResult<surfaceVectorField>(),
                    conditionField
                );

            } else if(driver.resultIsTyp<surfaceScalarField>()) {
                manipulateSurface(
                    driver.getResult<surfaceScalarField>(),
                    conditionField
                );
            } else if(driver.resultIsTyp<surfaceTensorField>()) {
                manipulateSurface(
                    driver.getResult<surfaceTensorField>(),
                    conditionField
                );
            } else if(driver.resultIsTyp<surfaceSymmTensorField>()) {
                manipulateSurface(
                    driver.getResult<surfaceSymmTensorField>(),
                    conditionField
                );
            } else if(driver.resultIsTyp<surfaceSphericalTensorField>()) {
                manipulateSurface(
                    driver.getResult<surfaceSphericalTensorField>(),
                    conditionField
                );
            } else {
                WarningIn("Foam::manipulateField::execute()")
                    << "Expression '" << expression_
                        << "' evaluated to an unsupported type "
                        << driver.typ() << " that is incompatible with a mask defined on faces"
                        << endl;
            }
        } else if(driver.resultIsTyp<pointScalarField>(true)) {
            pointScalarField conditionField(driver.getResult<pointScalarField>());

            driver.parse(expression_);

            if(driver.resultIsTyp<pointVectorField>()) {
                manipulate(
                    driver.getResult<pointVectorField>(),
                    conditionField,
                    "points"
                );

            } else if(driver.resultIsTyp<pointScalarField>()) {
                manipulate(
                    driver.getResult<pointScalarField>(),
                    conditionField,
                    "points"
                );
            } else if(driver.resultIsTyp<pointTensorField>()) {
                manipulate(
                    driver.getResult<pointTensorField>(),
                    conditionField,
                    "points"
                );
            } else if(driver.resultIsTyp<pointSymmTensorField>()) {
                manipulate(
                    driver.getResult<pointSymmTensorField>(),
                    conditionField,
                    "points"
                );
            } else if(driver.resultIsTyp<pointSphericalTensorField>()) {
                manipulate(
                    driver.getResult<pointSphericalTensorField>(),
                    conditionField,
                    "points"
                );
            } else {
                WarningIn("Foam::manipulateField::execute()")
                    << "Expression '" << expression_
                        << "' evaluated to an unsupported type "
                        << driver.typ() << " that is incompatible with a mask defined on faces"
                        << endl;
            }
        } else {
                WarningIn("Foam::manipulateField::execute()")
                    << "Mask " << maskExpression_ << " evaluates to type "
                        << driver.typ() << " which is an unsupported field type for logical"
                        << endl;
        }
    }

    driver_->tryWrite();

#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
    return true;
#endif
}


void Foam::manipulateField::end()
{
    execute();
}

void Foam::manipulateField::execute()
{
}

// ************************************************************************* //
