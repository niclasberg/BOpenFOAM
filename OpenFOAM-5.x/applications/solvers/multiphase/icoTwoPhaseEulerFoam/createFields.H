
Info<< "Creating twoPhaseModel\n" << endl;

twoPhaseModel fluid(mesh);
PhaseModel & phase1 = fluid.phase1();
PhaseModel & phase2 = fluid.phase2();

// References to phase specific entiies
volScalarField & alpha1 = phase1.alpha();
volScalarField & alpha2 = phase2.alpha();
volVectorField & U1 = phase1.U();
volVectorField & U2 = phase2.U();
surfaceScalarField & phi1 = phase1.phi();
surfaceScalarField & phi2 = phase2.phi();
surfaceScalarField & alphaPhi1 = phase1.alphaPhi();
surfaceScalarField & alphaPhi2 = phase2.alphaPhi();

// References to mixture entities
surfaceScalarField & phi = fluid.phi();

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.U()
);

#include "readGravitationalAcceleration.H"

Info<< "Calculating field g.h\n" << endl;
volScalarField gh("gh", g & mesh.C());
surfaceScalarField ghf("ghf", g & mesh.Cf());

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + fluid.rho()*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);
mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());