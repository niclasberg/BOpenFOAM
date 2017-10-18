# BOpenFOAM

Modified version of OpenFOAM 5.x with added support for hematocrit dependent viscosity and platelet activation modelling.

The following solvers have been added:
* vfDependentViscosityFoam
* vfDependentViscosityLPTFoam

Each solver will soon have a tutorial case.

Some postprocessing tools have also been added:

## Compilation
Compile OpenFOAM
```
cd OpenFOAM
source etc/bashrc
./Allwmake
```

Compile swak4foam
```
cd swak4foam
./Allwmake
```
