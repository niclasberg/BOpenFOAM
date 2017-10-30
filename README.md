# BOpenFOAM

Modified version of OpenFOAM 5.x with added support for hematocrit dependent viscosity and platelet activation modelling.

The following solvers have been added:
* vfDependentViscosityFoam (see tutorials/multiphase/vfDependentViscosityFoam/ for example cases)
* vfDependentViscosityLPTFoam

Some postprocessing tools have also been added:
* wallStressGradients - Computes wall shear stress (WSS) and gradients from the results of vfDependentViscosityFoam/vfDependentViscosityLPTFoam solvers.
* steadyPlateletLPT - Platelet tracking and activation from a time-independent velocity field
* plateletLPT - Platelet tracking and activation from a time-dependent solution

## Compilation
Compile OpenFOAM
```
cd OpenFOAM
source etc/bashrc
./Allwmake
```

Compile swak4foam
```
cd swak4foam/maintainanceScripts/
./compileRequirements.sh
export PATH=`pwd`/privateRequirements/bin:$PATH
cd ..
ln -s swakConfiguration.automatic swakConfiguration
./Allwmake
```
