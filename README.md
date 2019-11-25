# BOpenFOAM

Modified version of OpenFOAM 5.x with added support for hematocrit dependent viscosity and platelet activation modelling.

The main contribution is the addition of the solver vfDependentViscosityFoam (see tutorials/multiphase/vfDependentViscosityFoam/ for example cases). Documentation about the implemented model and properties can be found at: https://www.mech.kth.se/~niber/docs/#/openfoam . In essence, the model corresponds to a one-fluid formulation with a hematocrit and shear rate dependent viscosity (multiple options are implemented such as Quemada, Walburn-Schneck and Casson). The hematocrit is solved through a convection-flux equation (the flux term can be modeled as Fickian, hematocrit dependent or with the Leighton & Acrivos model, taking shear-induced diffusion into account).

Some function objects and post processing tools have also been implemented:
* trackPlatelets - post processing program for tracking of platelets and the activation state (related to the probability of activation) along the pathlines. 
* plateletLPT - Function object for platelet tracking in runtime.
* vorticitySources - evaluates the terms of the vorticity transport equation over time.

## Compilation
Compile OpenFOAM
```
cd OpenFOAM-5.x
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
