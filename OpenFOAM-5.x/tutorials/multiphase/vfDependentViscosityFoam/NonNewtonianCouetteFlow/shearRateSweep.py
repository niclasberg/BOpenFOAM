import vtk
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa
import foamIO
import vtkio
import os
import matplotlib.pyplot as plt

def writeTransportDict(viscosityModel):
	with open('constant/transportProperties', 'w') as f:
		f.write('/*--------------------------------*- C++ -*----------------------------------*\\\n')
		f.write('| =========                 |                                                 |\n')
		f.write('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n')
		f.write('|  \\    /   O peration     | Version:  1.6                                   |\n')
		f.write('|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n')
		f.write('|    \\/     M anipulation  |                                                 |\n')
		f.write('\*---------------------------------------------------------------------------*/\n')
		f.write('FoamFile\n')
		f.write('{\n')
		f.write('    version     2.0;\n')
		f.write('    format      ascii;\n')
		f.write('    class       dictionary;\n')
		f.write('    location    "constant";\n')
		f.write('    object      transportProperties;\n')
		f.write('}\n')
		f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')
		f.write('Dab           Dab [0 2 -1 0 0 0 0] 0.91e-08;\n')
		f.write('alphatab           alphatab [0 0 0 0 0 0 0] 1.2;\n')
		f.write('phases (RBC plasma);\n')
		f.write('RBC\n')
		f.write('{\n')
		f.write('    rho            rho [ 1 -3 0 0 0 0 0 ] 1102;\n')
		f.write('}\n')
		f.write('plasma\n')
		f.write('{\n')
		f.write('    rho             rho [1 -3 0 0 0 0 0] 1025;\n')
		f.write('}\n')
		f.write('transportModel	' + viscosityModel + ';\n')
		f.write('nu              nu [0 2 -1 0 0 0 0] 3.3e-06;\n')
		f.write('NewtonianCCoeffs\n')
		f.write('{\n')
		f.write('	mu			mu [1 -1 -1 0 0 0 0] 3.5e-3;\n')
		f.write('}\n')
		f.write('WalburnSchneckCCoeffs\n')
		f.write('{\n')
		f.write('    C1	            C1 [ 1 -1 -1 0 0 0 0 ] 0.797e-03;\n')
		f.write('    TPMA	    	TPMA [ 0 0 0 0 0 0 0 ] 25;\n')
		f.write('    muMax           muMax [ 1 -1 -1 0 0 0 0 ] 39.05e-03;\n')
		f.write('    muMin           muMin [ 1 -1 -1 0 0 0 0 ] 3.024e-03;\n')
		f.write('}\n')
		f.write('QuemadaCCoeffs\n')
		f.write('{\n')
		f.write('    mup             mup [ 1 -1 -1 0 0 0 0 ] 1.32e-03;\n')
		f.write('    muMax           muMax [ 1 -1 -1 0 0 0 0 ] 66.26e-03;\n')
		f.write('}\n')
		f.write('CassonCCoeffs\n')
		f.write('{\n')
		f.write('    mup             mup [ 1 -1 -1 0 0 0 0 ] 1.32e-03;\n')
		f.write('    nuMax           nuMax [ 0 2 -1 0 0 0 0 ] 15.406e-06;\n')
		f.write('    A               A [ 0 0 0 0 0 0 0 ] 1.387133;\n')
		f.write('    B               B [ 0.5 -0.5 -1 0 0 0 0 ] 1.965353e-01;\n')
		f.write('}\n')
		f.write('// ************************************************************************* //\n')

rhoRbc = 1102
rhoPlasma = 1025
H = 0.45
rho = H*rhoRbc + (1. - H) * rhoPlasma

shears = np.array([1., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.])
nus = np.zeros_like(shears)

# Quemada viscosity
shearTheory = np.linspace(np.min(shears), np.max(shears), 100)
mup = 1.32e-3
k0 = np.exp(3.874 - 10.41*H + 13.8*H**2 - 6.738*H**3)
kinf = np.exp(1.3435 - 2.803*H + 2.711*H**2 - 0.6479*H**3)
sqrtGammaOverGammaC = np.sqrt(shearTheory / np.exp(-6.1508 + 27.923*H - 25.6*H**2 + 3.697*H**3))
nuQuemada = mup * np.power(1.0 - 0.5 * H * (k0 + kinf*sqrtGammaOverGammaC) / (1.0 + sqrtGammaOverGammaC), -2) / rho

# Walburn Schneck viscosity
C1 = 0.797e-03
TPMA = 25
HPercent = H * 100.
nuWalburnSchneck = C1*np.exp(0.0608*HPercent) * np.exp(14.585 * TPMA/ HPercent**2) * np.power(shearTheory, (-0.00499 * HPercent)) / rho

# Newtonian viscosity
nuNewtonian = 3.5e-3 / rho * np.ones(shearTheory.shape[0])

for viscosityModel, viscosityValues in zip(['NewtonianC', 'WalburnSchneckC', 'QuemadaC'], [nuNewtonian, nuWalburnSchneck, nuQuemada]):
	writeTransportDict(viscosityModel)
	for i, shear in enumerate(shears):
		Umax = shear * 0.01

		# Read foam case
		reader = foamIO.OpenFOAMReader('r.foam',
					patchArrays = ['internalMesh'],
					cellToPoint = False)
		reader.startReading()

		# Compute cell centers
		ccFilter = vtk.vtkCellCenters()
		ccFilter.SetInputData(vtkio.getBlockByName(reader.getDataSet(), 'internalMesh'))
		ccFilter.Update()

		ds = dsa.WrapDataObject(ccFilter.GetOutput())
		points = ds.Points

		# Create initial condition for velocity
		U = np.zeros_like(points)
		U[:, 0] = shear*points[:, 2]

		foamIO.writeFoamData('.', 'U', 
			dimension = foamIO.DimLength / foamIO.DimTime,
			fieldType = 'volVectorField',
			internalField = U,
			boundaryField = {'left': 'cyclic', 'right': 'cyclic', 'front': 'cyclic', 'back': 'cyclic', 'bottomWall': np.array([0, 0, 0]), 'topWall': np.array([Umax, 0, 0])})

		# Run simulation
		os.system('vfDependentViscosityFoam')

		# Read viscosity data
		reader = foamIO.OpenFOAMReader('r.foam',
				patchArrays = ['internalMesh'],
				cellToPoint = False,
				cellArrays = ['nu'])
		reader.readIteration(1)
		nu = dsa.WrapDataObject(vtkio.getBlockByName(reader.getDataSet(), 'internalMesh')).CellData['nu']

		nus[i] = np.mean(nu)

	plt.clf()
	plt.plot(shears, nus, 'o', label = 'Simulation')
	plt.plot(shearTheory, viscosityValues, label = 'Theory')
	plt.xlabel('Shear rate [1/s]')
	plt.ylim([0, 2e-5])
	plt.ylabel('Kinematic viscosity [m^2/s]')
	plt.legend(loc = 'upper right')
	plt.savefig(viscosityModel + '_validation.png', bbox_inches='tight')
