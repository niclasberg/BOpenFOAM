import numpy as np
import vtk

# Convenience wrapper around vtkOpenFOAMReader
# Allows all reader options to be set in the constructor,
# and simplifies iteration over a timeseries.
class OpenFOAMReader:
	_currentIteration = -1
	_isInit = False
	_finishedReading = False

	def __init__(self, fileName, **kwargs):
		self.reader = vtk.vtkOpenFOAMReader()
		self.reader.SetFileName(fileName)
		self.reader.UpdateInformation()
		self.timeValues = self.reader.GetTimeValues()
		self.Nt = self.timeValues.GetNumberOfTuples()

		# Default options (read everything)
		self.reader.DecomposePolyhedraOn()
		self.reader.CacheMeshOn()
		self.reader.EnableAllCellArrays()
		self.reader.EnableAllPointArrays()
		self.reader.CreateCellToPointOff()
		self.reader.EnableAllPatchArrays()
		self.reader.EnableAllLagrangianArrays()

		# Read options
		for k,v in kwargs.items():
			print k, ' => ', v
			if k == 'cellArrays':		
				self.reader.DisableAllCellArrays()
				if v:
					for cellArray in v:
						self.reader.SetCellArrayStatus(cellArray, 1)
			elif k == 'patchArrays':
				self.reader.DisableAllPatchArrays()
				if v:
					for patchArray in v:
						self.reader.SetPatchArrayStatus(patchArray, 1)
			elif k == 'pointArrays':
				self.reader.DisableAllPointArrays()
				if v:
					for pointArray in v:
						self.reader.SetPointArrayStatus(pointArray, 1)
			elif k == 'lagrangianArrays':
				self.reader.DisableAllLagrangianArrays()
				if v:
					for lagrangianArray in v:
						self.reader.SetLagrangianArrayStatus(lagrangianArray, 1)
			elif k == 'decomposePolyhedra':
				if v:
					self.reader.DecomposePolyhedraOn()
				else:
					self.reader.DecomposePolyhedraOff()
			elif k == 'cellToPoint':
				if v:
					self.reader.CreateCellToPointOn()
				else:
					self.reader.CreateCellToPointOff()		
			elif k == 'cacheMesh':
				if v:
					self.reader.CacheMeshOn()
				else:
					self.reader.CacheMeshOff()
			else:
				raise ValueError('Invalid argument ' + k)

	def getDataSet(self):
		if not self._isInit:
			raise RuntimeError('The reader has not read anything yet!')
		if self.finishedReading():
			raise RuntimeError('The reader has no more data to read')
		return self.reader.GetOutput()

	def getPatchNames(self):
		return [self.reader.GetPatchArrayName(i) for i in range(self.reader.GetNumberOfPatchArrays())]

	def readIteration(self, iteration):
		if iteration < 0 or iteration >= self.Nt:
			raise ValueError('Iteration number out of range')
		self._currentIteration = iteration
		self.reader.SetTimeValue(self.timeValues.GetTuple(iteration)[0])
		self.reader.Modified()
		self.reader.Update()
		self._isInit = True

	def readTime(self, t):
		for i in range(self.Nt):
			if self.timeValues.GetTuple(i)[0] == t:
				it = i
				break
		else:
			raise RuntimeError('The time ' + str(t) + ' did not exist in the dataset')
		self.readIteration(i)

	def currentTime(self):
		return self.timeValues.GetTuple(self._currentIteration)[0]

	def currentIteration(self):
		return self._currentIteration

	def finishedReading(self):
		return self._finishedReading

	def startReading(self):
		self._finishedReading = False
		self.readIteration(0)

	def skipAndRead(self, n):
		if not self._isInit:
			raise RuntimeError('The reader has not been intialized, call startReading() before readNext()')
		if (self._currentIteration+n) >= self.Nt:
			self._finishedReading = True
		else:
			self.readIteration(self._currentIteration+n)

	def readNext(self):
		self.skipAndRead(1)

def _writeFoamHeader(f, fieldType, location, fieldName):
	f.write('/*--------------------------------*- C++ -*----------------------------------*\ \n')
	f.write('| =========                 |                                                 |\n')
	f.write('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n')
	f.write('|  \\    /   O peration     | Version:  2.3.x                                 |\n')
	f.write('|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n')
	f.write('|    \\/     M anipulation  |                                                 |\n')
	f.write('\*---------------------------------------------------------------------------*/\n')
	f.write('FoamFile\n')
	f.write('{\n')
	f.write('	version     2.0;\n')
	f.write('	format      ascii;\n')
	f.write('	class       '+fieldType+';\n')
	f.write('	location    "'+location+'";\n')
	f.write('	object      ' + fieldName + ';\n')
	f.write('}\n')
	f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')
	

def _writeVector(f, vec):
	f.write('(')
	for i in range(vec.size):
		if i != 0:
			f.write(' ')
		f.write(str(vec[i]))
	f.write(')')	

def _writeField(f, fieldData, fieldType):
	if fieldType == 'volScalarField':
		if fieldData.ndim == 0:
			# Uniform data
			f.write('uniform ' + str(fieldData))
		else:
			# Non-uniform data
			f.write('nonuniform List<scalar> ' + str(fieldData.shape[0]) + '\n(\n')
			for i in range(0, fieldData.shape[0]):
				f.write(str(fieldData[i]) + '\n')
			f.write(')')
	else:
		if fieldData.ndim == 1:
			# Uniform data
			f.write('uniform ')
			_writeVector(f, fieldData)
		else:
			# Non-uniform data
			f.write('nonuniform List<')
			if fieldType == 'volVectorField':
				f.write('vector')
			elif fieldType == 'volSymmTensorField':
				f.write('symmTensor')
			elif fieldType == 'volTensorField':
				f.write('tensor')
			else:
				raise RuntimeError('Unknown field type ' + fieldType)
			f.write('> ' + str(fieldData.shape[0]) + '\n(\n')
			for i in range(fieldData.shape[0]):
				_writeVector(f, fieldData[i, :])
				f.write('\n')
			f.write(')')
	f.write(';\n')

class _DataDimension:
	def __init__(self, dimList):
		self.dimList = np.array(dimList)

	def __mul__(self, rhs):
		return _DataDimension(self.dimList + rhs.dimList)

	def __div__(self, rhs):
		return _DataDimension(self.dimList - rhs.dimList)

	def __pow__(self, rhs):
		return _DataDimension(rhs*self.dimList)

	def __str__(self):
		return '[{} {} {} {} 0 0 0]'.format(self.dimList[0], self.dimList[1], self.dimList[2], self.dimList[3])

Dimensionless = _DataDimension([0, 0, 0, 0])
DimMass = _DataDimension([1, 0, 0, 0])
DimLength = _DataDimension([0, 1, 0, 0])
DimTime = _DataDimension([0, 0, 1, 0])
DimTemperature = _DataDimension([0, 0, 0, 1])

def _determineDataType(fields):
	maxSize = [0, 0]
	for field in fields:
		for axis in range(field.ndim):
			maxSize[axis] = max(field.shape[axis], maxSize[axis])
	if maxSize[1] != 0:
		if maxSize[1] == 1:
			return 'volScalarField'
		if maxSize[1] == 3:
			return 'volVectorField'
		if maxSize[1] == 6:
			return 'volSymmTensorField'
		if maxSize[1] == 9:
			return 'volTensorField'
	else:
		if maxSize[0] != 0:
			return 'volScalarField'
	raise RuntimeError('Unable to determine data type')

def writeFoamData(caseFolder, fieldName, **kwargs):
	internalField = kwargs.get('internalField')
	boundaryField = kwargs.get('boundaryField')
	dimension = kwargs.get('dimension', Dimensionless)
	time = kwargs.get('time', 0)
	fieldType = kwargs.get('fieldType')

	if int(time) == time:
		timeStr = str(int(time)) 
	else: 
		timeStr = str(time)

	# Determine datatype if not provided
	if fieldType == None:
		fieldType = _determineDataType([internalField] + boundaryField.values())
	
	# Verify that all the fields conform to the provided data type
	fieldTests = dict()
	fieldTests['volTensorField'] = [lambda field: ((field.ndim == 1 and field.shape[0] == 9) or (field.ndim == 2 and field.shape[1] == 9)), '(9, ) or (x, 9)']
	fieldTests['volScalarField'] = [lambda field: (field.ndim == 0 or field.ndim == 1), '() or (x,)']
	fieldTests['volVectorField'] = [lambda field: ((field.ndim == 1 and field.shape[0] == 3) or (field.ndim == 2 and field.shape[1] == 3)), '(x, 3) or (3,)']
	fieldTests['volSymmTensorField'] = [lambda field: ((field.ndim == 1 and field.shape[0] == 6) or (field.ndim == 2 and field.shape[1] == 6)), '(x, 6) or (6,)']

	fieldTest = fieldTests[fieldType][0]
	expectedShape = fieldTests[fieldType][1]
	
	for field, domainName in zip([internalField] + boundaryField.values(), ['internalField'] + boundaryField.keys()):
		if not isinstance(field, basestring):
			if not fieldTest(field):
				raise RuntimeError('Invalid ' + fieldType +' shape ' + str(field.shape) + ' in ' + domainName + ' (expected ' + expectedShape + ')')

	with open(caseFolder + '/'+timeStr+'/' + fieldName, 'w') as f:
		_writeFoamHeader(f, fieldType, timeStr, fieldName)
		f.write('dimensions ' + str(dimension) + ';\n')
		
		# Write the internal field
		if internalField != None:
			f.write('internalField ')
			_writeField(f, internalField, fieldType)

		# Write boundaries
		if boundaryField != None:
			f.write('boundaryField\n')
			f.write('{\n')
			for boundaryName, data in boundaryField.iteritems():
				f.write('	'+ boundaryName+'\n')
				f.write('	{\n')
				if isinstance(data, basestring):
					f.write('		type ' + data + ';\n')
				else:
					f.write('		type fixedValue;\n')
					f.write('		value ')
					_writeField(f, data, fieldType)
				f.write('	}\n')
			f.write('}\n')
