import numpy as np
import vtk
import vtk.numpy_interface
import os
from vtk.numpy_interface import dataset_adapter as dsa

def _verifyPolyData(dataSet):
	if not dataSet.IsA('vtkPolyData'):
		raise RuntimeError('The dataset is not a vtkPolyData')

def _verifyUnstructuredGrid(dataSet):
	if not dataSet.IsA('vtkUnstructuredGrid'):
		raise RuntimeError('The dataset is not a vtkUnstructuredGrid')

def writeDataSet(dataSet, fileName):
	filePrefix, fileExtension = os.path.splitext(fileName)
	if fileExtension == '.vtp':
		_verifyPolyData(dataSet)
		writer = vtk.vtkXMLPolyDataWriter()
	elif fileExtension == '.vtu':
		_verifyUnstructuredGrid(dataSet)
		writer = vtk.vtkXMLUnstructuredGridWriter()
	elif fileExtension == '.vtk':
		if dataSet.IsA('vtkUnstructuredGrid'):
			writer = vtk.vtkUnstructuredGridWriter()
		elif dataSet.IsA('vtkPolyData'):
			writer = vtk.vtkPolyDataWriter()
		else:
			raise RuntimeError('Incompatible data type ' + dataSet.GetClassName() + ' to format .vtk')
		writer.SetFileTypeToBinary()
	elif fileExtension == '.stl':
		_verifyPolyData(dataSet)
		writer = vtk.vtkSTLWriter()
	else:
		raise RuntimeError('Unknown file extension', fileExtension)
	writer.SetInputData(dataSet)
	writer.SetFileName(fileName)
	writer.Write()

def readDataSet(fileName):
	if not os.path.isfile(fileName):
		raise RuntimeError('The file', fileName, 'did not exist')

	filePrefix, fileExtension = os.path.splitext(fileName)
	if fileExtension == '.vtp':
		reader = vtk.vtkXMLPolyDataReader()
		reader.SetFileName(fileName)
	elif fileExtension == '.vtu':
		reader = vtk.vtkXMLUnstructuredGridReader()
		reader.SetFileName(fileName)
	elif fileExtension == '.stl':
		reader = vtk.vtkSTLReader()
		reader.SetFileName(fileName)
	elif fileExtension == '.case':
		reader = vtk.vtkEnSightGoldBinaryReader()
		reader.SetCaseFileName(fileName)
	else:
		raise RuntimeError('Unknown file extension', fileExtension)
	reader.Update()
	return reader.GetOutput()

def createPolyData(points):
	pd = vtk.vtkPolyData()
	pts = vtk.vtkPoints()
	pts.SetNumberOfPoints(len(points))
	for i, pt in enumerate(points):
		pts.SetPoint(i, pt[0], pt[1], pt[2])
	pd.SetPoints(pts)
	return pd

def addPolyLineToDataSet(polyData, points):
	origNumPts = polyData.GetNumberOfPoints()

	# Verify that the polyDataSet has a point array
	if origNumPts == 0:
		polyData.SetPoints(vtk.vtkPoints())
		polyData.SetLines(vtk.vtkCellArray())

	cell = vtk.vtkPolyLine()
	cell.GetPointIds().SetNumberOfIds(points.shape[0])
	for i, pt in enumerate(points):
		# Insert point
		polyData.GetPoints().InsertNextPoint(pt[0], pt[1], pt[2])
		# Associate point id to line
		cell.GetPointIds().SetId(i, origNumPts+i)
	lines = polyData.GetLines().InsertNextCell(cell)

def _insertIntoVTKArray(array, value):
	if value.ndim == 0:
		array.InsertNextValue(value.take(0))
	elif array.GetNumberOfComponents() > 1 and value.ndim == 1:
		array.InsertNextTuple(value.tolist())
	else:
		for val in value:
 			_insertIntoVTKArray(array, val)

def _appendData(dataContainer, arrayName, valuesToInsert):
	values = np.array(valuesToInsert)
	if not dataContainer.HasArray(arrayName):
		if values.dtype == 'int64':
			array = vtk.vtkLongArray()
		elif values.dtype == 'int32':
			array = vtk.vtkIntArray()
		elif values.dtype == 'float64':
			array = vtk.vtkDoubleArray()
		elif values.dtype == 'float32':
			array = vtk.vtkFloatArray()
		else:
			raise RuntimeError('Unsupported data type ' + str(values.dtype))

		array.SetName(arrayName)

		# Determine number of components
		if values.ndim == 1 or values.ndim == 0:
			array.SetNumberOfComponents(1)
		elif values.ndim > 1:
			array.SetNumberOfComponents(values.shape[1])
		else:
			raise ValueError('Invalid value array dimension ' + str(values.shape))

		dataContainer.AddArray(array)
	else:
		array = dataContainer.GetArray(arrayName)

	_insertIntoVTKArray(array, values)

def appendCellData(dataSet, arrayName, values):
	_appendData(dataSet.GetCellData(), arrayName, values)

def appendPointData(dataSet, arrayName, values):
	_appendData(dataSet.GetPointData(), arrayName, values)

def createPolyLine(points):
	pd = createPolyData(points)

	# Create cells
	cell = vtk.vtkPolyLine()
	cell.GetPointIds().SetNumberOfIds(pd.GetNumberOfPoints())
	
	for i in range(pd.GetNumberOfPoints()):
		cell.GetPointIds().SetId(i, i)

	cells = vtk.vtkCellArray()
	cells.InsertNextCell(cell)
	pd.SetLines(cells)

	return pd
	

def getBlockByName(dataSet, name):
	if not dataSet.IsA('vtkCompositeDataSet'):
		raise RuntimeError('Cannot get a block from a non-composite dataset')
	d = _getBlockByName(dataSet, name)
	if d == None:
		raise RuntimeError('No block found in the dataset with the name ' + name)
	return d

def _getBlockByName(dataSet, name):
	it = dataSet.NewIterator()
	while not it.IsDoneWithTraversal():
		blockName = it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME())
		currentBlock = it.GetCurrentDataObject()
		if blockName.rstrip() == name:
			return currentBlock
		if it.GetCurrentDataObject().IsA('vtkCompositeDataSet'):
			d = _getBlockByName(currentBlock, name)
			if d != None:
				return d
		it.GoToNextItem()
	return None

def createFolder(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
