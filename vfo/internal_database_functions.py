import numpy as np
import os
import warnings
import vfo.classTags as classTags
import openseespy.opensees as ops



def _getNodesandElements():
	"""
	This function returns the nodes and elments for an active model, in a 
	standardized format. The OpenSees model must be active in order for the 
	function to work.
    
	Returns
	-------
	nodes : 2dArray
		An array of all nodes in the model.
		Returns nodes in the shape:
		[Nodes, 3] in 2d and [Nodes, 4]
		For each node the information is tored as follows:
		[NodeID, x, y] or [NodeID, x, y, z]
	elements : Array 
		An list of all elements in. Each entry in the list is it's own'
		[element1, element2,...],   element1 = [element#, node1, node2,...]
	"""
   
	# Get nodes and elements
	nodeList = ops.getNodeTags()
	eleList = ops.getEleTags()   
    
	# Check Number of dimensions and intialize variables
	ndm = len(ops.nodeCoord(nodeList[0]))
	Nnodes = len(nodeList)
	Nele = len(eleList)
	nodes = np.zeros([Nnodes, ndm + 1])
	eleClassTags = np.zeros([Nele, 2])
    
	# Get Node list
	for ii, node in enumerate(nodeList):
		nodes[ii,0] = node
		nodes[ii,1:] = ops.nodeCoord(nodeList[ii])           
    
	Nele = len(eleList)
	elements = [None]*Nele
	
    
	# Generate the element list by looping through all emenemts
	for ii, ele in enumerate(eleList):
		tempNodes = ops.eleNodes(ele)
	
		tempNnodes = len(tempNodes)
		tempEle = np.zeros(tempNnodes + 1)
        
		tempEle[0] = int(ele)
		tempEle[1:] = tempNodes
        
		elements[ii] = tempEle     

	# Generate element class tags by looping through all elements
	for ii, ele in enumerate(eleList):
		eleClassTags[ii,0] = ele
		eleClassTags[ii,1] = ops.getEleClassTags(ele)[0]
    
	return nodes, elements, eleClassTags

def _saveNodesandElements(ModelName):
	"""   
	This file saves the node and element information for the structure. 
	For each node information is saved in the following format:
		Nodes:    [NodeID, xcord, ycord] or [NodeID, xcord, ycord, zcord]
    
	For elements, the element is saved with the element connectivity. 
	A different file is created for each type of element
	each possible element type.
		Elements: [EleID, eleNode1, eleNode2, ... , eleNodeN]
	Parameters
	----------
	nodeName : str, optional
		The name of the file to be saved. The default is 'Nodes'.
    eleName : str, optional
		The name of the . The default is 'Elements'.
	delim : str, optional
		The delimeter for the output file. The default is ','.
	fmt : str, optional
		the format of the file to be saved in. The default is '%.5e'.
	"""
    

	# Consider making these optional arguements
	nodeName = 'Nodes'
	eleName = 'Elements'
	eleClassName = 'EleClassTags'
	delim = ' '
	fmt = '%.5e'
	ftype = '.out'
    
	ODBdir = ModelName+"_ODB"		# ODB Dir name

	# Read noades and elements
	nodes, elements, eleClassTags = _getNodesandElements()

	# Sort through the element arrays
	ele2Node = np.array([ele for ele in elements if len(ele) == 3])
	ele3Node = np.array([ele for ele in elements if len(ele) == 4])
	ele4Node = np.array([ele for ele in elements if len(ele) == 5])
	ele8Node = np.array([ele for ele in elements if len(ele) == 9])
    
	nodeFile = os.path.join(ODBdir, nodeName + ftype)
    
	ele2File = os.path.join(ODBdir, eleName + "_2Node" + ftype)
	ele3File = os.path.join(ODBdir, eleName + "_3Node" + ftype)
	ele4File = os.path.join(ODBdir, eleName + "_4Node"  + ftype)
	ele8File = os.path.join(ODBdir, eleName + "_8Node"  + ftype)
	
	eleClassTagsFile = os.path.join(ODBdir, eleClassName + ftype)

	# SaveNodes
	np.savetxt(nodeFile, nodes, delimiter = delim, fmt = fmt)
    
	# Save element arrays
	np.savetxt(ele2File, ele2Node, delimiter = delim, fmt = fmt)
	np.savetxt(ele3File, ele3Node, delimiter = delim, fmt = fmt)
	np.savetxt(ele4File, ele4Node, delimiter = delim, fmt = fmt)
	np.savetxt(ele8File, ele8Node, delimiter = delim, fmt = fmt)
	
	# Save Element Class Tags
	np.savetxt(eleClassTagsFile, eleClassTags, delimiter = delim, fmt = fmt)
	



def _readNodesandElements(ModelName):
	"""   
	This function reads input node/element information, assuming it is in the 
	standard format. 
	If outputDir == False, the base directory will be used.    
    
	Parameters
    ----------
    nodeName : str, optional
        The base name for the node file. It will be appended to include
        the file type. The default is 'Nodes.out'.
    eleName : str, optional
        The base nae for the element files. The default is 'Elements.out'.
    delim : str, optional
        The delimiter for files to be read. The default is ','.
    dtype : TYPE, optional
        The data type to read in. The default is 'float32'.
    Returns
    -------
    nodes : Array
        An output vector in standard format
    elements : List
        An output Element vector in standard format.
        elements = [ele1, ele2,..., elen], 
        ele1 = [element, node 1, node 2, ... , node n]
    """

	# Consider making these optional arguements
	nodeName = 'Nodes'
	eleName = 'Elements'
	eleClassName = 'EleClassTags'
	delim = ' '
	dtype ='float32' 
	ftype = '.out'
        
	ODBdir = ModelName+"_ODB"		# ODB Dir name
    
	# Check if output database exists
	if not os.path.exists(ODBdir):
		print('No directory found for nodes and elements')
        
	# Generate the file names
	nodeFile = os.path.join(ODBdir, nodeName + ftype)
	ele2File = os.path.join(ODBdir, eleName + "_2Node" + ftype)
	ele3File = os.path.join(ODBdir, eleName + "_3Node" + ftype)
	ele4File = os.path.join(ODBdir, eleName + "_4Node"  + ftype)
	ele8File = os.path.join(ODBdir, eleName + "_8Node"  + ftype)     
       
	eleFileNames = [ele2File, ele3File, ele4File, ele8File]    
	
	eleClassTagsFile = os.path.join(ODBdir, eleClassName + ftype)
    
	## Load Node information
	try:
		nodes = np.loadtxt(nodeFile, dtype, delimiter = delim, unpack=False)
	except:
		print("Reading node data from a OpenSees Tcl model")
		nodes = np.transpose(np.loadtxt(nodeFile, dtype=float, delimiter=None, converters=None, unpack=True))
			
	# Populate an array with the input element information
	TempEle = [[]]*4
    
	# Check if the file exists, read it if it does. Ignore warnings if the files are empty
	for ii, FileName in enumerate(eleFileNames):
		if os.path.isfile(FileName):
			with warnings.catch_warnings():
				warnings.simplefilter("ignore")
				try:
					TempEle[ii] = np.loadtxt(FileName, dtype,  delimiter = delim, skiprows=0, ndmin=2, unpack=False)
				except:
					print("Reading element data from a OpenSees Tcl model")
					TempEle[ii] = np.transpose(np.loadtxt(FileName, dtype=float, delimiter=None,  skiprows=0, ndmin=2,converters=None, unpack=True))

	# define the final element array
	elements = [*TempEle[0],*TempEle[1],*TempEle[2],*TempEle[3]]

	## Load Element Class Tags information
	try:
		eleClassTags = np.loadtxt(eleClassTagsFile, dtype, delimiter = delim, unpack=False)
	except:
		print("Reading element ClassTag data from a OpenSees Tcl model")
		eleClassTags = np.loadtxt(eleClassTagsFile, dtype=float, delimiter=None, unpack=False)
		# print(eleClassTags)
	
	# Check if any files were read
	if elements is []:
		raise Exception('No element information files were found!')

	return nodes, elements, eleClassTags
	

################ Element lists for stress and strain ###########

def _getEleListforStress():
	
	# Sort through element array and create groups for stress and strain recorder 
	# Since this saveNodesandElements command is called in createODB before recorders, these lists can be saved in the ODB and also used for recorders. 
	# return shell_4N4GP_eleList , etc.
	## FUTURE: Also save element information for beam-column elements for element force and fiber recorder
	
	eleList = ops.getEleTags()  
	
	eleList_shell_4N4GP = [ele for ele in eleList if ops.getEleClassTags(ele)[0] in classTags.shell_4N4GP_EleTags]
	eleList_shell_3N3GP = [ele for ele in eleList if ops.getEleClassTags(ele)[0] in classTags.shell_3N3GP_EleTags]
	eleList_tetra_4N4GP = [ele for ele in eleList if ops.getEleClassTags(ele)[0] in classTags.tetra_4N4GP_EleTags]
	eleList_brick_8N8GP = [ele for ele in eleList if ops.getEleClassTags(ele)[0] in classTags.brick_8N8GP_EleTags]
	
	return eleList_shell_4N4GP, eleList_shell_3N3GP, eleList_tetra_4N4GP, eleList_brick_8N8GP
	
################ ModeShapes #############################

def _getModeShapeData(modeNumber):
	
	# Get nodes and elements
	nodeList = ops.getNodeTags()
    
	# Check Number of dimensions and intialize variables
	ndm = len(ops.nodeCoord(nodeList[0]))
	Nnodes = len(nodeList)
	nodes_modeshape = np.zeros([Nnodes, ndm + 1])
    
	for ii, node in enumerate(nodeList):
		nodes_modeshape[ii,0] = node
		tempData = ops.nodeEigenvector(nodeList[ii], modeNumber)
		nodes_modeshape[ii,1:] = tempData[0:ndm]

	return nodes_modeshape
	
	
def _saveModeShapeData(ModelName,modeNumber):
    
	nodes_modeshape = _getModeShapeData(modeNumber)
	
	# Consider making these optional arguements
	modeName = "ModeShape"
	delim = ' '
	fmt = '%.5e'
	ftype = '.out'
    
	ODBdir = ModelName+"_ODB"		# ODB Dir name
	ModeShapeDir = os.path.join(ODBdir,"ModeShapes")
	modeFile = os.path.join(ModeShapeDir, modeName+str(modeNumber)+ftype)
	
	## ModeShapeDir is a default name
	np.savetxt(modeFile, nodes_modeshape, delimiter = delim, fmt = fmt)    
	
	
def _readModeShapeData(ModelName,modeNumber):

	# Consider making these optional arguements
	modeName = "ModeShape"
	delim = ' '
	fmt = '%.5e'
	dtype ='float32'
	ftype = '.out'
        
	ODBdir = ModelName+"_ODB"		# ODB Dir name
	ModeShapeDir = os.path.join(ODBdir,"ModeShapes")
	
    # Check if output database exists
	if not os.path.exists(ModeShapeDir):
		print('Error: No directory found for modeshapes. Use recordODB() command to save modeshapes.')

	modeFile = os.path.join(ModeShapeDir, modeName+str(modeNumber)+ftype)
	modeTFile = os.path.join(ModeShapeDir, "ModalPeriods.out")
    
	## Read modal period data to display
	periods = np.loadtxt(modeTFile, dtype, delimiter = delim, unpack=False)
	
	## Load Node information
	try:
		nodes_modeshape = np.loadtxt(modeFile, dtype, delimiter = delim, unpack=False)
	except:
		print("Reading modeshape data from a OpenSees Tcl model")
		nodes_modeshape = np.transpose(np.loadtxt(modeFile, dtype=float, delimiter=None, converters=None, unpack=True))

	return nodes_modeshape, periods


############## Node Displacement Data ######################################

def _readNodeDispData(ModelName,LoadCaseName):
	
	ODBdir = ModelName+"_ODB"		# ODB Dir name
	LoadCaseDir = os.path.join(ODBdir, LoadCaseName)
	
	# Get number of nodes in the model to set a node displacement array
	nodes,elements, eleClassTags = _readNodesandElements(ModelName)
	Nnodes = len(nodes)
	ndm = len(nodes[0,1:])
	
	NodeDispFile = os.path.join(LoadCaseDir,"NodeDisp_All.out")
	Disp = np.transpose(np.loadtxt(NodeDispFile, dtype=float, delimiter=None, converters=None, unpack=True))
	
	timeSteps = Disp[:,0]
	Ntime = len(Disp[:,0])

	tempDisp = np.zeros([Ntime,Nnodes,ndm])
	tempDisp[:,:,0] = Disp[:,1::ndm]
	tempDisp[:,:,1] = Disp[:,2::ndm]
	
	if ndm == 3:
		tempDisp[:,:,2] = Disp[:,3::ndm]
		
	nodes_displacement = tempDisp
	
	return timeSteps, nodes_displacement


#### Read fibre data

def _readFiberData2D(ModelName, LoadCaseName, eleNumber, sectionNumber):
    
	# Consider making these optional arguements
	FibreName = "FiberData"
	delim = ' '
	# fmt = '%.5e'
	dtype ='float32'
	ftype = '.out'    
    
	ODBdir = ModelName+"_ODB"		# ODB Dir name
	FibreFileName = FibreName  + '_ele_' + str(eleNumber) + '_section_' + str(sectionNumber) + ftype
	FiberDir = os.path.join(ODBdir, LoadCaseName, FibreFileName)
	# Check if output database exists
	if not os.path.exists(FiberDir):
		print('Error: No file for Fiber data. Use saveFiberData2D() to create a recorder.')    
    
	FiberData = np.loadtxt(FiberDir, dtype=dtype, delimiter=delim)
	timeSteps = FiberData[:,0]
	FiberData = FiberData[:,1:]

	return timeSteps, FiberData



def _getMaxNodeDisp(timeSteps, displacement_nodeArray):
	
	"""
	This internal function calculates the absolute maximum deformation of nodes in a model. 
	The output is used in creating animation with contours.
	
	Do not use scale factor for these calculations. Use "scale" option from PyVista.
	
	"""
	print("calculating maximum node deflections")
	
	jj=0
	
	maxDispArray = np.empty([len(timeSteps), 7])  # Each row of this array will store [jj, min_X, max_x, min_y, max_y, min_z, max_z]
	
	
	for jj, tstep in enumerate(timeSteps):
		x_max = np.max(displacement_nodeArray[jj,:,0])
		x_min = np.min(displacement_nodeArray[jj,:,0])
		y_max = np.max(displacement_nodeArray[jj,:,1])
		y_min = np.min(displacement_nodeArray[jj,:,1])
		z_max = np.max(displacement_nodeArray[jj,:,2])
		z_min = np.min(displacement_nodeArray[jj,:,2])
		
		thisArray = np.array([jj, x_max, x_min, y_max, y_min, z_max, z_min])
		maxDispArray[jj,:] = thisArray
		
		# print(jj)

	return maxDispArray




def _saveMonitorElementData(monitorEleType, monitorOutput, GroupMonitorDir, monitorEleTags, deltaT, dofList_ele):

	monitorGroupArray = np.asarray(monitorEleTags)
	
	if monitorEleType == "TwoNodeLink" or monitorEleType == "twoNodeLink":
		if monitorOutput=="deformation" or monitorOutput=="deformations":
			"""
			3D ELements: 1,2,3,4,5,6;
			2D Elements: 1,2,3
			"""
			MonitorDefFile = os.path.join(GroupMonitorDir,"MonitorGroup_Deformation.out")
			MonitorForceFile = os.path.join(GroupMonitorDir,"MonitorGroup_Force.out")				# Needed to plot element hysteresis
			MonitorEleFile = os.path.join(GroupMonitorDir,"MonitorGroup_Tags_Deformation.out")
			MonitorInfoFile = os.path.join(GroupMonitorDir,"MonitorGroup_Info_Deformation.out")
		ops.recorder('Element', '-file', MonitorDefFile,  '-time', '-dT', deltaT, '-ele', *monitorEleTags, '-dof',*dofList_ele, monitorOutput) # Two node link element needs dofs input to record
		ops.recorder('Element', '-file', MonitorForceFile,  '-time', '-dT', deltaT, '-ele', *monitorEleTags, '-dof',*dofList_ele, 'localForce') # 6 columns for each column
		np.savetxt(MonitorEleFile, monitorGroupArray, delimiter = ' ', fmt = '%.5e')			# This file will be read by plotting functions
		np.savetxt(MonitorInfoFile, np.asarray(dofList_ele), delimiter = ' ', fmt = '%.5e')			# This file will be read by plotting functions

	elif monitorEleType in ["forceBeamColumn","beamWithHinges","dispBeamColumn","nonlinearBeamColumn"]:
		if monitorOutput in ["chordRotation", "chordDeformations"]:
			"""
			3D ELements: eps, thetaZ_1, thetaZ_2, thetaY_1, thetaY_2, thetaX;
			2D Elements: eps, theta_1, theta_2
			"""
			chordOutput = [1,2,3,4,5,6]  
			MonitorDefFile = os.path.join(GroupMonitorDir,"MonitorGroup_Rotation.out")
			MonitorEleFile = os.path.join(GroupMonitorDir,"MonitorGroup_Tags_Rotation.out")
			MonitorInfoFile = os.path.join(GroupMonitorDir,"MonitorGroup_Info_Rotation.out")
			ops.recorder('Element', '-file', MonitorDefFile,  '-time', '-dT', deltaT, '-ele', *monitorEleTags, monitorOutput)  # Records chord output
			np.savetxt(MonitorEleFile, monitorGroupArray, delimiter = ' ', fmt = '%.5e')			# This file will be read by plotting functions
			np.savetxt(MonitorInfoFile, np.asarray(dofList_ele), delimiter = ' ', fmt = '%.5e')			# This file will record number of outputs 
		elif monitorOutput == "deformation" or monitorOutput=="deformations":
			"""
			3D ELements: AxialStrain, CurvY, CurvZ, TorsionStrain
			2D Elements: AxialStrain, Curvature
			"""
			# Records axial strain, curvature
			curvatureOutput = [1,2,3,4,5,6]  #
			MonitorInfoFile = os.path.join(GroupMonitorDir,"MonitorGroup_Info_SecDeformation.out")
			MonitorEleFile = os.path.join(GroupMonitorDir,"MonitorGroup_Tags_SecDeformation.out")
			for SectionNum in range(1,numSections+1):
				MonitorDefFile = os.path.join(GroupMonitorDir,"MonitorGroup_Deformation_section"+str(SectionNum)+".out")
				ops.recorder('Element', '-file', MonitorDefFile,  '-time', '-dT', deltaT, '-ele', *monitorEleTags, 'section',str(SectionNum), monitorOutput)  # Records chord output
			np.savetxt(MonitorEleFile, monitorGroupArray, delimiter = ' ', fmt = '%.5e')			# This file will be read by plotting functions
			np.savetxt(MonitorInfoFile, np.asarray(dofList_ele), delimiter = ' ', fmt = '%.5e')			# This file will record number of outputs 
		else:
		#### FUTURE: Add recorder for nonlinear beam column element sections; hinge rotation for beam with hinges
			pass
	else:
		pass
					
					
def _readMonitorElementData(monitorOutput,GroupMonitorDir):

	if monitorOutput=="deformation" or monitorOutput=="deformations":
		MonitorDefFile = os.path.join(GroupMonitorDir,"MonitorGroup_Deformation.out")
		MonitorForceFile = os.path.join(GroupMonitorDir,"MonitorGroup_Force.out")				# Needed to plot element hysteresis
		MonitorEleFile = os.path.join(GroupMonitorDir,"MonitorGroup_Tags_Deformation.out")
		MonitorInfoFile = os.path.join(GroupMonitorDir,"MonitorGroup_Info_Deformation.out")
		
	elif monitorOutput in ["chordRotation", "chordDeformations"]:
		"""
		3D ELements: eps, thetaZ_1, thetaZ_2, thetaY_1, thetaY_2, thetaX;
		2D Elements: eps, theta_1, theta_2
		"""
		chordOutput = [1,2,3,4,5,6]  
		MonitorDefFile = os.path.join(GroupMonitorDir,"MonitorGroup_Rotation.out")
		MonitorEleFile = os.path.join(GroupMonitorDir,"MonitorGroup_Tags_Rotation.out")
		MonitorInfoFile = os.path.join(GroupMonitorDir,"MonitorGroup_Info_Rotation.out")
		
	else:
		print("The output quantity is not recognized. Use 'deformation' or 'chordRotation' based on the user input in createODB command earlier.")
		
		
	# ops.recorder('Element', '-file', MonitorDefFile,  '-time', '-dT', deltaT, '-ele', *monitorEleTags, '-dof',*dofList, monitorOutput)
	MonitorEleDef = np.transpose(np.loadtxt(MonitorDefFile, dtype=float, delimiter=None, converters=None, unpack=True))		# Returns the element deformations
	MonitorEleForce = np.transpose(np.loadtxt(MonitorForceFile, dtype=float, delimiter=None, converters=None, unpack=True))		# Returns the element localForce
	MonitorEleTags = np.transpose(np.loadtxt(MonitorEleFile, dtype=float, delimiter=None, converters=None, unpack=True))	# Returns the element tags
	MonitorEleInfo = np.transpose(np.loadtxt(MonitorInfoFile, dtype=float, delimiter=None, converters=None, unpack=True))  # Returns array of DOFs recorded
		
	return MonitorEleDef, MonitorEleForce, MonitorEleTags, MonitorEleInfo
	
	

def _elementMonitorCheck(eleTag, dof, monitorOutput, limitStates, limStateColors, MonitorEleInfo, MonitorEleTags, MonitorEleDef, tStep):
	eleMonitorColor = "blue"
	dofLength = len(MonitorEleInfo)
	
	if monitorOutput in ["deformation", "deformations"]:
	
		iMonitor, = np.where(MonitorEleInfo == float(dof))				# Get the DOF column number to read for each element
		jMonitor, = np.where(MonitorEleTags == float(eleTag))			# Get the element number
		readCol = 1 + iMonitor + dofLength*jMonitor						# First column is time or load factor
		# print(iMonitor, MonitorEleTags, float(eleTag), readCol)
		for kk in range(1,int(tStep)):
			# print(tStep, kk)
			if abs(MonitorEleDef[kk, readCol]) > abs(MonitorEleDef[kk-1, readCol]):
				# print(abs(MonitorEleDef[kk, readCol]))
				for lim in range(0,4):
					if limitStates[lim] == 0:
						eleMonitorColor = "solid"
					elif abs(MonitorEleDef[kk, readCol]) >= limitStates[lim]:
						eleMonitorColor = limStateColors[lim]
					else:
						pass
			else:
				pass
							
		return eleMonitorColor
	
	elif monitorOutput in ["chordRotation", "chordDeformations"]:
	
		iMonitor, = np.where(MonitorEleInfo == float(dof))				# Get the DOF column number to read for each element
		jMonitor, = np.where(MonitorEleTags == float(eleTag))			# Get the element number
		readCol = 1 + iMonitor + dofLength*jMonitor						# First column is time or load factor
		
		pass
	

def _elementHysteresis(eleTag, dof, MonitorEleInfo, MonitorEleTags, MonitorEleDef, MonitorEleForce):

	dofLength = len(MonitorEleInfo)
	
	iMonitor, = np.where(MonitorEleInfo == float(dof))				# Get the DOF column number to read for each element
	jMonitor, = np.where(MonitorEleTags == float(eleTag))			# Get the element number
	readCol = 1 + iMonitor + dofLength*jMonitor						# First column is time or load factor
	readForceCol = 1 + iMonitor + (1*dofLength)*jMonitor			# First column is time or load factor, 6 columns for each elements
	# print(iMonitor, MonitorEleTags, float(eleTag), readCol)
	
	eleDeformation = MonitorEleDef[:, readCol]
	eleForce = MonitorEleForce[:, readForceCol]
		
	return eleDeformation, eleForce

						
def _readShellStressData(ModelName,LoadCaseName):
	
	"""
	"""	
	
	ODBdir = ModelName+"_ODB"		# ODB Dir name
	LoadCaseDir = os.path.join(ODBdir, LoadCaseName)
	
	# Get number of nodes in the model to set a node displacement array
	nodes,elements, eleClassTags = _readNodesandElements(ModelName)
	Nnodes = len(nodes)
	ndm = len(nodes[0,1:])
	
	# Read shell element tags saved in ODB
	
	Elements_shell_4N4GP = os.path.join(ODBdir, "Elements_shell_4N4GP.out")
	
	if len(Elements_shell_4N4GP)>0:
		"""
		Read 4Node 4GP shell data if exists
		extrapolate stress from GP to nodes
		Calculate average stress at each node from elements
		Save node stress data for plotting
		"""
		
		EleStressFile_shell_4N4GP = os.path.join(LoadCaseDir,"EleStress_shell_4N4GP.out")
		GPstress_shell_4N4GP = np.transpose(np.loadtxt(EleStressFile_shell_4N4GP, dtype=float, delimiter=None, converters=None, unpack=True))
		
		t_row = 0 # get this dynamically from tstep input 
		ii=0
		for ii, ele in enumerate(Elements_shell_4N4GP):
			##### initialize a 2D array to store GP stress data for matrix calcs
			_thisGPstress = np.zeros([4,32])
			_thisGPstress = np.array([[GPstress_shell_4N4GP[t_row, 1:9+ii]],
										[GPstress_shell_4N4GP[t_row, 10+ii:18+ii]],
										[GPstress_shell_4N4GP[t_row, 19+ii:27+ii]],
										[GPstress_shell_4N4GP[t_row, 28+ii:36+ii]]])		
	
	
	
	# each element has 3 stress vectors
	# create a [Nnodes x Nele] array to store node stress from each element 
	# Each row is for a node. Store stress from element j at node i in the array location ij
	
	nodeStressArray = np.zeros([Nnodes,Nele])
	
	
	
	timeSteps = Disp[:,0]
	Ntime = len(Disp[:,0])

	tempDisp = np.zeros([Ntime,Nnodes,ndm])
	tempDisp[:,:,0] = Disp[:,1::ndm]
	tempDisp[:,:,1] = Disp[:,2::ndm]
	
	if ndm == 3:
		tempDisp[:,:,2] = Disp[:,3::ndm]
		
	nodes_displacement = tempDisp
	
	return timeSteps, nodes_displacement
