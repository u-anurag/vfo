##########################################################################################
## *Only user functions are to be included in this file. For helper/internal functions	##
## go to internal_database_functions.py and internal_plotting_functions.                ##
## *As of now, this procedure does not work for 9 and 20 node brick elements, and       ##
## tetrahedron elements.																##
##																						##
##																						##
## Created By - Anurag Upadhyay, University of Utah. https://github.com/u-anurag		##
##            									                                        ##
## 																						##
##########################################################################################

# Check if the script is executed on Jupyter Notebook Ipython. 
# If yes, force inline, interactive backend for Matplotlib.

import sys
import os
import matplotlib
import math

from math import asin
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

for line in range(0,len(sys.argv)):
    if "ipykernel_launcher.py" in sys.argv[line]:
        # matplotlib.use('nbagg')
        pv.set_jupyter_backend('panel')
        break
    else:
        pass

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import matplotlib.animation as animation
from matplotlib.widgets import Slider


#### CHANGE THESE BEFORE COMMITTING, call them before calling openseespy here ####
import vfo.classTags as classTags
import vfo.internal_database_functions as idbf
import vfo.internal_plotting_functions as ipltf
import openseespy.opensees as ops

#####################################################################
####    All the plotting related definitions start here.
####
#####################################################################

def createODB(model=None,loadcase=None, Nmodes=0, deltaT=0.0):
	
	"""
	This function creates a directory to save all the output data.

	Command: createODB(model="ModelName",<loadcase="LoadCase Name">, <Nmodes=Nmodes>, <recorders=*recorder>)
	
	PARAMETERS
	-----------
	
	model    : (string) 
		Name of the model. The main output folder will be named "ModelName_ODB" in the current directory.
		
	loadcase : (string), Optional. 
		Name of the load case forder to be created inside the ModelName_ODB folder. If not provided,
					no load case data will be read.
					
	Nmodes	 : (int) Optional 
		key argument to save modeshape data. Default is 0, no modeshape data is saved.
	
	deltaT	 : (float) Optional 
		time interval for recording. will record when next step is deltaT greater than last recorder step. 
					(default: records at every time step)
	
	recorders	: (string) Optional
		A list of additional quantities a users would like to record in the output database.
					The arguments for these additional inputs match the standard OpenSees arguments to avoid any confusion.
					'localForce','basicDeformation', 'plasticDeformation','stresses','strains'
					The recorders for node displacement and reactions are saved by default to help plot the deformed shape.
	
	Example: createODB(model="TwoSpanBridge", loadcase="Pushover", Nmodes=3, recorders=['stresses', 'strains'])
	
	Future: The integrationPoints output works only for nonlinear beam column elements. If a model has a combination 
			of elastic and nonlienar elements, we need to create a method distinguish. 
	
	"""
	
	ModelName = model
	ODBdir = ModelName+"_ODB"		# ODB Dir name
	if not os.path.exists(ODBdir):
			os.makedirs(ODBdir)

	nodeList = ops.getNodeTags()
	eleList = ops.getEleTags()
	
	## get individual lists of elements for stress and strain recorder 
	# for ele in eleList:
		
	
	if len(ops.nodeCoord(nodeList[0])) == 2:
		dofList_node = [1, 2]
		dofList_ele  = [1, 2, 3]				# By default records all 3 dofs
	if len(ops.nodeCoord(nodeList[0])) == 3:
		dofList_node = [1, 2, 3]
		dofList_ele  = [1, 2, 3, 4, 5, 6]		# By default records all 6 dofs. The dof recorder in 3D doesn't work with random dof input.
	
	# Save node and element data in the main Output folder
	idbf._saveNodesandElements(ModelName)
	
	#########################
	## Create mode shape dir
	#########################
	if Nmodes > 0:
		ModeShapeDir = os.path.join(ODBdir,"ModeShapes")
		if not os.path.exists(ModeShapeDir):
			os.makedirs(ModeShapeDir)
			
		## Run eigen analysis internally and get information to print
		Tarray = np.zeros([1,Nmodes])  # To save all the periods of vibration
		ops.wipeAnalysis()
		eigenVal = ops.eigen(Nmodes+1)
	
		for m in range(1,Nmodes+1):
			Tarray[0,m-1]=4*asin(1.0)/(eigenVal[m-1])**0.5
		
		modeTFile = os.path.join(ModeShapeDir, "ModalPeriods.out")
		np.savetxt(modeTFile, Tarray, delimiter = ' ', fmt = '%.5e')   
		
		### Save mode shape data
		for i in range(1,Nmodes+1):
			idbf._saveModeShapeData(ModelName,i)
		
		ops.wipeAnalysis()
		
	# Define standard outout filenames			
		
	if loadcase is not None:
		LoadCaseName = loadcase
		LoadCaseDir = os.path.join(ODBdir,LoadCaseName)

		if not os.path.exists(LoadCaseDir):
			os.makedirs(LoadCaseDir)
			
		eleList_shell_4N4GP, eleList_shell_3N3GP, eleList_tetra_4N4GP, eleList_brick_8N8GP = idbf._getEleListforStress()
		
		delim = ' '
		dtype ='float32' 
		ftype = '.out'
		fmt = '%.5e'
		
		ele_shell_4N4GP_File = os.path.join(ODBdir, 'Elements' + "_shell_4N4GP" + ftype)
		ele_shell_3N3GP_File = os.path.join(ODBdir, 'Elements' + "_shell_3N3GP" + ftype)
		ele_tetra_4N4GP_File = os.path.join(ODBdir, 'Elements' + "_tetra_4N4GP" + ftype)
		ele_brick_8N8GP_File = os.path.join(ODBdir, 'Elements' + "_brick_8N8GP" + ftype)
		
		# Save element arrays for stress and recorder
		np.savetxt(ele_shell_4N4GP_File, eleList_shell_4N4GP, delimiter = delim, fmt = fmt)
		np.savetxt(ele_shell_3N3GP_File, eleList_shell_3N3GP, delimiter = delim, fmt = fmt)
		np.savetxt(ele_tetra_4N4GP_File, eleList_tetra_4N4GP, delimiter = delim, fmt = fmt)
		np.savetxt(ele_brick_8N8GP_File, eleList_brick_8N8GP, delimiter = delim, fmt = fmt)
	
		
		NodeDispFile = os.path.join(LoadCaseDir,"NodeDisp_All.out")
		EleForceFile = os.path.join(LoadCaseDir,"EleForce_All.out")
		ReactionFile = os.path.join(LoadCaseDir,"Reaction_All.out")
		EleStressFile_shell_4N4GP = os.path.join(LoadCaseDir,"EleStress_shell_4N4GP.out")
		EleStrainFile_shell_4N4GP = os.path.join(LoadCaseDir,"EleStrain_shell_4N4GP.out")
		EleStressFile_shell_3N3GP = os.path.join(LoadCaseDir,"EleStress_shell_3N3GP.out")
		EleStrainFile_shell_3N3GP = os.path.join(LoadCaseDir,"EleStrain_shell_3N3GP.out")
		EleStressFile_tetra_4N4GP = os.path.join(LoadCaseDir,"EleStress_tetra_4N4GP.out")
		EleStrainFile_tetra_4N4GP = os.path.join(LoadCaseDir,"EleStrain_tetra_4N4GP.out")
		EleStressFile_brick_8N8GP = os.path.join(LoadCaseDir,"EleStress_brick_8N8GP.out")
		EleStrainFile_brick_8N8GP = os.path.join(LoadCaseDir,"EleStrain_brick_8N8GP.out")
		EleBasicDefFile = os.path.join(LoadCaseDir,"EleBasicDef_All.out")
		ElePlasticDefFile = os.path.join(LoadCaseDir,"ElePlasticDef_All.out")
# 		EleIntPointsFile = os.path.join(LoadCaseDir,"EleIntPoints_All.out")	

		# Save recorders in the ODB folder
		ops.recorder('Node', '-file', NodeDispFile,  '-time', '-dT', deltaT, '-node', *nodeList, '-dof',*dofList_node, 'disp')
		ops.recorder('Node', '-file', ReactionFile,  '-time', '-dT', deltaT, '-node', *nodeList, '-dof',*dofList_node, 'reaction')
		
		if len(eleList_shell_4N4GP)>0:
			ops.recorder('Element', '-file', EleStressFile_shell_4N4GP,  '-time', '-dT', deltaT, '-ele', *eleList_shell_4N4GP,  'stresses')
			ops.recorder('Element', '-file', EleStrainFile_shell_4N4GP,  '-time', '-dT', deltaT, '-ele', *eleList_shell_4N4GP,  'strains')
		
		if len(eleList_shell_3N3GP)>0:
			ops.recorder('Element', '-file', EleStressFile_shell_3N3GP,  '-time', '-dT', deltaT, '-ele', *eleList_shell_3N3GP, 'stresses')
			ops.recorder('Element', '-file', EleStrainFile_shell_3N3GP,  '-time', '-dT', deltaT, '-ele', *eleList_shell_3N3GP, 'strains')
		
		if len(eleList_tetra_4N4GP)>0:
			ops.recorder('Element', '-file', EleStressFile_tetra_4N4GP,  '-time', '-dT', deltaT, '-ele', *eleList_tetra_4N4GP, 'stresses')
			ops.recorder('Element', '-file', EleStrainFile_tetra_4N4GP,  '-time', '-dT', deltaT, '-ele', *eleList_tetra_4N4GP, 'strains')
		
		if len(eleList_brick_8N8GP)>0:
			ops.recorder('Element', '-file', EleStressFile_brick_8N8GP,  '-time', '-dT', deltaT, '-ele', *eleList_brick_8N8GP, 'stresses')
			ops.recorder('Element', '-file', EleStrainFile_brick_8N8GP,  '-time', '-dT', deltaT, '-ele', *eleList_brick_8N8GP, 'strains')
		

	else:
		print(">>> VFO Warning: No 'loadcase' argument was provided to save output.<<<")
		
		
def readODB(*argv):
	
	"""
	This function reads saved data from a directory.
	
	Command: readODB("ModelName",<"LoadCase Name">)
	
	ModelName    : (string) Name of the model. The main output folder will be named "ModelName_ODB" in the current directory.
	LoadCase Name: (string), Optional. Name of the load case forder to be created inside the ModelName_ODB folder. If not provided,
					no load case data will be read.
    
	"""
    
	ModelName = argv[0]
	ODBdir = ModelName+"_ODB"		# ODB Dir name

	# Read node and element data in the main Output folder
	nodes, elements, eleClassTags = idbf._readNodesandElements(ModelName)
	
	if len(argv)>=2:
		LoadCaseName = argv[1]
		LoadCaseDir = os.path.join(ODBdir, LoadCaseName)

		if not os.path.exists(LoadCaseDir):
			print("No database found")
		
		# Define standard outout filenames
		NodeDispFile = os.path.join(LoadCaseDir,"NodeDisp_All.out")
		EleForceFile = os.path.join(LoadCaseDir,"EleForce_All.out")
		ReactionFile = os.path.join(LoadCaseDir,"Reaction_All.out")
		
		# Read recorders in the ODB folder
		# FUTURE: Gives warning if the files are empty. Create a procedure to check if files are empty.
		NodeDisp = np.loadtxt(NodeDispFile,delimiter=' ')
		EleForce = np.loadtxt(EleForceFile,delimiter=' ')   
		Reaction = np.loadtxt(ReactionFile,delimiter=' ')
		
		return nodes, elements, NodeDisp, Reaction, EleForce
	
	else:
		return nodes, elements




		
"""
Change PyVista theme
https://docs.pyvista.org/examples/02-plot/themes.html
"""

pv.set_plot_theme("document")
# pv.set_jupyter_backend('panel')
#pv.set_jupyter_backend('none')



def _get_modelDisplay(nodeArray, elementArray, eleClassTags):
    
        
    vertices = nodeArray[:,1:]     # point coordinates
    nodeTags = nodeArray[:,0].astype(int) 
	
    eleTagCoord = np.empty(0, dtype=int)    # To store coordinates for element tags
    
    def _getNodeIndex(nodeTag):
        nodeIndex, = np.where(nodeTags[:] == int(nodeTag))
        return int(nodeIndex)

    
    def _eleClassTag(eleTag):
        """
        Returns element class tag, MAY BE MOVE ALL INTERNAL DEF TO A SEPARATE FILE
        """
        try:
            i, = np.where(eleClassTags[:,0] == int(eleTag))
            return eleClassTags[i,1]
        except:
            return eleClassTags[1]
    
    
    
    def _getPV_2NodeELements(elementArray):
        lineElements = np.empty(0, dtype=int)  # Initiate line element array
        for ii in elementArray:
            # print(ii)
            if len(ii)==3:
                # print("ok")
                tempArray = np.array([2, _getNodeIndex(ii[1]), _getNodeIndex(ii[2])])
                # print(tempArray)
			
                lineElements = np.hstack((lineElements,tempArray))

        return lineElements.astype(int)

    
    
    def _getPV_surfaces(elementArray, eleClassTags):
        """
        ALL POSSIBLE SURFACES SHOULD BE ADDED TO ONE np.hstack() FOR PLOTTING

        """
        faces = np.empty(0, dtype=int)
        
        for ii in elementArray:
            """
            _eleClassTag procedure is ALSO defined in main vfo file
            """
            eleTag = int(ii[0])
            classTag = _eleClassTag(ii[0])
            # print(classTag, type(classTag))
            tempArray = np.empty(0, dtype=int)

            if _eleClassTag(ii[0]) in classTags.triNodeEleTags:
                # print("triNodeEleTags")
                tempArray = np.array([3, _getNodeIndex(ii[1]), _getNodeIndex(ii[2]), _getNodeIndex(ii[3])])

            if _eleClassTag(ii[0]) in classTags.fourNodeEleTags:
                tempArray = np.array([4, _getNodeIndex(ii[1]), _getNodeIndex(ii[2]), _getNodeIndex(ii[3]), _getNodeIndex(ii[4])])

            if _eleClassTag(ii[0]) in classTags.MVLEMEleTags:
                #print("mvlem3d")
                tempArray = np.array([4, _getNodeIndex(ii[1]), _getNodeIndex(ii[2]), _getNodeIndex(ii[4]), _getNodeIndex(ii[3])])

            if _eleClassTag(ii[0]) in classTags.tetEleTags:
                tempArray = np.array([3, _getNodeIndex(ii[1]), _getNodeIndex(ii[2]), _getNodeIndex(ii[3]),
                                      3, _getNodeIndex(ii[1]), _getNodeIndex(ii[2]), _getNodeIndex(ii[4]),
                                      3, _getNodeIndex(ii[1]), _getNodeIndex(ii[3]), _getNodeIndex(ii[4]),
                                      3, _getNodeIndex(ii[2]), _getNodeIndex(ii[3]), _getNodeIndex(ii[4])])

            if _eleClassTag(ii[0]) in classTags.eightNodeBrickEleTags:
                tempArray = np.array([4, _getNodeIndex(ii[1]), _getNodeIndex(ii[2]), _getNodeIndex(ii[3]), _getNodeIndex(ii[4]),
                                      4, _getNodeIndex(ii[5]), _getNodeIndex(ii[6]), _getNodeIndex(ii[7]), _getNodeIndex(ii[8]),
                                      4, _getNodeIndex(ii[1]), _getNodeIndex(ii[2]), _getNodeIndex(ii[6]), _getNodeIndex(ii[5]),
                                      4, _getNodeIndex(ii[4]), _getNodeIndex(ii[3]), _getNodeIndex(ii[7]), _getNodeIndex(ii[8]),
                                      4, _getNodeIndex(ii[1]), _getNodeIndex(ii[4]), _getNodeIndex(ii[8]), _getNodeIndex(ii[5]),   
                                      4, _getNodeIndex(ii[2]), _getNodeIndex(ii[3]), _getNodeIndex(ii[7]), _getNodeIndex(ii[6]),
                                     ])    


            faces = np.hstack((faces,tempArray))
            # print(faces)

        return faces

    
    lines = _getPV_2NodeELements(elementArray)
    surf = _getPV_surfaces(elementArray, eleClassTags)
    
    mesh = pv.PolyData(vertices, surf)
    mesh_lines = pv.PolyData(vertices)
    mesh_lines.lines = lines
	
    return mesh, mesh_lines, vertices, nodeTags

def _roundNearestFive(num):
    nearestFive = np.around(np.floor(num)/5, decimals=0)*5
    # print(np.floor(num), nearestFive)
    if abs(np.floor(num))> abs(nearestFive):
        return (abs(nearestFive)+5)*num/abs(num)
    else:
        return abs(nearestFive)*num/abs(num)
		
def _nodecoords(nodetag, nodeArray):
	"""
	Returns an array of node coordinates: works like nodeCoord() in opensees.
	"""
	i, = np.where(nodeArray[:,0] == float(nodetag))
	return nodeArray[int(i),1:]



def _getEleTagCoord(elementArray, nodeArray):
	"""
	coordinate for element tags
	"""
	
	eleTagCoord = np.empty([len(elementArray), 4])  # Each row of this array will store [eleTag, x, y, z]
	# eleTagCoord = np.empty(0, dtype=int)    # To store coordinates for element tags
	count_ii = 0
	for ii in elementArray:
		tempCoordArray = np.empty([len(ii)-1,3])  # Array of shape [numNodes x 3] to store x, y, z coordinate of each node in this element 
		count_jj = 0
		# testArray = ii
		# print(testArray)
		for jj in ii[1:]:
			# print("jj :", jj)
			# print("jj type :", type(jj))
			tempCoordArray[count_jj] = _nodecoords(jj, nodeArray)
			count_jj+=1
			
		# print(tempCoordArray)
		thisTag = ii[0]
		# aveX = np.average(tempCoordArray[:,0])
		# aveY = np.average(tempCoordArray[:,1])
		# aveZ = np.average(tempCoordArray[:,2])
		# print(thisTag, aveX, aveY, aveZ)
		
		eleTagCoord[count_ii] = np.array([ii[0], np.average(tempCoordArray[:,0]), np.average(tempCoordArray[:,1]), np.average(tempCoordArray[:,2])])
		count_ii+=1
	
	# print(eleTagCoord)
		
	return eleTagCoord
	

def plot_model(model="none",show_nodes="no",show_nodetags="no",show_eletags="no",font_size=10,setview="3D",elementgroups=None,line_width=1, filename=None):
    
	"""
	a function to plot the model of structure.
	
	Command:
	
		plot_model(model="none",show_nodes="no",show_nodetags="no",show_eletags="no",font_size=10,setview="3D",elementgroups=None,line_width=1)
	
	INPUT:
		model		: (string)
			name of the model output database as used in createODB() function. If no name is provided, the function tries to get data from the active model. 
							Default is "none".
							
		show_nodes  : (string)
			Renders nodes as spheres. Default is "no".
							
		show_nodetags	: (string)
			Displays node tags if "yes". Default is "no".
			
		show_eletags	: (string)
			Displays element tags if "yes". Default is "no".
			
		font_size	: (int)
			Size of tag font. Default is 10.
			
		setview		: (str)
			sets the camera view to predefined angles. Valid entries are "xy","yz","xz","3D", or a list with [x,y,z] unit vector. Default is "3D".
			
		elementgroups : (list)
			List of lists of elements of groups and respective colors.
			
		line_width  : (float)
			Line width of the rendered beam-column elements.
			
		filename	: (str)
			Filename to save an image of the deformed shape.
		
		
	"""
		

	# Check if their is an output database or not.
	if model == "none":
		print("No Model_ODB specified, trying to get data from the active model.")
		try:
			this_nodeArray, elementArray, eleClassTags = idbf._getNodesandElements()
		except:
			raise Exception(">>>> VFO ERROR: No Model_ODB specified. No active model found.")
	else:
		print("Reading data from the "+model+"_ODB.")
		try:
			this_nodeArray, elementArray, eleClassTags = idbf._readNodesandElements(model) ## Add eleClassTags to the Tcl data
			# print(np.shape(nodeArray))
		except:
			raise Exception(">>>> VFO ERROR: No Model_ODB found. Exiting 'plot_model()' function now.<<<<")

	# print("this_nodeArray", this_nodeArray)
	# print("shape this_nodeArray", np.shape(this_nodeArray))
	# print("elementArray", elementArray)
	# print("eleClassTags", eleClassTags)
	
	## Check if the model is 2D or 3D
	nodeArray = np.zeros([len(this_nodeArray[:,0]), 4])
	ndm = len(this_nodeArray[0,:]) -1
	if ndm == 2:
		for ii in range(0,len(this_nodeArray[:,0])):
			nodeArray[ii,0:3] = this_nodeArray[ii,:]
	else:
		nodeArray = this_nodeArray
			
			
	
	mesh_original, mesh_lines_original, vertices, nodeTags = _get_modelDisplay(nodeArray, elementArray, eleClassTags)

	pl = pv.Plotter()
	pl.show(interactive_update=True)
	
	if show_nodes == "yes":
		point_size = 5.0
		spheres=True
	else:
		point_size = 0.0
		spheres=False
		
	pl.add_mesh(mesh_original, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, color="green", opacity=1.0, render_lines_as_tubes=False, line_width=1)  # "show_edges=True" will NOT work with PlotterITK
	pl.add_mesh(mesh_lines_original, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, color="green", opacity=1.0, render_lines_as_tubes=False, line_width=line_width)  # "show_edges=True" will NOT work with PlotterITK
	
	eleTagCoord = _getEleTagCoord(elementArray, nodeArray) 
	
	## Get nodes, eleClassTags for the elements in each group
	if elementgroups is not None:
		for gg, thisgroup in enumerate(elementgroups[0]):
			thisgroup_elementArray = [None]*len(thisgroup)
			thisgroup_color = elementgroups[1][gg]
			
			for ii, ele in enumerate(thisgroup):
				i, = np.where(eleTagCoord[:,0] == float(ele))
				thisgroup_elementArray[ii] = elementArray[int(i)]
				
			mesh_thisgroup, mesh_lines_thisgroup, vertices_thisgroup, nodeTags_thisgroup = _get_modelDisplay(nodeArray, thisgroup_elementArray, eleClassTags)
			pl.add_mesh(mesh_thisgroup, render_points_as_spheres=False, point_size=point_size, show_edges=False, color=thisgroup_color, opacity=1.0, render_lines_as_tubes=False, line_width=1)  # "show_edges=True" will NOT work with PlotterITK
			pl.add_mesh(mesh_lines_thisgroup, render_points_as_spheres=False, point_size=point_size, show_edges=False, color=thisgroup_color, opacity=1.0, render_lines_as_tubes=True, line_width=line_width)  # "show_edges=True" will NOT work with PlotterITK
		
		
		
		
	if setview in ["XY","xy"]:
		pl.view_xy()
	elif setview in ["YZ","yz"]:
		pl.view_yz()
	elif setview in ["XZ","xz"]:
		pl.view_xz()
	elif setview == "3D":
		pl.view_isometric()
	else:
		pl.set_viewup(setview)

	## Overwrite setview option to "xy" if the model is 2D
	if ndm == 2:
		pl.view_xy()
	else:
		pass
	 
	pl.add_text('VFO - Visualization for OpenSees', position='lower_left', color='green', font_size=6)
	
	if show_nodetags=="yes":
		pl.add_point_labels(vertices, nodeTags.astype(int), point_size=1,render_points_as_spheres=True, font_size=font_size ,shape_color="white", shape_opacity=0.2, render=True,always_visible=True)

	if show_eletags=="yes":
		# print(eleTagCoord[:,1:], "and", eleTagCoord[:,0])
		pl.add_point_labels(eleTagCoord[:,1:], eleTagCoord[:,0].astype(int), font_size=font_size,shape_color="gray", shape_opacity=0.2, render=True,always_visible=True)

	if filename is not None:
		pl.screenshot(filename+".png") 
		
	pl.show()




def plot_deformedshape(model="none",loadcase="none",scale=10,tstep=-1,overlap="no",contour="none",setview="3D", line_width=1, contourlimits=None, filename=None):
    
	"""
	a function to plot deformed shape of structure.
	
	Command:
	
		plot_deformedshape(model="none",loadcase="none",scale=10,tstep=-1,overlap="no",contour="none",setview="3D", line_width=1, contourlimits=None, filename=None)
	
	INPUT:
		model		: (string)
			Name of the model output database as used in createODB() function. Default is "none".
			
		loadcase	: (string)
			Name of the load case to save output data.
			
		scale		: (int)
			Scale factor to be applied to the deformed shape. Default is 10.
			
		tstep		: (float)
			Time step or analysis step at which the deformed shape is to be plotted.
			
		contour		: (str)
			Contours of displacement in x, y, or z. Default is "none".
			
		contourlimits: (list)
			A list of minimum and maximum limits of the displacement contours.
		
		setview		: (str)
			sets the camera view to predefined angles. Valid entries are "xy","yz","xz","3D", or a list with [x,y,z] unit vector. Default is "3D".
			
		line_width  : (int)
			Line thickness for the beam-column elements. Default is 1.
			
		filename	: (str)
			Filename to save an image of the deformed shape.
		
	"""
	print("Reading structure data from "+model+"_ODB.")
	try:
		this_nodeArray, elementArray, eleClassTags = idbf._readNodesandElements(model) ## Add eleClassTags to the Tcl data
		# print(np.shape(this_nodeArray))
	except:
		raise Exception(">>>> VFO ERROR: No Model_ODB found. No active model found. Exiting 'plot_modeshape()' now.<<<<")
        
	print("Reading loadcase '"+loadcase+"' data from "+model+"_ODB.")
	
    
	# timeSteps, this_displacement_nodeArray = idbf._readNodeDispData(model,loadcase) ## Add eleClassTags to the Tcl data


	try:
		timeSteps, this_displacement_nodeArray, singleStep = idbf._readNodeDispData(model,loadcase) ## Add eleClassTags to the Tcl data
		# print('timeSteps = ', timeSteps)
		# print(np.shape(this_displacement_nodeArray[0,:,:]))
		# print("this_displacement_nodeArray = ", this_displacement_nodeArray)
		# print(this_displacement_nodeArray[0,:,:])
	except:
		raise Exception(">>>> VFO ERROR: No node displacements from loadcase was found in "+model+"_ODB. Exiting 'plot_deformedshape()' now. <<<<")
		
	## Check if the model is 2D or 3D
	nodeArray = np.zeros([len(this_nodeArray[:,0]), 4])
	displacement_nodeArray = np.zeros([np.shape(this_displacement_nodeArray)[0],np.shape(this_displacement_nodeArray)[1], 3])
	ndm = len(this_nodeArray[0,:]) -1
	if ndm == 2:
		for ii in range(0,len(this_nodeArray[:,0])):
			nodeArray[ii,0:3] = this_nodeArray[ii,:]
			
		for ii in range(0,np.shape(this_displacement_nodeArray)[0]):
			for jj in range(0,np.shape(this_displacement_nodeArray)[1]):
				# print("this_displacement_nodeArray ", this_displacement_nodeArray[ii,jj,:])
				displacement_nodeArray[ii,jj,0:2] = this_displacement_nodeArray[ii,jj,:]
	else:
		nodeArray = this_nodeArray
		displacement_nodeArray = this_displacement_nodeArray
		
	# print("displacement_nodeArray = ", displacement_nodeArray)
	
	pl = pv.Plotter()		
	pl.show(interactive_update=True)


	if tstep == -1:
		jjj = len(timeSteps)-1    #  May 20, 2023    len(timeSteps)-1
		print("Final deformed shape")
	else:
		jjj = (np.abs(timeSteps - tstep)).argmin()			# index closest to the time step requested.
		print("Deformation at time: " + str(round(timeSteps[jjj], 2)))
				
				
				
	def _get_deformed_mesh(tstep):
	
		if tstep == -1:
			jj = len(timeSteps)-1    #  May 2    len(timeSteps)-1
			printLine = "Final deformed shape"
		else:
			jj = (np.abs(timeSteps - tstep)).argmin()			# index closest to the time step requested.
			if timeSteps[-1] < tstep:
				print(">>> VFO Warning: plot_deformedshape: Time-Step has exceeded maximum analysis time step. <<<")
			printLine = "Deformation at time: " + str(round(timeSteps[jj], 2))
						
		DeflectedNodeCoordArray = nodeArray[:,1:]+ scale*displacement_nodeArray[jj,:,:]
		
		DeflectedNodeArray = np.hstack((nodeArray[:,0].reshape(len(nodeArray[:,0]),1),DeflectedNodeCoordArray))
		
		mesh_deflected, mesh_lines_deflected, vertices, nodeTags = _get_modelDisplay(DeflectedNodeArray, elementArray, eleClassTags)
		
		# try:
			# pl.remove_actor("thisMesh")
		# except:
			# pass
			
		point_size = 0.0
		spheres=False
		
		Clim = contourlimits
		
		if contour == "none":
			pl.add_mesh(mesh_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=1, color="green", opacity=1.0,name="thisMesh1")  # "show_edges=True" will NOT work with PlotterITK
			pl.add_mesh(mesh_lines_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=line_width, color="green", opacity=1.0,name="thisMesh2")  # "show_edges=True" will NOT work with PlotterITK
		elif contour in ["x","X"]:
			x = displacement_nodeArray[jj,:,0]
			pl.add_mesh(mesh_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=1, scalars=x, clim=Clim, opacity=1.0,name="thisMesh1")  # "show_edges=True" will NOT work with PlotterITK
			pl.add_mesh(mesh_lines_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=line_width, scalars=x, clim=Clim, opacity=1.0,name="thisMesh2")  # "show_edges=True" will NOT work with PlotterITK
		elif contour in ["y","Y"]:
			y = displacement_nodeArray[jj,:,1]
			pl.add_mesh(mesh_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=1, scalars=y, clim=Clim, opacity=1.0,name="thisMesh1")  # "show_edges=True" will NOT work with PlotterITK
			pl.add_mesh(mesh_lines_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=line_width, scalars=y, clim=Clim, opacity=1.0,name="thisMesh2")  # "show_edges=True" will NOT work with PlotterITK
		elif contour in ["z","Z"]:
			z = displacement_nodeArray[jj,:,2]
			pl.add_mesh(mesh_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=1, scalars=z, clim=Clim, opacity=1.0,name="thisMesh1")  # "show_edges=True" will NOT work with PlotterITK
			pl.add_mesh(mesh_lines_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=line_width, scalars=z, clim=Clim, opacity=1.0,name="thisMesh2")  # "show_edges=True" will NOT work with PlotterITK
			
			# ALTERNATIVELY, USE https://docs.pyvista.org/api/plotting/_autosummary/pyvista.Plotter.view_xy.html
    
		if setview in ["XY","xy"]:
			pl.view_xy()
		elif setview in ["YZ","yz"]:
			pl.view_yz()
		elif setview in ["XZ","xz"]:
			pl.view_xz()
		elif setview == "3D":
			pl.view_isometric()
		else:
			pl.set_viewup(setview)
		
		## Overwrite setview option to "xy" if the model is 2D
		if ndm == 2:
			pl.view_xy()
		else:
			pass

	pl.add_text('VFO - Visualization for OpenSees', position='lower_left', color='green', font_size=6)
	
	if singleStep is False:
		slider = pl.add_slider_widget(_get_deformed_mesh, [0, timeSteps[-1]], value=timeSteps[jjj], title="time steps", title_opacity=0.5, title_color="red", fmt="%0.2f", title_height=0.008,)
	else:
		_get_deformed_mesh(-1)
		
	if filename is not None:
		pl.screenshot(filename+".png") 
	
	# if overlap=="yes":
		# mesh_original, vertices, nodeTags = _get_modelDisplay(nodeArray, elementArray, eleClassTags)
		# pl.add_mesh(mesh_original, show_edges=True, color="gray", opacity=0.3)  # "show_edges=True" will NOT work with PlotterITK
    
	

	point_size = 1
	#font_size = 10
    
	# if nodes == "yes":
		# point_size = 8
		# pl.add_mesh(vertices, color='black', point_size=point_size, render_points_as_spheres=True)
        
	# if node_tags=="yes":
		# pl.add_point_labels(vertices, nodeTags, point_size=1,render_points_as_spheres=True, font_size=font_size ,shape_color="white", shape_opacity=0.2, render=True,always_visible=True)

	pl.show()


def plot_modeshape(model="none",modenumber=1,scale=10,overlap="yes",contour="none",setview="3D", line_width=1, contourlimits=None, filename=None):
    
	"""
	
	Command:
	
		plot_modeshape(model="none",modenumber=1,scale=10,overlap="no",contour="none",setview="3D", contourlimits=None, filename=None)
	
	INPUT:
		model		: (string)
			Name of the model output database as used in createODB() function. Default is "none".
			
		loadcase	: (string)
			Name of the load case to save output data.
			
		scale		: (int)
			Scale factor to be applied to the deformed shape. Default is 10.
		
		contour		: (str)
			Contours of displacement in x, y, or z. Default is "none".
			
		contourlimits: (list)
			A list of minimum and maximum limits of the displacement contours.
		
		setview		: (str)
			sets the camera view to predefined angles. Valid entries are "xy","yz","xz","3D", or a list with [x,y,z] unit vector. Default is "3D".
			
		line_width  : (int)
			Line thickness for the beam-column elements. Default is 1.
			
		filename	: (str)
			Filename to save an image of the mode shape.
		
	"""
	
	modeNumber = modenumber
	if model == "none":
		print("No Model_ODB specified to plot modeshapes")
		ops.wipeAnalysis()
		eigenVal = ops.eigen(modeNumber+1)
		Tn=4*asin(1.0)/(eigenVal[modeNumber-1])**0.5
		this_nodeArray, elementArray, eleClassTags = idbf._getNodesandElements()
		Mode_nodeArray = idbf._getModeShapeData(modeNumber)		# DOES NOT GIVE MODAL PERIOD
		ops.wipeAnalysis()
	else:
		print("Reading modeshape data from "+str(model)+"_ODB")
		this_nodeArray, elementArray, eleClassTags = idbf._readNodesandElements(model)
		Mode_nodeArray, Periods = idbf._readModeShapeData(model,modeNumber)
		Tn = Periods[modeNumber-1]
				

	## Check if the model is 2D or 3D
	nodeArray = np.zeros([len(this_nodeArray[:,0]), 4])
	ndm = len(this_nodeArray[0,:]) -1
	if ndm == 2:
		for ii in range(0,len(this_nodeArray[:,0])):
			nodeArray[ii,0:3] = this_nodeArray[ii,:]
	else:
		nodeArray = this_nodeArray
		
	pl = pv.Plotter()		
	pl.show(interactive_update=True)
	
	
	# DeflectedNodeCoordArray = nodeArray[:,1:]+ scale*displacement_nodeArray[jj,:,:]
	DeflectedNodeCoordArray = nodeArray[:,1:]+ scale*Mode_nodeArray[:,1:]
	
	DeflectedNodeArray = np.hstack((nodeArray[:,0].reshape(len(nodeArray[:,0]),1),DeflectedNodeCoordArray[:,0:]))
	# DeflectedNodeArray = Mode_nodeArray
	# print("nodeArray is ", nodeArray)
	# print("Mode_nodeArray ", Mode_nodeArray)
	# print("DeflectedNodeCoordArray is ", DeflectedNodeCoordArray)
	# print("DeflectedNodeArray is ", DeflectedNodeArray)
	
	# mesh_deflected, vertices, nodeTags = _get_modelDisplay(DeflectedNodeArray, elementArray, eleClassTags)
	mesh_deflected, mesh_lines_deflected, vertices, nodeTags = _get_modelDisplay(DeflectedNodeArray, elementArray, eleClassTags)
	
	point_size = 0.0
	spheres=False
		
	Clim = contourlimits
	
	if contour == "none":
		pl.add_mesh(mesh_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=1, color="green", opacity=1.0,name="thisMesh1")  # "show_edges=True" will NOT work with PlotterITK
		pl.add_mesh(mesh_lines_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=line_width, color="green", opacity=1.0,name="thisMesh2")  # "show_edges=True" will NOT work with PlotterITK
	elif contour in ["x","X"]:
		x = Mode_nodeArray[:,1]
		pl.add_mesh(mesh_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=1, scalars=x, clim=Clim, opacity=1.0,name="thisMesh1")  # "show_edges=True" will NOT work with PlotterITK
		pl.add_mesh(mesh_lines_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=line_width, scalars=x, clim=Clim, opacity=1.0,name="thisMesh2")  # "show_edges=True" will NOT work with PlotterITK
	elif contour in ["y","Y"]:
		y = Mode_nodeArray[:,2]
		pl.add_mesh(mesh_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=1, scalars=y, clim=Clim, opacity=1.0,name="thisMesh1")  # "show_edges=True" will NOT work with PlotterITK
		pl.add_mesh(mesh_lines_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=line_width, scalars=y, clim=Clim, opacity=1.0,name="thisMesh2")  # "show_edges=True" will NOT work with PlotterITK
	elif contour in ["z","Z"]:
		z = Mode_nodeArray[:,3]
		pl.add_mesh(mesh_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=1, scalars=z, clim=Clim, opacity=1.0,name="thisMesh1")  # "show_edges=True" will NOT work with PlotterITK
		pl.add_mesh(mesh_lines_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=line_width, scalars=z, clim=Clim, opacity=1.0,name="thisMesh2")  # "show_edges=True" will NOT work with PlotterITK
			
		# ALTERNATIVELY, USE https://docs.pyvista.org/api/plotting/_autosummary/pyvista.Plotter.view_xy.html

	if setview in ["XY","xy"]:
		pl.view_xy()
	elif setview in ["YZ","yz"]:
		pl.view_yz()
	elif setview in ["XZ","xz"]:
		pl.view_xz()
	elif setview == "3D":
		pl.view_isometric()
	else:
		pl.set_viewup(setview)
	
	## Overwrite setview option to "xy" if the model is 2D
	if ndm == 2:
		pl.view_xy()
	else:
		pass

	pl.add_text("Mode = "+str(modenumber), color="black", font_size=12)
	pl.add_text('VFO - Visualization for OpenSees', position='lower_left', color='green', font_size=6)
		
	if filename is not None:
		pl.screenshot(filename+".png") 
	
	point_size = 1
	#font_size = 10
    
	pl.show()


def animate_deformedshape(model="none",loadcase="none",scale=10,speedup=1,overlap="yes",contour="none",setview="3D",line_width=1,node_for_th=None, node_dof=1, moviename=None,gifname=None):
    
	"""
	a function to plot mode shape of model.
	
	Command:
	
	animate_deformedshape(model="none",loadcase="none",scale=10,speedup=1,overlap="yes",contour="none",setview="3D",line_width=1,node_for_th=None, node_dof=1, moviename=None,gifname=None)
    	
	INPUT:
		model		: (string)
			Name of the model output database as used in createODB() function. Default is "none".
			
		loadcase	: (string)
			Name of the load case to save output data.
			
		scale		: (int)
			Scale factor to be applied to the deformed shape. Default is 10.
		
		contour		: (str)
			Contours of displacement in x, y, or z. Default is "none".
			
		contourlimits: (list)
			A list of minimum and maximum limits of the displacement contours.
		
		setview		: (str)
			sets the camera view to predefined angles. Valid entries are "xy","yz","xz","3D", or a list with [x,y,z] unit vector. Default is "3D".
			
		line_width  : (int)
			Line thickness for the beam-column elements. Default is 1.
			
		node_for_th : (int)
			Node ID to display displacement time-history.
		
		node_dof	: (int)
			Degree-of-freedom to display displacement time-history of node_for_th.
			
		moviename	: (str)
			Filename to save animation as a movie in .mp4 format.
			
		gifname	: (str)
			Filename to save animation as a movie in .gif format.
		
	"""
	print("Reading structure data from "+model+"_ODB.")
	try:
		this_nodeArray, elementArray, eleClassTags = idbf._readNodesandElements(model) ## Add eleClassTags to the Tcl data
		# print(np.shape(this_nodeArray))
	except:
		raise Exception(">>>> VFO ERROR: No Model_ODB found. No active model found. Exiting 'plot_modeshape()' now.<<<<")
        
	print("Reading loadcase '"+loadcase+"' data from "+model+"_ODB.")
	try:
		timeSteps, this_displacement_nodeArray, singleStep = idbf._readNodeDispData(model,loadcase) ## Add eleClassTags to the Tcl data
		# print('timeSteps')
		# print(np.shape(displacement_nodeArray[0,:,:]))
	except:
		raise Exception(">>>> VFO ERROR: No node displacements from loadcase was found in "+model+"_ODB. Exiting 'plot_modeshape()' now. <<<<")
		
	
	## Check if the model is 2D or 3D
	nodeArray = np.zeros([len(this_nodeArray[:,0]), 4])
	displacement_nodeArray = np.zeros([np.shape(this_displacement_nodeArray)[0],np.shape(this_displacement_nodeArray)[1], 3])
	ndm = len(this_nodeArray[0,:]) -1
	if ndm == 2:
		for ii in range(0,len(this_nodeArray[:,0])):
			nodeArray[ii,0:3] = this_nodeArray[ii,:]
			
		for ii in range(0,np.shape(this_displacement_nodeArray)[0]):
			for jj in range(0,np.shape(this_displacement_nodeArray)[1]):
				# print("this_displacement_nodeArray ", this_displacement_nodeArray[ii,jj,:])
				displacement_nodeArray[ii,jj,0:2] = this_displacement_nodeArray[ii,jj,:]
	else:
		nodeArray = this_nodeArray
		displacement_nodeArray = this_displacement_nodeArray
		
		

	pl = pv.Plotter(window_size=[960, 528])
	pl.show(interactive_update=True)
	# pl.show(auto_close=False)
	
	
	#### Additional chart to overlay  ######S
	if node_for_th is not None:
		
		i, = np.where(nodeArray[:,0] == float(node_for_th))
		th_nodeDisp = displacement_nodeArray[:,int(i),node_dof]
	
		th_chart = pv.Chart2D(size=(0.30, 0.30), loc=(0.06, 0.1), x_label="Time (s)", y_label="")
		th_chart.border_width=0
		th_chart.border_color = 'w'
		th_line = th_chart.line(timeSteps,th_nodeDisp)
		th_chart.x_range = (0, _roundNearestFive(max(timeSteps)))
		th_chart.y_range = (_roundNearestFive(min(th_nodeDisp)), _roundNearestFive(max(th_nodeDisp)))
		th_chart.x_axis.tick_count = 5
		th_chart.y_axis.tick_count = 5
		th_chart.title = 'Node '+str(node_for_th)+' disp in dof '+str(node_dof)
		# th_chart.y_range = (-10, 15)
		# print("chart limits= ", np.floor(min(th_nodeDisp)), np.ceil(max(th_nodeDisp)), max(timeSteps))
		th_chart.background_color = (1.0, 1.0, 1.0, 0.4)
		pl.add_chart(th_chart)
	################################
	
	
	
	
	if moviename is not None:
		pl.open_movie(moviename+".mp4")
		# pl.write_frame()
		
	if gifname is not None:
		pl.open_gif(gifname+".gif")
	
	# mesh_original, vertices, nodeTags = _get_modelDisplay(nodeArray, elementArray, eleClassTags)
	mesh_original, mesh_lines_original, vertices, nodeTags = _get_modelDisplay(nodeArray, elementArray, eleClassTags)

	
	point_size = 0.0
	spheres=False
	
	# Add original mesh to set the view first
	# thisMesh1 = pl.add_mesh(mesh_original, show_edges=True, color="gray", opacity=0.3, name="thisMesh")  # "show_edges=True" will NOT work with PlotterITK
	# thisMesh2 = pl.add_mesh(mesh_original, show_edges=True, color="gray", opacity=0.3, name="thisMesh")  # "show_edges=True" will NOT work with PlotterITK
	
	thisMesh1 = pl.add_mesh(mesh_original, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, color="gray", opacity=0.3, render_lines_as_tubes=False, line_width=1)  # "show_edges=True" will NOT work with PlotterITK
	thisMesh2 = pl.add_mesh(mesh_lines_original, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, color="gray", opacity=0.3, render_lines_as_tubes=False, line_width=line_width)  # "show_edges=True" will NOT work with PlotterITK
	
	
	text1 = pl.add_text("")
	
	# ALTERNATIVELY, USE https://docs.pyvista.org/api/plotting/_autosummary/pyvista.Plotter.view_xy.html
	
	if setview in ["XY","xy"]:
		pl.view_xy()
	elif setview in ["YZ","yz"]:
		pl.view_yz()
	elif setview in ["XZ","xz"]:
		pl.view_xz()
	elif setview == "3D":
		pl.view_isometric()
		pos=pl.camera.position
		pl.set_position(pos)
	else:
		pl.set_viewup(setview)
		
	## Overwrite setview option to "xy" if the model is 2D
	if ndm == 2:
		pl.view_xy()
	else:
		pass
		
	point_size = 1
		
	maxDispArray = idbf._getMaxNodeDisp(timeSteps, displacement_nodeArray)
	# print(maxDispArray)
	# index_x_max = np.argmax(maxDispArray[:,1])
	# index_x_min = np.argmax(abs(maxDispArray[:,1]))
	# index_y_max = np.argmax(maxDispArray[:,3])
	
	if contour in ["x","X"]:
		_thisScalar=0
		Clim = [0,max(max(maxDispArray[:,1]),max(abs(maxDispArray[:,2])))]
	elif contour in ["y","Y"]:
		_thisScalar=1
		Clim = [0,max(max(maxDispArray[:,3]),max(abs(maxDispArray[:,4])))]
	elif contour in ["z","Z"]:
		_thisScalar=2
		Clim = [0,max(max(maxDispArray[:,5]),max(abs(maxDispArray[:,6])))]
		
	
	jj=0
	
	# Animation
	print("Creating animation gif/movie file. On-screen updating may seem choppy/slow. Please don't close the window.")
	
	for jj in range(0,len(timeSteps),math.ceil(speedup)):
		tstep = timeSteps[jj]
		DeflectedNodeCoordArray = nodeArray[:,1:]+ scale*displacement_nodeArray[jj,:,:]
		DeflectedNodeArray = np.hstack((nodeArray[:,0].reshape(len(nodeArray[:,0]),1),DeflectedNodeCoordArray))
		# print(jj)
		
		mesh_deflected, mesh_lines_deflected, vertices, nodeTags = _get_modelDisplay(DeflectedNodeArray, elementArray, eleClassTags)
		
		if contour == "none":
			Scalar = None
			Clim = None
		else:
			Scalar = np.absolute(displacement_nodeArray[jj,:,_thisScalar])
			
		# pl.clear()
		pl.remove_actor(thisMesh1)
		pl.remove_actor(thisMesh2)
		pl.remove_actor(text1)
		# thisMesh = pl.add_mesh(mesh_deflected, show_edges=True, color="green", opacity=1.0, scalars=Scalar, clim=Clim, show_scalar_bar=True)  # "show_edges=True" will NOT work with PlotterITK		
		thisMesh1=pl.add_mesh(mesh_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=1, scalars=Scalar, clim=Clim, opacity=1.0,show_scalar_bar=True,name="thisMesh1")  # "show_edges=True" will NOT work with PlotterITK
		thisMesh2=pl.add_mesh(mesh_lines_deflected, show_edges=True, render_points_as_spheres=spheres, point_size=point_size, line_width=line_width, scalars=Scalar, clim=Clim, opacity=1.0,name="thisMesh2")  # "show_edges=True" will NOT work with PlotterITK
	
	
		text1 = pl.add_text("Time = "+str(tstep), color="black", font_size=12)
		pl.add_text('VFO - Visualization for OpenSees', position='lower_left', color='green', font_size=6)
		
		try:
			th_line.update(timeSteps[0:jj,], displacement_nodeArray[0:jj,50,1])
		except:
			pass
		
		if moviename or gifname is not None:
			pl.write_frame()

		pl.update()
			
	print("Movie/Gif file is saved.")
	


def saveFiberData2D(ModelName, LoadCaseName, eleNumber, sectionNumber = 1, deltaT = 0.0, ZLE = False):
    """
    Model : string
        The name of the input model database.    
    LoadCase : string
        The name of the input loadcase.    
    element : int
        The input element to be recorded
    section : int
        The section in the input element to be recorded.
    deltaT : float, optional
        The time step to be plotted. The program will find the closed time 
        step to the input value. The default is -1.    
    """
    
    #TODO Allow for inputing more than one element/section?
    
	# Consider making these optional arguements
    FibreName = "FiberData"
    ftype = '.out'
    
    ODBdir = ModelName+"_ODB"		# ODB Dir name
    FibreFileName = FibreName  + '_ele_' + str(eleNumber) + '_section_' + str(sectionNumber) + ftype
    FiberDir = os.path.join(ODBdir, LoadCaseName, FibreFileName)
	
    if ZLE == True:
        ops.recorder('Element' , '-file', FiberDir, '-time', '-dT', deltaT, '-ele', eleNumber, 'section', 'fiberData')
    else:
        ops.recorder('Element' , '-file', FiberDir, '-time', '-dT', deltaT, '-ele', eleNumber, 'section', str(sectionNumber), 'fiberData')
        

### All the plotting related definitions start here.

ele_style = {'color':'black', 'linewidth':1, 'linestyle':'-'} # elements
node_style = {'color':'black', 'marker':'o', 'facecolor':'black','linewidth':0.}
node_style_animation = {'color':'black', 'marker':'o','markersize':2., 'linewidth':0.} 

node_text_style = {'fontsize':8, 'fontweight':'regular', 'color':'green'} 
ele_text_style = {'fontsize':8, 'fontweight':'bold', 'color':'darkred'} 

WireEle_style = {'color':'black', 'linewidth':1, 'linestyle':':'} # elements
Eig_style = {'color':'red', 'linewidth':1, 'linestyle':'-'} # elements
	













def plot_fiberResponse2D(Model, LoadCase, element, section, LocalAxis = 'y', InputType = 'stress', tstep = -1):
    """
    

    Parameters
    ----------
    Model : string
        The name of the input model database.    
    LoadCase : string
        The name of the input loadcase.    
    element : int
        The input element to be plotted
    section : TYPE
        The section in the input element to be plotted.
    LocalAxis : TYPE, optional
        The local axis to be plotted on the figures x axis. 
        The default is 'y', 'z' is also possible.
    InputType : TYPE, optional
        The quantity to be plotted. The default is 'stress', 'strain' is 
        also possible
    tstep : TYPE, optional
        The time step to be plotted. The program will find the closed time 
        step to the input value. The default is -1.

    """
    
    
    
    # Catch invalid input types
    if InputType not in ['stress', 'strain']:
        raise Exception('Invalid input type. Valid Entries are "stress" and "strain"')
    
    # Catch invalid Direction types
    if LocalAxis not in ['z', 'y']:
        raise Exception('Invalid LocalAxis type. Valid Entries are "z" and "y"')
        
    # get the local axis and response labes and indexes
    axisIndex, axisXlabel = ipltf._getAxisInfo(LocalAxis)
    responseIndex, axisYlabel = ipltf._getResponseInfo(InputType)
    
    
    timeSteps, fiberData  = idbf._readFiberData2D(Model, LoadCase, element, section)
    
    # find the appropriate time step
    if tstep == -1:
        LoadStep = -1
        printLine = "Final deformed shape"
    else:
        LoadStep = (np.abs(timeSteps - tstep)).argmin()			# index closest to the time step requested.
        if timeSteps[-1] < tstep:
            print("XX Warining: Time-Step has exceeded maximum analysis time step XX")
        printLine = 'Fibre '+  InputType + " at time: " + str(round(timeSteps[LoadStep], 2))
            

    fiberYPosition = fiberData[LoadStep,axisIndex::5]
    fiberResponse  = fiberData[LoadStep, responseIndex::5]
    
    # Sort indexes so they appear in an appropraiate location
    sortedIndexes = np.argsort(fiberYPosition)
    fibrePositionSorted = fiberYPosition[sortedIndexes]
    fibreResponseSorted = fiberResponse[sortedIndexes]
    
    
    fig, ax = plt.subplots()
    Xline = ax.plot([fibrePositionSorted[0],fibrePositionSorted[-1]],[0, 0], c ='black', linewidth = 0.5)
    line = ax.plot(fibrePositionSorted, fibreResponseSorted)
    
    xyinput = np.array([fibrePositionSorted,fibreResponseSorted]).T
    ipltf._setStandardViewport(fig, ax, xyinput, 2)
    
    ax.set_ylabel(axisYlabel)  
    ax.set_xlabel(axisXlabel)    
        
    plt.show()
    return fig, ax
    

def animate_fiberResponse2D(Model, LoadCase, element, section, LocalAxis = 'y', InputType = 'stress', skipStart = 0, 
                            skipEnd = 0, rFactor=1, outputFrames=0, fps = 24, Xbound = [], Ybound = []):
    """
    Parameters
    ----------
    Model : string
        The name of the input model database.    
    LoadCase : string
        The name of the input loadcase.    
    element : int
        The input element to be plotted
    section : TYPE
        The section in the input element to be plotted.
    LocalAxis : string, optional
        The local axis to be plotted on the figures x axis. 
        The default is 'y', 'z' is also possible.
    InputType : string, optional
        The quantity 
    skipStart : int, optional
        If specified, this many datapoints will be skipped from the analysis
        data set, before reductions.
        The default is 0, or no reduction
    skipEnd : int, optional
        If specified, this many frames will be skipped at the end of the 
        analysis dataset, before reduction. The default is 0, or no reduction.
    rFactor : int, optional
        If specified, only every "x" frames will be plotted. e.g. x = 2, every 
        other frame is shown.
        The default is 1.
    outputFrames : int, optional
        The number of frames to be included after all other reductions. If the
        reduced number of frames is less than this value, no change is made.
        The default is 0.
    fps : int, optional
        Number of animation frames to be displayed per second. The default is 24.
    Xbound : [xmin, xmax], optional
        The domain of the chart. The default is 1.1 the max and min values.
    Ybound : [ymin, ymax], optional
        The range of the chart. The default is 1.1 the max and min values.

    
    """
    
    # Catch invalid input types
    if InputType not in ['stress', 'strain']:
        raise Exception('Invalid input type. Valid Entries are "stress" and "strain"')
    
    # Catch invalid Direction types
    if LocalAxis not in ['z', 'y']:
        raise Exception('Invalid LocalAxis type. Valid Entries are "z" and "y"')
    
    # get the local axis and response labes and indexes
    axisIndex, axisXlabel = ipltf._getAxisInfo(LocalAxis)
    responseIndex, axisYlabel = ipltf._getResponseInfo(InputType)
    
    timeSteps, fiberData  = idbf._readFiberData2D(Model, LoadCase, element, section)
                
    # Get the desired iformation.
    fiberYPosition = fiberData[:,axisIndex::5]
    fiberResponse  = fiberData[:, responseIndex::5]
    
    # Sort indexes so they appear in an appropraiate location.
    sortedIndexes = np.argsort(fiberYPosition[0,:])
    fibrePositionSorted = fiberYPosition[:,sortedIndexes]
    fibreResponseSorted = fiberResponse[:,sortedIndexes]    
    
    # Get the bounds on the final output data.
    xmin, xmax, ymin, ymax = ipltf._getFiberBounds(fibrePositionSorted, fibreResponseSorted, Xbound, Ybound)
    
    xinputs, yinputs = ipltf._skipFiberData(fibrePositionSorted, fibreResponseSorted, fiberYPosition,  skipStart, skipEnd)

    # Reduce the data if the user specifies
    if rFactor != 1:
        xinputs = xinputs[::rFactor, :]
        yinputs = yinputs[::rFactor, :]    
        timeSteps = timeSteps[::rFactor]
    
    # If the Frames isn't specified, use the length of the reduced vector.
    if outputFrames == 0:
        outputFrames = len(xinputs[:, 0])
    else:
        outputFrames = min(outputFrames,len(xinputs[:, 0]))
    
    # Get the final output frames. X doesn't change
    xinputs = xinputs[:outputFrames, :]
    yinputs = yinputs[:outputFrames, :]    
    xinput = xinputs[0,:]
    
    # Initialize the plot
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=.15) # Add extra space bellow graph
    
    line, = ax.plot(xinput, yinputs[0,:])
    Xline = ax.plot([fibrePositionSorted[0,0],fibrePositionSorted[0,-1]], [0, 0], c ='black', linewidth = 0.5)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
        
    ax.set_ylabel(axisYlabel)  
    ax.set_xlabel(axisXlabel)    

    Frames = np.arange(0, outputFrames)
    FrameStart = int(Frames[0])
    FrameEnd = int(Frames[-1])
    
    # Slider Location and size relative to plot
    # [x, y, xsize, ysize]
    axSlider = plt.axes([0.25, .03, 0.50, 0.02])
    plotSlider = Slider(axSlider, 'Time', timeSteps[FrameStart], timeSteps[FrameEnd], valinit=timeSteps[FrameStart])

    # Animation controls
    global is_paused
    is_paused = False # True if user has taken control of the animation   
    
    def on_click(event):
        # Check where the click happened
        (xm,ym),(xM,yM) = plotSlider.label.clipbox.get_points()
        if xm < event.x < xM and ym < event.y < yM:
            # Event happened within the slider, ignore since it is handled in update_slider
            return
        else:
            # Toggle on/off based on click
            global is_paused
            if is_paused == True:
                is_paused=False
            elif is_paused == False:
                is_paused=True       
    
    # Define the update function
    def update_line_slider(Time):
        global is_paused
        is_paused=True

        TimeStep = (np.abs(timeSteps - Time)).argmin()
        # Get the current data        
        y = yinputs[TimeStep,:]
        
        # Update the background line
        line.set_data(xinput, y)
        fig.canvas.draw_idle()    
        
        return line,
    
    
    def update_plot(ii):
    
        # If the control is manual, we don't change the plot    
        global is_paused
        if is_paused:
            return line,
       
        # Find the close timeStep and plot that
        # CurrentFrame = int(np.floor(plotSlider.val))

        # Find the close timeStep and plot that
        CurrentTime = plotSlider.val
        CurrentFrame = (np.abs(timeSteps - CurrentTime)).argmin()

        CurrentFrame += 1
        if CurrentFrame >= FrameEnd:
            CurrentFrame = FrameStart
        
        # Update the slider
        plotSlider.set_val(timeSteps[CurrentFrame])        
        
        # Update the slider
        is_paused = False # the above line called update_slider, so we need to reset this
        return line,  
    
    
    plotSlider.on_changed(update_line_slider)
    
    # assign click control
    fig.canvas.mpl_connect('button_press_event', on_click)    
    
    interval = 1000/fps
    
    line_ani = animation.FuncAnimation(fig, update_plot, outputFrames, 
                                       # fargs=(xinput, yinputs, line), 
                                       interval=interval)
									   
    plt.show()
    return line_ani



def plot_eleHysteresis(Model="none", LoadCase="none", element=[], monitorGroupName="none", dof=1, limitStates=[]):
	
	"""
	plots the output of specific elements.
	
	Model: str
		Name of the Model as defined in createODB() command.
		
	LoadCase: str
		Name of the LoadCase as defined in createODB() command.
		
	element: list
		List of elements whose output is to be displayed. Should be a part of monitorGroupName group.
		
	monitorGroupName: str
		Name of the monitorGroupName as defined in createODB() command.
		
	dof: int
		Degree-of-freedom to be displayed.
		
	limitStates: list (Optional)
		Limit states to be displayed.
	
	"""

	if Model == "none" or LoadCase=="none":
		print("Command should be plot_eleHysteresis(Model='modelname',loadCase='loadcase',element=int(eleTag), monitorGroupName='monitorGroupName', dof=int(dof), <limitStates=[]>")
		raise Exception("No output database specified in plot_eleHysteresis() command. Exiting now.")
	
	if element==[]:
		print("Command should be plot_eleHysteresis(Model='modelname',loadCase='loadcase',element=int(eleTag), monitorGroupName='monitorGroupName', dof=int(dof), <limitStates=[]>")
		raise Exception('No element tags provided to plot the data in plot_eleHysteresis() command. Exiting now.')
		# print("Plotting structure hysteresis from load case ")	
		
	if monitorGroupName=="none":
		print("Command should be plot_eleHysteresis(Model='modelname',loadCase='loadcase',element=int(eleTag), monitorGroupName='monitorGroupName', dof=int(dof), <limitStates=[]>")
		raise Exception('No monitorGroupName provided in plot_eleHysteresis() command. Exiting now.')
		
	else:
		print("Reading data from "+str(Model)+"_ODB/"+LoadCase)
		nodeArray, elementArray, eleClassTags = idbf._readNodesandElements(Model)
		timeSteps, Disp_nodeArray = idbf._readNodeDispData(Model,LoadCase)
		
	
		
	#### Read the monitoring element deformation data
	
	if monitorGroupName !="none":
		monitorOutput = "deformations"
		ODBdir = Model+"_ODB"
		LoadCaseDir = os.path.join(ODBdir, LoadCase)
		GroupMonitorDir = os.path.join(LoadCaseDir,monitorGroupName)
		if not os.path.exists(GroupMonitorDir):
			raise Exception("No Group monitor data saved for "+str(monitorGroupName)+". Check inputs. Or, use createODB command to save group data.")
			
		MonitorEleDef, MonitorEleForce, MonitorEleTags, MonitorEleInfo = idbf._readMonitorElementData(monitorOutput,GroupMonitorDir)
	else:
		pass
		
		
	### Initiate a matrix to save specific element data
	Nele = len(element)
	eleDeformation = np.array([Nele, len(MonitorEleDef)])
	eleForce = np.array([Nele, len(MonitorEleDef)])
	NNele = 1
	displayEleTags = []
	
	for ele in element:
		if ele not in MonitorEleTags:
			print('element '+str(ele)+' was not found in the group '+str(monitorGroupName)+'. Use createODB() command to add this element to the group')
		else:
			displayEleTags.append(ele)
			NNele=+1
			
	### Initialize the figures

	yLabelText = "DOF "+str(dof)+" Force"
	xLabelText = "DOF "+str(dof)+" Def"

	for ele in displayEleTags:
		fig, ax = plt.subplots()
		eleDeformation, eleForce = idbf._elementHysteresis(ele, dof, MonitorEleInfo, MonitorEleTags, MonitorEleDef, MonitorEleForce)
		ax.plot(eleDeformation,eleForce, color ='black', linewidth = 1.0)
		ax.set_ylabel(yLabelText)  
		ax.set_xlabel(xLabelText)    
		ax.set_title("Element "+str(ele)+" Hysteresis - DOF "+str(dof)+", LC:"+str(LoadCase))  
		ax.grid(linestyle='dotted')		
		ax.axhline(y=0, color='k', linewidth = 0.5)
		ax.axvline(x=0, color='k', linewidth = 0.5)		
		
	plt.show()
	return fig, ax
