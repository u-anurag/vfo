

####################################################
## 
## Test Examples
##
####################################################

import numpy as np
import openseespy.opensees as ops
import vfo.vfo as vfo

import math

def portalframe2d():

	print("portal frame 2D test")

	ops.wipe()

	# set modelbuilder
	ops.model('basic', '-ndm', 2, '-ndf', 3)

	############################################
	### Units and Constants  ###################
	############################################

	inch = 1;
	kip = 1;
	sec = 1;

	# Dependent units
	sq_in = inch*inch;
	ksi = kip/sq_in;
	ft = 12*inch;

	# Constants
	g = 386.2*inch/(sec*sec);
	pi = math.acos(-1);

	#######################################
	##### Dimensions 
	#######################################

	# Dimensions Input
	H_story=10.0*ft;
	W_bayX=16.0*ft;
	W_bayY_ab=5.0*ft+10.0*inch;
	W_bayY_bc=8.0*ft+4.0*inch;
	W_bayY_cd=5.0*ft+10.0*inch;

	# Calculated dimensions
	W_structure=W_bayY_ab+W_bayY_bc+W_bayY_cd;

	# ###########################
	# ##### Nodes
	# ###########################

	# Create All main nodes
	ops.node(1, 0.0, 0.0)
	ops.node(2, W_bayX, 0.0)

	ops.node(11, 0.0, H_story)
	ops.node(12, W_bayX, H_story)

	# ###############
	#  Constraints
	# ###############

	ops.fix(1, 1, 1, 1)
	ops.fix(2, 1, 1, 1)

	# #######################
	# ### Elements 
	# #######################

	ColTransfTag=1
	BeamTranfTag=2

	ops.geomTransf('Linear', ColTransfTag)
	ops.geomTransf('Linear', BeamTranfTag)

	# Assign Elements  ##############
	# ## Add non-linear column elements
	ops.element('elasticBeamColumn', 1, 1, 11, 20., 1000., 1000., ColTransfTag, '-mass', 0.0)
	ops.element('elasticBeamColumn', 2, 2, 12, 20., 1000., 1000., ColTransfTag, '-mass', 0.0)

	#  ### Add linear main beam elements, along x-axis
	ops.element('elasticBeamColumn', 101, 11, 12, 20., 1000., 1000., BeamTranfTag, '-mass', 0.0)

	# Visualize the model
	return  vfo.plot_model()
	

def portalframe3d():
	print("portal frame 3D test")

	ops.wipe()

	# set modelbuilder
	ops.model('basic', '-ndm', 3, '-ndf', 6)

	############################################
	### Units and Constants  ###################
	############################################

	inch = 1;
	kip = 1;
	sec = 1;

	# Dependent units
	sq_in = inch*inch;
	ksi = kip/sq_in;
	ft = 12*inch;

	# Constants
	g = 386.2*inch/(sec*sec);
	pi = math.acos(-1);

	#######################################
	##### Dimensions 
	#######################################

	# Dimensions Input
	H_story=10.0*ft;
	W_bayX=16.0*ft;
	W_bayY_ab=5.0*ft+10.0*inch;
	W_bayY_bc=8.0*ft+4.0*inch;
	W_bayY_cd=5.0*ft+10.0*inch;

	# Calculated dimensions
	W_structure=W_bayY_ab+W_bayY_bc+W_bayY_cd;

	# ###########################
	# ##### Nodes
	# ###########################

	# Create All main nodes
	ops.node(1, 0.0, 0.0, 0.0)
	ops.node(2, W_bayX, 0.0, 0.0)

	ops.node(11, 0.0, H_story, 0.0)
	ops.node(12, W_bayX, H_story, 0.0)

	# ###############
	#  Constraints
	# ###############

	ops.fix(1, 1, 1, 1, 1, 1, 1)
	ops.fix(2, 1, 1, 1, 1, 1, 1)

	# #######################
	# ### Elements 
	# #######################

	ColTransfTag=1
	BeamTranfTag=2

	ops.geomTransf('Linear', ColTransfTag, 1, 0, 0)
	ops.geomTransf('Linear', BeamTranfTag, 0, 0, 1)

	# Assign Elements  ##############
	# ## Add non-linear column elements
	ops.element('elasticBeamColumn', 1, 1, 11, 20., 1000., 1000., 1000., 1000., 1000., ColTransfTag, '-mass', 0.0)
	ops.element('elasticBeamColumn', 2, 2, 12, 20., 1000., 1000., 1000., 1000., 1000.,ColTransfTag, '-mass', 0.0)

	#  ### Add linear main beam elements, along x-axis
	ops.element('elasticBeamColumn', 101, 11, 12, 20., 1000., 1000., 1000., 1000., 1000.,BeamTranfTag, '-mass', 0.0)

	# Visualize the model
	return vfo.plot_model()
	

def tri3d():
	print("portal frame 3D test")

	ops.wipe()

	# set modelbuilder
	ops.model('basic', '-ndm', 3, '-ndf', 3)

	# Create All main nodes
	ops.node(1, 0.0, 0.0, 0.0)
	ops.node(2, 10.0, 0.0, 0.0)
	ops.node(3, 0.0, 10.0, 0.0)

	# ###############
	#  Constraints
	# ###############

	ops.fix(1, 1, 1, 1)
	ops.fix(2, 1, 1, 1)

	ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25, 6.75) 

	# Assign Elements  ##############
	# ## Add non-linear column elements
	ops.element('Tri31', 1, 1,2,3, 0.1, 'PlaneStress', 1)
	
	# Visualize the model
	return vfo.plot_model()
	
	
def quad2d():
	print("quad3d 3D test")

	ops.wipe()

	# set modelbuilder
	ops.model('basic', '-ndm', 2, '-ndf', 2)
	
	# ###########################
	# ##### Nodes
	# ###########################

	# Create All main nodes
	ops.node(1, 0.0, 0.0)
	ops.node(2, 50.0, 0.0)
	ops.node(3, 0.0, 30.0)
	ops.node(4, 50.0, 30.0)
	
	# ops.fix(1, 1, 1, 1, 1, 1, 1)
	# ops.fix(2, 1, 1, 1, 1, 1, 1)

	ops.fix(1, 1, 1)
	ops.fix(2, 1, 1)

	ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25, 6.75) 

	# Assign Elements  ##############
	# ## Add non-linear column elements
	ops.element('quad', 1, 1,2,4,3, 0.1, 'PlaneStress', 1)
	
	# vfo.createODB("Test_mvlem3d", "none", 0)
	# Visualize the model
	return vfo.plot_model()
	
	
	
def mvlem_3d():
	print("MVLEM 3D test")

	ops.wipe()

	# set modelbuilder
	ops.model('basic', '-ndm', 3, '-ndf', 6)
	
	# ###########################
	# ##### Nodes
	# ###########################

	# Create All main nodes
	ops.node(1, 0.0, 0.0, 0.0)
	ops.node(2, 50.0, 0.0, 0.0)
	ops.node(3, 0.0, 0.0, 30.0)
	ops.node(4, 50.0, 0.0, 30.0)
	ops.node(5, 0.0, 0.0, 60.0)
	ops.node(6, 50.0, 0.0, 60.0)

	ops.fix(1, 1, 1, 1, 1, 1, 1)
	ops.fix(2, 1, 1, 1, 1, 1, 1)

	ops.uniaxialMaterial('Elastic', 1, 314705)
	# Concrete Materials
	ops.uniaxialMaterial('Concrete02', 201, -7.934, -0.0023, 0, -0.01, 0.079, 0.356292015066294, 253.858060734734)

	# Steel Materials
	ops.uniaxialMaterial('SteelMPF', 301, 68.313, 68.313, 27847.246, 0.0055, 0.0055, 20, 0.925, 0.15)

	# Assign Elements  ##############
	ops.element('MVLEM_3D', 1,1,2,4,3, 10, '-thick', *[3.937, 3.937, 3.937, 3.937, 3.937, 3.937, 3.937, 3.937, 3.937, 3.937], 
												'-width', *[5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0], 
												'-rho', *[0.0387, 0.0387, 0.00342, 0.00342, 0.00342, 0.00342, 0.00342, 0.00342, 0.0226, 0.0226], 
												'-matConcrete', *[201, 201, 201, 201, 201, 201, 201, 201, 201, 201], 
												'-matSteel', *[301, 301, 301, 301, 301, 301, 301, 301, 301, 301], '-matShear', 1)
	ops.element('MVLEM_3D', 2,3,4,6,5, 10, '-thick', *[3.937, 3.937, 3.937, 3.937, 3.937, 3.937, 3.937, 3.937, 3.937, 3.937], 
												'-width', *[5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0], 
												'-rho', *[0.0387, 0.0387, 0.00342, 0.00342, 0.00342, 0.00342, 0.00342, 0.00342, 0.0226, 0.0226], 
												'-matConcrete', *[201, 201, 201, 201, 201, 201, 201, 201, 201, 201], 
												'-matSteel', *[301, 301, 301, 301, 301, 301, 301, 301, 301, 301], '-matShear', 1)

	# vfo.createODB("Test_mvlem3d", "none", 0)
	# Visualize the model
	return vfo.plot_model()
	
def tetra():
	print("Tetrahedral test")

	ops.wipe()

	# set modelbuilder
	ops.model('basic', '-ndm', 3, '-ndf', 3)

	# Create All main nodes
	ops.node(1, 0.0, 0.0, 0.0)
	ops.node(2, 10.0, 0.0, 0.0)
	ops.node(3, 5.0, 10.0, 0.0)
	ops.node(4, 5.0, 5.0, 10.0)

	# ###############
	#  Constraints
	# ###############

	ops.fix(1, 1, 1, 1)
	ops.fix(2, 1, 1, 1)

	ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25, 6.75) 

	# Assign Elements  ##############
	# ## Add non-linear column elements
	ops.element('FourNodeTetrahedron', 1, 1,2,3,4, 1)
	
	# Visualize the model
	return vfo.plot_model()
	

fig1 = portalframe2d()
fig2 = portalframe3d()
fig3 = quad2d()
fig3 = mvlem_3d()
fog4 = tri3d()
fog5 = tetra()
