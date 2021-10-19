"""
    A code to show a live model rendering of an OpenSeesPy modeling script. Works only for 2D beam-column elements for now.
	
	Use: 
		import vfo.live_model_view as opspre
		
		opspre.live_view("ModelFileName.py")
		
    Copyright 2021 @ vfo
"""

import numpy as np
import os
import sys
import string
import matplotlib.pyplot as pl
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

def live_view(filename):

	refresh_rate = 1

	def genTempFile(filename):
		# List of words to ignore while parsing
		ignoreWords = ['openseespy.opensees','chdir','dir','makedirs','path','model(','wipe','wipeAnalysis','geomTransf(','nodeCoord(','plot_modeshape','eigen','print','recorder',
							'recorder(','timeSeries(','pattern(','eleLoad(','constraints(','numberer','system(',
							'test(','algorithm(','integrator(','analysis(','analyze','loadConst(','timeSeries(',
							'rayleigh','nodeDisp','betaKcomm','matplotlib','numpy','scipy','os','math','Pushover',
							'beamIntegration(','equalDOF(','load(','getLoadFactor(','plt.','nodeReaction(',
							'uniaxialMaterial(','section(','eleResponse']
		selectWords = ['node(','mass(','fix(','element(']
		
		ops = ""	# Initialize a string to record how the opensees is imported
		ndm = 2		# Initialize model space 2D or 3D
		iVar = 0
		with open (filename, 'r') as script, open('temp_model.py', 'w') as f_out:
			# CHECK if you need it here. You may not need it once using the procedure as a function
			f_out.write("import numpy as np \n")
				
			for line in script:
				checkLine = line.lstrip().replace(" ","").replace("\t","")
				checkLine1 = line.lstrip()
				if '\t' in line:
					line1=line.split('\t')
				else:
					line1=line.split(' ')
				
				if len(checkLine) !=0 and checkLine[0] == "#":
					pass
				elif 'openseespy.opensees' in line:
					if line.split()[3] != "*":
						ops = line.split()[3]
						line = line.replace((ops+"."),"")
						checkLine = line.lstrip().replace(" ","").replace("\t","")
					else:
						pass
					f_out.write("from openseespy.opensees import * \n")
				else:
					pass
				
				if checkLine[:6] == 'model(':
					if checkLine.split('(')[1].split(',')[2] == str(2):
						ndm = 2
						# print("2D Model")
						f_out.write("ndm = 2\n")
						f_out.write("constraints = np.zeros([1,4]) # initialize the constaints matrix \n")
						f_out.write("nodeMass = np.zeros([1,4]) # initialize the nodal mass matrix \n")
					else:
						# print("3D Model")
						ndm = 3
						f_out.write("ndm = 3\n")
						f_out.write("constraints = np.zeros([1,7]) # initialize the matrix \n")
						f_out.write("nodeMass = np.zeros([1,7]) # initialize the nodal mass matrix \n")

					f_out.write("nodes = np.zeros([1,ndm+1]) # initialize the matrix \n")
					f_out.write("elements = np.zeros([1,4]) # initialize the matrix \n")

				
				# print(line)
				
				if any(x in line for x in ignoreWords):
					pass
				elif len(checkLine) !=0 and checkLine[0] == "#":
					pass
				elif 'node(' in line:
					# newLine = line.replace("node(","nodes=np.vstack((nodes,[").replace(")","]))")
					#print(checkLine)
					nodeLine = checkLine.split("(")[1].split(")")[0].replace("\n","").split(",")
					preLine = line.split('node(')[0]
					# print(preLine)
					# print(nodeLine)
					if ndm == 2:
						newLine = str(preLine)+"nodes=np.vstack((nodes,["+str(nodeLine[0])+","+str(nodeLine[1])+","+str(nodeLine[2])+"]))\n"
					else:
						newLine = str(preLine)+"nodes=np.vstack((nodes,["+str(nodeLine[0])+","+str(nodeLine[1])+","+str(nodeLine[2])+","+str(nodeLine[3])+"]))\n"
					f_out.write(newLine)
				elif 'mass(' in line:
					newLine = line.replace("mass(","nodeMass=np.vstack((nodeMass,[").replace(")","]))")
					f_out.write(newLine)
				elif 'fix(' in line:
					newLine = line.replace("fix(","constraints=np.vstack((constraints,[").replace(")","]))")
					f_out.write(newLine)
				elif 'element(' in line:
					eleLine = checkLine.split("(")[1].split(",")
					preLine = line.split('element(')[0]
					newLine = str(preLine)+"elements=np.vstack((elements,["+str(eleLine[0])+","+str(eleLine[1])+","+str(eleLine[2])+","+str(eleLine[3])+"]))\n"
					f_out.write(newLine)
					# create a subclass for NL elements to get section tag
				else:
					f_out.write(line)

			f_out.write("np.savetxt('temp_nodes.txt',nodes,fmt='%s') \n")
			f_out.write("np.savetxt('temp_elements.txt',elements,fmt='%s') \n")
			
		f_out.close()

	def update_plot2D(frame, filename):
	
		lines=[]
	
		ele_style = {'color':'black', 'linewidth':1, 'linestyle':'-'} # elements
		node_style = {'color':'black', 'marker':'o', 'facecolor':'black'}
		node_text_style = {'fontsize':6, 'fontweight':'regular', 'color':'green'} 
		ele_text_style = {'fontsize':6, 'fontweight':'bold', 'color':'darkred'}
	
		genTempFile(filename)
		exec(open('temp_model.py').read())
		Nodes = np.loadtxt('temp_nodes.txt', delimiter=' ', unpack=True)
		Elements = np.loadtxt('temp_elements.txt', dtype=str, delimiter=' ', unpack=True)
		
		if frame > 1:
			xlim_low, xlim_high = ax.get_xlim()
			ylim_low, ylim_high = ax.get_ylim()
		else:
			pass
		
		
		# function to get node coordinates	
		def nodeCoord(node):
			i=1
			while Nodes[0,i] != node:
				i += 1
			return Nodes[1,i], Nodes[2,i]
		
		ax.clear()
		ax.scatter(Nodes[1,:],Nodes[2,:])
		
		for node in Nodes[0,1:]:
			# print(node)
			ax.text(nodeCoord(node)[0]*1.02, nodeCoord(node)[1]*1.02, str(int(node)),**node_text_style) #label nodes		
					
						
		Nele=len(Elements[0,:])
		x = []
		y = []
		for i in range(1,Nele):
			EleType=Elements[0,i]
			Element=Elements[1,i]
			iNode = nodeCoord(int(Elements[2,i]))
			jNode = nodeCoord(int(Elements[3,i]))
			
			x.append(iNode[0])  # list of x coordinates to define plot view area
			y.append(iNode[1])	# list of y coordinates to define plot view area
			
			ax.plot((iNode[0], jNode[0]), (iNode[1], jNode[1]),marker='', **ele_style)
			line = [[iNode[0], jNode[0]], [iNode[1], jNode[1]]]
			ax.text((iNode[0]+jNode[0])/2, (iNode[1]+jNode[1])/2, str(Element), **ele_text_style) #label elements
			lines.append(line)
			
		if frame > 1:
			ax.set_xlim(xlim_low, xlim_high)
			ax.set_ylim(ylim_low, ylim_high)
			# print(xlim_low, xlim_high)
		else:
			pass
		

	fig = pl.figure(figsize=(8, 8))
	ax = fig.add_subplot(1,1,1)
	ani = animation.FuncAnimation(fig, update_plot2D, interval=refresh_rate*1000, blit=False, fargs=(filename,))
	
	pl.show()
