#!/usr/bin/env python
#coding: utf-8

import os
import os.path
import sys
import numpy as np
from numpy import *
import mmap
import glob


baseDir = os.getcwd()

filename = str(baseDir)+'/ns-3d_RefBl3/amrs/Channel_3d_d/'

dim = 3
cellsPerBlock = 32
#SizeOneBlock

contFilesamrs = 0

for amrFile in glob.glob(os.path.join(filename, '*.amr')):
	#print amrFile
	path_amr = os.path.join(filename, amrFile)

	f_in = open(path_amr, 'r') #opens a file in the reading mode

	in_lines = f_in.readlines()           #reads it line by line
	out=[]
	contLines = 0

	for line in in_lines:
	    list_values = np.array(line.split()) #separate elements by the spaces, returning a list with the numbers as array
	    #points = list_values.astype(np.float) #converts them to floats

	    if (contLines==0):
	    	points = list_values.astype(np.float) #converts them to floats
	        sizeX = points[1]-points[0]
	        sizeY = points[3]-points[2]

	        sizeOneBlock = min(sizeX,sizeY)

	        if (dim == 3):
	        	sizeZ = points[5]-points[4]
	        	sizeOneBlock = min(sizeOneBlock,sizeZ)

	    if (contLines == 1):
	    	points = list_values.astype(np.int) #converts them to floats

	    if (contLines == 2):
	    	#points = list_values #converts them to floats
	    	#points = list_values.astype(np.float) #converts them to floats
	    	points = [None, None, None, None]
	    	print points
	    	ds = float(sizeOneBlock/cellsPerBlock)
	       	points[0] = float(ds)
	       	points[0] = float(ds)
	       	points[1] = float(ds)
	       	if (dim == 2):
	       		points[2] = '1'

	       	if (dim == 3):
	       		points[2] = float(ds)
	       		points[3] = '1'
	       	print points
	       	#points.compressed()


	    if (contLines == 3):
	    	points = list_values.astype(np.int) #converts them to floats
	        points[3] = int(sizeX/ds)
	        points[4] = int(sizeY/ds)
	        if (dim == 3):
	        	points[5] = int(sizeZ/ds)
	        #points = list_values.astype(np.int) #converts them to floats

	    out.append(points)         #stores the numbers in a list, where each list corresponds to a lines' content
	    contLines = contLines+1
	f_in.close()	                        #closes the file
	contFilesamrs = contFilesamrs + 1

	#------ Writing refined domain
	f_out=open(amrFile, "w")     #opens amrs files in the writing mode
	for cur_list in out:
	    for i in cur_list:
	        f_out.write(str(i)+" ")    #writes each number, plus a space
	    f_out.write("\n")               #writes a newline
	f_out.close()                       #closes the file

#	f = open(amrFile, 'r+')
#	f.read()
#	for line in f:
#		line = line.replace("1.0", "1")
#		f.write(line + '\n')
#	f.close()
