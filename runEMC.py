#!/usr/bin/env python
"""
Usage instructions


"""

import os, sys, re
import glob as G
import numpy as N
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as PL


class runEMC(object):
	def __init__(self):
		if os.path.isfile("contrast.dat"):
			tmp = N.fromfile("contrast.dat", sep=" ")
			self.l = int(tmp[0])
			self.numImgs = int(tmp[1])
			self.contrast = (tmp[2:]).reshape(self.l, self.l)
		else:
			print("Warning: contrast.dat is missing. You'll have to make one.")
			self.l = 50
			self.contrast = N.random.rand(self.l, self.l)

	def compile(self):	
		print("Compiling make_data.c")
		try:
			os.system("gcc -O3 make_data.c -lm -o make_data")
		except:
			print("Error: make_data.c was not be compiled.")
		print("Compiling EMC.c")
		try:
			os.system("gcc -O3 EMC.c -lm -o EMC")
		except:
			print("Error: EMC.c could not be compiled.")

	def removeBinaryObjects(self):
		print("Removing binary: make_data")
		try:
			os.remove("make_data")
		except:
			print("Error: could not remove binary, make_data.")
		print("Removing binary: EMC")
		try:
			os.remove("EMC")
		except:
			print("Could not remove binary, EMC.")

	def makeData(self, numData, sigNoise, bgNoise, hitRate):
		"""
		The number of conformations (-c flag) always defaults to here. 
		Reconstructing with multiple conformations is possible, but breaks the examples in runEMC.py.  
		"""
		cmd = "./make_data -d %d -s %lf -n %lf -h %lf -c 1"%(numData, sigNoise, bgNoise, hitRate)	
		print("Making data with the following command:\n%s"%cmd)
		os.system(cmd)
		print("data.dat created.")
		
	def EMC(self, numIter):
		cmd = "./EMC -n %d"%(numIter)	
		print("Running EMC with the following command:\n%s"%cmd)
		os.system(cmd)
		print("Iterations completed.")
		
	def viewResults(self):
		

	def viewData(self):
		pass

	def runExample1(self):
		pass

	def runExample2(self):
		pass

	def archiveResult(self):
		#Prompt if ASCII output should be deleted
		pass

	def loadArchivedResults(self):
		#Prompt that local ASCII files will be overwritten
		pass

	def computeStats(self):
		pass

	def createNewImage(self):
		try:
			import PIL
		except:
			print("Unable to import Python Imaging Library. Creation of new Image aborted.")



