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
		self.contrastFile = "contrast.dat"

	def compile(self)	
		print("Compiling make_data.c")
		try:
			os.system("gcc -03 make_data.c -lm -o make_data")
		except:
			print("Error: make_data.c was not be compiled.")
		print("Compiling EMC.c")
		try:
			os.system("gcc -03 EMC.c -lm -o EMC")
		except:
			print("Error: EMC.c could not be compiled.")

	def makeData(self):
		pass

	def EMC(self):
		pass

	def viewResults(self):
		pass

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



