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
		
	def viewLog(self):
		f = open("EMC.log", 'r')
		lines = f.readlines()
		f.close()
		m = re.compile('iteration \= (?P<it>\d+)\s+error \= (?P<err>[0-9.]+)\s+m_info \= (?P<mi>[0-9.]+)\s+iter_time \= (?P<t>[0-9.]+)')
		self.log = []
		for ll in lines:
			g = m.search(ll)
			if g != None:
				self.log.append([float(g.group('it')), float(g.group('err'))])
		self.log = N.array(self.log)
		plt.figure()
		plt.xlabel("iterations")
		plt.ylabel("change in model")
		plt.plot(self.log[:,0], self.log[:,1])
		plt.show()

	def viewRecon(self, recon="latest"):
		if recon == "latest":
			fn = G.glob("recon*.dat")[-1]
			print("Opening %s"%(fn))
			d = N.fromfile(fn, sep=" ")
		else:
			d = N.fromfile(recon, sep=" ")

		self.fg = d[:self.l*self.l].reshape(self.l,-1)
		self.bg = d[self.l*self.l:].reshape(self.l,-1)
		plt.figure()
		plt.subplot(1,3,1)
		plt.title("true signal")
		plt.imshow(self.contrast)

		plt.subplot(1,3,2)
		plt.title("recon.\nbackg.+sig.")
		plt.imshow(self.fg)
		
		plt.subplot(1,3,3)
		plt.title("recon.\nbackg.+sig.")
		plt.imshow(self.bg)
		
		plt.show()

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


if __name__ == "__main__":
	emc = runEMC()

