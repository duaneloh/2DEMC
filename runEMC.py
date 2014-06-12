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
		#If contrast.dat is present, use it to initialize reconstruction parameters.
		if os.path.isfile("contrast.dat"):
			tmp = N.fromfile("contrast.dat", sep=" ")
			self.l = int(tmp[0])
			self.numConf = int(tmp[1])
			self.contrast = (tmp[2:]).reshape(self.l, self.l)
		else:
			print("Warning: contrast.dat is missing. You'll have to make one.")
			self.l = 50
			self.contrast = N.random.rand(self.l, self.l)
			self.numConf = 1
		#Default signal parameters
		self.numPatterns = 10000
		self.sigLvl = 10.
		#Default has no blanks
		self.hasBackground = False
		self.hitRate = 1.
		self.bgLvl = 0.

		#Default number of iterations
		self.numIter = 100

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

	def removeBinaries(self):
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

	def recompile(self):
		self.removeBinaries()
		self.compile()

	def deleteIntermediates(self):
		files = ["EMC.log", "data.dat", "background.dat", "hidden_variables.dat", "hitFractions.dat"]
		files += G.glob("recon*.dat")
		files += G.glob("cond_prob*.dat")
		for f in files:
			if os.path.isfile(f):
				os.remove(f)

	def makeData(self):
		"""
		The number of conformations (-c flag) always defaults to here. 
		Reconstructing with multiple conformations is possible, but breaks the examples in runEMC.py.  
		"""
		cmd = "./make_data -d %d -s %lf -n %lf -h %lf -c 1"%(self.numPatterns, self.sigLvl, self.bgLvl, self.hitRate)	
		print("Making data with the following command:\n%s"%cmd)
		os.system(cmd)
		print("data.dat created.")
		
	def sampleData(self):
		#Check if data.dat is present
		#Read and store the sparse data format as a dictionary
		#View random images in this.
		pass
	
	def sparseToDense(self, imgSlice):
		#Convert a slice of the sparse data for the format into a dense array
		#This function should return an array
		pass

	def recon(self):
		cmd = "./EMC -n %d"%(self.numIter)	
		if self.hasBackground:
			cmd += " -b"
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
			fn = "recon%03d.dat"%(self.numIter)
			print("Opening %s"%(fn))
			d = N.fromfile(fn, sep=" ")
		else:
			d = N.fromfile(recon, sep=" ")
		if self.hasBackground:
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
		else:
			self.fg = d[:self.l*self.l].reshape(self.l,-1)
			plt.figure()
			plt.subplot(1,2,1)
			plt.title("true signal")
			plt.imshow(self.contrast)
			plt.subplot(1,2,2)
			plt.title("recon. sig.")
			plt.imshow(self.fg)
			plt.show()

	def viewCondProb(self):
		f = open("hitFractions.dat", "r")
		d = [float(i) for i in (f.readline()).split("\t")]
		[self.hitPh, self.hitPhSq, self.blankPh, self.blankPhSq, self.recHitRate, self.recFalseHitRate, self.recBlankRate, self.recFalseBlankRate] = d
		f.close()
		print("-"*50)
		print("Reconstructed hit rate:			%lf"%self.recHitRate)
		print("Actual hit rate:			%lf"%self.hitRate)
		print("Reconstructed hit purity:		%lf"%(1.-self.recFalseHitRate - self.recFalseBlankRate))
		print("Reconstructed false hit rate:		%lf"%self.recFalseHitRate)
		print("Reconstructed false blank rate:		%lf"%self.recFalseBlankRate)
		print("-"*50)

	def runDefaultCase(self):
		print("."*50)
		print("Running default case.")
		print("\tNo background noise (and hitRate=1).")
		print("\tNumber of patterns (hits): %d"%self.numPatterns)
		print("\tAverage photons per hit: %lf"%self.sigLvl)
		print("\tNumber of reconstruction iterations: %d"%self.numIter)
		print("."*50)
		self.makeData()
		print("."*50)
		self.recon()
		self.viewLog()
		self.viewRecon()

	def runExample(self, numPatterns=10000, numIter=100,sigLvl=10., bgLvl=10., hitRate=0.5,):
		print("."*50)
		print("Running case with background.")
		self.hasBackground = True
		self.bgLvl = bgLvl
		self.numPatterns = numPatterns
		self.sigLvl = sigLvl
		self.hitRate = hitRate
		self.numIter = numIter

		print("\tHit rate: %lf."%self.hitRate)
		print("\tBackground noise level: %lf photons per blank."%self.bgLvl)
		print("\tNumber of patterns (hits+blanks): %d"%self.numPatterns)
		print("\tAverage photons per hit: %lf"%self.sigLvl)
		print("\tNumber of reconstruction iterations: %d"%self.numIter)
		print("."*50)
		self.makeData()
		print("."*50)
		self.recon()
		self.viewCondProb()
		self.viewLog()
		self.viewRecon()

if __name__ == "__main__":
	emc = runEMC()

