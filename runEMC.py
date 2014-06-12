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
import struct

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
		self.sparseData = {}
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
		
	def sampleData(self, readFrac=0.1):
		#Check if data.dat is present
		#Read and store the sparse data format as a dictionary
		if not os.path.isfile("data.dat"):
			print("data.dat is not present. Skipping file read.")
			return 1
		#Start to parse file
		f = open("data.dat", "rb")
		inD = f.read(32)
		(self.num_data,) = struct.unpack('i', inD[:4])
		(self.present_conf,) = struct.unpack('i', inD[4:8])
		(self.mean_total_photons,) = struct.unpack('d', inD[8:16])
		(self.len_,) = struct.unpack('i', inD[16:20])
		(self.tomo_len,) = struct.unpack('i', inD[20:24])
		(self.m_info_given_choice,) = struct.unpack('d', inD[24:32])
		#Print header for sanity check
		print("num_data: %d"%self.num_data)
		print("present_conf: %d"%self.present_conf)
		print("mean_total_photons: %lf"%self.mean_total_photons)
		print("len: %d"%self.len_)
		print("tomo_len: %d"%self.tomo_len)
		print("m_info_given_choice: %lf"%self.m_info_given_choice)
		#Flush out existing sparseData dictionary.
		self.sparseData = {}
		#Read only a fraction of images. 
		for dNum in range(int(readFrac*self.num_data)):
			self.sparseData[dNum] = {}
			#Read the one photons
			bNumOnes = f.read(2)
			(numOnes,) = struct.unpack("H", bNumOnes)
			bOnes = f.read(2*numOnes)
			self.sparseData[dNum]['o'] = []
			for o in range(numOnes):
				(oLoc,) = struct.unpack("H", bOnes[2*o:2*(o+1)]) 
				(self.sparseData[dNum]['o']).append(oLoc)
			#Read the multiple photons
			bNumMultis = f.read(2)
			(numMultis,) = struct.unpack("H", bNumMultis)
			bMultis = f.read(2*numMultis)
			self.sparseData[dNum]['m'] = []
			for m in range(numMultis/2):
				(mLoc,) = struct.unpack("H", bMultis[4*m:4*m+2])
				(mPh,) = struct.unpack("H", bMultis[4*m+2:4*m+4])
				(self.sparseData[dNum]['m']).append([mLoc, mPh])
		f.close()

	def sparseToDense(self, imgSlice):
		#Convert a slice of the sparse data for the format into a dense array
		#This function should return an array
		if len(self.sparseData) == 0:
			res = self.sampleData()
			if res is 1:
				return 
		arr = N.zeros(self.len_*self.len_)
		cD = self.sparseData[imgSlice]
		for o in cD['o']:
			arr[o] = 1.
		for m in cD['m']:
			arr[m[0]] = arr[m[1]]
		return arr.reshape(self.len_,-1)

	def showNineRandomData(self):
		if len(self.sparseData) == 0:
			res = self.sampleData()
			if res is 1:
				return 
		keys = self.sparseData.keys()
		from random import shuffle
		shuffle(keys)
		fig,ax = plt.subplots(3,3, figsize=(9,10))
		ccmap = PL.cm.get_cmap('bone', 2)
		for nn,kk in enumerate(keys[:9]):
			arr = self.sparseToDense(kk)
			plt.subplot(3,3,nn+1)
			plt.imshow(arr, cmap=ccmap, vmax=2)	
			plt.title("Img %d"%kk)
		ax2 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
		cb = matplotlib.colorbar.ColorbarBase(ax2, cmap=ccmap, ticks=range(3), format='%1i')
		plt.show()

	def showAverageData(self):
		if len(self.sparseData) == 0:
			res = self.sampleData()
			if res is 1:
				return 
		keys = self.sparseData.keys()
		hostArr = N.zeros((self.len_,self.len_))
		for nn,kk in enumerate(keys):
			arr = self.sparseToDense(kk)
			hostArr += arr
		hostArr /= 1.*len(keys)
		fig = plt.figure(figsize=(6,5))
		plt.imshow(hostArr, cmap='bone')	
		plt.title("%d-image average data"%len(keys))
		ax2 = fig.add_axes([0.90, 0.1, 0.03, 0.8])
		plt.title("photons")
		cb = matplotlib.colorbar.ColorbarBase(ax2, cmap='bone', ticks=N.arange(0,1,0.1))
		plt.show()

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
		plt.figure(figsize=(9,5))
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
			plt.figure(figsize=(9,3))
			plt.subplot(1,3,1)
			plt.title("true signal")
			plt.imshow(self.contrast, cmap='bone')
			plt.subplot(1,3,2)
			plt.title("recon.\nbackg.+sig.")
			plt.imshow(self.fg, cmap='bone')
			plt.subplot(1,3,3)
			plt.title("recon.\nbackg.+sig.")
			plt.imshow(self.bg, cmap='bone')
			plt.show()
		else:
			self.fg = d[:self.l*self.l].reshape(self.l,-1)
			plt.figure(figsize=(6,3))
			plt.subplot(1,2,1)
			plt.title("true signal")
			plt.imshow(self.contrast, cmap='bone')
			plt.subplot(1,2,2)
			plt.title("recon. sig.")
			plt.imshow(self.fg, cmap='bone')
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
		self.showAverageData()
		self.showNineRandomData()

	def runExample(self, numPatterns=10000, numIter=100,sigLvl=10., bgLvl=10., hitRate=0.5, hasBackground=True):
		print("."*50)
		print("Running case with background.")
		self.hasBackground = hasBackground
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
		self.showAverageData()
		self.showNineRandomData()

if __name__ == "__main__":
	emc = runEMC()
	if not (os.path.isfile("make_data") and os.path.isfile("EMC")):
		emc.recompile()

	try:
		if sys.argv[1] == 'd':
			print(("="*100+"\n")*3)
			print("Starting with the default background-free case (Hit rate = 100%)")
			emc.runDefaultCase()
			emc.deleteIntermediates()
			print(("="*100+"\n")*3)
		
			raw_input("Press any key to continue.\Next: background included but too few dataset.")
			emc.runExample()
			emc.deleteIntermediates()
		 	print(("="*100+"\n")*3) 
		
			raw_input("Press any key to continue.\Next: background included now with 10 times more datasets.")
			emc.runExample(numPatterns=100000)
			emc.deleteIntermediates()
		 	print(("="*100+"\n")*3) 
	except:
		pass
	
