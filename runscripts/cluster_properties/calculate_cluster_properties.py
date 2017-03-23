"""
- calculates cluster props.

- to be used with special version of order parameter
  prints out all the clusters.
  found in the same folder.

- needs a configuration for seed.

- needs a list of seed atoms.

- needs a list of seed surface atoms.

- requires tistools module
 
"""

import numpy as np
import sys
import os
from  tistools_helpers import tistools_helpers
import logging
import subprocess as sub

#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("analysis.log")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False 

#create helpers class
helpers = tistools_helpers.tistools_helpers()

#function to read dump files
def read_alles(filename,filetype="dump"):

	if (filetype=="dump"):
		#ninenumber of lines are not required
		#after that column 0,3,4,5 to be read.
		count = 0
		data = []
		for line in open(filename,'r'):
			data.append(line)

		boxsizelist = []
		natoms = int(data[3])
		atoms = np.empty([natoms,5])
		i = 0
		for line in data:
			if (count==5) or (count==6) or (count==7):
				raw = line.split()
				boxsizelist.append(float(raw[0]))
				boxsizelist.append(float(raw[1]))

			elif (count>8):
				raw = line.split()
				atoms[i][0] = int(raw[0])
				atoms[i][1] = float(raw[3])
				atoms[i][2] = float(raw[4])
				atoms[i][3] = float(raw[5])
				atoms[i][4] = 99999.00
				#atoms[i][4] = False
				i+=1
			count+=1

		#print atoms
		#print boxsizelist
		return atoms


#class for seed
class Seed(object):
	
	def __init__(self,seedfileaddress):
		#self.atoms=np.empty([natoms,5])
		self.seedfile=seedfileaddress
	
	#read the positions of seed atoms
	def ReadSeed(self):
		self.seedids = np.loadtxt(self.seedfile,unpack=True,dtype=int)
		self.atoms=np.empty([len(self.seedids),5])

	
	#populate the positions of seed atoms
	def PopulateSeed(self,filename,read=True):
		if read==True:
			atoms = read_alles(filename)
		else:
			atoms = filename

		#add the atoms and positions
		k = 0
		for i in range(len(atoms)):
			if atoms[i][0] in self.seedids:
				self.atoms[k][0]=atoms[i][0]
				self.atoms[k][1]=atoms[i][1]
				self.atoms[k][2]=atoms[i][2]
				self.atoms[k][3]=atoms[i][3]
				k+=1

	
	def CalculateDistances(self,otheratoms):
		for atom in self.atoms:
			dist = []
			for oatom in otheratoms.atoms:
				#print 'seed'
				#print seedatom[0]
				a = oatom[1]
				b = oatom[2]
				c = oatom[3]
				distance = np.sqrt((a-atom[1])**2 + (b-atom[2])**2 + (c-atom[3])**2 )
				dist.append(distance)
			mindist=min(dist)
			atom[4]=mindist


	def EvaluateMinDistance(self):
		
		min_d = 100000.00 
		for atomito in self.atoms:
			#a=1
			#print type(atomito[4])
			#print type(min_d)
			if atomito[4]<min_d:
				min_d=atomito[4]
		return min_d

	
	def CalculateInclusion(self,otheratoms):

		commonatoms = list(set(self.seedids).intersection(otheratoms.seedids))
		nocommonatoms = len(commonatoms)
		seedsize = len(self.seedids)
		seedinclusion = float(nocommonatoms)/float(seedsize)

		return nocommonatoms,seedinclusion



if __name__=="__main__":

	traj = sys.argv[1]
	outputfile = sys.argv[2]
	tmpname = sys..argv[3]
	binary = sys.argv[4]
	seedfileaddress = sys.argv[5]
	seedsurfaceaddress = sys.argv[6]
	gzip = sys.argv[7]

	if gzip=='True':
		gzip=True
	else:
		gzip=False



	#first create the seed class
	 
	seed = Seed(seedfileaddress)
	seed.ReadSeed()
	seeddump = 'conf.dump'
	seed.PopulateSeed(seeddump)


	#first create the seed surface class
	 
	seedsurface = Seed(seedsurfaceaddress)
	seedsurface.ReadSeed()
	seeddump = 'conf.dump'
	seedsurface.PopulateSeed(seeddump)

	#now we can start the reading of data
	#data = helpers.separate_traj(traj)
	data = helpers.combine_paths_return(traj,gzip=gzip)
	#needs path number argument, change it later
	#data = helpers.combine_paths_return(traj)
	fout = open(outputfile,'w')

	#write header
	fout.write(('# %10s  %10s  %10s  %10s %15s %10s %15s %15s\n')%('nucsize','mindist','seedincluster','inclusion','seedsurfaceincluster','inclusion','seedinsurface','seedinsurfacep'))

    	for i in range(len(data)):

        	tmpfile = tmpname+str(i)+".dat"
        	outfile = open(tmpfile,'w') 
        	for j in range(len(data[i])):
            		outfile.write(data[i][j])
            	timestep = data[i][1]
        	outfile.flush()
        	outfile.close()
       
        	cmd = []
        	cmd.append(binary)
        	cmd.append(tmpfile)
        	proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
        	out,err = proc.communicate(input="")
        	proc.wait()
        	nucsize = out
        	
	        #read the atoms in
	        atoms = read_alles(tmpfile)
	        
	        #first create the surface class
		surfacefileaddress = "surface.dat" 
		surface = Seed(surfacefileaddress)
		surface.ReadSeed()
		surface.PopulateSeed(atoms,read=False)


		#first create the surface class
		clusterfileaddress = "cluster.dat" 
		cluster = Seed(clusterfileaddress)
		cluster.ReadSeed()
		cluster.PopulateSeed(atoms,read=False)

	    	#calculate distance from surface atoms to seed
	    	surface.CalculateDistances(seed)
	    	min_dist = surface.EvaluateMinDistance()


	    	#calculate inclusion of seed in cluster
	    	numberofatoms,percent = seed.CalculateInclusion(cluster)

	    	#calculate inclusion of seed surface in cluster
	    	surfacenumberofatoms,surfacepercent = seedsurface.CalculateInclusion(cluster)

	    	#seed in surface atoms
	    	seedinsurface = len(list(set(surface.seedids).intersection(seed.seedids)))
	    	seedinsurfacep = float(seedinsurface)/float(len(seed.seedids))


		fout.write(('%10s  %10.4f  %10d  %10.4f %15d  %10.4f  %15d  %15.4f\n')%(nucsize.strip(),min_dist,numberofatoms,percent,surfacenumberofatoms,surfacepercent,seedinsurface,seedinsurfacep))
		
	fout.close()


	