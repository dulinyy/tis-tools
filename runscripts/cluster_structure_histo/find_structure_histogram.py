"""
- calculates cluster props.

- to be used with special version of order parameter
  prints the distribution
  found in the same folder.

- needs a configuration for seed.

- needs a list of seed atoms.

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
		self.exists=True
	
	#read the positions of seed atoms
	def ReadSeed(self,seedignore=False,seedclass=None):
		
		self.seedids = []
		if os.path.isfile(self.seedfile) and os.path.getsize(self.seedfile) > 0:
			for line in open(self.seedfile):
				self.seedids.append(int(line.strip()))
			if seedignore==True:
				dummy = [x for x in self.seedids if x not in seedclass.seedids]
				self.seedids = dummy
				


			self.atoms=np.empty([len(self.seedids),5])
		else:
			self.exists=False

	
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

class Histogram(object):
	
	def __init__(self,histomin,histomax,histobins):
		self.histomax = histomax
		self.histomin = histomin
		self.histobins = histobins
		self.histo = np.zeros(histobins)
		self.histox = np.linspace(self.histomin,self.histomax,self.histobins)

	def addAtomtoHisto(self,atom,addvalue):
		distance = atom[4]
		value = (float(self.histobins)*(float(distance) - float(self.histomin)))/(float(self.histomax-self.histomin))
		value = int(round(value))
		self.histo[value]+=addvalue

def AssignHistograms(atomsclass,histogram,nucsize):
	for atom in atomsclass.atoms:
		addvalue = 1.00/float(nucsize)
		histogram.addAtomtoHisto(atom,addvalue)
	return histogram


if __name__=="__main__":

	traj = sys.argv[2]
	outputfile = sys.argv[4]
	tmpname = sys.argv[3]
	binary = sys.argv[1]
	gzip = sys.argv[5]	
	
	deletefiles=True
	seedfileaddress = '/home/users/menonsqr/SeedFCC19/seed.dat'
	seeddump = '/home/users/menonsqr/SeedFCC19/conf.dump'
	

	if gzip=='True':
		gzip=True
	else:
		gzip=False

	
	#first create the seed class 
	seed = Seed(seedfileaddress)
	seed.ReadSeed()
	seed.PopulateSeed(seeddump)

	#set up histogram variable
	histomin = 0.0
	histomax = 60.00
	histobins = 6000

	#set up histograms
	bcchisto = Histogram(histomin,histomax,histobins)
	fcchisto = Histogram(histomin,histomax,histobins)
	hcphisto = Histogram(histomin,histomax,histobins)
	udfhisto = Histogram(histomin,histomax,histobins)

	sbcchisto = Histogram(histomin,histomax,histobins)
	sfcchisto = Histogram(histomin,histomax,histobins)
	shcphisto = Histogram(histomin,histomax,histobins)
	sudfhisto = Histogram(histomin,histomax,histobins)	
	#now we can start the reading of data
	#data = helpers.separate_traj(traj)
	data = helpers.combine_paths_return(traj,gzip=gzip)
	if deletefiles:
		forpath = os.path.join(traj,'forward')
        	bakpath = os.path.join(traj,'backward')
        	os.system(('rm -rf %s')%(forpath))
        	os.system(('rm -rf %s')%(bakpath))

	#needs path number argument, change it later
	#data = helpers.combine_paths_return(traj)



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
        	distribution = out
        	out = out.split()

        	#read in the output
        	nucsize = out[0]
        	fbcc = out[1]
        	ffcc = out[2]
        	fhcp = out[3]
        	fudf = out[4]
        	
	        #read the atoms in
	        atoms = read_alles(tmpfile)
	        os.system(('rm %s')% tmpfile) 
	        #first create the surface class
		surfacefileaddress = "surface.dat" 
		surface = Seed(surfacefileaddress)
		surface.ReadSeed()
		if surface.exists:
			surface.PopulateSeed(atoms,read=False)

	        #first create the surface class
		bccfileaddress = "bccid.dat" 
		bcc = Seed(bccfileaddress)
		bcc.ReadSeed(seedignore=True,seedclass=seed)
		if bcc.exists:
			bcc.PopulateSeed(atoms,read=False)

	        #first create the fcc class
		fccfileaddress = "fccid.dat" 
		fcc = Seed(fccfileaddress)
		fcc.ReadSeed(seedignore=True,seedclass=seed)
		if fcc.exists:
			fcc.PopulateSeed(atoms,read=False)

	        #first create the hcp class
		hcpfileaddress = "hcpid.dat" 
		hcp = Seed(hcpfileaddress)
		hcp.ReadSeed(seedignore=True,seedclass=seed)
		if hcp.exists:
			hcp.PopulateSeed(atoms,read=False)

	        #first create the udf class
		udffileaddress = "udfid.dat" 
		udf = Seed(udffileaddress)
		udf.ReadSeed(seedignore=True,seedclass=seed)
		if udf.exists:
			udf.PopulateSeed(atoms,read=False)


	    	#calculate distance from surface atoms to seed
	    	if bcc.exists:
	    		bcc.CalculateDistances(surface)
	    		bcchisto = AssignHistograms(bcc,bcchisto,nucsize)
	    		bcc.CalculateDistances(seed)
	    		sbcchisto = AssignHistograms(bcc,sbcchisto,nucsize)
	    	if fcc.exists:
	    		#print 'boolo'
	    		fcc.CalculateDistances(surface)
	    		#print fcchits
	    		fcchisto = AssignHistograms(fcc,fcchisto,nucsize)
	    		#print fcc
	    		fcc.CalculateDistances(seed)
	    		#print fcc
	    		sfcchisto = AssignHistograms(fcc,sfcchisto,nucsize)
	    		#print fcc
	    	if hcp.exists:
	    		hcp.CalculateDistances(surface)
	    		hcphisto = AssignHistograms(hcp,hcphisto,nucsize)
	    		hcp.CalculateDistances(seed)
	    		shcphisto = AssignHistograms(hcp,shcphisto,nucsize)
	    	if udf.exists:
	    		udf.CalculateDistances(surface)
	    		udfhisto = AssignHistograms(udf,udfhisto,nucsize)
	    		udf.CalculateDistances(seed)
	    		sudfhisto = AssignHistograms(udf,sudfhisto,nucsize)

	#write out data
	#print fcchisto
	histo_surface = np.column_stack((bcchisto.histox,bcchisto.histo,fcchisto.histo,hcphisto.histo,udfhisto.histo,sbcchisto.histo,sfcchisto.histo,shcphisto.histo,sudfhisto.histo))

	np.savetxt(outputfile,histo_surface)


		
		
		
		






		
