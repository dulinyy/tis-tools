"""
- calculates time evolution of all clusters.

- to be used with special version of order parameter

- prints out all the clusters along with various properties.

- requires seed list: give the list:

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
		self.seedids = np.loadtxt(seedfileaddress,unpack=True,dtype=int)
		self.atoms=np.empty([len(self.seedids),4])

	#populate the positions of seed atoms
	def PopulateSeed(self,filename):
		atoms = read_alles(filename)
		#add the atoms and positions
		k = 0
		for i in range(len(atoms)):
			if atoms[i][0] in self.seedids:
				self.atoms[k][0]=atoms[i][0]
				self.atoms[k][1]=atoms[i][1]
				self.atoms[k][2]=atoms[i][2]
				self.atoms[k][3]=atoms[i][3]
				k+=1


class Cluster(object):
	#will contain the atoms in each timestep
	

	def __init__(self,id):
		self.id = id
		self.atoms = []
		self.min_distance = []
		self.time_evolution = []
		self.lasttimestepatomids = []
		self.atomlist=[]
		self.seedinclusion = []
		self.dead=True

	
	def ClusterAtTime(self,clusterids,atoms):
		self.atomlist = []
		dummyatom = np.empty([5])
		self.lasttimestepatomids = []
		#clusterids = set(clusterids)
		#print clusterids
		for i in range(len(atoms)):
			dummyatom = np.empty([5])
			if atoms[i][0] in clusterids:
				dummyatom[0]=atoms[i][0]
				dummyatom[1]=atoms[i][1]
				dummyatom[2]=atoms[i][2]
				dummyatom[3]=atoms[i][3]
				
				self.atomlist.append(dummyatom)
			#print(len(self.atomlist))
		self.atoms.append(self.atomlist)

	def DistanceToSeed(self,seed):
		for atom in self.atomlist:
			#print 'main'
			#print atom[0]
			dist = []
			for seedatom in seed.atoms:
				#print 'seed'
				#print seedatom[0]
				a = seedatom[1]
				b = seedatom[2]
				c = seedatom[3]
				distance = np.sqrt((a-atom[1])**2 + (b-atom[2])**2 + (c-atom[3])**2 )
				dist.append(distance)
			mindist=min(dist)
			atom[4]=mindist
			dummyagain = self.atomlist
		#self.atomlist.append(dummyagain)
		self.time_evolution.append(len(dummyagain))

	def EvaluateMinDistance(self):
		min_d = 100000.00 
		for atomito in self.atomlist:
			a=1
			#print type(atomito[4])
			#print type(min_d)
			if atomito[4]<min_d:
				min_d=atomito[4]
		self.min_distance.append(min_d)

	def FillLastTimestepAtoms(self):
		for dummyatom in self.atomlist:
			self.lasttimestepatomids.append(dummyatom[0])


def ProcessBinaryOutput():
	clusterids,clusters = np.loadtxt("clusterlist.dat",unpack=True,dtype=int)
	#print clusters
	unique_clusters = np.unique(clusters)
	#print unique_clusters

	atomsincluster = [[] for i in range(len(unique_clusters))]
	#print atomsincluster
	#print len(unique_clusters)

	for i in range(len(unique_clusters)):
		for j in range(len(clusterids)):
			if clusters[j]==unique_clusters[i]:
				atomsincluster[i].append(clusterids[j])

	#print atomsincluster
	return atomsincluster

def DecideWhichCluster(bigclusterlist,atomsincluster,atoms,seed,idstatus):
	#first check for existing clusters:
	#for each new cluster
	#print atomsincluster
	for newcluster in atomsincluster:
		#print 'cluster'
		#print newcluster
		assigned=False
		#print seed.seedids
		seedincluster = list(set(newcluster).intersection(seed.seedids))
		#seedincluster = [x for x in newcluster and seed.seedids]
		#print seedincluster
		seedsize = len(seed.seedids)
		seedinclusion = float(len(seedincluster))/float(seedsize)
		#set all clusters as dead

		#for existing clusters
		for cluster in bigclusterlist:
			commonids = list(set(newcluster).intersection(cluster.lasttimestepatomids))
			#commonids = [x for x in newcluster and cluster.lasttimestepatomids]
			if len(commonids) > 3:
				#print 'coo'
				assigned=True
				#now assign this guy to this cluster
				cluster.ClusterAtTime(newcluster,atoms)
				cluster.DistanceToSeed(seed)
				cluster.EvaluateMinDistance()
				cluster.seedinclusion.append(seedinclusion)
				#reverse death of clusters that propogate
				cluster.dead=False
				break
		if assigned==False:
			#print "fila"
			#print idstatus
			neucluster = Cluster(idstatus)
			neucluster.ClusterAtTime(newcluster,atoms)
			neucluster.DistanceToSeed(seed)
			neucluster.EvaluateMinDistance()
			neucluster.seedinclusion.append(seedinclusion)
			#set the new cluster as alive
			neucluster.dead=False
			bigclusterlist.append(neucluster)
			#print bigclusterlist
			idstatus+=1

	for cluster in list(bigclusterlist):
		if cluster.dead==True:
			bigclusterlist.remove(cluster)
		else:
			cluster.FillLastTimestepAtoms()

	for cluster in bigclusterlist:
		cluster.dead=True

	return bigclusterlist,idstatus


if __name__=="__main__":


	traj = sys.argv[1]
	outputfile = sys.argv[2]
	tmpname = sys.argv[3]
	binary = sys.argv[4]
	seedfileaddress = sys.argv[5]
	gzip = sys.argv[7]

	if gzip=='True':
		gzip=True
	else:
		gzip=False

	

	#first create the seed class
	seed = Seed(seedfileaddress)
	seed.ReadSeed()
	#give one snapshot with the seed positions. Set it globally so that rereading can be avoided.
	seeddump = 'conf.dump'
	seed.PopulateSeed(seeddump)


	#now we can start the reading of data
	#data = helpers.separate_traj(traj)
	data = helpers.combine_paths_return(path,gzip=gzip)
	#needs path number argument, change it later
	#data = helpers.combine_paths_return(traj)
	fout = open(outputfile,'w')

	#tmpfilelist = []
	qtraj = []
	idstatus = 0
	bigclusterlist = []

	
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
        	qd = out
        	print i
        	#os.system(("rm %s")% tmpfile)
        
        	#opfile.write(("%s")%(qd))
        	#opfile.flush()
        	qtraj.append(qd)
   			
	        #now process the time slice.
	        #gives the output array; len of array number of clusters
	        #each slice has the ids of atoms in that
	        atomsincluster = ProcessBinaryOutput()
	        #read the time slice again
	        atoms = read_alles(tmpfile)
	        #now assign cluster
	        bigclusterlist,idstatus = DecideWhichCluster(bigclusterlist,atomsincluster,atoms,seed,idstatus)

	        #print len(bigclusterlist)
	        #print atomsincluster
	    	#now remove tmpfile
	    	os.system(("rm %s")% tmpfile)
	    	fout.write('TIMESTEP\n')
	    	fout.write(('%s')%timestep)
	    	fout.write(('%-10s %-10s %-10s %-10s\n')%('CLUSTERID','NUCSIZE','MINDIST','INCLUSION'))
		#writing the output
		for k in range(len(bigclusterlist)):
			fout.write(('%10d  %10d  %10.4f  %10.4f\n')%(bigclusterlist[k].id,bigclusterlist[k].time_evolution[-1],bigclusterlist[k].min_distance[-1],bigclusterlist[k].seedinclusion[-1]))
		fout.write('ENDOFTIMESTEP\n')	

	fout.close()


	
