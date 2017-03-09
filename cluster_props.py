"""
calculates cluster props. Edit it to make it nicer and parallelise alles
"""

import numpy as np
import sys
import os
import tistools as tt
#import matplotlib.pyplot as plt



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
		return atoms,boxsizelist


def calc_distance(atoms,centreid):
	
	#provide centre atom here
	for i in range(len(atoms)):
		if atoms[i][0]==centreid:
			centre = i
	
	a = atoms[centre][1]
	b = atoms[centre][2]
	c = atoms[centre][3]
	#print a,b,c
	
	#print centre


	for i in range(len(atoms)):
		#print "=================="
		#print i
		#print (a-atoms[i][1])**2
		#print (b-atoms[i][2])**2
		#print (c-atoms[i][3])**2
		distance = np.sqrt((a-atoms[i][1])**2 + (b-atoms[i][2])**2 + (c-atoms[i][3])**2 )
		atoms[i][4] = round(distance,2)

	return atoms


if __name__=="__main__":


	traj = sys.argv[1]
	centreid = int(sys.argv[2])
	outputfile = sys.argv[3]
	tmpname = sys.argv[4]
	binary = sys.argv[5]
	seedlist = np.loadtxt("/home/users/menonsqr/bin/seed.dat",unpack=True)


	#data = tt.combine_paths_return(traj)
	data = tt.separate_traj(traj)
	fout = open(outputfile,'w')

	tmpfilelist = []
	qtraj = []

    	for i in range(len(data)):
        	tmpfile = tmpname+str(i)+".dat"
        	outfile = open(tmpfile,'w') 
        	for j in range(len(data[i])):
            		outfile.write(data[i][j])
        		outfile.flush()
        	outfile.close()
        	tmpfilelist.append(tmpfile)	

	for tmp in tmpfilelist:
		#generate the ids of atoms in cluster
		os.system(("%s %s")%(binary,tmp))
		
	        #open that, calculate the minum distance
	        clusterlist = np.loadtxt("cluster.dat",unpack=True)
	        atoms,boxsizelist = read_alles(tmp)
	        atoms = calc_distance(atoms,centreid)
	
		sortedids = np.argsort(atoms[:,0])
		
		seedcount = len(seedlist)
		seedinclustercount = 0
	
		min_dist_sorted_args = np.argsort(atoms[:,4])

		found = False
		for i in range(len(min_dist_sorted_args)):
			if atoms[min_dist_sorted_args[i]][0] in clusterlist:
				mindist = atoms[min_dist_sorted_args[i]][4]
				found = True
				break

		if found==False:
			mindist = 999999999
			print "min dist not found"

		#find the percentage of atoms in the cluster
		for atom in seedlist:
			if atom in clusterlist:
				seedinclustercount+=1

		perc = float(seedinclustercount)/float(seedcount)

		#now the centre of mass code
		#
		for i in range(len(atoms)):
			if atoms[i][0]==centreid:
				centre = i
	
		a = atoms[centre][1]
		b = atoms[centre][2]
		c = atoms[centre][3]

		distx = 0.00
		disty = 0.00
		distz = 0.00
		distcount = 0
		
		for i in range(len(clusterlist)):
			for j in range(len(atoms)):
				if clusterlist[i]==atoms[j][0]:
					distx+=atoms[j][1]
					disty+=atoms[j][2]
					distz+=atoms[j][3]
					distcount+=1

		distx/=float(distcount)
		disty/=float(distcount)
		distz/=float(distcount)

		#print distx
		#print disty
		#print distz

		dcom = ((a-distx)**2+(b-disty)**2+(c-distz)**2)**0.5

		fout.write(("%f %f %f\n")%(mindist,perc,dcom))
		os.remove(tmp)

	fout.close()





