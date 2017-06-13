#!/usr/bin/python
"""
Set of functions written for analysis of data after a tis run using
TPS wrapper.
These functions are called in different scripts to perform the different
functions.
"""
import os
import sys
import subprocess as sub
import numpy as np
import time
import logging
import gzip
import shutil

#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("analysis.log")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False 

#class for seed
class Seed(object):
    
    def __init__(self,seedfileaddress):
        #self.atoms=np.empty([natoms,5])
        self.seedfile=seedfileaddress
        self.exists=True
    
    def ReadAlles(self,filename):
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


    #read the positions of seed atoms
    def ReadSeed(self,read=True,idlist=None):
        
        self.seedids = []
        if read==True:
                if os.path.isfile(self.seedfile) and os.path.getsize(self.seedfile) > 0:
                        for line in open(self.seedfile):
                                self.seedids.append(int(line.strip()))
                        self.atoms=np.empty([len(self.seedids),5])
                else:
                        self.exists=False
        else:
                self.seedids=idlist
                if len(self.seedids)>0:
                        self.atoms=np.empty([len(self.seedids),5])
                else:
                        self.exists=False

    
    #populate the positions of seed atoms
    def PopulateSeed(self,filename,read=True):
        if read==True:
            atoms = self.ReadAlles(filename)
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
        loo = []
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
	    #print mindist
	    #print (mindist<1e-5)
	    if mindist<1e-5:
		#print "oh"
		mindist=0.00
            atom[4]=mindist
	    loo.append(mindist)
	return loo