    
import os
import sys
import subprocess as sub
import numpy as np
import time
import logging
import tistools_helpers.tistools_helpers as tistools_helpers

#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("analysis.log")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False 



histomax = 19.0
histomin = 0.0
histobins = 19
dhistomax = 20.0
dhistomin = 0.0
dhistobins = 40
maxconfs=10
tstcluster = 404
#create helpers class
helpers = tistools_helpers.tistools_helpers()


#class for histogram
class Histogram(object):
    
    def __init__(self,histomin,histomax,histobins):
        self.histomax = histomax
        self.histomin = histomin
        self.histobins = histobins
        self.histo = np.zeros(histobins)
        self.histox = np.linspace(self.histomin,self.histomax,self.histobins)

    def addAtomtoHisto(self,nucsize,addvalue):
	#print nucsize
        value = int(self.histobins*(float(nucsize) - float(self.histomin))/(float(self.histomax-self.histomin)))
	#print value
	#print value
	if value<len(self.histo):
        	self.histo[value]+=addvalue
	else:
		print "weird value"
		print value
		print nucsize
    	return value
	
    def getBoxX(self,hbox):
        x = (float(hbox)*float(self.histomax-self.histomin) )/float(self.histobins) + float(self.histomin)
        return x


#main function that is to be called
def MakeStructureHistogram(pathtype,manual=False):
    """
    special function to make histograms
    hardcoded. Remove at some point.
    """
    total=0
    #set up histograms
    seed_nuc = Histogram(histomin,histomax,histobins)
    seed_dist = Histogram(dhistomin,dhistomax,dhistobins)
    count=np.zeros(histobins)
    dcount=np.zeros(dhistobins)

    if manual==False:
        interfacelist = helpers.generate_intflist()
    else:
        interfacelist = helpers.read_intflist()
    

    for interface in interfacelist:
        if total>maxconfs:
                break
        interface = interface.strip()
        intfpath = os.path.join(os.getcwd(),"tis","la",interface)
        intfpath = intfpath.strip()
        pathpath = os.path.join(intfpath,pathtype+".dat")
        pathpath = pathpath.strip()
        pathlist = []
        filenamelist = []

        #we get the list of all paths that needs to be analysed
        for path in open(pathpath,'r'):
            pathlist.append(path)

        #may the analysis start
        for path in pathlist:
            if total>maxconfs:
                break
            path = path.strip()
            pathpath= os.path.join(intfpath,path)
            identifier = interface+path

            histofile = os.path.join(pathpath,(identifier+'.clu.dat'))
	    
            if os.path.exists(histofile):
                nucsize,mindist,seedincluster,a,b,c,seedinsurface,d = np.loadtxt(histofile,comments='#',unpack=True)
		#print nucsize
                nucsize = nucsize.astype(int)
		#print nucsize
                seedincluster = seedincluster.astype(int)
                seedinsurface = seedinsurface.astype(int)
            else:
                continue
            #loooping over each slice in the trajectory
            for i in range(len(nucsize)):
                if (nucsize[i] <= tstcluster+3) and (nucsize[i] >= tstcluster-3):
                        if total>maxconfs:
                                break
                        value = seed_nuc.addAtomtoHisto(seedincluster[i],1)
                        if value<len(count):
                                count[value]+=1
                        value = seed_dist.addAtomtoHisto(mindist[i],1)
                        if value<len(dcount):
                                dcount[value]+=1
                        total+=1


    
    #normalise the histograms
    #histogram x values
    histox = np.zeros(len(seed_nuc.histox))
    dhistox = np.zeros(len(seed_dist.histox))

    for i in range(len(seed_nuc.histox)):
        
        if count[i]>0:
		seed_nuc.histo[i]/=float(count[i])
        
        histox[i] = seed_nuc.getBoxX(i)

    for i in range(len(seed_dist.histox)):
        
        if dcount[i]>0:
                seed_dist.histo[i]/=float(dcount[i])

        
        dhistox[i] = seed_dist.getBoxX(i)

    histo_nuc = np.column_stack((histox,seed_nuc.histo))
    histo_dist = np.column_stack((dhistox,seed_dist.histo))
  
    np.savetxt('averaged_cluster_nuc.dat',histo_nuc)
    np.savetxt('averaged_cluster_dist.dat',histo_dist)
  
 
if __name__=='__main__':

    MakeStructureHistogram('AB',manual=False)


            

            
            
