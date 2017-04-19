    
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



histomax = 700.0
histomin = 0.0
histobins = 100
maxconfs=1000
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

        value = int(self.histobins*(float(nucsize) - float(self.histomin))/(float(self.histomax-self.histomin)))
	#print value
	if value<len(self.histo):
        	self.histo[value]+=addvalue
	else:
		print "weird value"
		print value
		print nucsize
    
    def getBoxX(self,hbox):
        x = (float(hbox)*float(self.histomax-self.histomin) )/float(self.histobins) + float(self.histomin)
        return x


#main function that is to be called
def MakeStructureHistogram(pathtype,manual=False):
    """
    special function to make histograms
    hardcoded. Remove at some point.
    """
    #set up histograms
    seed_nuc = Histogram(histomin,histomax,histobins)
    seed_sur = Histogram(histomin,histomax,histobins)
    seed_dist = Histogram(histomin,histomax,histobins)
    count=0

    if manual==False:
        interfacelist = helpers.generate_intflist()
    else:
        interfacelist = helpers.read_intflist()
    

    for interface in interfacelist:
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
            path = path.strip()
            pathpath= os.path.join(intfpath,path)
            identifier = interface+path

            histofile = os.path.join(pathpath,(identifier+'.clu.dat'))
	    
            if os.path.exists(histofile):
                nucsize,mindist,seedincluster,seedinsurface = np.loadtxt(histofile,comments='#',usecols=(1,2,3,7),unpack=True)
                nucsize = nucsize.astype(int)
                seedincluster = seedincluster.astype(int)
                seedinsurface = seedinsurface.astype(int)
            else:
                continue
            #loooping over each slice in the trajectory
            for i in range(len(nucsize)):
                seed_nuc.addAtomtoHisto(nucsize[i],seedincluster[i])
                seed_sur.addAtomtoHisto(nucsize[i],seedinsurface[i])
                seed_dist.addAtomtoHisto(nucsize[i],mindist[i])
                count+=1
    
    #normalise the histograms
    #histogram x values
    histox = np.zeros(len(seed_nuc.histox))

    for i in range(len(seed_nuc.histox)):
        
        #if sum_sur>0:
	seed_nuc.histo[i]/=float(count)
        seed_sur.histo[i]/=float(count)
        seed_dist.histo[i]/=float(count)
        
        histox[i] = seed_nuc.getBoxX(i)

    histo_nuc = np.column_stack((histox,seed_nuc,seed_sur,seed_dist))
  
    np.savetxt('averaged_cluster_histo.dat',histo_sur)
  
    print count
 
if __name__=='__main__':

    MakeStructureHistogram('AB',manual=False)


            

            
            
