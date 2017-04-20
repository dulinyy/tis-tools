    
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

#seedfile
seedfileaddress = '/home/users/menonsqr/SeedFCC43/seed.dat'

#create helpers class
helpers = tistools_helpers.tistools_helpers()

#histogram variables
histomin = 0.0
histomax = 800.00
histobins = 200

#class for histogram
class Histogram(object):
    
    def __init__(self,histomin,histomax,histobins):
        self.histomax = histomax
        self.histomin = histomin
        self.histobins = histobins
        self.histo = np.zeros(histobins)
        self.histox = self.getBoxX(range(histobins))


    def addAtomtoHisto(self,putvalue,addvalue):
        #print distance
        #print atom[0]
        value = int(self.histobins*(float(putvalue) - float(self.histomin))/(float(self.histomax-self.histomin)))
        #print value
        if value<len(self.histo):
                self.histo[value]+=addvalue
        else:
                print "weird value"
                print value
                print distance
    
    def getBoxX(self,hbox):
        x = (float(hbox)*float(self.histomax-self.histomin) )/float(self.histobins) + float(self.histomin)
        return x


#function to assign histograms
def AssignHistograms(atomsclass,histogram,nucsize):
    counter=0
    for atom in atomsclass.atoms:
        counter+=1
        addvalue = 1.00/float(nucsize)
        histogram.addAtomtoHisto(atom,addvalue)
        #print atom[4]
    #print counter
    #print nucsize
    return histogram


#main function that is to be called
def MakeStructureHistogram(pathtype,manual=False,gzip=False):
    """
    special function to make histograms
    hardcoded. Remove at some point.
    """

    #read seedids in an array
    seedids = np.loadtxt(seedfileaddress,dtype=int,unpack=True)
    
    #set up histograms
    #for normal without seed
    nbcc = Histogram(histomin,histomax,histobins)
    nfcc = Histogram(histomin,histomax,histobins)
    nhcp = Histogram(histomin,histomax,histobins)
    nudf = Histogram(histomin,histomax,histobins)
    #for core without seed
    cbcc = Histogram(histomin,histomax,histobins)
    cfcc = Histogram(histomin,histomax,histobins)
    chcp = Histogram(histomin,histomax,histobins)
    cudf = Histogram(histomin,histomax,histobins)
    #for surface without seed
    sbcc = Histogram(histomin,histomax,histobins)
    sfcc = Histogram(histomin,histomax,histobins)
    shcp = Histogram(histomin,histomax,histobins)
    sudf = Histogram(histomin,histomax,histobins)

    #for normalizing histograms later
    norm = np.zeros(histobins)

    #read interfaces
    if manual==False:
        interfacelist = helpers.generate_intflist()
    else:
        interfacelist = helpers.read_intflist()
    
    #set up files for writing output

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
            #we are in the folder now

            histofile = os.path.join(pathpath,(identifier+'.histo.list'))
            histodataslices = []
            histodata = []
            count=0
	    #print "ll"
            #print histofile
	    if os.path.exists(histofile):
		#print histofile
		#print "kkll"
                for line in open(histofile,'r'):
                        histodata.append(line.strip())
                        count+=1
                        if count==12:
                                histodataslices.append(histodata)
                                histodata = []
                                count =0
            else:
                continue
            #loooping over each slice in the trajectory
            for i in range(len(histodataslices)):
		#print snapshots
                bccids = map(int,histodataslices[i][3].split())
                fccids = map(int,histodataslices[i][5].split())
                hcpids = map(int,histodataslices[i][7].split())
                udfids = map(int,histodataslices[i][9].split())
                surids = map(int,histodataslices[i][11].split())

                nucsize = len(bccids)+len(fccids)+len(hcpids)+len(udfids)
 
                #delete the seed particles from the lists
                bccids = [x for x in bccids if x not in seedids]
                fccids = [x for x in fccids if x not in seedids]
                hcpids = [x for x in hcpids if x not in seedids]
                udfids = [x for x in udfids if x not in seedids]

                #calculate the fractions now
                clustersize =  len(bccids)+len(fccids)+len(hcpids)+len(udfids)

                if clustersize>0:
                        fnbcc = float(len(bccids))/float(clustersize)
                        fnfcc = float(len(fccids))/float(clustersize)
                        fncp = float(len(hcpids))/float(clustersize)
                        fnudf = float(len(udfids))/float(clustersize)
                else:
                        fnbcc = 0
                        fnfcc = 0
                        fnhcp = 0
                        fnudf = 0

                #add to histo
                value = nbcc.addAtomtoHisto(nucsize,fnbcc)
                value = nfcc.addAtomtoHisto(nucsize,fnfcc)
                value = nhcp.addAtomtoHisto(nucsize,fnhcp)
                value = nudf.addAtomtoHisto(nucsize,fnudf)

                
                #now for core
                cbccids = [x for x in bccids if x not in surids]
                cfccids = [x for x in fccids if x not in surids]
                chcpids = [x for x in hcpids if x not in surids]
                cudfids = [x for x in udfids if x not in surids]
        
                #calculate the fractions now
                coresize =  len(cbccids)+len(cfccids)+len(chcpids)+len(cudfids)

                if coresize>0:
                        fcbcc = float(len(cbccids))/float(coresize)
                        fcfcc = float(len(cfccids))/float(coresize)
                        fchcp = float(len(chcpids))/float(coresize)
                        fcudf = float(len(cudfids))/float(coresize)
                else:
                        fcbcc = 0
                        fcfcc = 0
                        fchcp = 0
                        fcudf = 0

                #add to histo
                value = cbcc.addAtomtoHisto(nucsize,fcbcc)
                value = cfcc.addAtomtoHisto(nucsize,fcfcc)
                value = chcp.addAtomtoHisto(nucsize,fchcp)
                value = cudf.addAtomtoHisto(nucsize,fcudf)

                #for surface
                sbccids = [x for x in bccids if x in surids]
                sfccids = [x for x in fccids if x in surids]
                shcpids = [x for x in hcpids if x in surids]
                sudfids = [x for x in udfids if x in surids]

                #calculate the fractions now
                surfacesize =  len(sbccids)+len(sfccids)+len(shcpids)+len(sudfids)

                if surfacesize>0:
                        fsbcc = float(len(sbccids))/float(surfacesize)
                        fsfcc = float(len(sfccids))/float(surfacesize)
                        fshcp = float(len(shcpids))/float(surfacesize)
                        fsudf = float(len(sudfids))/float(surfacesize)
                else:
                        fsbcc = 0
                        fsfcc =0
                        fshcp = 0
                        fsudf = 0

		
                #add to histo
                value = sbcc.addAtomtoHisto(nucsize,fsbcc)
                value = sfcc.addAtomtoHisto(nucsize,fsfcc)
                value = shcp.addAtomtoHisto(nucsize,fshcp)
                value = sudf.addAtomtoHisto(nucsize,fsudf)

                if value < len(norm):
                        norm[value]+=1



    #normalise histograms
    for i in range(histobins):
        if norm[i]>0:
                nbcc[i]/=float(norm[i])
                nfcc[i]/=float(norm[i])
                nhcp[i]/=float(norm[i])
                nudf[i]/=float(norm[i])

                cbcc[i]/=float(norm[i])
                cfcc[i]/=float(norm[i])
                chcp[i]/=float(norm[i])
                cudf[i]/=float(norm[i])

                sbcc[i]/=float(norm[i])
                sfcc[i]/=float(norm[i])
                shcp[i]/=float(norm[i])
                sudf[i]/=float(norm[i])

    savefile1 = 'averaged_histo_structure_normal.dat'
    savefile2 = 'averaged_histo_structure_core.dat'
    savefile3 = 'averaged_histo_structure_surface.dat'

    #stack the histos
    nhisto = np.column_stack((nbcc.histox,nbcc.histo,nfcc.histo,nhcp.histo,nudf.histo))
    chisto = np.column_stack((cbcc.histox,cbcc.histo,cfcc.histo,chcp.histo,cudf.histo))
    shisto = np.column_stack((sbcc.histox,sbcc.histo,sfcc.histo,shcp.histo,sudf.histo))

    np.savetxt(savefile1,nhisto)
    np.savetxt(savefile2,chisto)
    np.savetxt(savefile3,shisto)


if __name__=='__main__':

    MakeStructureHistogram('AB',manual=False,gzip=True)


            

            
            
