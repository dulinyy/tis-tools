    
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
seedfileaddress = '/home/users/menonsqr/SeedFCC27/seed.dat'
orgfiles = '/home/users/menonsqr/storage/NS27/tis_run'
#create helpers class
helpers = tistools_helpers.tistools_helpers()



#main function that is to be called
def MakeStructureHistogram(pathtype,manual=False,gzip=False):
    """
    special function to make histograms
    hardcoded. Remove at some point.
    """

    seedids = np.loadtxt(seedfileaddress,dtype=int,unpack=True)


    #read interfaces
    if manual==False:
        interfacelist = helpers.generate_intflist()
    else:
        interfacelist = helpers.read_intflist()
    
    #set up files for writing output

    for interface in interfacelist:

        interface = interface.strip()
        intfpath = os.path.join(os.getcwd(),"tis","la",interface)
	orgintfpath = os.path.join(orgfiles,"tis","la",interface)
        intfpath = intfpath.strip()
	orgintfpath = orgintfpath.strip()
        pathpath = os.path.join(intfpath,pathtype+".dat")
	orgpathpath = os.path.join(orgintfpath,pathtype+".dat")
        pathpath = pathpath.strip()
	orgpathpath = orgpathpath.strip()
        pathlist = []
        filenamelist = []

        #we get the list of all paths that needs to be analysed
        for path in open(orgpathpath,'r'):
            pathlist.append(path)

        #may the analysis start
        for path in pathlist:
	    nopath=False
            path = path.strip()
            pathpath= os.path.join(intfpath,path)
            identifier = interface+path
            #we are in the folder now
#	    nwritefile = os.path.join(pathpath,(identifier+'.opd.normal'))
#            cwritefile = os.path.join(pathpath,(identifier+'.opd.core'))
#	    swritefile = os.path.join(pathpath,(identifier+'.opd.surface'))
            histofile = os.path.join(pathpath,(identifier+'histo.list'))
#	    nfout = open(nwritefile,'w')
#	    cfout = open(cwritefile,'w')
#	    sfout = open(swritefile,'w')	


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
		nopath=True
                continue

            nwritefile = os.path.join(pathpath,(identifier+'.opd.normal'))
            cwritefile = os.path.join(pathpath,(identifier+'.opd.core'))
            swritefile = os.path.join(pathpath,(identifier+'.opd.surface'))
#            histofile = os.path.join(pathpath,(identifier+'.histo.list'))
            nfout = open(nwritefile,'w')
            cfout = open(cwritefile,'w')
            sfout = open(swritefile,'w')

            #loooping over each slice in the trajectory
            for i in range(len(histodataslices)):
		#print snapshots
                bccids = map(int,histodataslices[i][3].split())
                fccids = map(int,histodataslices[i][5].split())
                hcpids = map(int,histodataslices[i][7].split())
                udfids = map(int,histodataslices[i][9].split())
                surids = map(int,histodataslices[i][11].split())

                nucsize = len(bccids)+len(fccids)+len(hcpids)+len(udfids)
		#print nucsize 
                #delete the seed particles from the lists
                bccids = [x for x in bccids if x not in seedids]
                fccids = [x for x in fccids if x not in seedids]
                hcpids = [x for x in hcpids if x not in seedids]
                udfids = [x for x in udfids if x not in seedids]

                #calculate the fractions now
                clustersize =  len(bccids)+len(fccids)+len(hcpids)+len(udfids)
		#print clustersize
                if clustersize>0:
                        fnbcc = float(len(bccids))/float(clustersize)
                        fnfcc = float(len(fccids))/float(clustersize)
                        fnhcp = float(len(hcpids))/float(clustersize)
                        fnudf = float(len(udfids))/float(clustersize)
                else:
                        fnbcc = 0
                        fnfcc = 0
                        fnhcp = 0
                        fnudf = 0

                #add to hist
		#print str(nucsize)
		#print str(fnbcc)
		#print str(fnfcc)
		#print str(fnhcp)
		#print str(fnudf)
		nfout.write(("%d %.4f %.4f %.4f %.4f\n")%(nucsize,fnbcc,fnfcc,fnhcp,fnudf))
                
                #now for core
                cbccids = [x for x in bccids if x not in surids]
                cfccids = [x for x in fccids if x not in surids]
                chcpids = [x for x in hcpids if x not in surids]
                cudfids = [x for x in udfids if x not in surids]
        
                #calculate the fractions now
                coresize =  len(cbccids)+len(cfccids)+len(chcpids)+len(cudfids)
		#print coresize
		#print cudfids
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
		cfout.write(("%d %.4f %.4f %.4f %.4f\n")%(nucsize,fcbcc,fcfcc,fchcp,fcudf))

                #for surface
                sbccids = [x for x in bccids if x in surids]
                sfccids = [x for x in fccids if x in surids]
                shcpids = [x for x in hcpids if x in surids]
                sudfids = [x for x in udfids if x in surids]

                #calculate the fractions now
                surfacesize =  len(sbccids)+len(sfccids)+len(shcpids)+len(sudfids)
		#print surfacesize
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
		sfout.write(("%d %.4f %.4f %.4f %.4f\n")%(nucsize,fsbcc,fsfcc,fshcp,fsudf))
	    nfout.close()
	    sfout.close()
	    cfout.close()




if __name__=='__main__':

    MakeStructureHistogram('AB',manual=True,gzip=True)


            

            
            
