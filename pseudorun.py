    
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

#create helpers class
helpers = tistools_helpers.tistools_helpers()

"""
Combine the outut of the previous function from different files into an averaged data

This is hard coded.
 """

def combine_averages(pathtype,bintype,manual=False):
    """
    Combine the outut of the previous function from different files into an averaged data

    This is hard coded.
    """
    #read stable states
    sstateA,sstateB = helpers.read_op()
    
    #read interfaces
    if manual==False:
        interfacelist = helpers.generate_intflist()
    else:
        interfacelist = helpers.read_intflist()
    
    if bintype=='struct':
            
        oparray = np.array(range(int(sstateA),int(sstateB+1)))
        bccavg = np.zeros(len(oparray))
        fccavg = np.zeros(len(oparray))
        hcpavg = np.zeros(len(oparray))
        udfavg = np.zeros(len(oparray))
        count = np.zeros(len(oparray))

    elif bintype=='cluster':

        oparray = np.array(range(int(sstateA),int(sstateB+1)))
        min_dist = np.zeros(len(oparray))
        numberofatoms = np.zeros(len(oparray))
        percent = np.zeros(len(oparray))
        surfacenumberofatoms = np.zeros(len(oparray))
        surfacepercent = np.zeros(len(oparray))
        seedinsurface = np.zeros(len(oparray))
        seedinsurfacep = np.zeros(len(oparray))
        count = np.zeros(len(oparray))
    
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

        #setup the tables and columns
        ##this part is hardcoded now.
        #has to be set up individually everytime
        #module for structure dstributions
        if bintype=='struct':
            #may the analysis start
            for path in pathlist:
                path = path.strip()
                #points to specific path folder
                pathpath = os.path.join(intfpath,path)
                #generate a random identifier
                indentifier = interface+path
                #generate a tempname
                filedummy = indentifier+'.opd.dat'
                #a dummy file
                filename = os.path.join(pathpath,filedummy)

                #next hardcoded part
		if os.path.exists(filename):
                	op,bcc,fcc,hcp,udf = np.loadtxt(filename,unpack=True)
		else:
			continue
                op = op.astype(int)
                #print fcc
                for i in range(len(op)):
                    for j in range(len(oparray)):
                            if op[i]==oparray[j]:
                                    #print op[i]
                                    bccavg[j]+=bcc[i]
                                    #print bccavg[j]
                                    fccavg[j]+=fcc[i]
                                    hcpavg[j]+=hcp[i]
                                    udfavg[j]+=udf[i]
                                    count[j]+=1


        elif bintype=='cluster':

            #may the analysis start
            for path in pathlist:
                path = path.strip()
                #points to specific path folder
                pathpath = os.path.join(intfpath,path)
                #generate a random identifier
                indentifier = interface+path
                #generate a tempname
                filedummy = indentifier+'.clu.dat'
                filename = os.path.join(pathpath,filedummy)
                op,smin_dist,snumberofatoms,spercent,ssurfacenumberofatoms,ssurfacepercent,sseedinsurface,sseedinsurfacep = np.loadtxt(filename,unpack=True,comments='#') 
                op = op.astype(int)
                snumberofatoms = snumberofatoms.astype(int)
                ssurfacenumberofatoms = ssurfacenumberofatoms.astype(int)
                sseedinsurface = sseedinsurface.astype(int)

                for i in range(len(op)):
                    for j in range(len(oparray)):
                            if op[i]==oparray[j]:
                                min_dist[j] += smin_dist[i]
                                numberofatoms[j] += snumberofatoms[i]
                                percent[j] += spercent[i]
                                surfacenumberofatoms[j] += ssurfacenumberofatoms[i]
                                surfacepercent[j] += ssurfacepercent[i]
                                seedinsurface[j] += sseedinsurface[i]
                                seedinsurfacep[j] += sseedinsurfacep[i]
                                count += 1

        
    if bintype=='struct':

        for i in range(len(oparray)):
                if count[i]!=0:
                    bccavg[i]/=float(count[i])
                    fccavg[i]/=float(count[i])
                    hcpavg[i]/=float(count[i])
                    udfavg[i]/=float(count[i])
            #print bccavg
        X = np.column_stack((oparray,bccavg,fccavg,hcpavg,udfavg))
        np.savetxt("averaged_data_struct.dat",X)
            
    elif bintype=='cluster':
        
        for i in range(len(oparray)):
                if count[i]!=0:
                    min_dist[i]/=count[i]
                    numberofatoms[i]/=count[i]
                    percent[i]/=count[i]
                    surfacenumberofatoms[i]/=count[i]
                    surfacepercent[i]/=count[i]
                    seedinsurface[i]/=count[i]
                    seedinsurfacep[i]/=count[i]

        X = np.column_stack((oparray,min_dist,numberofatoms,percent,surfacenumberofatoms,surfacepercent,seedinsurface,seedinsurfacep))
        np.savetxt("averaged_data_cluster.dat",X)

if __name__=='__main__':
        
        pathtype = 'AB'
        bintype = 'struct'
        combine_averages(pathtype,bintype,manual=False)
