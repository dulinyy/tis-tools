    
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


#workdir
workdir = '/home/users/menonsqr/storage/20UC_TIS/tis_run'
seedfileaddress = '/home/users/menonsqr/SeedFCC19/seed.dat'
tstcluster = 400
histomax = 60.0
histomin = 0.0
histobins = 6000

#create helpers class
helpers = tistools_helpers.tistools_helpers()


#class for seed
class Seed(object):
    
    def __init__(self,seedfileaddress):
        #self.atoms=np.empty([natoms,5])
        self.seedfile=seedfileaddress
        self.exists=True
    
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


#class for histogram
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


#function to assign histograms
def AssignHistograms(atomsclass,histogram):
    for atom in atomsclass.atoms:
        addvalue = 1.00
        histogram.addAtomtoHisto(atom,addvalue)
    return histogram


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


#main function that is to be called
def MakeStructureHistogram(pathtype,manual=False,gzip=False):
    """
    special function to make histograms
    hardcoded. Remove at some point.
    """
    tmpfile = 'my_tmp'
    #set up histograms
    bcc_sur = Histogram(histomin,histomax,histobins)
    fcc_sur = Histogram(histomin,histomax,histobins)
    hcp_sur = Histogram(histomin,histomax,histobins)
    udf_sur = Histogram(histomin,histomax,histobins)

    bcc_see = Histogram(histomin,histomax,histobins)
    fcc_see = Histogram(histomin,histomax,histobins)
    hcp_see = Histogram(histomin,histomax,histobins)
    udf_see = Histogram(histomin,histomax,histobins)
    
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
            indentifier = interface+path
            #we are in the folder now
            #we have to read the actual trajectory
            actualtraj = os.path.join(workdir,'tis','la',interface,path)
            data = helpers.combine_paths_return(actualtraj,gzip=gzip)
            #we have the data on standby
            #time to read the output raw data histo file.
            
            histofile = os.path.join(pathpath,(identifier+'histo.dat'))
            histodataslices = []
            histodata = []
            count=0
            if os.path.exists(histofile):
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
                bccids = map(int,histodataslices[i][3].split())
                fccids = map(int,histodataslices[i][5].split())
                hcpids = map(int,histodataslices[i][7].split())
                udfids = map(int,histodataslices[i][9].split())
                surids = map(int,histodataslices[i][11].split())
                nucsize = len(bccids)+len(fccids)+len(hcpids)+len(udfids) 
                #check if the guy should be part of histo,  and which histo
                if (nucsize <= tstcluster+3) and (nucsize >= tstcluster-3):
                        #he is part of tst cluster
                        #write the data down to a tempfile

                        #write the slice
                        outfile = open(tmpfile,'w')
                        for j in range(len(dataslices[i])):
                                outfile.write(dataslices[i][j])
                        outfile.flush()
                        outfile.close()

                        #read the slice
                        atoms = read_alles(tmpfile)
                        os.system(('rm %s')% tmpfile) 

                        #set up the seed classes
                        seed = Seed(seedfileaddress)
                        seed.ReadSeed()
                        seed.PopulateSeed(atoms,read=False)

                        #delete the seed particles from the lists
                        bccids = [x for x in bccids if x not in seed.seedids]
                        fccids = [x for x in bccids if x not in seed.seedids]
                        hcpids = [x for x in bccids if x not in seed.seedids]
                        udfids = [x for x in bccids if x not in seed.seedids]

                        #set up surface class
                        surface = Seed('dummy')
                        surface.ReadSeed(read=False,idlist=surids)
                        if surface.exists:
                                surface.PopulateSeed(atoms,read=False)

                        #set up BCC class
                        bcc = Seed('dummy')
                        bcc.ReadSeed(read=False,idlist=bccids)
                        if bcc.exists:
                                bcc.PopulateSeed(atoms,read=False)

                        #set up FCC class
                        fcc = Seed('dummy')
                        fcc.ReadSeed(read=False,idlist=fccids)
                        if fcc.exists:
                                fcc.PopulateSeed(atoms,read=False)

                        #set up HCP class
                        hcp = Seed('dummy')
                        hcp.ReadSeed(read=False,idlist=hcpids)
                        if hcp.exists:
                                hcp.PopulateSeed(atoms,read=False)

                        #set up UDF class
                        udf = Seed('dummy')
                        udf.ReadSeed(read=False,idlist=udfids)
                        if udf.exists:
                                udf.PopulateSeed(atoms,read=False)

                        if bcc.exists:
                                bcc.CalculateDistances(surface)
                                bcc_sur = AssignHistograms(bcc,bcc_sur)
                                bcc.CalculateDistances(seed)
                                bcc_see = AssignHistograms(bcc,bcc_see)
                        if fcc.exists:
                                fcc.CalculateDistances(surface)
                                fcc_sur = AssignHistograms(fcc,fcc_sur)
                                fcc.CalculateDistances(seed)
                                fcc_see = AssignHistograms(fcc,fcc_see)
                        if hcp.exists:
                                hcp.CalculateDistances(surface)
                                hcp_sur = AssignHistograms(hcp,hcp_sur)
                                hcp.CalculateDistances(seed)
                                hcp_see = AssignHistograms(hcp,hcp_see)
                        if udf.exists:
                                udf.CalculateDistances(surface)
                                udf_sur = AssignHistograms(udf,udf_sur)
                                udf.CalculateDistances(seed)
                                udf_see = AssignHistograms(udf,udf_see)

    #normalise the histograms
    for i in range(len(bcc_sur.histox)):
        sum_sur = bcc_sur[i]+fcc_sur[i]+hcp_sur[i]+udf_sur[i]
        sum_see = bcc_see[i]+fcc_see[i]+hcp_see[i]+udf_see[i]
        bcc_sur[i]/=sum_sur
        fcc_sur[i]/=sum_sur
        hcp_sur[i]/=sum_sur
        udf_sur[i]/=sum_sur
        bcc_see[i]/=sum_see
        fcc_sur[i]/=sum_see
        hcp_sur[i]/=sum_see
        udf_sur[i]/=sum_see

    histo_sur = np.column_stack((bcc_sur.histox,bcc_sur.histo,fcc_sur.histo,hcp_sur.histo,udf_sur.histo))
    histo_see = np.column_stack((bcc_see.histox,bcc_see.histo,fcc_see.histo,hcp_see.histo,udf_see.histo))

    np.savetxt('averaged_histo_surface',histo_sur)
    np.savetxt('averaged_histo_seed',histo_see)


if __name__=='__main__':

        MakeStructureHistogram('AB',manual=True,gzip=True):


            

            
            