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

histomax = 44.0
histomin = 0.0
histobins = 20
maxconfs=1000000000000000000000000000000000000000
tstcluster = 180
#create helpers class
helpers = tistools_helpers.tistools_helpers()
workdir = '/home/users/menonsqr/storage/20UC_TIS/tis_run'
seedfileaddress = '/home/users/menonsqr/SeedFCC19/seed.dat'

#read seedids
##################################################################################################
seedonlyids = np.loadtxt(seedfileaddress,dtype=int,unpack=True)
wfile = os.path.join(os.getcwd(),'seedinsurface.dat')
foutt = open(wfile,'w')
#class for histogram
class Histogram(object):
    
    def __init__(self,histomin,histomax,histobins):
        self.histomax = histomax
        self.histomin = histomin
        self.histobins = histobins
        self.histo = np.zeros(histobins)
        self.histox = np.linspace(self.histomin,self.histomax,self.histobins)
	self.count = np.zeros(histobins)

    def addAtomtoHisto(self,atom,addvalue):
        distance = atom[4]
	#print distance
	#print atom[0]
        value = int(self.histobins*(float(distance) - float(self.histomin))/(float(self.histomax-self.histomin)))
	#print value
	if value<len(self.histo):
        	self.histo[value]+=addvalue
		self.count[value]+=1
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
	#print nucsize
	#print addvalue
        histogram.addAtomtoHisto(atom,addvalue)
	#print atom[4]
    #print counter
    #print nucsize
    return histogram



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
	    #print mindist
	    #print (mindist<1e-5)
	    if mindist<1e-5:
		#print "oh"
		mindist=0.00
            atom[4]=mindist




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
    snapshots=1
    
    #set up histograms
    seed_dist = Histogram(histomin,histomax,histobins)
    count=np.zeros(histobins)
    dcount=np.zeros(dhistobins)


    if manual==False:
        interfacelist = helpers.generate_intflist()
    else:
        interfacelist = helpers.read_intflist()
    

    for interface in interfacelist:
	if snapshots>maxconfs:
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
            if snapshots>maxconfs:
		break
            path = path.strip()
            pathpath= os.path.join(intfpath,path)
            identifier = interface+path
            #we are in the folder now
            #we have to read the actual trajectory
            actualtraj = os.path.join(workdir,'tis','la',interface,path)
            data = helpers.combine_paths_return(actualtraj,gzip=gzip)
            #we have the data on standby
            #time to read the output raw data histo file.
            
            histofile = os.path.join(pathpath,(identifier+'histo.list'))
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

                ########################################################################
                ########################################################################
                ########################################################################
		nucids = bccids + fccids + hcpids + udfids
		seedinnuc = [x for x in seedonlyids if x in nucids]
		seedinsurf = [x for x in seedinnuc if x in surids]
		seedfrac = float(len(seedinsurf))/float(len(seedonlyids))
		foutt.write(("%d %f\n")%(nucsize,seedfrac))

		########################################################################
		########################################################################
	"""	
		#print fccids
		#print nucsize 
                #check if the guy should be part of histo,  and which histo
                if (nucsize <= tstcluster+3) and (nucsize >= tstcluster-3):
                        print "in"
			#he is part of tst cluster
                        #write the data down to a tempfile
                        #write the slice
                        outfile = open(tmpfile,'w')
                        for j in range(len(data[i])):
                                outfile.write(data[i][j])
                        outfile.flush()
                        outfile.close()

                        #read the slice
                        atoms = read_alles(tmpfile)
                        os.system(('rm %s')% tmpfile) 

                        #set up the seed classes
                        seed = Seed(seedfileaddress)
                        seed.ReadSeed()
                        seed.PopulateSeed(atoms,read=False)
			
                        #set up surface class
                        surface = Seed('dummy')
                        surface.ReadSeed(read=False,idlist=surids)
                        if surface.exists:
                                surface.PopulateSeed(atoms,read=False)
			

                        if surface.exists:
				#print "UDF"
                                seed.CalculateDistances(surface)
				#print "UDF"
                                seed_dist = AssignHistograms(seed,seed_dist,1.00)

			snapshots+=1
			print snapshots
			if snapshots>maxconfs:
				break
    
    #normalise the histograms
    #histogram x values
    histox_sur = np.zeros(len(seed_dist.histox))
    totvalue = np.sum(seed_dist.histo)

    for i in range(len(seed_dist.histox)):
        histox_sur[i] = seed_dist.getBoxX(i)
        seed_dist.histo[i] = float(seed_dist.histo[i])/float(totvalue)
    
    histo_sur = np.column_stack((histox_sur,seed_dist.histo))
    

    savefile1 = 'averaged_histo_cluster_'+str(tstcluster)+'.dat'

    np.savetxt(savefile1,histo_sur)

    #print snapshots
    #extdist = sum(extdist,[])
    #extdist2 = sum(extdist2,[])
    #print max(extdist)
    #print min(extdist)
    #print max(extdist2)
    #print min(extdist2)
    """
    ##################################################################
    #
    foutt.close()
#########################################################################
if __name__=='__main__':

    MakeStructureHistogram('AB',manual=True,gzip=True)




