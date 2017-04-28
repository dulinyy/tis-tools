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
workdir = '/home/users/menonsqr/NS19/tis_run'
seedfileaddress = '/home/users/menonsqr/SeedFCC43/seed.dat'
binary = 'orderparameter/main'
tstcluster = 200

maxconfs=10
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
                        self.atoms=np.empty([len(self.seedids),7])
                else:
                        self.exists=False
        else:
                self.seedids=idlist
                if len(self.seedids)>0:
                        self.atoms=np.empty([len(self.seedids),7])
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
                self.atoms[k][4]=atoms[i][4]
                self.atoms[k][5]=atoms[i][5]
                self.atoms[k][6]=atoms[i][6]
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
        #atoms values are as follows
        # 0 : id
        # 1,2,3 : x,y,z
        # 4 : whichever distance value
        # 5,6 : avg q4 and q6 respectively
        # #total of seven parameters
        atoms = np.empty([natoms,7])
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
                atoms[i][5] = 99999.00
                atoms[i][6] = 99999.00
                #atoms[i][4] = False
                i+=1
            count+=1

        idn, aq4,aq6 = np.loadtxt('result.dat',unpack=True)
        idn = idn.astype(int)

        for i in range(len(idn)):
                for j in range(len(atoms)):
                        if idn[i]==atoms[j][0]:
                                atoms[j][5] = aq4[i]
                                atoms[j][6] = aq6[i]
                                break

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
    distance1 = []
    distance2 = []
    distance3 = []
    distance4 = []
    distance5 = []
    distance6 = []
    
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
		#print fccids
		#print nucsize 
                #check if the guy should be part of histo,  and which histo
                if (nucsize <= tstcluster+3) and (nucsize >= tstcluster-3):
                        outfile = open(tmpfile,'w')
                        for j in range(len(data[i])):
                                outfile.write(data[i][j])
                        outfile.flush()
                        outfile.close()


                        #apply order parameter and read histo stuff
                        cmd = [binary,tmpfile]
                        proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
                        out,err = proc.communicate()
                        proc.wait()

                        #read the slice
                        #modify read alles to read in the q4 q6 too.- done
                        atoms = read_alles(tmpfile)
                        os.system(('rm %s')% tmpfile) 

                        #set up the seed classes
                        seed = Seed(seedfileaddress)
                        seed.ReadSeed()
                        seed.PopulateSeed(atoms,read=False)

                        #delete the seed particles from the lists
                        bccids = [x for x in bccids if x not in seed.seedids]
                        fccids = [x for x in fccids if x not in seed.seedids]
                        hcpids = [x for x in hcpids if x not in seed.seedids]
                        udfids = [x for x in udfids if x not in seed.seedids]
                        
                        #set up surface class
                        surface = Seed('dummy')
                        surface.ReadSeed(read=False,idlist=surids)
                        if surface.exists:
                                surface.PopulateSeed(atoms,read=False)

                        #find udf ids in surface
                        udfsurids = [x for x in udfids if x in surids]
                        udfcoreids = [x for x in udfids if x not in surids]

                        #set up UDF class
                        udfsur = Seed('dummy')
                        udfsur.ReadSeed(read=False,idlist=udfsurids)
                        if udfsur.exists:
                                udfsur.PopulateSeed(atoms,read=False)


                        udfcore = Seed('dummy')
                        udfcore.ReadSeed(read=False,idlist=udfcoreids)
                        if udfcore.exists:
                                udfcore.PopulateSeed(atoms,read=False)
                        
                        #seeds are populated. Now find distance of each atom to the surface.
                        udfcore.CalculateDistances(surface)

                        #now add the points to the arrays.
                        for atomcito in udfcore.atoms:
                                if atomcito[4]<=1.0:
                                        distance1.append([atomcito[5],atomcito[6]])
                                elif atomcito[4]<=2.0:
                                        distance2.append([atomcito[5],atomcito[6]])
                                elif atomcito[4]<=3.0:
                                        distance3.append([atomcito[5],atomcito[6]])                                
                                elif atomcito[4]<=4.0:
                                        distance4.append([atomcito[5],atomcito[6]])
                                elif atomcito[4]<=5.0:
                                        distance5.append([atomcito[5],atomcito[6]])
                                elif atomcito[4]<=6.0:
                                        distance6.append([atomcito[5],atomcito[6]])
                                else:
                                        print "jsjsjsj"

    #write out the files
    fout = open('distance1.dat','w')
    for i in range(len(distance1)):
        fout.write(("%f %f\n")%(distance1[i][0],distance1[i][1]))
    fout.close()

    fout = open('distance2.dat','w')
    for i in range(len(distance2)):
        fout.write(("%f %f\n")%(distance2[i][0],distance2[i][1]))
    fout.close()

    fout = open('distance3.dat','w')
    for i in range(len(distance3)):
        fout.write(("%f %f\n")%(distance3[i][0],distance3[i][1]))
    fout.close()

    fout = open('distance4.dat','w')
    for i in range(len(distance4)):
        fout.write(("%f %f\n")%(distance4[i][0],distance4[i][1]))
    fout.close()

    fout = open('distance5.dat','w')
    for i in range(len(distance5)):
        fout.write(("%f %f\n")%(distance5[i][0],distance5[i][1]))
    fout.close()

    fout = open('distance1.dat','w')
    for i in range(len(distance6)):
        fout.write(("%f %f\n")%(distance6[i][0],distance6[i][1]))
    fout.close()


if __name__=='__main__':
    MakeStructureHistogram('AB',manual=False,gzip=True)
                        
                