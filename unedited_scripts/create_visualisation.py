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
workdir = '/home/users/menonsqr/storage/HCP19/tis_run'
seedfileaddress = '/home/users/menonsqr/SeedFCC19/seed.dat'


#create helpers class
helpers = tistools_helpers.tistools_helpers()


def VisualiseTrajectory(interface,pathno):
	seedids = np.loadtxt(seedfileaddress,dtype=int,unpack=True)
	intfpath = os.path.join(os.getcwd(),"tis","la",interface)
	pathpath= os.path.join(intfpath,pathno)
	identifier = interface+pathno
	histofile = os.path.join(pathpath,(identifier+'.histo.list'))
        actualtraj = os.path.join(workdir,'tis','la',interface,path)
        data = helpers.combine_paths_return(actualtraj,gzip=gzip)

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

        for i in range(len(histodataslices)):
		#print snapshots
                bccids = map(int,histodataslices[i][3].split())
                fccids = map(int,histodataslices[i][5].split())
                hcpids = map(int,histodataslices[i][7].split())
                udfids = map(int,histodataslices[i][9].split())
                surids = map(int,histodataslices[i][11].split())

                bccids = [x for x in bccids if x not in seedids]
                fccids = [x for x in fccids if x not in seedids]
                hcpids = [x for x in hcpids if x not in seedids]
                udfids = [x for x in udfids if x not in seedids]  
                clusterids = bccids + fccids + hcpids + udfids
                coreids = [x for x in clusterids if x not in surids ]

                itervar = 0
                for line in data[i]:
			if (itervar==8):
				line = line + " structure cluster seed core\n"

			elif (itervar>8):
				sline = line.split()
				ident = int(sline[0])
				
				if ident in bccids:
					structure = " 1"
					cluster = " 1"

				elif ident in fccids:
					structure = " 2"
					cluster = " 1"

				elif ident in hcpids:
					structure = " 3"
					cluster = " 1"

				elif ident in udfids:
					structure = " 4"
					cluster = " 1"
					
				else:
					structure = " 0"
					cluster = " 0"

				if ident in seedids:
					seed = " 1"
				else:
					seed = " 0"

				if ident in coreids:
					core = " 1"
				else:
					core = " 0"

				line = line + structure +cluster+seed+core
			itervar+=1

	#write out the file
	outfile = os.path.join(pathpath,'traj.mod.dat')
	fout = open(outfile,'w')
	for datacito in data:
		fout.write(datacito)
	fout.close()

if __name__=="__main__":
	interface = 'f13'
	pathno = ''
	VisualiseTrajectory(interface,pathno)
