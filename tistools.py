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
import tistools_helpers.tistools_helpers as tistools_helpers
import tistools_helpers.atoms as atomsclass
import tistools_helpers.histogram as histogramclass
import multiprocessing as mp

#SRM:set up logger for the general error messages,
logger = logging.getLogger(__name__)
handler = logging.FileHandler("analysis.log")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False 

#create helpers class
helpers = tistools_helpers.tistools_helpers()

#function to zip all paths
def zip_paths(start,stop,manual=False):
    if manual==False:
        interfacelist = helpers.generate_intflist()
    else:
        interfacelist = helpers.read_intflist()
    helpers.zip_all_paths(start,stop,interfacelist)

#calculates the binary value over a trajectory and writes to a file of choice
def calc_trajectory(binary,traj,tmpname,filename,gzip=False,writetofile=True):
    """
    calculates the binary value over a trajectory and writes to a file of choice
    """
    datacmb = helpers.combine_paths_return(traj,gzip=gzip)
    opfile = open(filename,'w')
    tmpfilelist = []
    qtraj = []
    for i in range(len(datacmb)):
        tmpfile = tmpname+str(i)+".dat"
        outfile = open(tmpfile,'w') 
        for j in range(len(datacmb[i])):
            outfile.write(datacmb[i][j])
        outfile.flush()
        outfile.close()
        tmpfilelist.append(tmpfile)
        
    i = 0
    for tmp in tmpfilelist: 
        cmd = []
        cmd.append(binary)
        cmd.append(tmp)
        #logger.info("command created")
        #logger.info(cmd)
        proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
        out,err = proc.communicate(input="")
        proc.wait()
        qd = out
        os.system(("rm %s")% tmp)
        if writetofile==True:
            opfile.write(("%s")%(qd))
            opfile.flush()
        qtraj.append(qd)
        i+=1

    opfile.close()
    return qtraj




def find_paths(start,stop):
    """
    Goes through every path between start and stop variable and assigns them as AA,AB,BA,BB or UN type.
    UN paths are defective paths.
    """
    logger.info("starting to sort paths")
    sstateA,sstateB = helpers.read_op()
    int_list = helpers.generate_intflist()
    pathno_list = helpers.generate_pathlist(start,stop)
    ABcount=0

    for interface in int_list:
        abfile_path = os.path.join(os.getcwd(),'tis/la',interface,'AB.dat')
        aafile_path = os.path.join(os.getcwd(),'tis/la',interface,'AA.dat')
        bbfile_path = os.path.join(os.getcwd(),'tis/la',interface,'BB.dat')
        bafile_path = os.path.join(os.getcwd(),'tis/la',interface,'BA.dat')

        unfile_path = os.path.join(os.getcwd(),'tis/la',interface,'UN.dat')

        fAB = open(abfile_path,'w')
        fAA = open(aafile_path,'w')
        fBB = open(bbfile_path,'w')
        fBA = open(bafile_path,'w')
        fUN = open(unfile_path,'w')


        for path in pathno_list:
            file_name = 'trajectory.00.'+path+'.dat'
            read_file = os.path.join(os.getcwd(),'tis/la',interface,path,file_name)
            sA,sB = helpers.read_fline(read_file)
	    sA,sB = helpers.extract_opval(sA,sB)
            path_type = helpers.check_type(sA,sB,sstateA,sstateB)
            if path_type=='AB':
                fAB.write(("%s\n")% path)
                ABcount+=1
            elif path_type=='BA':
                fBA.write(("%s\n")% path)
            elif path_type=='AA':
                fAA.write(("%s\n")% path)
            elif path_type=='BB':
                fBB.write(("%s\n")% path)
            else:
                fUN.write(("%s\n")% path)

        fAB.close()
        fBA.close()
        fAA.close()
        fBB.close()
        fUN.close()
    logger.info("path sorting completed")
    logger.info(("%d AB paths found.")%ABcount)

#average the properties of the path. Both vulcan and non vulcan here
#interfacelist = list of interfaces for which to be calculated
##trial one did not work. Changing to text file based approach

def average_trajectory(binary,pathtype,vulcan=False,jobs=50,pythonscript=None,manual=False,cores=1,gzip=False,queue='serial',extension='.opd.dat'):
    """
    Run an binary on the selected type of paths and gather the output into text files.
    Outputs are now with an opd.dat extension. This can be changed to allow for custom names.

    If vulcan is set to True, it runs jobs on vulcan based on the total number of jobs allowed.

    If manual is set to True, interfaces are read from read_interfaces.txt from the sim directory.

    Option to be added later.
    Arglist option:
    should have all the arguments required for the corresponding python script
    should be a dictionary

    """
    logger.info(('binary selected: %s')%binary)
    logger.info(('pathtype selected: %s')%pathtype)
    logger.info(('vulcan set to: %s')%str(vulcan))
    logger.info(('pythonscript selected: %s')%pythonscript)
    logger.info(('manual reading of interfaces set to: %s')%str(manual))
    logger.info(('zipped file support set to: %s')%str(gzip))
    logger.info(('cores selected: %d')%cores)

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
            #points to specific path folder
            pathpath = os.path.join(intfpath,path)
            #combine the paths and return
            #datacmb = combine_paths_return(pathpath,gzip)
            #generate a random identifier
            indentifier = interface+path
            #generate a tempname
            tmpname = os.path.join(pathpath,indentifier)
            filedummy = indentifier+extension
            #a dummy file
            filename = os.path.join(pathpath,filedummy)
            #pass it on
            if vulcan==False:
                qtrajs = calc_trajectory(binary,pathpath,tmpname,filename,gzip=gzip,writetofile=True)
                #add the filename to the list that has to read later
            else:
                #scriptpath,jobname,pythonname,argarray

                argarray = [binary,pathpath,tmpname,filename,gzip]
                scriptname = os.path.join(pathpath,'subscript.job')
                os.system(("cp %s %s")%(pythonscript,pathpath))
	  	jobname = interface+str(path)
                helpers.create_vulcan_script(scriptname,jobname,pythonscript,cores,queue,argarray)
                currentpath = os.getcwd()
                os.chdir(pathpath)
                logger.info(('deploying job in %s')%pathpath)
                
                
                runstatus=False
                while(runstatus==False):
                    no_jobs = helpers.monitor_jobs()
                    if no_jobs < jobs:
                        helpers.run_job(scriptname)
                        runstatus=True
                    else:
                        time.sleep(60)
                logger.info('job deployed')

                os.chdir(currentpath)

def average_trajectory_storage(binary,pathtype,vulcan=False,jobs=50,pythonscript=None,manual=False,cores=1,gzip=False,queue='serial',extension='.opd.dat',folder='FCC19'):
    """
    Run an binary on the selected type of paths and gather the output into text files.
    Outputs are now with an opd.dat extension. This can be changed to allow for custom names.

    If vulcan is set to True, it runs jobs on vulcan based on the total number of jobs allowed.

    If manual is set to True, interfaces are read from read_interfaces.txt from the sim directory.

    Option to be added later.
    Arglist option:
    should have all the arguments required for the corresponding python script
    should be a dictionary

    """

    #hardcoded part
    workdir = folder 
    pathrun=0
    if not os.path.exists(workdir):
	os.mkdir(workdir)

    logger.info(('binary selected: %s')%binary)
    logger.info(('pathtype selected: %s')%pathtype)
    logger.info(('vulcan set to: %s')%str(vulcan))
    logger.info(('pythonscript selected: %s')%pythonscript)
    logger.info(('manual reading of interfaces set to: %s')%str(manual))
    logger.info(('zipped file support set to: %s')%str(gzip))
    logger.info(('cores selected: %d')%cores)

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

        #make a folder for the interface in the workdir
        workinfpath = os.path.join(workdir,interface)
	os.mkdir(workinfpath)

        #we get the list of all paths that needs to be analysed
        for path in open(pathpath,'r'):
            pathlist.append(path)

        #may the analysis start
        for path in pathlist:
            path = path.strip()
            #points to specific path folder
            #copy the path now
            #helpers.copy_path(pathpath,workinfpath)
	    pathpath= os.path.join(intfpath,path)
	    os.system(('cp -rf %s %s')%(pathpath,workinfpath))
            #reset pathpath
            pathpath=os.path.join(workinfpath,path)
            #combine the paths and return
            #datacmb = combine_paths_return(pathpath,gzip)
            #generate a random identifier
            indentifier = interface+path


            #now copy the file
            #make a dir inside main dir
            
            #generate a tempname
            tmpname = os.path.join(pathpath,indentifier)
            filedummy = indentifier+extension
            #a dummy file
            filename = os.path.join(pathpath,filedummy)
            #pass it on
            if vulcan==False:
                qtrajs = calc_trajectory(binary,pathpath,tmpname,filename,gzip=gzip,writetofile=True)
                #add the filename to the list that has to read later
            else:
                #scriptpath,jobname,pythonname,argarray

                argarray = [binary,pathpath,tmpname,filename,gzip]
                scriptname = os.path.join(pathpath,'subscript.job')
                os.system(("cp %s %s")%(pythonscript,pathpath))
                jobname = interface+str(path)
                helpers.create_vulcan_script(scriptname,jobname,pythonscript,cores,queue,argarray)
                currentpath = os.getcwd()
                os.chdir(pathpath)
                logger.info(('deploying job in %s')%pathpath)
                
                
                #runstatus=False
                #while(runstatus==False):
                #    no_jobs = helpers.monitor_jobs()
                #    if no_jobs < jobs:
                helpers.run_job(scriptname)
                #        runstatus=True
                #    else:
                #        time.sleep(60)
                logger.info('job deployed')
		pathrun+=1
		if pathrun>=300:
			raise SystemExit()
                os.chdir(currentpath)
            


def combine_averages(pathtype,columns=4,extension='.opd.dat',manual=False):
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
    
    oparray = np.array(range(int(sstateA),int(sstateB+1)))

    count = np.zeros(len(oparray))
    read_cols = [np.zeros(len(oparray)) for x in range(columns)]
    
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
        #may the analysis start
        for path in pathlist:
                path = path.strip()
                #points to specific path folder
                pathpath = os.path.join(intfpath,path)
                #generate a random identifier
                indentifier = interface+path
                #generate a tempname
                filedummy = indentifier+extension
                #a dummy file
                filename = os.path.join(pathpath,filedummy)

                #next hardcoded part
                if os.path.isfile(filename) and os.path.getsize(filename) > 0:
		#if os.path.exists(filename):
                        ops = np.loadtxt(filename,unpack=True)
                else:
                        continue
                
                #assuming first one is order parameter
                #ops[0] = ops[0].astype(int)
                
                #print fcc
                for i in range(len(ops[0])):
                    for j in range(len(oparray)):
                            if ops[0][i]==oparray[j]:
                                    for k in range(len(read_cols)):
                                        read_cols[k][j]+=ops[k+1][i]
                                    count[j]+=1

    #implement different methods here

    for i in range(len(oparray)):
        if count[i]!=0:
                for k in range(len(read_cols)):
                        read_cols[k][i]/=float(count[i])   
        
    #write out data
    fout = open('averaged_data.dat','w')
    for i in range(len(oparray)):
        fout.write(("%d ")%oparray[i])
        for k in range(len(read_cols)):
                fout.write((" %f")%read_cols[k][i])
        fout.write('\n')
    fout.close()



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



"""
Create a function to create histograms
"""
#create histograms
def setup_histogram(nhistograms,histomin,histomax,histobins):

        histolist = []
        for i in range(nhistograms):
                hist = histogramclass.Histogram(histomin,histomax,histobins)
                histolist.append(hist)


        for histo in histolist:
                for i in range(len(histo.histox)):
                        histo.histox[i] = histo.getBoxX(i)

        return histolist



#function to wrap histograms
def create_histogram(file,histograms,binary,addvalue=1.):
        
        atoms = read_alles(file)
                                
        #apply the op
        cmd = [binary,file]
        proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
        out,err = proc.communicate(input="")
        proc.wait()

        #read the histofile created
        surfacefileaddress = "surface.dat" 
        surface = atomsclass.Seed(surfacefileaddress)
        surface.ReadSeed()
        if surface.exists:
                surface.PopulateSeed(atoms,read=False)

        #first create the surface class
        bccfileaddress = "bccid.dat" 
        bcc = atomsclass.Seed(bccfileaddress)
        bcc.ReadSeed()
        if bcc.exists:
                bcc.PopulateSeed(atoms,read=False)
                bcc.CalculateDistances(surface)
                for atomito in bcc.atoms:
                        histograms[0].addAtomtoHisto(atomito,addvalue)

        #first create the fcc class
        fccfileaddress = "fccid.dat" 
        fcc = atomsclass.Seed(fccfileaddress)
        fcc.ReadSeed()
        if fcc.exists:
                fcc.PopulateSeed(atoms,read=False)
                fcc.CalculateDistances(surface)
                for atomito in fcc.atoms:
                        histograms[1].addAtomtoHisto(atomito,addvalue)
                                
        #first create the hcp class
        hcpfileaddress = "hcpid.dat" 
        hcp = atomsclass.Seed(hcpfileaddress)
        hcp.ReadSeed()
        if hcp.exists:
                hcp.PopulateSeed(atoms,read=False)
                hcp.CalculateDistances(surface)
                for atomito in hcp.atoms:
                        histograms[2].addAtomtoHisto(atomito,addvalue)

        #first create the udf class
        udffileaddress = "udfid.dat" 
        udf = atomsclass.Seed(udffileaddress)
        udf.ReadSeed()
        if udf.exists:
                udf.PopulateSeed(atoms,read=False)
                udf.CalculateDistances(surface)
                for atomito in udf.atoms:
                        histograms[3].addAtomtoHisto(atomito,addvalue)


        return histograms



#normalise and write histograms
def normalise_histograms(histograms,normalise='perbin',outputfile='histo.dat'):
        if normalise=='perbin':
                #normalise values per bin per different histograms
                for i in range(len(histograms[0].histo)):
                        binsum = 0
                        for histo in histograms:
                                binsum += histo.histo[i]
                        for histo in histograms:
                                if binsum!=0:
                                        histo.histo[i]/=float(binsum)

        elif normalise=='peraverage':
                #normalise with the total value
                for i in range(len(histograms[0].histo)):
                        for histo in histograms:
                                if histo.count[i]!=0:
                                        histo.histo[i]/=float(histo.count[i])

        elif normalise=='perhisto':
                for histo in histograms:
                        histosum = np.sum(histo.histo)
                        for i in range(len(histo.histo)):
                                if histosum!=0:
                                        histo.histo[i]/=float(histosum)

        #print out the data
        #hardcoded at this point
        histo_output = np.column_stack((histograms[0].histox,histograms[0].histo,histograms[1].histo,histograms[2].histo,histograms[3].histo))
        np.savetxt(outputfile,histo_output)


     

def eval_op_parallel(pythonscript,filename,outfilename,binary,queue='serial',cores=8,gzip=True,vulcan=True):
        """
        Op evaluation code for running on vulcan
        """

        if vulcan==True:
                scriptpath = os.path.join(os.getcwd(),'subscript.job')
                jobname = 'opeval'
                argarray = [binary,filename,outfilename,vulcan,gzip]
                helpers.create_vulcan_script(scriptpath,jobname,pythonscript,cores,queue,argarray)
                helpers.run_job(scriptpath)

        else:
                os.system(("python %s %s %s %s %s %s")%(pythonscript,binary,filename,outfilename,str(vulcan),str(gzip)))














