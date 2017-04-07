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

def average_trajectory_storage(binary,pathtype,vulcan=False,jobs=50,pythonscript=None,manual=False,cores=1,gzip=False,queue='serial',extension='.opd.dat'):
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
    workdir = '/home/users/menonsqr/temp' 
    
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
            

            
def combine_averages(pathtype,bintype,manual=False,histomin=0,histomax=0,histobins=0):
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
    
    elif bintype=='structhisto':
        histox = np.linspace(histomin,histomax,histobins)
        surface_bcc = np.zeros(len(histox))
        surface_fcc = np.zeros(len(histox))
        surface_hcp = np.zeros(len(histox))
        surface_udf = np.zeros(len(histox))
        seed_bcc = np.zeros(len(histox))
        seed_fcc = np.zeros(len(histox))
        seed_hcp = np.zeros(len(histox))
        seed_udf = np.zeros(len(histox))
        count=0

    #struct histo doesnt need anything here
    
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

        elif bintype=='structhisto':

            #may the analysis start
            for path in pathlist:
                path = path.strip()
                #points to specific path folder
                pathpath = os.path.join(intfpath,path)
                #generate a random identifier
                indentifier = interface+path
                #generate a tempname
                filedummy = indentifier+'.histo.dat'
                filename = os.path.join(pathpath,filedummy)
		if os.path.exists(filename):
                	dist,surbcc,surfcc,surhcp,surudf,seebcc,seefcc,seehcp,seeudf = np.loadtxt(filename,unpack=True) 
                else:
			continue
                surface_bcc+=surbcc
                surface_fcc+=surfcc
                surface_hcp+=surhcp
                surface_udf+=surudf
                seed_bcc+=seebcc
                seed_fcc+=seefcc
                seed_hcp+=seehcp
                seed_udf+=seeudf
                count+=1

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

    elif bintype=='structhisto':
        
        surface_bcc/=float(count)
        surface_fcc/=float(count)
        surface_hcp/=float(count)
        surface_udf/=float(count)
        seed_bcc/=float(count)
        seed_fcc/=float(count)
        seed_hcp/=float(count)
        seed_udf/=float(count)

        X = np.column_stack((histox,surface_bcc,surface_fcc,surface_hcp,surface_udf,seed_bcc,seed_fcc,seed_hcp,seed_udf))
        np.savetxt("averaged_structure_histo.dat",X)


def average_trajectory_md(binary,folderlistfile,vulcan=False,jobs=50,pythonscript=None):
    """
    Does the same function as average_trajectory function, but for normal md paths (AB)

    folder names are read from a file.

    DO NOT USE
    """

    folderlist = []
    for line in open(folderlistfile,'r'):
        line.strip()
        folderlist.append(line)


    for folder in folderlist:
        pathpath = os.path.join(os.getcwd(),folder)
        pathpath = pathpath.strip()
        traj = os.path.join(pathpath,'traj.dat')
        #data = separate_traj(traj)
        if vulcan==False:
            qtrajs = calc_trajectory(binary,pathpath,tmpname,filename,writetofile=True)
            #add the filename to the list that has to read later
        else:
            #scriptpath,jobname,pythonname,argarray
            tmpname = os.path.join(pathpath,'traj.temp')
            filename = os.path.join(pathpath,'traj.opd.dat')
            argarray = [binary,traj,tmpname,filename]
            scriptname = os.path.join(pathpath,'subscript.job')
            os.system(("cp %s %s")%(pythonscript,pathpath))
            jobname = 'fd'+folder
            create_vulcan_script(scriptname,jobname,pythonscript,argarray)
            currentpath = os.getcwd()
            os.chdir(pathpath)
	    print pathpath
            runstatus=False
            while(runstatus==False):
                no_jobs = monitor_jobs()
                if no_jobs < jobs:
                    run_job(scriptname)
                    runstatus=True
                else:
                    time.sleep(60)
            
            os.chdir(currentpath)


def average_cluster(binary,pathtype,vulcan=False,jobs=50,pythonscript=None,manual=False,queue='serial'):
    """
    Run an binary on the selected type of paths and gather the output into text files.
    Outputs are now with an opd.dat extension. This can be changed to allow for custom names.

    If vulcan is set to True, it runs jobs on vulcan based on the total number of jobs allowed.

    If manual is set to True, interfaces are read from read_interfaces.txt from the sim directory.

    Used for calculating cluster props. Not completely done.

    """
    
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
            #datacmb = combine_paths_return(pathpath)
            #generate a random identifier
            indentifier = interface+path
            #generate a tempname
            tmpname = os.path.join(pathpath,indentifier)
            filedummy = indentifier+'.clu.dat'
            #a dummy file
            filename = os.path.join(pathpath,filedummy)
            #pass it on
            if vulcan==False:
                os.system(("python  %s %s 4393 %s %s %s")%(pythonscript,pathpath,filename,tmpname,binary))
            else:
                #scriptpath,jobname,pythonname,argarray
                argarray = [pathpath,'4393',filename,tmpname,binary]
                scriptname = os.path.join(pathpath,'subscript.job')
                os.system(("cp %s %s")%(pythonscript,pathpath))
                jobname = interface+path
                helpers.create_vulcan_script(scriptname,jobname,pythonscript,queue,argarray)
                currentpath = os.getcwd()
                os.chdir(pathpath)
                print pathpath
                runstatus=False
                while(runstatus==False):
                        no_jobs = monitor_jobs()
                        if no_jobs < jobs:
                                run_job(scriptname)
                                runstatus=True
                        else:
                        	time.sleep(60)
            
                os.chdir(currentpath)








