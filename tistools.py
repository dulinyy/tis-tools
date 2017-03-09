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


#function to read first and last line
def read_fline(file):
    """
    Function which reads the first and last line of a file
    Argument : filename of the file to be read
    Return   : the first and last line in that order.
    """
    with open(file, "rb") as f:
        first = f.readline()      # Read the first line.
        f.seek(-2, 2)             # Jump to the second last byte.
        while f.read(1) != b"\n": # Until EOL is found...
            f.seek(-2, 1)         # ...jump back the read byte plus one more.
        last = f.readline()
    return first,last

#extract order parameter value
def extract_opval(first,last):
    """
    Function which reads the order parameter value from the line of a trajectory.00.0000XXX.dat file
    Argument : first and last lines of the above file
    Return   : the order parameter minimum and maximum value.
    """
    rawf = first.split()
    rawl = last.split()
    return rawf[1],rawl[1]

#to check the type of the path
def check_type(sA,sB,sstateA,sstateB):
    """
    Function which checks the type of path, if it is Ab,BA,AA,BB or UN(which is error paths)
    Argument  : first two values are the beginning and end of the path from traj file, next two values are the order parameter limits.
    Return    : Path type string 
    """
    #AB path condition
    sA = float(sA)
    sB = float(sB)
    sstateA = float(sstateA)
    sstateB = float(sstateB)
    if ((sA<=sstateA) and (sB>=sstateB)):
        path_type = "AB"
    elif ((sA>=sstateB) and (sB<=sstateA)):
        path_type = "BA"
    elif ((sA<=sstateA) and (sB<=sstateA)):
        path_type = "AA"
    elif ((sA>=sstateB) and (sB>=sstateB)):
        path_type = "BB"
    else:
        path_type = "UN"
    #print path_type
    return path_type


#generate a list of path numbers first
def generate_pathlist(start,stop):
    """
    """
    pathno_list = []
    for i in range(start,stop+1):
        pathno_list.append(str(("%07d")%i))
    #okay lets read the stable states   
    return pathno_list


#lets read the stable states
def read_op():
    im_here = os.getcwd()
    for line in open(os.path.join(im_here,"options","orderparameters.txt")):
        raw = line.split()
        if (raw[0]=="custom"):
            sstateA = float(raw[3])
            sstateB = float(raw[4])
            break
    return sstateA,sstateB



#okay. The stable states are in.
#time to read the interfaces list
def generate_intflist():
    int_type = []
    int_val = []
    int_list = []
    for line in open(os.path.join(os.getcwd(),"options","interfaces.txt")):
        raw = line.split()
        int_type.append(raw[0])
        int_val.append(("%02d")%int(raw[1]))

    for i in range(len(int_type)):
        int_list.append(int_type[i]+str(int_val[i]))

    return int_list


def separate_traj(traj,writetofile=False):
    """
    Seperate the trajectory and return in a data file
    data array consits of sub array of slices
    """
    infile = open(traj,'r')
    data = []
    datasliced = []
    for line in infile:
        data.append(line)

    nlines = len(data)
    natoms = int(data[3])
    nblock = natoms+9
    nslices = nlines/nblock

    for j in range(nslices):
        start = j*nblock
        end = (j+1)*nblock
        dummy = []
        for i in range(start,end):
            dummy.append(data[i])
        datasliced.append(dummy)

    if writetofile==True:
        filenamelist = []
        for i in range(len(datasliced)):
            splitfile = os.path.join(traj,str(i))
            fout = open(splitfile,'w')
            fout.write(datasliced[i])
            fout.close()
        return datasliced,filenamelist
    else:
        return datasliced


def combine_paths_return_md(pathno):
    """
    Combines forward and backward part in right order
    and returns the combined data string.
    """
    print pathno
    fwdpath = os.path.join(pathno,"forward","traj.dat")
    bkdpath = os.path.join(pathno,"backward","traj.dat")
    datafwd = separate_traj(fwdpath)
    databkd = separate_traj(bkdpath)
    datacmb = []

    firstslice = True
    for data in databkd[::-1]:
        datacmb.append(data)
    for data in datafwd[1:]:
        datacmb.append(data)
    return datacmb


def combine_paths_return(pathno):
    """
    Combines forward and backward part in right order
    and returns the combined data string.
    """
    print pathno
    fwdpath = os.path.join(pathno,"forward","traj.dat")
    bkdpath = os.path.join(pathno,"backward","traj.dat")
    datafwd = separate_traj(fwdpath)
    databkd = separate_traj(bkdpath)
    datacmb = []

    firstslice = True
    for data in databkd[::-1]:
        datacmb.append(data)
    for data in datafwd[1:]:
        datacmb.append(data)
    return datacmb


def combine_paths_write(pathno,writefile):
    """
    Combine the backward and forward part of the trajectory
    remember to delete the combined paths after use!
    might take a lot of memory.
    """
    print pathno
    outfile = open(writefile,'w')
    datacmb = combine_paths_return(pathno)
    for i in range(len(datacmb)):
        for j in range(len(datacmb[i])):
            outfile.write(datacmb[i][j])

    outfile.close()


#calculates the binary value over a trajectory and writes to a file of choice
def calc_trajectory(binary,traj,tmpname,filename,writetofile=True):
    """
    calculates the binary value over a trajectory and writes to a file of choice
    """
    datacmb = combine_paths_return(traj)
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



#make a sumbmission script for vulcan
def create_vulcan_script(scriptpath,jobname,pythonname,argarray):
    """
    Creates a submission script for vulcan
    """

    of = open(scriptpath,"w")
    of.write("#!/bin/bash\n")
    of.write(("#$ -N %s\n")% (jobname))
    of.write("#$ -S /bin/bash\n")
    of.write("#$ -r n\n")
    of.write("#$ -cwd\n")
    of.write("#$ -l h_rt=05:59:00\n")
    of.write("#$ -l qname=serial.q\n")
    of.write("#$ -j y\n")
    of.write("#$ -R y\n")
    of.write("source $HOME/.bashrc\n")
    of.write("module load numpy/1.7.0b2\n")
    of.write("hostname\n")
    #modify from here

    of.write(("python %s ")% (pythonname))
    for i in range(len(argarray)):
        of.write(("%s ")% argarray[i])
    of.write("\n")

    of.close()

#code to run qsub on vulcan
def run_job(scriptpath):
    """
    code to run qsub on vulcan
    """
    os.system(("qsub %s")% scriptpath)

#monitor jobs on vulcan
def monitor_jobs():
    """
    monitor jobs on vulcan
    """
    os.system("qstat > qstat.dat")
    i=0
    for line in open("qstat.dat",'r'):
        i+=1
    no_jobs = i-2
    return no_jobs

#test binary
def test_binary(binary):
    """
    code to test the kind of outputs a binary can produce. Not used for now.
    """
    pathtotest = os.path.join(os.getcwd(),"tis","standardfiles","conf.dump")
    cmd = []
    cmd.append(binary)
    cmd.append(pathtotest)
    proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
    out,err = proc.communicate(input="")
    proc.wait()
    qd = out.split()
    return len(qd)

#make path type lists
#
def find_paths(start,stop):
    """
    Goes through every path between start and stop variable and assigns them as AA,AB,BA,BB or UN type.
    UN paths are defective paths.
    """
    sstateA,sstateB = read_op()
    int_list = generate_intflist()
    pathno_list = generate_pathlist(start,stop)


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
            sA,sB = read_fline(read_file)
	    sA,sB = extract_opval(sA,sB)
            path_type = check_type(sA,sB,sstateA,sstateB)
            if path_type=='AB':
                fAB.write(("%s\n")% path)
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


#average the properties of the path. Both vulcan and non vulcan here
#interfacelist = list of interfaces for which to be calculated
##trial one did not work. Changing to text file based approach

def average_trajectory(binary,pathtype,vulcan=False,jobs=50,pythonscript=None,manual=False):
    """
    Run an binary on the selected type of paths and gather the output into text files.
    Outputs are now with an opd.dat extension. This can be changed to allow for custom names.

    If vulcan is set to True, it runs jobs on vulcan based on the total number of jobs allowed.

    If manual is set to True, interfaces are read from read_interfaces.txt from the sim directory.

    """
    
    if manual==False:
        interfacelist = generate_intflist()
    else:
        interfacelist = []
        for line in open('read_interfaces.txt','r'):
            interfacelist.append(line)

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
            datacmb = combine_paths_return(pathpath)
            #generate a random identifier
            indentifier = interface+path
            #generate a tempname
            tmpname = os.path.join(pathpath,indentifier)
            filedummy = indentifier+'.opd.dat'
            #a dummy file
            filename = os.path.join(pathpath,filedummy)
            #pass it on
            if vulcan==False:
                qtrajs = calc_trajectory(binary,pathpath,tmpname,filename,writetofile=True)
                #add the filename to the list that has to read later
            else:
                #scriptpath,jobname,pythonname,argarray

                argarray = [binary,pathpath,tmpname,filename]
                scriptname = os.path.join(pathpath,'subscript.job')
                os.system(("cp %s %s")%(pythonscript,pathpath))
	  	jobname = interface+str(path)
                create_vulcan_script(scriptname,jobname,pythonscript,argarray)
                currentpath = os.getcwd()
                os.chdir(pathpath)

                
                
                runstatus=False
                while(runstatus==False):
                    no_jobs = monitor_jobs()
                    if no_jobs < jobs:
                        run_job(scriptname)
                        runstatus=True
                    else:
                        time.sleep(60)
                os.chdir(currentpath)

            
def combine_averages(pathtype,bintype,manual=True):
    """
    Combine the outut of the previous function from different files into an averaged data

    This is hard coded.
    """
    #read stable states
    sstateA,sstateB = read_op()
    
    #read interfaces
    if manual==False:
        interfacelist = generate_intflist()
    else:
        interfacelist = []
        for line in open('read_interfaces.txt','r'):
            interfacelist.append(line)
    
    if bintype=='struct':
            
        oparray = np.array(range(int(sstateA),int(sstateB+1)))
        bccavg = np.zeros(len(oparray))
        fccavg = np.zeros(len(oparray))
        hcpavg = np.zeros(len(oparray))
        udfavg = np.zeros(len(oparray))
        count = np.zeros(len(oparray))

    elif bintype=='cluster':

        oparray = np.array(range(int(sstateA),int(sstateB+1)))
        bccavg = np.zeros(len(oparray))
        fccavg = np.zeros(len(oparray))
        hcpavg = np.zeros(len(oparray))
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
                op,bcc,fcc,hcp,udf = np.loadtxt(filename,unpack=True)
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
		filedummy2 = indentifier+'.opd.dat'
                #a dummy file
                filename = os.path.join(pathpath,filedummy)
		filename2 = os.path.join(pathpath,filedummy2)
                #next hardcoded part
                bcc,fcc,hcp = np.loadtxt(filename,unpack=True)
		op,a,b,c,f = np.loadtxt(filename2,unpack=True)
                op = op.astype(int)

                for i in range(len(op)):
                    for j in range(len(oparray)):
                            if op[i]==oparray[j]:
                                    bccavg[j]+=bcc[i]
                                    fccavg[j]+=fcc[i]
                                    hcpavg[j]+=hcp[i]
                                    count[j]+=1
        
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
                    bccavg[i]/=count[i]
                    fccavg[i]/=count[i]
                    hcpavg[i]/=count[i]

        X = np.column_stack((oparray,bccavg,fccavg,hcpavg))
        np.savetxt("averaged_data_cluster.dat",X)

            



def average_trajectory_md(binary,folderlistfile,vulcan=False,jobs=50,pythonscript=None):
    """
    Does the same function as average_trajectory function, but for normal md paths (AB)

    folder names are read from a file.
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


def average_cluster(binary,pathtype,vulcan=False,jobs=50,pythonscript=None,manual=False):
    """
    Run an binary on the selected type of paths and gather the output into text files.
    Outputs are now with an opd.dat extension. This can be changed to allow for custom names.

    If vulcan is set to True, it runs jobs on vulcan based on the total number of jobs allowed.

    If manual is set to True, interfaces are read from read_interfaces.txt from the sim directory.

    """
    
    if manual==False:
        interfacelist = generate_intflist()
    else:
        interfacelist = []
        for line in open('read_interfaces.txt','r'):
            interfacelist.append(line)

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








