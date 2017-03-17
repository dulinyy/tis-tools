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
import gzip
import shutil

#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("analysis.log")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False 


class tistools_helpers(object):
    def __init__(self):
        logger.info("module helpers created")

#---------------------------------------------------------------------------------------------------------
#                        LEVEL ONE HELPERS : STANDALONE FUNCTIONS
#---------------------------------------------------------------------------------------------------------
    #function to read first and last line
    def read_fline(self,file):
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
    def extract_opval(self,first,last):
        """
        Function which reads the order parameter value from the line of a trajectory.00.0000XXX.dat file
        Argument : first and last lines of the above file
        Return   : the order parameter minimum and maximum value.
        """
        rawf = first.split()
        rawl = last.split()
        return rawf[1],rawl[1]

    def readgzlines(self,fname):
        """
        Function to read zip files (.gz)
        argument :filename
        """
        f = sub.Popen(['zcat', fname], stdout=sub.PIPE)
        for line in f.stdout:
                yield line

    def converttogz(self,fname):
        """
        function to zip files
        """
        gzipfile = fname + '.gz'
	f_out = gzip.open(gzipfile, 'wb')
        with open(fname, 'rb') as f_in:
            shutil.copyfileobj(f_in, f_out)
	f_out.close()
        os.remove(fname)

    #to check the type of the path
    def check_type(self,sA,sB,sstateA,sstateB):
        """
        Function which checks the type of path, if it is Ab,BA,AA,BB or UN(which is error paths)
        Argument  : first two values are the beginning and end of the path from traj file, next two values are the order parameter limits.
        Return    : Path type string 
        """
        #AB path condition
        sA = float(sA)
        sAB = float(sB)
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
    def generate_pathlist(self,start,stop):
        """
        Helper function to generate a list of paths
        Argument  : fstart and stop values of path numbers
        Return    : Path no list
        """
        pathno_list = []
        for i in range(start,stop+1):
            pathno_list.append(str(("%07d")%i))
        #okay lets read the stable states   
        return pathno_list


    #lets read the stable states
    def read_op(self):
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
    def generate_intflist(self):
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

    #function to read interface list
    #should be provided in read_interfaces.dat
    def read_intflist(self):
        interfacelist=[]
        for line in open('read_interfaces.dat','r'):
            interfacelist.append(line)
        return interfacelist


    def separate_traj(self,traj,gzip=False,writetofile=False):
        """
        Seperate the trajectory and return in a data file
        data array consits of sub array of slices
        zip support added
        """
        data = []
        datasliced = []
        if gzip==False:
            infile = open(traj,'r')
            for line in infile:
                data.append(line)
        elif gzip==True:
            for line in self.readgzlines(traj):
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

    #make a sumbmission script for vulcan
    def create_vulcan_script(self,scriptpath,jobname,pythonname,cores,queue,argarray):
        """
        Creates a submission script for vulcan
        multi core support added
        """
        queue = queue + '.q'
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
        of.write("#$ -pe smp " + str(cores) + "\n")
        of.write("#$ -P ams.p\n")
        of.write("source $HOME/.bashrc\n")
        of.write("module load numpy/1.7.0b2\n")
        of.write("module load intel/016.0.047\n")
        of.write("hostname\n")
        #modify from here

        of.write(("python %s ")% (pythonname))
        for i in range(len(argarray)):
            of.write(("%s ")% argarray[i])
        of.write("\n")

        of.close()

    #code to run qsub on vulcan
    def run_job(self,scriptpath):
        """
        code to run qsub on vulcan
        """
        os.system(("qsub %s")% scriptpath)

    #monitor jobs on vulcan
    def monitor_jobs(self):
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
    def test_binary(self,binary):
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

#----------------------------------------------------------------------------------------------------------------
#                                       LEVEL TWO HELPERS
#----------------------------------------------------------------------------------------------------------------
    def combine_paths_return(self,pathno,gzip=False):
        """
        Combines forward and backward part in right order
        and returns the combined data string.
        zip support added
        """
        print pathno
        fwdpath = os.path.join(pathno,"forward","traj.dat")
        bkdpath = os.path.join(pathno,"backward","traj.dat")
        datafwd = self.separate_traj(fwdpath,gzip=gzip)
        databkd = self.separate_traj(bkdpath,gzip=gzip)
        datacmb = []

        firstslice = True
        for data in databkd[::-1]:
            datacmb.append(data)
        for data in datafwd[1:]:
            datacmb.append(data)
        return datacmb


    def combine_paths_write(self,pathno,writefile,gzip=False):
        """
        Combine the backward and forward part of the trajectory
        remember to delete the combined paths after use!
        might take a lot of memory.
        """
        print pathno
        outfile = open(writefile,'w')
        datacmb = self.combine_paths_return(pathno,gzip=gzip)
        for i in range(len(datacmb)):
            for j in range(len(datacmb[i])):
                outfile.write(datacmb[i][j])

        outfile.close()


    def zip_all_paths(self,start,stop,int_list):
        """
        Zip all traj files
        """
        #int_list = self.generate_intflist()
        pathno_list = self.generate_pathlist(start,stop)
        for interface in int_list:
            interface = interface.strip()
            for path in pathno_list:
                path = path.strip()
                fwd_traj = os.path.join(os.getcwd(),'tis/la',interface,path,'forward','traj.dat')
                bkd_traj = os.path.join(os.getcwd(),'tis/la',interface,path,'backward','traj.dat')
                if os.path.exists(fwd_traj):
                        self.converttogz(fwd_traj)
                else:
                        logger.error(('%s not found')%fwd_traj)
                if os.path.exists(bkd_traj):
                        self.converttogz(bkd_traj)
                else:
                        logger.error(('%s not found')%bkd_traj)
            logger.info(("%s completed.")%interface)




