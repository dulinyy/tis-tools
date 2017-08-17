"""
Sub file used for op calculation; to be used for tis paths
Have to overwrite to make it nicer
"""

import os
import sys
import subprocess as sub
import numpy as np
import tistools as tt
from  tistools_helpers import tistools_helpers
import logging

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


#calculates the binary value over a trajectory and writes to a file of choice
def calc_trajectory(binary,traj,tmpname,filename):
    datacmb = helpers.separate_traj(traj)
    opfile = open(filename,'w')
    #tmpfilelist = []
    qtraj = []
    for i in range(len(datacmb)):
        tmpfile = tmpname+str(i)+".dat"
        outfile = open(tmpfile,'w') 
        for j in range(len(datacmb[i])):
            outfile.write(datacmb[i][j])
        outfile.flush()
        outfile.close()
        #tmpfilelist.append(tmpfile)
        
    #i = 0
    #for tmp in tmpfilelist: 
        cmd = []
        cmd.append(binary)
        cmd.append(tmpfile)
        proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
        out,err = proc.communicate(input="")
        proc.wait()
        qd = out
        os.system(("rm %s")% tmpfile)
        
        opfile.write(("%s")%(qd))
        opfile.flush()
        qtraj.append(qd)
        #i+=1

    opfile.close()

if __name__=="__main__":

    binary = sys.argv[1]
    traj = sys.argv[2]
    tmpname = sys.argv[3]
    filename = sys.argv[4]
    
    calc_trajectory(binary,traj,tmpname,filename)

