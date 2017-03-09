import os
import sys
import subprocess as sub
import numpy as np
import tistools as tt

#calculates the binary value over a trajectory and writes to a file of choice
def calc_trajectory(binary,traj,tmpname,filename):
    datacmb = tt.separate_traj(traj)
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
        
        opfile.write(("%s")%(qd))
        opfile.flush()
        qtraj.append(qd)
        i+=1

    opfile.close()

if __name__=="__main__":

    binary = sys.argv[1]
    traj = sys.argv[2]
    tmpname = sys.argv[3]
    filename = sys.argv[4]

    calc_trajectory(binary,traj,tmpname,filename)

