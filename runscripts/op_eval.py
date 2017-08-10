"""
Parallelised op evaluation script. To be used with tis tools.
"""

import os
import sys
import time
import shutil
import subprocess as sub
import multiprocessing as mp
import copy_reg
import types

import tistools_helpers.tistools_helpers as tistools_helpers

#create helpers class
helpers = tistools_helpers.tistools_helpers()


#stuff for multiprocessing
def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

def return_result(result):
    return result

copy_reg.pickle(types.MethodType, _pickle_method)


def getQTrajectoryParallel(data,instance,binary):
        
        qtraj = []
        
        for j in range(len(data)):
                
                tmpfile = os.path.join(os.getcwd(),"my_tmp_conf" + str(instance))
                
                outfile = open(tmpfile, "w")
                
                for i in range(len(data[j])):
                        outfile.write(data[j][i])
                outfile.flush()
                outfile.close()
                
                cmd = []
                cmd.append(binary)
                cmd.append(tmpfile)
                
                proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
                out,err = proc.communicate(input="")
                proc.wait()
                qd = out          #JR print
                
                qtraj.append(qd)

        os.remove(tmpfile)

        return instance,qtraj


def getFullTrajectoryParallel(filename,outfilename,binary,vulcan,gzip):
        
        #pick number of cores
        try:
                cores = int(os.environ['NSLOTS'])
        except:
                cores = 8
        
        if vulcan == False:
                cores = mp.cpu_count() - 1
        
        
        results = []
        kernel_list = [ 0 for x in range(cores)]
        
        data = []
        
        if os.path.exists(filename):
                datasliced = helpers.separate_traj(filename,gzip)
        
        else:
                print "file not found"


        if len(datasliced) > 0:
                
                per_core = len(datasliced)/cores
                #print per_core
                #JR: in case nslices < ncore
                if (per_core == 0):
                        kernel_list = [ 0 for x in range(len(datasliced))]

                for j in range(len(kernel_list)):
                        kernel_list[j] = per_core
                for j in range(len(datasliced)%cores):
                        kernel_list[j] += 1
                
                #kernel_list contains no. of slices for each core
                #assigen start and end for the loop on each core in kernel_slices
                start = 0
                end = 0
                kernel_slices = []
                
                for j in range(len(kernel_list)):
                        end += kernel_list[j]
                        kernel_slices.append([start,end])
                        start += kernel_list[j]
                #print kernel_slices
                #JR distribute data per core
                data_per_core = [ [] for x in range(len(kernel_list))]
                
                for k in range(len(kernel_list)):
                        start = kernel_slices[k][0]
                        end = kernel_slices[k][1]

                        for j in range(start,end):
                                data_per_core[k].append(datasliced[j])

                #now launch process
                pool = mp.Pool(processes=len(kernel_list))
                
                for x in range(len(kernel_list)):
                        results.append(pool.apply_async(getQTrajectoryParallel, args=(data_per_core[x],x,binary,)))
                
                pool.close()
                pool.join()
                
                output = [p.get() for p in results]
                
                #sort the results
                output.sort()
                
                #trim the unnecessary variable
                output = [out[1] for out in output]
                
                #join output
                output = sum(output,[])

                fout = open(outfilename,'w')

                for i in range(len(output)):
                        fout.write(("%s")%(str(output[i])))

                fout.close()


if __name__=="__main__":
        
        binary = sys.argv[1]
        filename = sys.argv[2]
        outfilename = sys.argv[3]
        vulcan = sys.argv[4]
        gzip = sys.argv[5]

        if gzip=='False':
                gzip=False
        elif gzip=='True':
                gzip=True

        if vulcan=='False':
                vulcan=False
        elif vulcan=='True':
                vulcan=True

        getFullTrajectoryParallel(filename,outfilename,binary,vulcan,gzip)



