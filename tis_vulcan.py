#!/usr/bin/env python

'''
Created May 28, 2013

@author: Jutta Rogal
'''

binary = '/home/users/menonsqr/bin/struct'

import os
import sys
import subprocess as sub

def getQTrajectory(data,tmpfile):
	nlines = len(data)
	natoms = int(data[3])
	nblock = natoms+9
	nslices = nlines/nblock

	if nlines%nblock != 0:
		print 'Error in getQTrajectory'
		print 'nlines%nblock != 0'
		print 'Exiting program...'
		sys.exit(1)

	of = open("opdata.dat","w")
	qtraj=[]
	print '#No. of slices = ',nslices
	for j in range(nslices):
		start = j*nblock
		end = (j+1)*nblock
		outfile = open(tmpfile,'w')
		for i in range(start,end):
			outfile.write(data[i])
		outfile.flush()
		outfile.close()
		cmd = []
		cmd.append(binary)
		cmd.append(tmpfile)
		proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
		out,err = proc.communicate(input="")
		proc.wait() 
		#q = float(out)
		qtraj.append(out)

	os.remove(tmpfile)
	return qtraj




if __name__ == '__main__':
	
	
	if len(sys.argv) < 2:
		print "Call with one argument, argument is the input file"
		sys.exit(1)

	infile = open(sys.argv[1],'r')
	tmpfile = sys.argv[1]+'.tmp'
	data = []
	
	for line in infile:
		data.append(line)

	if len(data) > 4:
		nlines = len(data)
		natoms = int(data[3])
		nblock = natoms+9
		nslices = nlines/nblock
		if nlines%nblock != 0:
			print 'trajectory not fully written in workdir, removing unfinished slice...'
			newlines = nblock*nslices
			data_tmp = []
			for i in range(newlines):
				data_tmp.append(data[i])

			data = []
			data = data_tmp

			nlines = len(data)
			natoms = int(data[3])
			nblock = natoms+9
			nslices = nlines/nblock

			if nlines%nblock != 0:
				print 'Cannot trajectory file with proper number of lines'
				print 'Exiting program...'
				sys.exit(1)
	
	if len(data) > 0:
		qtraj = []
		qtraj = getQTrajectory(data,binary,tmpfile)
		outfile = sys.argv[1]+'.opd'
		fout = open(outfile,'w')
		for i in range(len(qtraj)):
			fout.write(qtraj[i])
			fout.write('\n')
		fout.close()


