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


#class for histogram
class Histogram(object):
    
    def __init__(self,histomin,histomax,histobins):
        self.histomax = histomax
        self.histomin = histomin
        self.histobins = histobins
        self.histo = np.zeros(histobins)
        self.histox = np.zeros(histobins)
	self.count = np.zeros(histobins)
    
    def addAtomtoHisto(self,atom,addvalue):
        distance = atom[4]
	#print distance
	#print atom[0]
        value = int(self.histobins*(float(distance) - float(self.histomin))/(float(self.histomax-self.histomin)))
	#print value
	if value<len(self.histo):
        	self.histo[value]+=addvalue
		self.count[value]+=1
	else:
		print "weird value"
		print value
		print distance
    
    def getBoxX(self,hbox):
        x = (float(hbox)*float(self.histomax-self.histomin) )/float(self.histobins) + float(self.histomin)
        return x
