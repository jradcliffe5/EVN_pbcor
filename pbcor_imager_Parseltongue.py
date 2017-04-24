#!/usr/bin/env ParselTongue
import AIPS
from Wizardry.AIPSData import AIPSUVData as wizAIPSUVData
from Wizardry.AIPSData import AIPSTableRow
import sys, operator, pickle, os
import os, re, time, datetime, sys, math, fnmatch
from os.path import join, getsize
from datetime import date
from collections import deque
import Utilities
#from multiprocessing import Process	# No longer needed now SERPent is parallel
#from multiprocessing import Pool
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import math, time, datetime
from numpy import *
import itertools
from time import gmtime, strftime, localtime
ti = time.time()    # To time the script

## Inputs
UVFITSFILESpath = '/scratch/users/radcliff/EG078B_pbcor/large_fits/'
TASAVfilespath = '/net/10.0.6.249/volume1/data/radcliff/EG078B/MSSC_PBCOR/TASAV/'
IMAGEFITSFILESout = '/net/10.0.6.249/volume1/data/radcliff/EG078B/MSSC_PBCOR/IM/'
AIPS.userno = 10
disk = 1

def get_tab(uvdata, table):
	# find the number of tables of a certain type
	ver = 0
	for i in range(len(uvdata.tables)) :
		if table in uvdata.tables[i][1] :
			ver = uvdata.tables[i][0]
	print "HIGHEST TABLE OF TYPE", table, "is", ver
	return ver



for file in os.listdir(UVFITSFILESpath):
    if file.endswith('.UV'):
	for i in os.listdir(TASAVfilespath):
		if i[:8] == file[:8]:
		        fitld = AIPSTask('FITLD')
		        fitld.datain = UVFITSFILESpath+file
		        fitld.digicor = -1
		        fitld.outname = file[:8]
		        fitld.outclass = 'UVDATA'
		        fitld.go()
	 		fitld = AIPSTask('FITLD')
		        fitld.datain = TASAVfilespath+file[:8]+'_PBCOR_TASAV.fits'
		        fitld.digicor = -1
		        fitld.outname = file[:8]
		        fitld.outclass = 'TASAV'
		        fitld.go()
		        uvdata = AIPSUVData(file[:8],'UVDATA',disk,1)
			tasavfile = AIPSUVData(file[:8],'TASAV',disk,1)
			tacop = AIPSTask('TACOP')
			tacop.inext = 'CL'
			tacop.inver = 0
			tacop.ncount = 1
			tacop.outvers = 0
			tacop.indata = tasavfile
			tacop.outdata = uvdata
			tacop.go()
		        imagr = AIPSTask('IMAGR')
		        imagr.indata = uvdata
		        imagr.docalib=2
		        imagr.doband=0
		        imagr.gainuse = get_tab(uvdata,'CL')
		        imagr.sources[1:] = str(file[:8]),''
		        imagr.outname = file[:8]+'_PB'
		        imagr.nchav = 32
		        imagr.niter = 200
		        imagr.imsize[1:] = 2048, 2048
		        imagr.cellsize[1:] = 0.001,0.001
			imagr.uvwtfn = 'NA'
			imagr.robust = 5
		        imagr.go()
		        imagedata = AIPSImage(file[:8]+'_PB','ICL001',1,1)
		        fittp = AIPSTask('FITTP')
		        fittp.indata = imagedata
		        fittp.dataout = IMAGEFITSFILESout+file[:8]+'_NA_PBCOR_IM.fits'
		        fittp.go()
		        uvdata.zap()
		        imagedata.zap()
			tasavfile.zap()
		        AIPSImage(file[:8]+'_PB','IBM001',1,1).zap()
