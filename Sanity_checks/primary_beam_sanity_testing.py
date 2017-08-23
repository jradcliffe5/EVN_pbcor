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

#### Inputs ####
AIPS.userno = 11
uvdata = AIPSUVData('HDFC0214','MULTI',1,1)
path = '/net/10.0.6.249/volume1/data/radcliff/EG078B/MSSC_PBCOR/Correct_model_multiple_pointings/PB_testing/'
################

for i in os.listdir('./flag_tables'):
    x = os.listdir('./flag_tables')
    uvdata.zap_table('FG',1)
    while i in x:
        x.remove(i)
    for j in x:
        uvflg = AIPSTask('UVFLG')
        uvflg.indata = uvdata
        uvflg.opcode = 'FLAG'
        uvflg.intext = '%sflag_tables/%s' % (path,j)
        uvflg.go()
    imagr = AIPSTask('IMAGR')
    imagr.indata = uvdata
    imagr.sources[1] = 'HDFC0214'
    imagr.docalib = 2
    imagr.cellsize[1:] = 0.001,0.001
    imagr.imsize[1:] = 512, 512
    imagr.nchav = 32
    imagr.outname = i
    imagr.niter = 300
    imagr.go()
