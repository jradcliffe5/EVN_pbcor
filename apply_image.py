from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
from code_mailer import headless
import numpy as np
import sys
import os

def get_tab(uvdata, table):
	# find the number of tables of a certain type
	ver = 0
	for i in range(len(uvdata.tables)) :
		if table in uvdata.tables[i][1] :
			ver = uvdata.tables[i][0]
	print "HIGHEST TABLE OF TYPE", table, "is", ver
	return ver

inputs = headless('inputs.txt')

AIPS.userno = int(inputs['AIPSnumber'])

x = np.load('aips_info.npy')
uvdata = AIPSUVData(str(x[0]),str(x[1]),int(x[2]),int(x[3]))

clcal = AIPSTask('CLCAL')
clcal.indata = uvdata
clcal.interpol = 'SELN'
clcal.snver = get_tab(uvdata,'SN')
clcal.invers = get_tab(uvdata,'SN')
clcal.gainver = get_tab(uvdata,'CL')
clcal.gainuse = get_tab(uvdata,'CL') + 1
clcal.go()

tasav = AIPSTask('TASAV')
tasav.indata = uvdata
tasav.outname = uvdata.name + 'PB'
tasav.go()
imagr = AIPSTask('IMAGR')
imagr.indata = uvdata
imagr.docalib=2
imagr.doband=0
imagr.flagver = 0
imagr.gainuse = get_tab(uvdata,'CL')
imagr.sources[1:] = str(uvdata.name),''
imagr.outname = uvdata.name+'PB'
imagr.nchav = 32
imagr.niter = 1500
imagr.imsize[1:] = 4096,4096
imagr.cellsize[1:] = 0.0001,0.0001
imagr.go()
imagr.uvwtfn = 'NA'
imagr.cellsize[1:]=0.001,0.001
imagr.go()

imagedata = AIPSImage(uvdata.name+'PB','ICL001',1,1)
imagedata2 = AIPSImage(uvdata.name+'PB','ICL001',1,2)
fittp = AIPSTask('FITTP')
fittp.indata = imagedata
fittp.dataout = 'PWD:%s_PBCOR_IM.fits' % uvdata.name
fittp.go()
fittp.indata = imagedata2
fittp.dataout = 'PWD:%s_PBCOR_NA_IM.fits' % uvdata.name
fittp.go()
fittp.indata = AIPSUVData(uvdata.name+'PB','TASAV',1,1)
fittp.dataout = 'PWD:%s_PBCOR_TASAV.fits' % uvdata.name
fittp.go()
AIPSCat().zap()
