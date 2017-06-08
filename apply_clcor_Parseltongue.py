#!/usr/bin/env ParselTongue
import AIPS
from Wizardry.AIPSData import AIPSUVData as wizAIPSUVData
from Wizardry.AIPSData import AIPSTableRow
from AIPSTask import AIPSTask
import sys, operator, pickle, os

## Inputs
UVFITSFILESpath = '/net/10.0.6.249/volume1/data/radcliff/EG078B/MSSC_PBCOR/Scripts/FITS/UV/'
IMAGEFITSFILESout = '/net/10.0.6.249/volume1/data/radcliff/EG078B/MSSC_PBCOR/Scripts/FITS/IM/'
UVFITSFILESout = '/net/10.0.6.249/volume1/data/radcliff/EG078B/MSSC_PBCOR/Scripts/FITS/TASAV/'
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

CLCOR_params = pickle.load(open("CLCOR_params.pckl", "rb" )) ## pickle file made by primary beam sim

text = open('CLCOR_params.txt','w')
text.close()

print os.listdir(UVFITSFILESpath)
print CLCOR_params
for file in os.listdir(UVFITSFILESpath):
    if file.endswith('.UV'):
        for i in range(len(CLCOR_params)):
            if CLCOR_params[i][0] == file[:8]:
                fitld = AIPSTask('FITLD')
                fitld.datain = UVFITSFILESpath+file
                fitld.digicor = -1
                fitld.outname = file[:8]
                fitld.outclass = 'UVDATA'
                fitld.go()
                uvdata = AIPSUVData(file[:8],'UVDATA',disk,1)
                for j in range(len(CLCOR_params[i][3])):
                    clcor = AIPSTask('CLCOR')
                    clcor.indata = uvdata
                    clcor.sources[1] ='*'
                    clcor.antennas[1] = j+1
                    clcor.opcode = 'GAIN'
                    clcor.clcorprm[1] = float(CLCOR_params[i][3][j])
                    clcor.gainver = get_tab(uvdata,'CL')
                    if j == 0:
                        clcor.gainuse = get_tab(uvdata,'CL') + 1
                    else:
                        clcor.gainuse = get_tab(uvdata,'CL')
                    clcor.go()
		f = open('CLCOR_params.txt','a')
		f.write(str(CLCOR_params[i][0])+' '+str(CLCOR_params[i][3])+'\n')
		f.close()
		tasav = AIPSTask('TASAV')
		tasav.indata = uvdata
		tasav.outname = file[:8] + 'PB'
		tasav.go()
                imagr = AIPSTask('IMAGR')
                imagr.indata = uvdata
                imagr.docalib=2
                imagr.doband=0
                imagr.gainuse = get_tab(uvdata,'CL')
                imagr.sources[1:] = str(file[:8]),''
                print imagr.sources
                imagr.outname = file[:8]+'PB'
                imagr.nchav = 32
                imagr.niter = 200
                imagr.imsize[1:] = 1024, 1024
                imagr.cellsize[1:] = 0.001,0.001
                imagr.go()
                imagedata = AIPSImage(file[:8]+'PB','ICL001',1,1)
                fittp = AIPSTask('FITTP')
                fittp.indata = imagedata
                fittp.dataout = IMAGEFITSFILESout+file[:8]+'_PBCOR_IM.fits'
                fittp.go()
                fittp.indata = AIPSUVData(file[:8]+'PB','TASAV',1,1)
                fittp.dataout = UVFITSFILESout+file[:8]+'_PBCOR_TASAV.fits'
                fittp.go()
                uvdata.zap()
                imagedata.zap()
                AIPSImage(file[:8]+'PB','IBM001',1,1).zap()
		AIPSUVData(file[:8]+'PB','TASAV',1,1).zap()
