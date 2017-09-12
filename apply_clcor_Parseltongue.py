#!/usr/bin/env ParselTongue
import AIPS
from Wizardry.AIPSData import AIPSUVData as wizAIPSUVData
from Wizardry.AIPSData import AIPSTableRow
from AIPSTask import AIPSTask, AIPSList
from AIPSData import  AIPSImage
import sys, operator, pickle, os

## Inputs ##
UVFITSFILESpath = './'
IMAGEFITSFILESout = './'
UVFITSFILESout = './'
FLAGFILESin = './flag_tables/'
AIPS.userno = 9
disk = 1
pointing_names = ['HDFN','P1','P2','P3','P4']
#############

def parse_flag_times(flagfile,FLAGFILESin):
	file=open(FLAGFILESin+flagfile,'r')
	row = file.readlines()
	dict_return = {}
	times = []
	try:
		for line in row:
			if line.startswith('!!'):
				Pointing_name = line.rstrip('\n').split('!! ')[1]
			elif line.startswith('TIMERANG'):
				times = times +[map(int,line.split('ANTENNAS', 1)[0].split('=',1)[1].replace(" ", "").split(','))]
		dict_return[Pointing_name] = times
	except UnboundLocalError:
		print 'Need correct flag format'
		sys.exit()

	return dict_return

def get_tab(uvdata, table):
	# find the number of tables of a certain type
	ver = 0
	for i in range(len(uvdata.tables)) :
		if table in uvdata.tables[i][1] :
			ver = uvdata.tables[i][0]
	print "HIGHEST TABLE OF TYPE", table, "is", ver
	return ver

def input_files(pickle_file):
	params = pickle.load(open("%s" % pickle_file, "rb" )) ## pickle file made by primary beam sim
	text = open('%s.txt' % pickle_file,'w')
	text.close()
	return params

#### These are the correction factors
central_pointing_params = input_files('central_pointing_params.pckl')
outside_pointing_params = input_files('outside_pointing_params.pckl')

###################################################################
### This function parses the timeranges from the user made fits files
### The !! in the file is needed to match with the name of the pointing centers
### It creates a dictionary that can get a list of timeranges using timerange[pointing name]
timerange = {}
for file in os.listdir(FLAGFILESin):
	for i in pointing_names:
		if i in file:
			timerange.update(parse_flag_times(file,FLAGFILESin))
###################################################################


########### Firstly get a list of all the UV data files to derive corrections for
for file in os.listdir(UVFITSFILESpath):
	####### Check if each of these files end in .UV... just out of courtesy
    if file.endswith('.UV'):
		### Now need to firstly correct for the central pointing

		for i in range(len(central_pointing_params)):
			if central_pointing_params[i][0] == file[:8]:
				### Slice the correction factors based upon source name
				correction_factor = central_pointing_params[i]
				print correction_factor
				print 'Correcting %s' % file[:8]

				### Load the data into an AIPS directory
				fitld = AIPSTask('FITLD')
				fitld.datain = UVFITSFILESpath+file
				fitld.digicor = -1
				fitld.outname = file[:8]
				fitld.outclass = 'UVDATA'
				fitld.go()

				uvdata = wizAIPSUVData(file[:8],'UVDATA',disk,1)

				### Apply the corrections for the central telescopes
				i = 0
				for j in range(len(correction_factor[3])):
					i = i + 1
					clcor = AIPSTask('CLCOR')
					clcor.indata = uvdata
					clcor.sources[1] ='*'
					print float(correction_factor[3][j].keys()[0])
					clcor.antennas[1] = int(correction_factor[3][j].keys()[0])
					clcor.opcode = 'GAIN'
					print float(correction_factor[3][j].values()[0])
					clcor.clcorprm[1] = float(correction_factor[3][j].values()[0])
					clcor.gainver = get_tab(uvdata,'CL')
					if i == 1:
						clcor.gainuse = get_tab(uvdata,'CL') + 1
					else:
						clcor.gainuse = get_tab(uvdata,'CL')
					clcor.go()
				f = open('central_pointing_params.pckl.txt','a')
				f.write(str(correction_factor[0])+' '+str(correction_factor[3])+'\n')
				f.close()

				### Now apply the corrections to the other telescopes
		for pointing in outside_pointing_params.keys():
			for i in range(len(outside_pointing_params[pointing])):
				if outside_pointing_params[pointing][i][0] == file[:8]:
					### Slice the correction factors based upon source name
					correction_factor = outside_pointing_params[pointing][i]
					pointing_timeranges = timerange[pointing]
					print pointing
					print pointing_timeranges
					print correction_factor
					print 'Correcting %s' % file[:8]

					uvdata = wizAIPSUVData(file[:8],'UVDATA',disk,1)

					### Apply the corrections for the central telescopes
					i = 0
					for k in pointing_timeranges:
						i = i + 1
						for j in range(len(correction_factor[3])):
							i = i + 1
							clcor = AIPSTask('CLCOR')
							clcor.indata = uvdata
							clcor.sources[1] ='*'
							print float(correction_factor[3][j].keys()[0])
							clcor.antennas[1] = int(correction_factor[3][j].keys()[0])
							clcor.opcode = 'GAIN'
							clcor.timerang = AIPSList(k)
							print float(correction_factor[3][j].values()[0])
							clcor.clcorprm[1] = float(correction_factor[3][j].values()[0])
							clcor.gainver = get_tab(uvdata,'CL')
							if i == 1:
								clcor.gainuse = get_tab(uvdata,'CL') + 1
							else:
								clcor.gainuse = get_tab(uvdata,'CL')
							clcor.go()
					f = open('central_pointing_params.pckl.txt','a')
					f.write(str(correction_factor[0])+' '+str(correction_factor[3])+'\n')
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
		imagr.uvwtfn = 'NA'
		imagr.go()
		imagedata = AIPSImage(file[:8]+'PB','ICL001',1,1)
		imagedata2 = AIPSImage(file[:8]+'PB','ICL001',1,2)
		fittp = AIPSTask('FITTP')
		fittp.indata = imagedata
		fittp.dataout = IMAGEFITSFILESout+file[:8]+'_PBCOR_IM.fits'
		fittp.go()
		fittp.indata = imagedata2
		fittp.dataout = IMAGEFITSFILESout+file[:8]+'_PBCOR_NA_IM.fits'
		fittp.go()
		fittp.indata = wizAIPSUVData(file[:8]+'PB','TASAV',1,1)
		fittp.dataout = UVFITSFILESout+file[:8]+'_PBCOR_TASAV.fits'
		fittp.go()
		uvdata.zap()
		imagedata.zap()
		imagedata2.zap()
		AIPSImage(file[:8]+'PB','IBM001',1,1).zap()
		AIPSImage(file[:8]+'PB','IBM001',1,2).zap()
		wizAIPSUVData(file[:8]+'PB','TASAV',1,1).zap()
