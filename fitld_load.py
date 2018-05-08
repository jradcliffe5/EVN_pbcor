from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
from code_mailer import headless
import numpy as np
import sys
import os

inputs = headless('inputs.txt')
AIPS.userno = int(inputs['AIPSnumber'])
indisk = int(inputs['AIPSdisk'])
x = sys.argv[1]

print sys.argv[1]
fitld = AIPSTask('FITLD')
fitld.datain = 'PWD:%s' % sys.argv[1]
fitld.outname = '%s' % x[:8]
fitld.outclass = 'PBCOR'
fitld.outseq = 0
fitld.outdisk = indisk
fitld.digicor = -1
fitld.go()

os.system('rm aips_info.npy')
np.save('aips_info.npy',[x[:8],'PBCOR',1,indisk])
