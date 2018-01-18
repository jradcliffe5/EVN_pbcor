import os, sys

### Defaults
import os, re, sys, exceptions
### Plotter
import matplotlib.pyplot as plt
### Numerics
import numpy as np
from datetime import datetime
startTime = datetime.now()
### For the emailer and logger
from code_mailer import setup_logging_to_file, log_exception, gmail_emailer, headless
import platform
import logging

### Setup logger
log_name = "%s.log" % os.path.basename(__file__).split('.py')[0]
setup_logging_to_file(log_name)
logging.info('Beginning %s' % os.path.basename(__file__))
print 'Beginning %s' % os.path.basename(__file__)

### Email credentials
email_creds = headless('mail_passwords.txt')
user = email_creds['username']
pwd = email_creds['pwd']

inputs = headless('inputs.txt')
vexfile = str(inputs['vex_file'])
uvdir = str(inputs['UV_dir'])
parseltongue = str(inputs['parseltongue'])
aipsno = str(inputs['AIPSnumber'])
imdir = str(inputs['IM_output_dir'])
tasavdir = str(inputs['TASAV_output_dir'])

try:
    ### First make the pbcor file
    if (vexfile.endswith('.vex') or vexfile.endswith('.vax') or vexfile.endswith('.vix')):
        os.system('rm pbcor.txt')
        os.system('python create_table.py %s pbcor.txt' % vexfile)
    ### Now apply the table to every dataset
    for i in os.listdir(uvdir):
        if i.endswith('.UV'):
            ### Copy data
            logging.info('Copying %s to ./' % i)
            logging.info('rsync -ar --progress %s%s ./' % (uvdir,i))
            os.system('rsync -ar --progress %s%s ./' % (uvdir,i))

            ### Load into AIPS
            logging.info('Loading %s to AIPS number %s' % (i,aipsno))
            logging.info('%s fitld_load.py %s' % (parseltongue,i))
            os.system('%s fitld_load.py %s' % (parseltongue,i))

            ### Run applycal
            aips_info = np.load('aips_info.npy')
            logging.info('Running apply_cal.py on %s %s %s %s (AIPSno %s)' % (aips_info[0],aips_info[1],aips_info[2],aips_info[3],aipsno))
            os.system('%s apply_cal.py pbcor.txt %s %s %s %s %s' % (parseltongue,aipsno,aips_info[0],aips_info[1],aips_info[2],aips_info[3]))

            ## Run CLCAL and make IMAGES + TASAV
            logging.info('Running imaging + calibration on %s %s %s %s (AIPSno %s)' % (aips_info[0],aips_info[1],aips_info[2],aips_info[3],aipsno))
            os.system('%s apply_image.py' % (parseltongue))

            logging.info('Moving products to desired directories')
            os.system('rm %s' % i)
            if os.path.exists(imdir) == False:
                os.system('mkdir %s' % imdir)
            if os.path.exists(tasavdir) == False:
                os.system('mkdir %s' % tasavdir)
            os.system('mv *IM.fits %s' % imdir)
            os.system('mv *TASAV.fits %s' % tasavdir)
    
    gmail_emailer(user=user,pwd=pwd,recipient='j.f.radcliffe@rug.nl',subject='CODE %s RUN SUCCESSFULLY - %s' % (os.path.basename(__file__),platform.node()),body='The code %s has run successfully for %s.\n\n Please see %s:%s  for the results.\n\n The logger output of %s is as follows: %s' % (os.path.basename(__file__),datetime.now() - startTime, platform.node(),os.path.dirname(os.path.realpath(__file__)), log_name, open(log_name,'r').read()))
    
except exceptions.Exception as e:
    log_exception(e)
    
    gmail_emailer(user=user,pwd=pwd,recipient='j.f.radcliffe@rug.nl',subject='CODE %s FAILED - %s' % (os.path.basename(__file__),platform.node()),body='The code %s has failed.\n\n Please see %s:%s for the details.\n\n The logger and associated error messages of %s is as follows:\n\n %s' % (os.path.basename(__file__), platform.node(),os.path.dirname(os.path.realpath(__file__)), log_name, open(log_name,'r').read()))
    
