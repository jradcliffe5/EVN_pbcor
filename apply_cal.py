#!/usr/bin/env ParselTongue
import AIPS
from Wizardry.AIPSData import AIPSUVData as wizAIPSUVData
from Wizardry.AIPSData import AIPSTableRow
import sys, operator, pdb
import vex_time, datetime, math

def parse_config(filename):
  config = open(filename, 'r')
  line = config.readline()
  FIND_BLOCK, READ_PCENTERS_MODE, READ_PCENTERS_POINTINGS, READ_PCENTERS_SOURCES  = range(4)
  READ_PCENTERS_STATIONS, READ_PCENTERS_CHANS, READ_SCANS, READ_NEW_SCAN = range(4,8)
  state = FIND_BLOCK
  line_nr = 1
  reason = ''
  pcenters = None
  scans = None
  error = False
  while line != "":
    sline = line.partition('*')[0].upper().split()
    if sline == []:
      pass
    elif state == FIND_BLOCK:
      if sline[0] == "$PHASE_CENTERS":
        pcenters = {}
        state = READ_PCENTERS_MODE
      elif sline[0] == '$SCANS':
        scans = []
        state = READ_SCANS
      else:
        reason = 'expected either $PHASE_CENTERS or $SCANS'
        error = True
    elif state == READ_PCENTERS_MODE:
      if sline[0] == 'DEF' and sline[1] == 'MODE':
        mode = sline[3]
        pcenters[mode] = {}
        state = READ_PCENTERS_POINTINGS
      elif sline[0] == '$SCANS':
        scans = []
        state = READ_SCANS
      else:
        if scans == None:
          reason = "Expected new mode or EOF"
        else:
          reason = "Expected new scan or $PHASE_CENTERS block"
        error = True
    elif state == READ_PCENTERS_POINTINGS:
      if sline[0] == 'DEF' and sline[1] == 'POINTING':
        pointing = sline[3]
        pcenters[mode][pointing] = {}
        state = READ_PCENTERS_SOURCES
      elif sline[0] == 'ENDDEF':
        state = READ_PCENTERS_MODE 
      else:
        reason = 'Expected new pointing or enddef.'
        error = True
    elif state == READ_PCENTERS_SOURCES:
      if sline[0] == 'DEF' and sline[1] == 'SOURCE':
        source = sline[3]
        pcenters[mode][pointing][source] = {}
        state = READ_PCENTERS_STATIONS
      elif sline[0] == 'ENDDEF':
        state = READ_PCENTERS_POINTINGS
      else:
        reason = 'expected new source or enddef.'
        error = True
    elif state == READ_PCENTERS_STATIONS:
      if sline[0] == 'DEF' and sline[1] == 'STATION':
        station = sline[3]
        pcenters[mode][pointing][source][station] = {'freqs':[], 'values' : []}
        state = READ_PCENTERS_CHANS
      elif sline[0] == 'ENDDEF':
        state = READ_PCENTERS_SOURCES
      else:
        reason = 'expected new station or enddef.'
        error = True
    elif state == READ_PCENTERS_CHANS:
      if sline[0] == 'ENDDEF':
        state = READ_PCENTERS_STATIONS
      elif len(sline) == 2:
        f = int(sline[0])
        c = float(sline[1])        
        pcenters[mode][pointing][source][station]['freqs'].append(f)
        if abs(c) < 1e-6:
          pcenters[mode][pointing][source][station]['values'].append(1.0)
        else:
          pcenters[mode][pointing][source][station]['values'].append(1.0 / abs(c))
      else:
        print 'len = ', len(sline), ', isint=', isinstance(sline[0], int), ', isfloat=',isinstance(sline[0], float)
        reason = 'Error parsing channel'
        error = True
    elif state == READ_SCANS:
      if sline[0] == 'DEF' and sline[1] == 'SCAN':
        new_scan = {}
        state = READ_NEW_SCAN
      elif sline[0] == '$PHASE_CENTERS':
        pcenters = {}
        state = READ_PCENTERS_MODE
      else:
        if pcenters == None:
          reason = "Expected new scan or EOF"
        else:
          reason = "Expected new scan or $PHASE_CENTERS block"
        error = True
    elif state == READ_NEW_SCAN:
      if sline[0] == 'ENDDEF':
        scans.append(new_scan)
        state = READ_SCANS
      elif sline[0] == 'START' or sline[0] == 'END':
        new_scan[sline[0]] = vex_time.get_time(sline[2].lower())
      elif sline[0] == 'POINTING':
        stripped = line.partition('*')[0].upper()
        start = stripped.index(':')
        stations = stripped[start+1:].replace(',', ' ').split()
        source = sline[2]
        try:
          new_scan['POINTING'][source] = stations
        except KeyError:
          new_scan['POINTING'] = {source : stations}
      else:
        new_scan[sline[0]] = sline[2]
    if error:
      print 'Error in line ', line_nr, ' : ', reason
      print 'Offending line : ', line
      exit(1)
    line_nr += 1
    line = config.readline()
  return pcenters, scans

def get_sources(uvdata):
  su = uvdata.table('SU', 0)
  return dict((s['id__no'], s['source'].strip()) for s in su)

def get_scan(scans, curscan, time):
  if (curscan != None) and (curscan['START'] <= time) and (curscan['END'] > time):
    return curscan
 
  for scan in scans:
    if (scan['START'] <= time) and (scan['END'] > time):
      return scan
  print 'Error : time ',time, ' not found in scans'
  sys.exit(1)

def get_pointing(scan, ant):
  for src, stations in scan['POINTING'].iteritems():
    if ant in stations:
      return src
  print 'Error : station ', ant, ' not in scan start = ', scan['START'], ' end = ', scan['END']

def apply_cal(uvdata, pcenters, scans, sources):
  oldcl = uvdata.table('CL', 0)
  newcl = uvdata.attach_table('CL', 0, no_term=oldcl.keywords['NO_TERM'])
  newcl.keywords['NO_ANT'] = oldcl.keywords['NO_ANT']
  i = uvdata.header['ctype'].index('FREQ')
  f0 = uvdata.header.crval[i]
  start_obs = [int(t) for t in uvdata.header.date_obs.split('-')]
  t0 = datetime.datetime(start_obs[0], start_obs[1], start_obs[2])
  scan = None
  for row in oldcl:
    dt = round(row.time*24*60*60)
    time = t0 + datetime.timedelta(seconds = dt)
    scan = get_scan(scans, scan, time)
    ant = uvdata.antennas[row.antenna_no -1]
    mode = scan['MODE']
    pointing = get_pointing(scan, ant)
    src = sources[row.source_id]
    print mode,pointing,src,ant
    try:
      ctab = pcenters[mode][pointing][src][ant]
      row.real1 = map(operator.mul, row.real1, ctab['values'])
      row.imag1 = map(operator.mul, row.imag1, ctab['values'])
      if len(uvdata.polarizations) == 2:
        row.real2 = map(operator.mul, row.real2, ctab['values'])
        row.imag2 = map(operator.mul, row.imag2, ctab['values'])
    except KeyError:
      # Stations not participating in scan can still have CL entries
      pass
    newcl.append(row)
  newcl.close()
  oldcl.close()

def create_sn(uvdata, pcenters, scans, sources):
  try:
    nx = uvdata.table('NX', 0)
  except IOError, msg:
    print 'Could not open NX table. Reason ', msg
    sys.exit(1)
  sn = uvdata.attach_table('SN', 0)
  start_obs = [int(t) for t in uvdata.header.date_obs.split('-')]
  t0 = datetime.datetime(start_obs[0], start_obs[1], start_obs[2])
  scan = None
  for row in nx:
    dt = round(row.time*24*60*60)
    time = t0 + datetime.timedelta(seconds = dt)
    print 'dt = ', dt, ', t0 = ', t0, ', time = ', time
    scan = get_scan(scans, scan, time)
    mode = scan['MODE']
    src = sources[row.source_id]
    tstart = row.time - row.time_interval/2 
    tend = row.time + row.time_interval/2 
    times = [tstart, tend]
    for atime in times:
      for pointing, stations in scan['POINTING'].iteritems():
        for s in stations:
          snr = AIPSTableRow(sn)
          ant = uvdata.antennas.index(s) + 1
          print mode,pointing,src,ant
          ctab = pcenters[mode][pointing][src][s]
          nif = len(ctab['freqs'])
          snr.antenna_no = ant
          snr.time = atime
          snr.freq_id = 1
          snr.real1 = ctab['values']
          snr.imag1 = [0.]*nif
          snr.weight_1 = [1.]*nif
          if len(uvdata.polarizations) == 2:
            snr.real2 = ctab['values']
            snr.imag2 = [0.]*nif
            snr.weight_2 = [1.]*nif
          sn.append(snr)
  nx.close()
  sn.close()

################### MAIN ####################
######
######
argv = sys.argv
if len(argv) != 7:
  print 'Usage : ', argv[0], ' <pbcal file> <AIPS USER> <AIPS name> <AIPS class> <AIPS disk> <AIPS seq>'
  sys.exit(1)

pcenters, scans = parse_config(argv[1])
AIPS.userno = int(argv[2])
uvdata = wizAIPSUVData(argv[3], argv[4], int(argv[5]), int(argv[6]))
sources = get_sources(uvdata)
#apply_cal(uvdata, pcenters, scans, sources)
create_sn(uvdata, pcenters, scans, sources)
