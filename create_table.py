
# (c) 2012 k.keimpema keimpema@jive.nl
import os,sys
#sys.path.append("/usr/local/lib/python2.7/site-packages")
from pylab import *
from scipy.special import j1
from vex import Vex
import sys, pdb

def station_table():
  stations = {}
  # Effective appertuere size in meters for each station, currently only L band
  stations["ARECIBO"] = 206 # Effective, physical diameter = 305
  stations["BADARY"] = 32 / 1.05
  stations["CAMBG32M"] = 32 / 1.05
  stations["EFLSBERG"] = 78 # 100
  stations["HART"] = 26/1.05
  stations["IRBENE"] = 30 / 1.05
  stations["JODRELL1"] = 67. # 76
  stations["JODRELL2"] = 25/1.05
  stations["KUNMING"] = 40/1.05
  stations["MEDICINA"] = 25 #32/1.05
  stations["METSAHOV"] = 14/1.05
  stations["NOTO"] = 32 / 1.05
  stations["ONSALA60"] = 20/1.05
  stations["ONSALA85"] = 25/1.05
  stations["ONS_DBBC"] = stations["ONSALA85"]
  stations["SHANGHAI"] = 22.5 # 25
  stations["SVETLOE"] = 32/1.05
  stations["TORUN"] = 32/1.05
  stations["URUMQI"] = 25/1.05
  stations["WSTRBORK"] = 25 / 1.05
  stations["YEBES40M"] = 40 / 1.05
  stations["ZELENCHK"] = 32 / 1.05
  stations["GBT_VLBA"] = 100. / 1.05
  return stations

def vex_stations(vex, station_par):
  stations = {}
  count = 0
  for key, val in vex["STATION"].iteritems():
    ID = vex["SITE"][val["SITE"]]["site_ID"]
    name = vex["SITE"][val["SITE"]]["site_name"]
    stations[ID] = (count, station_par[name])
    count += 1
  return stations

def vex_frequencies(vex):
  frequencies = {}
  for key, val in vex['MODE'].iteritems():
    freqs = []
    fblock = val['FREQ'][0]
    for chan in vex['FREQ'][fblock].getall('chan_def'):
      freq = float(chan[1].split()[0]) * 1e6
      sb = -1 if chan[2] == 'L' else 1
      bw = float(chan[3].split()[0]) * 1e6
      mid = freq + sb*bw/2
      if mid not in freqs:
        freqs.append(mid)
    freqs.sort()
    frequencies[key] = freqs
  return frequencies

# This code is adapted from a script by Stephen Bourke
def calc_shift(ra_str, dec_str, ra0_str, dec0_str):
  h,m,s = (ra_str.find('h'), ra_str.find('m'), ra_str.find('s'))
  sra = [ra_str[:h], ra_str[h+1:m], ra_str[m+1:s]]
  d,q,qq = (dec_str.find("d"), dec_str.find("'"), dec_str.find('"'))
  sdec = [dec_str[:d], dec_str[d+1:q], dec_str[q+1:qq]]
  h,m,s = (ra0_str.find('h'), ra0_str.find('m'), ra0_str.find('s'))
  sra_0 = [ra0_str[:h], ra0_str[h+1:m], ra0_str[m+1:s]]
  d,q,qq = (dec0_str.find("d"), dec0_str.find("'"), dec0_str.find('"'))
  sdec_0 = [dec0_str[:d], dec0_str[d+1:q], dec0_str[q+1:qq]]

  # Convert to decimal arc-degrees
  ra = 15 * (float(sra[0]) + float(sra[1])/60 + float(sra[2])/3600)
  dec = float(sdec[0]) + float(sdec[1])/60 + float(sdec[2])/3600

  ra_0 = 15 * (float(sra_0[0]) + float(sra_0[1])/60 + float(sra_0[2])/3600)
  dec_0 = float(sdec_0[0]) + float(sdec_0[1])/60 + float(sdec_0[2])/3600

  # Distance calculation
  shift_x = cos(radians(dec_0)) * (ra - ra_0)
  shift_y = dec - dec_0
  return pi*hypot(shift_x, shift_y)/180.

def get_phase_centers(vex, stations, freqs):
  phase_centers={}
  for scan_name, scan in vex['SCHED'].iteritems():
    mode = scan['mode']
    if mode not in  phase_centers:
      phase_centers[mode]={}
    sources = scan.getall('source')
    pcenter = sources[0]
    if pcenter not in phase_centers[mode]:
      phase_centers[mode][pcenter] = {}
    pc = vex['SOURCE'][pcenter]
    for source in sources:
      if source not in phase_centers[mode][pcenter]:
        phase_centers[mode][pcenter][source] = {}
      src = vex['SOURCE'][source]
      delta = calc_shift(src['ra'], src['dec'], pc['ra'], pc['dec'])
      for station_line in scan.getall('station'):
        station = station_line[0]
        if station not in phase_centers[mode][pcenter][source]:
          f = [0]*len(freqs[mode])
          idx, D = stations[station]
          for i, freq in enumerate(freqs[mode]):
            k = 2 * pi * freq / 299792458
            if delta == 0:
              cfactor = 1.
            else:
              cfactor = 2*j1(k*(D/2.)*sin(delta))/(k*(D/2.)*sin(delta))
            print 'i=%d, f=%e, k=%e, D = %f, c= %e, delta=%e'%(i, freq,k,D,cfactor,delta)
            f[i] = (freq, cfactor)
          phase_centers[mode][pcenter][source][station] = f
  return phase_centers

def write_phase_centers(vex, phase_centers, outfile):
  outfile.write("$PHASE_CENTERS\n")
  for mode, pcenter_dict in phase_centers.iteritems():
    outfile.write("def mode = " + mode + '\n')
    for pcenter, src_dict in pcenter_dict.iteritems():
      outfile.write("  def pointing = " + pcenter + '\n')
      for src, station_dict in src_dict.iteritems():
        outfile.write("    def source = " + src + '\n')
        for station, channels in station_dict.iteritems():
          outfile.write("      def station = " + station + '\n')
          for chan in channels:
            outfile.write('         %d\t%f\n'%(int(chan[0]), chan[1]))
            #outfile.write(' '*8+`int(chan[0])` + '\t' + `chan[1]` + '\n')
          outfile.write('      enddef\n')
        outfile.write('    enddef\n')
      outfile.write('  enddef\n')
    outfile.write('enddef\n')

def write_scans(vex, outfile):
  outfile.write('$SCANS\n')
  for scan_name, scan in vex['SCHED'].iteritems():
    duration = scan_duration(scan)
    start = vtime(scan['start'])
    end = vtime_add(start, duration)
    outfile.write('def scan\n')
    outfile.write('  start = ' + scan['start'] + '\n')
    outfile.write('  end = ' + vtime_str(end) + '\n')
    stations = [s[0] for s in scan.getall('station')]
    outfile.write('  mode = ' + scan['mode']+'\n')
    outfile.write('  pointing = ' + scan['source'] + ' : ' + ", ".join(stations) + '\n')
    outfile.write('enddef\n')

def scan_duration(scan):
  stations = scan.getall('station')
  dur = 0
  for s in stations:
    end = s[2].index('s')
    dur = max(dur, int(s[2][:end]))
  return dur

def vtime(vstr):
  i = [vstr.index(c) for c in "ydhms"]
  y = int(vstr[:i[0]])
  d = int(vstr[i[0]+1:i[1]])
  h = int(vstr[i[1]+1:i[2]])
  m = int(vstr[i[2]+1:i[3]])
  s = int(vstr[i[3]+1:i[4]])
  return [y,d,h,m,s]

def vtime_add(vt, sec):
  days_per_year = 366 if vt[0]%4 == 0 else 365
  s = vt[4] + int(sec)
  m = vt[3] + s / 60
  h = vt[2] + m / 60
  d = vt[1] + h / 24
  y = vt[0] + d / days_per_year
  return [y, d%days_per_year, h%24, m%60, s%60]

def vtime_str(vt):
  return `vt[0]`+'y'+`vt[1]`+ 'd' + `vt[2]` + 'h' +`vt[3]` + 'm' + `vt[4]` + 's'

############################### MAIN #########################################
##############
#############
if len(sys.argv) != 3:
  print "Usage : ", sys.argv[0], " <vex file> <output file>"
  sys.exit(1)

vex = Vex(sys.argv[1])
try:
  outfile = open(sys.argv[2], "w")
except IOError:
  print "Error opening file : ", sys.argv[2]
  sys.exit(1)

station_par = station_table()
stations = vex_stations(vex, station_par)
frequencies = vex_frequencies(vex)
phase_centers = get_phase_centers(vex, stations, frequencies)
write_phase_centers(vex, phase_centers, outfile)
write_scans(vex, outfile)
