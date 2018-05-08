def parse_flag_times(flagfile):
    file=open(flagfile,'r')
    row = file.readlines()
    dict_return = {}
    times = []
    for line in row:
        if line.startswith('!!'):
            Pointing_name = line.rstrip('\n').split('!! ')[1]
            print Pointing_name
        elif line.startswith('TIMERANG'):
            times = times +[map(int,line.split('ANTENNAS', 1)[0].split('=',1)[1].replace(" ", "").split(','))]
    dict_return[Pointing_name] = times
    return dict_return



print parse_flag_times('EF_P1.txt')
