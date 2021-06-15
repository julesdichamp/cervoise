import configparser
import os
import time

nbTime = 3 # max time_end
randTry = 2 # number of tries with a different seed_num
runningTime = []

parser = configparser.ConfigParser()
parser.read('config_VALIDATION_SIR_compartments_fullConnectedRegions.ini')

for k in range(2, nbTime):
    parser.set('model', 's_ini', str(18 * (10**k)))
    with open('config_VALIDATION_SIR_compartments_fullConnectedRegions.ini', 'w') as configfile:
        parser.write(configfile)
    
    sumTry = 0
    for l in range(randTry):
        parser.set('settings', 'seed_num', str(l))
        with open('config_VALIDATION_SIR_compartments_fullConnectedRegions.ini', 'w') as configfile:
            parser.write(configfile)
        tps1 = time.time()
        os.system("python cervoise.py config_VALIDATION_SIR_compartments_fullConnectedRegions.ini")
        tps2 = time.time()
        sumTry += tps2 - tps1
            
    runningTime.append(sumTry / randTry)
    
with open('runTime.txt', 'w') as filehandle:
    for listitem in runningTime:
        filehandle.write('%s\n' % str(listitem))
