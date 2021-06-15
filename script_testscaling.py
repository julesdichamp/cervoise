import configparser
import os
import matplotlib.pyplot as plt

nbTime = 5 # max time_end
randTry = 5 # number of tries with a different seed_num
d_tot = 100
nbTravels = []
dProd = []
timeTab = []

parser = configparser.ConfigParser()
parser.read('config_VALIDATION_scaling_of_hopping_rates_on_regions_metropole.ini')

for k in range(1, nbTime):
    parser.set('settings', 'time_end', str(k))
    with open('config_VALIDATION_scaling_of_hopping_rates_on_regions_metropole.ini', 'w') as configfile:
        parser.write(configfile)
    
    sumTry = 0
    for l in range(randTry):
        parser.set('settings', 'seed_num', str(l))
        with open('config_VALIDATION_scaling_of_hopping_rates_on_regions_metropole.ini', 'w') as configfile:
            parser.write(configfile)
        os.system("python cervoise.py config_VALIDATION_scaling_of_hopping_rates_on_regions_metropole.ini")
        with open('nbTravels.txt', 'r') as travels:
            sumTry += int(travels.read())
            
    nbTravels.append(sumTry / randTry)
    dProd.append(d_tot * k)
    timeTab.append(k)
    
plt.plot(timeTab, nbTravels, "ro", label="nb of travels")
plt.plot(timeTab, dProd, label="d_tot * t")
plt.legend()

plt.show()
