import configparser
import os
import matplotlib.pyplot as plt

nbBeta = 9 # max beta
randTry = 1 # number of tries with a different seed_num
beta_tab = []
rt_tab = []

parser = configparser.ConfigParser()
par = configparser.ConfigParser()
parser.read('config_VALIDATION_simple_return_to_home_on_regions_metropole.ini')
output_dir = parser.get("outputs", "outputs_dir")

for k in range(0, nbBeta):
    parser.set('model', 'beta', str(k))
    with open('config_VALIDATION_simple_return_to_home_on_regions_metropole.ini', 'w') as configfile:
        parser.write(configfile)
    
    sumTry = 0
    for l in range(randTry):
        parser.set('settings', 'seed_num', str(l))
        with open('config_VALIDATION_simple_return_to_home_on_regions_metropole.ini', 'w') as configfile:
            parser.write(configfile)
        os.system("python cervoise.py config_VALIDATION_simple_return_to_home_on_regions_metropole.ini")
        
        k_S = parser.getfloat("model", "k_S")
        S_ini = parser.getint("model", "S_ini")
        t_end = parser.getfloat("settings", "time_end")
        d_tot = parser.getint("model", "d_tot")
        
        string_for_savefig = '_k_S='+ str(k_S) +'_S_ini='+ str(S_ini) + '_beta='+ str(float(k)) + '_d_tot='+ str(d_tot) +'_timeend='+ str(t_end) +'_seed_num='+str(l)
        par.read(output_dir +"/"+ "Informations" +str(string_for_savefig) + '.ini')
        sumTry += par.getfloat("informations", "average_return_time2")
            
    rt_tab.append(sumTry / randTry)
    beta_tab.append(k)
    
plt.plot(beta_tab, rt_tab)
plt.xlabel("values of beta")
plt.ylabel("average return times")
plt.tight_layout()

plt.show()