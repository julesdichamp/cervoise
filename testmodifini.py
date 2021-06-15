import configparser

k = 6
parser = configparser.SafeConfigParser()
parser.set('settings', 'time_end', 'k')


with open('config_VALIDATION_scaling_of_hopping_rates_on_regions_metropole.ini', 'w') as configfile:
    parser.write(configfile)