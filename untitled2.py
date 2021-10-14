import configparser
import os
import time

config = configparser.ConfigParser()
config.read('config_VALIDATION_Superspreaders_full_agent_based.ini)
            
print("InitParams")
config = configparser.ConfigParser()
config.read(config_name)AllGraphs
params = {}
#  load input data infos
params["graph_name"] = config.get("inputs", "graph_name")
params["data_dir"] = config.get("inputs", "data_dir")
params["roads_name"] = config.get("inputs", "roads_name")
params["shapefile_name"] = config.get("inputs", "shapefile_name")
params["gareFile_name"] = config.get("inputs", "gareFile_name")
params["ageFile_name"] = config.get("inputs", "ageFile_name")
params['airplaneFile_name'] = config.get("inputs", "airplaneFile_name")
params['citiesFile_name'] = config.get("inputs", 'citiesFile_name')
print("InitParams")
config = configparser.ConfigParser()
config.read(config_name)
params = {}
#  load input data infos
params["graph_name"] = config.get("inputs", "graph_name")
params["data_dir"] = config.get("inputs", "data_dir")
params["roads_name"] = config.get("inputs", "roads_name")
params["shapefile_name"] = config.get("inputs", "shapefile_name")
params["gareFile_name"] = config.get("inputs", "gareFile_name")
params["ageFile_name"] = config.get("inputs", "ageFile_name")
params['airplaneFile_name'] = config.get("inputs", "airplaneFile_name")
params['citiesFile_name'] = config.get("inputs", 'citiesFile_name')
params['covid19DepartmentsDir_name'] = config.get("inputs", "covid19DepartmentsDir_name")
params['covid19RegionsDir_name'] = config.get("inputs", "covid19RegionsDir_name")
# get all file names we are interested in in this directory (do not keep polygons files)s
params['allCovid19DepartmentsFileNames'] = glob.glob(
params['covid19DepartmentsDir_name'] + '/*departements-france-*.json')
params['allCovid19RegionsFileNames'] = glob.glob(params['covid19RegionsDir_name'] + '/*regions-france-2020*.json')
params['departments_name'] = config.get("inputs", "departments_name")
params['citiesFile_name'] = config.get("inputs", "citiesFile_name")

# contains all data of getFranceGeomDpts(params) or getFranceGeomRegs(params)
# we can directly read it if it was saved so we do not compute it again
params['franceGeom_name'] = config.get("inputs", "franceGeom_name")
# location of saved file
# params['graphFile'] = config.get("inputs", 'graphFile')
# if not params['graphFile'].strip():
#     params['graphFile'] = 'franceGeomRegions.txt'

# load options infos
# params["genAndQuit"] = distutils.util.strtobool(config.get("graph", "genAndQuit"))  # quit after graph is generated and saved
# JD: it fails for me otherwise...
params["genAndQuit"] = strtobool(config.get("graph", "genAndQuit"))  # quit after graph is generated and saved
params["location"] = config.get("graph", "location")
# if params["location"] != "metropole":
# print('[ERROR]: location can only be metropole for now. Exit here.')
# exit(0)

params["graph_type"] = config.get("graph", "graph_type")  # what type of graph is is

# load model infos
params["model_name"] = config.get("model", "model_name")

params["k_SE_I"] = config.getfloat("model", "k_SE_I")
params["k_IR"] = config.getfloat("model", "k_IR")
params["k_SE_Iss"] = config.getfloat("model", "k_SE_Iss")
params["k_EI"] = config.getfloat("model", "k_EI")
params["S_ini"] = config.getint("model", "S_ini")
params["I_ini"] = config.getint("model", "I_ini")
params["R_ini"] = config.getint("model", "R_ini")
params["Sss_ini"] = config.getint("model", "Sss_ini")
params["Iss_ini"] = config.getint("model", "Iss_ini")
params["E_ini"] = config.getint("model", "E_ini")
params["Ess_ini"] = config.getint("model", "Ess_ini")
params["beta"] = config.getfloat("model", "beta")
params["betass"] = config.getfloat("model", "betass")
params["d_tot"] = config.getint("model", "d_tot")

print("InitGraph")
"""
initialize a graph of regions with positions, and ids per regions
in: params, dictionnary with all needed parameters
out: a graph of type networkx.Graph()
"""
### --- 1/ read data we need depending of what is written on params["graph_type"]
franceGeom = getFranceGeom.readFranceGeom(params)
G = nx.Graph()
params['total_pop'] = 0 # feed dictionnary with aggregated data from graph
params['total_covid_data'] = {} # feed dictionnary with aggregated data from graph

list_outremer = ['Mayotte', 'La Réunion', 'Guyane', 'Guadeloupe', 'Martinique']
    
### --- Construct graph
# same routine for regions and departments as it is handled by readFranceGeom(params)
if (params['graph_type'] == 'regions') or (params['graph_type'] == 'departments'):
    dic_total_population = {}
    dic_data_covid = {}
    dic_data_hosp = {}
    dic_positions = {}
    dic_names = {}
    dic_areas = {}        
# add nodes to the graph
    for key in franceGeom.keys():
        id_node = franceGeom[key]["id"]  # get ID attributed from readFranceGeom(params)
    
        if (params['location'] == 'metropole'):
            # only add departments/ergions which are in metropole
            if not(key in list_outremer):
                G.add_node(id_node)                    
                # add data on node_attributes, for now:
                # - Total population
                # - Confirmed cases (old source)
                # - Deaths (old source)
                # - Recovered (old source)
                # - Hospital data
                # - Position 
                # - Area
                dic_total_population[id_node] = franceGeom[key]['agesSex']['ftot'] + franceGeom[key]['agesSex']['htot']
                dic_data_covid[id_node] = franceGeom[key]['data_covid']
                dic_data_hosp[id_node] = franceGeom[key]['data_hosp']
                centroid = [franceGeom[key]['centroid'][0], franceGeom[key]['centroid'][1]]
                dic_positions[id_node] = centroid                    
                dic_names[id_node] = key
                dic_areas[id_node] = franceGeom[key]['area']                    
                # print("id_node, dic_names[id_node] = " , id_node, dic_names[id_node])
                params['total_pop'] += dic_total_population[id_node]
                for _date in dic_data_covid[id_node].keys():
                    if not(_date in params['total_covid_data']):
                        params['total_covid_data'][_date] = {}
                        params['total_covid_data'][_date]['Confirmed'] = 0
                        params['total_covid_data'][_date]['Deaths']    = 0
                        params['total_covid_data'][_date]['Recovered'] = 0
                    params['total_covid_data'][_date]['Confirmed'] += dic_data_covid[id_node][_date]['Confirmed']
                    params['total_covid_data'][_date]['Deaths']    += dic_data_covid[id_node][_date]['Deaths']
                    params['total_covid_data'][_date]['Recovered'] += dic_data_covid[id_node][_date]['Recovered']

    # add attributes with positions to graph
    nx.set_node_attributes(G, name='positions', values=dic_positions)
    nx.set_node_attributes(G, name='names', values=dic_names)
    nx.set_node_attributes(G, name='total_pop', values=dic_total_population)
    nx.set_node_attributes(G, name='data_covid', values=dic_data_covid)
    nx.set_node_attributes(G, name='data_hosp', values=dic_data_hosp)
    nx.set_node_attributes(G, name='areas', values=dic_areas)        
        
    if (params['location'] == 'metropole'):
        # add edges to the graph
        for key in franceGeom.keys():
            if not(key in list_outremer):
                id_node1 = franceGeom[key]["id"]  # get ID attributed from readFranceGeom(params)
                for keyneigh in franceGeom[key]['neighbours']:
                    if not(keyneigh in list_outremer):
                        id_node2 = franceGeom[keyneigh]["id"]  # get ID attributed from readFranceGeom(params)
                        edge12 = (id_node1, id_node2)
                        edge21 = (id_node2, id_node1)
                        if (not (edge12 in G.edges())) or (not (edge21 in G.edges())):
                            G.add_edge(*edge12)  # unpack edge tuple*
            # add edges based on airways connections
        for key in franceGeom.keys():
            if not(key in list_outremer):
                id_node1 = franceGeom[key]["id"]  # get ID attributed from readFranceGeom(params)
                for connection in franceGeom[key]['airways']:
                    keyneigh = connection[0]
                    if not(keyneigh in list_outremer):
                        pax = connection[1]  # not used yet
                        id_node2 = franceGeom[keyneigh]["id"]  # get ID attributed from readFranceGeom(params)
                        # make sure to not add loop on itself
                        if id_node1 != id_node2:
                            edge12 = (id_node1, id_node2)
                            edge21 = (id_node2, id_node1)
                            if (not (edge12 in G.edges())) or (not (edge21 in G.edges())):
                                    G.add_edge(*edge12)  # unpack edge tuple*



def plot_SuperspreadersFullAgentBased(params,outputs, G):
    """
    plot simple SEIR model and hopping between nodes of graph
    with an additionnal state of a superspreader Iss which has a much 
    faster rate
    with an agent based fashion where we track every individual
    depending values contained in params['plot_options'] we can either:
    WellMixed: plot the well mixed limit by summing all compartments
    AllNodes: plot the dynamic on every node separately
    AllGraphs: plot all graphs over all time steps
    FinalGraph: plot only the final graph at the final time step
    """
    print('plot_SuperspreadersFullAgentBased')
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    nbTimeSteps = len(outputs)
    mytime = np.zeros(nbTimeSteps)    
    nbStates = params['nbStates']
    nbReactions = params['nbReactions']
    x = np.zeros([nbTimeSteps, nbStates])
    dic_graphs = {}
    # orders of states : S, E, I, Iss, R
    for i in range(0,nbTimeSteps):
        mytime[i] = outputs[i][0]
        dic_x = outputs[i][1]        
        for n in dic_x.keys():
            for k in range(0,nbStates):
                x[i,k] += len(dic_x[n][k])
                
    total_infected = []   
    max_infected = 0
    index_max_infected = 0
    
    for k in range(nbTimeSteps):
        total_infected.append(x[k,4] + x[k,5])
        if(max_infected < total_infected[k]):
            max_infected = total_infected[k]
            index_max_infected = k
            
    time_max_infected = mytime[index_max_infected]
                
    print("max number of infected people =", max_infected)
    print("time of max infected =", time_max_infected)
    
    string_for_savefig = 'SuperspreadersFullAgentBased_graph_type_'+str(params['graph_type'])\
    +'_ratio_k='+str(params['k_SE_Iss'] / params['k_SE_I']) +'_ratio_beta='+ str(params['beta'] / params['betass'])+'_ratio_ss='+str(params['S_ini'] / params['Sss_ini'])
    
    with open(params["outputs_dir"]+"/"+ "Informations" +str(string_for_savefig) + '.ini', 'w') as f:
        f.write("[" + string_for_savefig + "]" + "\n" + "\n")
        f.write("[informations]" + "\n")
        f.write("max_infected = ")
        f.write(str(max_infected) + "\n")
        f.write("time_max_infected = ")
        f.write(str(time_max_infected))
    
    if ('WellMixed' in params['plot_options']):
        print("Running: plot option: 'WellMixed'")
        ### --- plot wellmixed model ( = sum over all compartments i.e. nodes of the graph)
        fig10 = plt.figure(figsize=(16, 10)) # in inches!
        plt.plot(mytime,x[:,0], 'orange', label='Susceptible')
        plt.plot(mytime,x[:,1], 'red', label='Susceptible superspreader')
        plt.plot(mytime,x[:,2], 'green', label='Exposed')
        plt.plot(mytime,x[:,3], 'yellow', label='Exposed superspreader')
        plt.plot(mytime,x[:,4], 'blue', label='Infected')
        plt.plot(mytime,x[:,5], 'black', label='Infected superspreader')
        plt.plot(mytime,x[:,6], 'pink', label='Recovered')
        plt.legend()
        plt.xlabel('Arbitrary time unit')
        plt.ylabel('Population numbers')
        # string_for_title='SEIIssROnGraphAgentBased, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
        # plt.title(string_for_title)
        # string_for_savefig = 'SEIIssROnGraphAgentBased_graph_type_'+str(params['graph_type'])+'_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
        string_for_savefig = 'SuperspreadersFullAgentBased_graph_type_'+str(params['graph_type'])\
          +'_k_SE_I='+str(params['k_SE_I']) +'_k_SE_Iss='+ str(params['k_SE_Iss'])+'_k_EI='+str(params['k_EI'])\
          +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini'])+ '_Sss_ini='+ str(params['Sss_ini']) +'_E_ini='+ str(params['E_ini'])\
          + '_Ess_ini='+ str(params['Ess_ini']) +'_I_ini='+str( params['I_ini']) + '_Iss_ini='+ str(params['Iss_ini']) + '_R_ini='+str( params['R_ini']) 
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
    
    if ('AllNodes' in params['plot_options']):
        print("Running: plot option: 'AllNodes'")
        ### --- loop over all saved data and plot SIRWellMixed on each region
        # use dic_x available in node attributes from graph to get 
        # number of nodes and number of reactions
        dic_x = nx.get_node_attributes(G,'x')        
        nbNodes = len(dic_x)
        nbStates = len(dic_x[0])
        nodes = G.nodes()
        all_nodes_all_x_over_time = {}
        for n in nodes:
            x = np.zeros([nbTimeSteps, nbStates])
            for i in range(0,nbTimeSteps):
                dic_x = outputs[i][1]            
                for j in range(0,nbStates):
                    x[i,j] = len(dic_x[n][j])
            all_nodes_all_x_over_time[n] = x
            fig10 = plt.figure(figsize=(16, 10)) # in inches!
            plt.plot(mytime,x[:,0], 'orange', label='Susceptible')
            plt.plot(mytime,x[:,1], 'red', label='Susceptible superspreader')
            plt.plot(mytime,x[:,2], 'green', label='Exposed')
            plt.plot(mytime,x[:,3], 'yellow', label='Exposed superspreader')
            plt.plot(mytime,x[:,4], 'blue', label='Infected')
            plt.plot(mytime,x[:,5], 'black', label='Infected superspreader')
            plt.plot(mytime,x[:,6], 'pink', label='Recovered')
            plt.legend()
            plt.xlabel('Arbitrary time unit')
            plt.ylabel('Population numbers')
            # string_for_title='SEIIssROnGraphAgentBased, region num = '+str(n)+', k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
            # plt.title(string_for_title)            
            string_for_savefig = 'SuperspreadersFullAgentBased_graph_type_'+str(params['graph_type'])+'_region_num='+str(n)\
            +'_k_SE_I='+str(params['k_SE_I']) +'_k_SE_Iss='+ str(params['k_SE_Iss'])+'_k_EI='+str(params['k_EI'])\
            +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini'])+ '_Sss_ini='+ str(params['Sss_ini']) +'_E_ini='+ str(params['E_ini'])\
            + '_Ess_ini='+ str(params['Ess_ini']) +'_I_ini='+str( params['I_ini']) + '_Iss_ini='+ str(params['Iss_ini']) + '_R_ini='+str( params['R_ini']) 
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
                            
    if ('AllGraphs' in params['plot_options']):
        print("Running: plot option: 'AllGraphs'")
        ### --- loop over all saved graph and plot at each saved time step    
        for j in range(0,nbTimeSteps):
            print("saving ",j," th time step = " ,mytime[j])        
            # dic_x = outputs[j][1]
            dic_x = copy.deepcopy(outputs[j][1])
            dic_positions = nx.get_node_attributes(G,'positions')        
            dic_S = {}
            dic_I = {}
            dic_R = {}
            for i in dic_x.keys():
                dic_S[i] = len(dic_x[i][0])
                dic_I[i] = len(dic_x[i][2])
                dic_R[i] = len(dic_x[i][4])
            
            list_S_values = [v for v in dic_S.values()]
            list_I_values = [v for v in dic_I.values()]
            list_R_values = [v for v in dic_R.values()]                        
            fig20 = plt.figure(figsize=(16, 10)) # in inches!  
            string_for_title='Susceptible population on Graph at time = '+str(round(mytime[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_S_values, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')    
            string_for_savefig = 'Susceptible_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
            plt.close()
            
            fig30 = plt.figure(figsize=(16, 10)) # in inches!
            string_for_title='Infected population on Graph at time = '+str(round(mytime[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_I_values, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')    
            string_for_savefig = 'Infected_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
            
            fig40 = plt.figure(figsize=(16, 10)) # in inches!            
            string_for_title='Recovered population on Graph at time = '+str(round(mytime[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_R_values, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')
            string_for_savefig = 'Recovered_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
            
    if ('FinalGraph' in params['plot_options']):
        print("Running: plot option: 'FinalGraph'")
        # ### --- plot graph with final species update == state in Graph
        dic_positions = nx.get_node_attributes(G,'positions')        
        dic_x = nx.get_node_attributes(G,'x')
        dic_S = {}
        dic_I = {}
        dic_R = {}
        for i in dic_x.keys():
            dic_S[i] = len(dic_x[i][0])
            dic_I[i] = len(dic_x[i][2])
            dic_R[i] = len(dic_x[i][4])
        
        list_S_values = [v for v in dic_S.values()]
        list_I_values = [v for v in dic_I.values()]
        list_R_values = [v for v in dic_R.values()]
        fig20 = plt.figure(figsize=(16, 10)) # in inches!  
        string_for_title='Susceptible population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_S_values, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')    
        string_for_savefig = 'Susceptible_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
        fig30 = plt.figure(figsize=(16, 10)) # in inches!
        string_for_title='Infected population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_I_values, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')    
        string_for_savefig = 'Infected_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
        fig40 = plt.figure(figsize=(16, 10)) # in inches!            
        string_for_title='Recovered population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_R_values, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')
        string_for_savefig = 'Recovered_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
        plt.close()

    if ('TrackSuperSpreaders' in params['plot_options']):
        print("Running: plot option: 'TrackSuperSpreaders'")
        dic_x_ini = nx.get_node_attributes(G,'x_ini')
        nbNodes = G.number_of_nodes()
        nbStates = params['nbStates']
        nodes = G.nodes()
        all_nodes_all_x_over_time = {}        
        S_ini = params["S_ini"]
        E_ini = params["E_ini"]
        I_ini = params["I_ini"]
        Iss_ini = params["Iss_ini"]
        R_ini = params["R_ini"]
        # because of conservation of mass, there are the same numeber 
        # of individuals all along the simulation
        N = S_ini + E_ini + I_ini + Iss_ini + R_ini 
        
        # --- first version to compute number of jumps: TOO SLOW
        #______________________________________________________________#
        
        if False:
            t0 = time.time()
            all_individuals_id = np.arange(0,N)
            
            dic_all_individuals_residence_nodes = {}
            for id_ in all_individuals_id:
                dic_all_individuals_residence_nodes[id_] = []
                for n in nodes:
                    for j in range(0,nbStates):
                        if (id_ in dic_x_ini[n][j]):
                            starting_node = n
                ref_node = starting_node
                dic_all_individuals_residence_nodes[id_].append(ref_node)            
                for i in range(0,nbTimeSteps):
                    dic_x = outputs[i][1]
                    for n in nodes:
                        for j in range(0,nbStates):
                            if (id_ in dic_x[n][j]):
                                new_node = n
                    if (new_node != ref_node):
                        print('individual ', id_, ' initially in node ',  starting_node, ' hoped from node ', ref_node, ' to node ' , new_node, ' at time ' , mytime[i])
                        ref_node = new_node
                        dic_all_individuals_residence_nodes[id_].append(ref_node)
                    
            number_of_jumps = []
            for id_ in dic_all_individuals_residence_nodes.keys():
                # -1 because the first node is the starting node 
                number_of_jumps.append(len(dic_all_individuals_residence_nodes[id_]) - 1)
            number_of_jumps = np.array(number_of_jumps)
            print(number_of_jumps)
            t1 = time.time()
            print("number of jumps ran for: ", t1 - t0, " seconds\n")
            # --- plot distribution of number of jumps per id
            plt.hist(number_of_jumps, bins=20, density=True, label='number of jumps distribution')
            plt.legend()
            plt.xlabel('Number of jumps')
            plt.ylabel('PDF')
            string_for_title='SEIIssROnGraphAgentBased, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
            plt.title(string_for_title)
            string_for_savefig = 'SEIIssROnGraphAgentBased_number_of_jumps_distribution_graph_type_'+str(params['graph_type'])+'_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
        #______________________________________________________________#
                
        t0 = time.time()
        dic_total_individual_jumps = {}        
        # create dictionnary of all individuals in every node (node is the key)
        dic_all_individuals_in_nodes_old = {}        
        dic_x_ini = nx.get_node_attributes(G,'x_ini')
        for n in nodes:
            all_individuals_in_node_old = [dic_x_ini[n][i][j] for i in range(0,len(dic_x_ini[n])) for j in range(0,len(dic_x_ini[n][i]))]
            dic_all_individuals_in_nodes_old[n] = all_individuals_in_node_old
        for i in range(0,nbTimeSteps): # start from second time step
            dic_x = outputs[i][1]
            individuals_that_jumped_at_this_time_step = []
            for n in nodes:
                all_individuals_in_node_old = dic_all_individuals_in_nodes_old[n]
                # concatenate lists of individuals
                all_individuals_in_node_new = [dic_x[n][i][j] for i in range(0,len(dic_x[n])) for j in range(0,len(dic_x[n][i]))]                
                # set is like a mathematical set
                # test if there was a jump                
                if (set(all_individuals_in_node_new) != set(all_individuals_in_node_old)):
                    if (len(set(all_individuals_in_node_new)) > len(set(all_individuals_in_node_old))):
                        individual_that_jumped = set(all_individuals_in_node_new) - set(all_individuals_in_node_old)                    
                    else:
                        individual_that_jumped = set(all_individuals_in_node_old) - set(all_individuals_in_node_new)                        
                    if (len(individual_that_jumped) != 1): 
                        print('two individuals jumped at the same time. ERROR')
                        print(individual_that_jumped)
                        exit(0)
                    else:
                        # transform set into integer                        
                        individual_that_jumped = [i for i in individual_that_jumped][0]                        
                        if (individual_that_jumped == 1818):
                            print('superspreader jumped!')
                        # initialize in case it is the first time we register a jump
                        if not(individual_that_jumped in dic_total_individual_jumps):
                            dic_total_individual_jumps[individual_that_jumped] = 0
                        
                        # prevent double jumps at the same time step:
                        # e.g. we would count that id 1 leaves node 1 and that it arrives in node 2
                        if not(individual_that_jumped in individuals_that_jumped_at_this_time_step):
                            dic_total_individual_jumps[individual_that_jumped] += 1
                            individuals_that_jumped_at_this_time_step.append(individual_that_jumped)
                            # print('individual ',individual_that_jumped, ' left/entered node ',n,' at time ' , mytime[i])
                        
                        # update new into old
                        all_individuals_in_node_old = all_individuals_in_node_new
                        dic_all_individuals_in_nodes_old[n] = all_individuals_in_node_old
                        
        number_of_jumps_alternative = np.zeros(N)
        for id_ in dic_total_individual_jumps.keys():
            number_of_jumps_alternative[id_] = dic_total_individual_jumps[id_]
        print('\n', number_of_jumps_alternative)
        t1 = time.time()        

        print('dic_total_individual_jumps=', dic_total_individual_jumps,'\n')
        # --- plot distribution of number of jumps per id
        plt.hist(number_of_jumps_alternative, bins=20, density=True, label='number of jumps distribution alternative')
        # plt.legend()
        plt.xlabel('Number of jumps')
        plt.ylabel('PDF')
        fig = plt.gcf()
        fig.tight_layout()
        string_for_savefig = 'SEIIssROnGraphAgentBased_number_of_jumps_distribution_alternatives_graph_type_'+str(params['graph_type'])\
            +'_k_SE_I='+str(params['k_SE_I'])+'_k_SE_Iss='+ str(params['k_SE_Iss'])+'_k_EI='+str(params['k_EI'])\
            +'_k_IR='+ str(params['k_IssR']) + '_S_ini='+ str(params['S_ini'])+ '_E_ini='+ str(params['E_ini'])\
            + '_I_ini='+str( params['I_ini']) + '_Iss_ini='+ str(params['Iss_ini']) + '_R_ini='+str( params['R_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
        # --- track numer of infection each individual does
        # save all id individuals that infected others        
        dic_number_of_infectious_per_individual = {}
        for i in range(0,nbTimeSteps): 
            # it means we saved the id of the infected that trigger the infection
            if (len(outputs[i]) == 3):
                infected_id = outputs[i][2]
                # intiialize in case it's the first time that this individual
                # infects someone 
                if not(infected_id in dic_number_of_infectious_per_individual):
                    dic_number_of_infectious_per_individual[infected_id] = 0
                dic_number_of_infectious_per_individual[infected_id] +=1
            
        number_of_infections_per_individual = np.zeros(len(dic_number_of_infectious_per_individual))
        ii=0
        for jj in dic_number_of_infectious_per_individual.keys():
            number_of_infections_per_individual[ii] = dic_number_of_infectious_per_individual[jj]
            ii+=1
        
        print('dic_number_of_infectious_per_individual=', dic_number_of_infectious_per_individual,'\n')
        # --- plot distribution of number of jumps per id
        plt.hist(number_of_infections_per_individual, bins=20, density=True, label='number of jumps distribution alternative')
        plt.xlabel('Number of infections per individual')
        plt.ylabel('PDF')
        fig = plt.gcf()
        fig.tight_layout()
        string_for_savefig = 'SEIIssROnGraphAgentBased_number_of_infections_per_individual_graph_type_'+str(params['graph_type'])\
              +'_k_SE_I='+str(params['k_SE_I']) +'_k_SE_Iss='+ str(params['k_SE_Iss'])+'_k_EI='+str(params['k_EI'])\
              +'_k_IR='+ str(params['k_IssR']) + '_S_ini='+ str(params['S_ini'])+ '_E_ini='+ str(params['E_ini'])\
              + '_I_ini='+str( params['I_ini']) + '_Iss_ini='+ str(params['Iss_ini']) + '_R_ini='+str( params['R_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
        # sorted_keys_dic_total_individual_jumps=sorted(dic_total_individual_jumps, key=dic_total_individual_jumps.get)
        # for _id in sorted_keys_dic_total_individual_jumps:
            # if _id in dic_number_of_infectious_per_individual:
                # num_infections_by_id = dic_number_of_infectious_per_individual[_id]
            # else:
                # num_infections_by_id = 0
            # print('id, total_jumps, number_of_infections=', _id, dic_total_individual_jumps[_id], num_infections_by_id)

        num_max_infections = int(max(number_of_infections_per_individual))
        num_jumps_vs_num_infections = np.zeros(num_max_infections+1)
        ii = 0
        # look for all number of infections from 0 to =max        
        for ii in range(0,num_max_infections+1):
            _count = 0
            for kk in dic_number_of_infectious_per_individual:
                # sorted_number_of_infections_per_individual[ii] = dic_number_of_infectious_per_individual[kk]
                # get total number of jumps of individuals who infected ii other individuals
                if (dic_number_of_infectious_per_individual[kk]==ii):
                    # check if the individual jumped                    
                    if (kk in dic_total_individual_jumps.keys()):                        
                        num_jumps_vs_num_infections[ii] += dic_total_individual_jumps[kk]
                    else:
                        num_jumps_vs_num_infections[ii] += 0
                    _count += 1
            
            if(_count >0):
                num_jumps_vs_num_infections[ii] = float(num_jumps_vs_num_infections[ii]) / _count
            else:
                num_jumps_vs_num_infections[ii] = None
            _count += 1
        
        print('\n num_jumps_vs_num_infections=',num_jumps_vs_num_infections,'\n')
        # plot number of infection vs mean total jumps
        plt.plot(np.arange(0,num_max_infections+1), num_jumps_vs_num_infections, marker='o', linestyle='')
        plt.xlabel('Number of infections per individual')
        plt.ylabel('Average total number of jumps')
        # plt.xticks(np.arange(0,num_max_infections+1,10)) # have a tick every 10 infected
        fig = plt.gcf()
        fig.tight_layout()
        string_for_savefig = 'SEIIssROnGraphAgentBased_number_of_jumps_vs_number_of_infections_graph_type_'+str(params['graph_type'])\
              +'_k_SE_I='+str(params['k_SE_I']) +'_k_SE_Iss='+ str(params['k_SE_Iss'])+'_k_EI='+str(params['k_EI'])\
              +'_k_IR='+ str(params['k_IssR']) + '_S_ini='+ str(params['S_ini'])+ '_E_ini='+ str(params['E_ini'])\
              + '_I_ini='+str( params['I_ini']) + '_Iss_ini='+ str(params['Iss_ini']) + '_R_ini='+str( params['R_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()