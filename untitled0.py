def dichotomie(t, v):
   a = 0
   b = len(t) - 1
   while 1 < b - a:
      m = (a + b) // 2
      if t[m] < v:
         a = m 
      else:
         b = m 
   if(abs(v - t[a]) <= abs(v - t[b])):
      return a
   else:
      return b
  
                for j in range(1, len(next_states_sum)):
                    if (next_states_sum[j - 1] < xi) and (xi < next_states_sum[j]):
                        id_reaction = j
                        
                a = 0
                b = len(next_states_sum) - 1
                while 1 < b - a:
                    m = (a + b) // 2
                    if next_states_sum[m] < xi:
                        a = m
                    else:
                        b = m 
                id_reaction = b
                
            else: # we set a hopping process
                # print("hopping")                
                kk = nbReactions*len(G.nodes())
                kk_init = kk
                for ni in G.nodes(): # loop over all nodes
                    for nj in G.neighbors(ni): # loop over all neighbors
                        if (ni != nj):
                            for i in range(0, nbStates): # loop over all states                                
                                if (kk==id_reaction):
                                    if (i > 0):
                                        print('a hopping happens from someone who IS NOT susceptible i = ' , i)                                    
                                    dic_x[ni][i] -= 1 # add one
                                    dic_x[nj][i] += 1 # remove one
                                    # update attributes right away
                                    nx.set_node_attributes(G, name='x', values=dic_x)
                                    if ((dic_x[ni][i] < 0) or (dic_x[nj][i] < 0)):
                                        print('negative hopping: id_reaction, ni, nj, i, kk, dic_x[ni][i], dic_x[nj][i]  = ' , id_reaction,  ni, nj, i, kk, dic_x[ni][i], dic_x[nj][i])
                                        print('dic_t_between[(',ni,',',nj,')][',i,' = ' , dic_t_between[(ni,nj)][i])
                                        print('dic_t_between[(',nj,',',ni,')][',i,' = ' , dic_t_between[(nj,ni)][i])
                                        print('dic_t_between = ' , dic_t_between)
                                        exit(-1)
                                kk+=1

            else: # we set a hopping process
                # print("hopping")                
                kk_init = nbReactions*len(G.nodes())
                                
                ni, nj, i = tab[id_reaction - kk_init][0], tab[id_reaction - kk_init][1], tab[id_reaction - kk_init][2]
                if(ni != nj):
                    if(i > 0):
                        print('a hopping happens from someone who IS NOT susceptible i = ' , i)   
                    dic_x[ni][i] -= 1
                    dic_x[nj][i] += 1
                    
                    nx.set_node_attributes(G, name='x', values=dic_x)
                    if ((dic_x[ni][i] < 0) or (dic_x[nj][i] < 0)):
                        print('negative hopping: id_reaction, ni, nj, i, kk, dic_x[ni][i], dic_x[nj][i]  = ' , id_reaction,  ni, nj, i, kk, dic_x[ni][i], dic_x[nj][i])
                        print('dic_t_between[(',ni,',',nj,')][',i,' = ' , dic_t_between[(ni,nj)][i])
                        print('dic_t_between[(',nj,',',ni,')][',i,' = ' , dic_t_between[(nj,ni)][i])
                        print('dic_t_between = ' , dic_t_between)
                        exit(-1)  
                        
            tab = [] # will allow us to make the hop more easily
                            tab.append((ni, nj, i))