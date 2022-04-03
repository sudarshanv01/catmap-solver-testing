
def rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,mpf,matrix,mpexp,mpsqrt,include_derivatives=True):
    
    def interaction_function(coverages,energies,epsilon,F,include_derivatives=True,include_integral=False):
    
        ##this dictionary is passed in via the "template" so that it can be compiled
        site_info_dict = {'s': [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1, {'slope': 1.0, 'cutoff': 0.66, 'smoothing': 0.05}], 's&h': [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1.0, {'slope': 1.0, 'cutoff': 2.0, 'smoothing': 0.05}], 'h': [[0], 1.0, {'slope': 1.0, 'cutoff': 2.0, 'smoothing': 0.05}]}
    
        N_ads = len(coverages)
        single_sites = [s for s in site_info_dict if '&' not in s]
        N_sites = len(single_sites)
    
        idx_lists = []
        f = [[0]*N_sites for x in range(N_sites)]
        df = [[0]*N_sites for x in range(N_sites)]
        d2f = [[0]*N_sites for x in range(N_sites)]
    
        ##form matrix of f, df, d2f for each site and cross-site
        for si,s in enumerate(single_sites):
            for qi,q in enumerate(single_sites):
                if s == q:
                    idxs,max_cvg,F_params = site_info_dict[s]
                    idx_lists.append(site_info_dict[s][0])
                    theta_tot = 2*sum([coverages[i] for i in idxs])
                else:
                    key1 = '&'.join([s,q])
                    key2 = '&'.join([q,s])
                    if key1 in site_info_dict:
                        idxs,max_cvg,F_params = site_info_dict[key1]
                    elif key2 in site_info_dict:
                        idxs,max_cvg,F_params = site_info_dict[key2]
                    else:
                        raise UserWarning(
                        ('No cross-site interactions'
                         ' specified for {s},{q}')
                          .format(**locals()))
                    theta_tot = sum([coverages[i] for i in idxs])
    
                fs,dfs,d2fs = F(theta_tot,**F_params)
                f[si][qi] = f[qi][si] = fs
                df[si][qi] = df[qi][si] = dfs
                d2f[si][qi] = d2f[qi][si] = d2fs
    
        #initiate terms for first derivative
        term_1 = [0]*N_ads
        term_2 = [0]*N_ads
        
        #initate intermediate quantities
        eps_theta_theta = [[0]*N_sites for x in range(N_sites)]
        eps_theta = [[0]*N_sites for x in range(N_ads)]
        site_lookup = [0]*N_ads
    
        #form site_lookup and eps_theta_theta matrix.
        #site_lookup used to avoid empty sums over sites
        #eps_theta_theta used for term_2 and jacobian
        for s in range(N_sites):
            for i in idx_lists[s]:
                if i in idx_lists[s]:
                    site_lookup[i] = s
                for q in range(N_sites):
                    for j in idx_lists[q]:
                        ep_t_t = epsilon[i*N_ads+j]*coverages[i]*coverages[j]
                        eps_theta_theta[s][q] += ep_t_t
    
        #form term_1 and eps_theta matrix (eps_theta needed for jacobian)
        for k in range(N_ads):
            q_k = site_lookup[k]
            for s in range(N_sites):
                for i in idx_lists[s]:
                    eps_theta[k][s] += epsilon[i*N_ads+k]*coverages[i]
                term_1[k] += ((f[s][q_k])**2)*eps_theta[k][s]
    
        #form term_2
        for k in range(N_ads):
            q_k = site_lookup[k]
            for s in range(N_sites):
                term_2[k] += 2*f[s][q_k]*df[s][q_k]*eps_theta_theta[s][q_k]
    
        #combine terms with constant energy to give E_diff
        E_diff = [a+b+c for a,b,c in zip(energies,term_1,term_2)]
    
        if include_derivatives:
            #compute the jacobian
            E_jacob = [[0]*N_ads for x in range(N_ads)]
    
            for k in range(N_ads):
                for l in range(N_ads):
                    s_l = site_lookup[l]
                    q_k = site_lookup[k]
                    f_sq = f[s_l][q_k]
                    df_sq = df[s_l][q_k]
                    #obtained from tensor calculus, checked against numerical derivative
                    E_jacob[l][k] += (     \
                             (f_sq**2)*epsilon[l*N_ads+k] +\
                             2*(f_sq*df_sq*(eps_theta[k][s_l]+eps_theta[l][q_k]) + \
                             eps_theta_theta[s_l][q_k]*((df_sq**2) + (f_sq*d2f[s_l][q_k]))))
                    
                    #term is only included if l and k are on the same site
                    if l in idx_lists[q_k]:
                        for s in range(N_sites):
                            E_jacob[l][k] += 2*(     \
                                 f[s][q_k]*df[s][q_k]*(eps_theta[k][s]+eps_theta[l][s]) + \
                                 eps_theta_theta[s][q_k]*((df[s][q_k]**2) + \
                                                           f[s][q_k]*d2f[s][q_k]))
        else:
            E_jacob = None
    
        if include_integral:
            E = 0
            for i in range(N_ads):
                E += energies[i]*coverages[i]
    
            for s in range(N_sites):
                for q in range(N_sites):
                    for i in idx_lists[s]:
                        for j in idx_lists[q]:
                            E += 0.5*(f[s][q]**2)*epsilon[i*N_ads+j]*coverages[i]*coverages[j]
        else:
            E = 0
    
    
        return E, E_diff, E_jacob
    

    kfs = []
    krs = []
    dEfs = []
    dErs = []

    n_adsorbates = 8
    n_transition_states = 7
    n_tot = n_adsorbates+n_transition_states
    # Account for the numbers solver where the extra coverage 
    # of the clean slab might be needed
    theta = theta[:n_adsorbates]
    if len(theta) == n_adsorbates:
        theta = list(theta) + [0]*n_transition_states #put transition-state coverages to 0
    elif len(theta) != n_adsorbates+n_transition_states:
        raise ValueError('Coverage vector was not correct length')
    energies = rxn_parameters[:n_tot]
    if len(rxn_parameters) == n_tot + n_tot**2:
        interaction_vector = rxn_parameters[-n_tot**2:]
    elif len(rxn_parameters) == n_tot:
        interaction_vector = [0]*n_tot**2
    else:
        raise ValueError('Length of reaction parameters is not correct. '+ str(rxn_parameters))

    G_int, Gf,dGs =  interaction_function(theta,energies,interaction_vector,F,include_derivatives=include_derivatives, include_integral=False)

    kB = mpf('0.00008617332478000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000029')
    h = mpf('0.000000000000004135667516000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000012')
    prefactor_list = [kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h]

    def element_wise_addition(lists):
        return [sum(L) for L in zip(*lists)]

    G_IS = [0]*9
    G_TS = [0]*9
    G_FS = [0]*9
    G_af = [0]*9
    G_ar = [0]*9
    dG_IS = [0]*9
    dG_TS = [0]*9
    dG_FS = [0]*9
    G_IS[0] = gas_energies[3] + site_energies[2] + site_energies[2]
    G_FS[0] = Gf[0] + Gf[0]
    dG_IS[0] = element_wise_addition([[0]*8])
    dG_FS[0] = element_wise_addition([dGs[0] , dGs[0]])
    G_TS[0] = max([G_IS[0],G_FS[0]])
    dG_TS[0] = None #determined later
    
    G_IS[1] = site_energies[0] + gas_energies[1]
    G_FS[1] = Gf[6]
    dG_IS[1] = element_wise_addition([[0]*8])
    dG_FS[1] = element_wise_addition([dGs[6]])
    G_TS[1] = max([G_IS[1],G_FS[1]])
    dG_TS[1] = None #determined later
    
    G_IS[2] = Gf[6] + Gf[0]
    G_FS[2] = Gf[4] + site_energies[2]
    dG_IS[2] = element_wise_addition([dGs[6] , dGs[0]])
    dG_FS[2] = element_wise_addition([dGs[4]])
    G_TS[2] = Gf[12] + site_energies[2]
    dG_TS[2] = element_wise_addition([dGs[12]])
    
    G_IS[3] = Gf[4] + Gf[0]
    G_FS[3] = Gf[3] + site_energies[2]
    dG_IS[3] = element_wise_addition([dGs[4] , dGs[0]])
    dG_FS[3] = element_wise_addition([dGs[3]])
    G_TS[3] = Gf[14] + site_energies[2]
    dG_TS[3] = element_wise_addition([dGs[14]])
    
    G_IS[4] = Gf[3] + site_energies[0]
    G_FS[4] = Gf[5] + Gf[7]
    dG_IS[4] = element_wise_addition([dGs[3]])
    dG_FS[4] = element_wise_addition([dGs[5] , dGs[7]])
    G_TS[4] = Gf[9] + site_energies[0]
    dG_TS[4] = element_wise_addition([dGs[9]])
    
    G_IS[5] = Gf[5] + Gf[0]
    G_FS[5] = Gf[1] + site_energies[2]
    dG_IS[5] = element_wise_addition([dGs[5] , dGs[0]])
    dG_FS[5] = element_wise_addition([dGs[1]])
    G_TS[5] = Gf[8] + site_energies[2]
    dG_TS[5] = element_wise_addition([dGs[8]])
    
    G_IS[6] = Gf[1] + Gf[0]
    G_FS[6] = Gf[2] + site_energies[2]
    dG_IS[6] = element_wise_addition([dGs[1] , dGs[0]])
    dG_FS[6] = element_wise_addition([dGs[2]])
    G_TS[6] = Gf[10] + site_energies[2]
    dG_TS[6] = element_wise_addition([dGs[10]])
    
    G_IS[7] = Gf[2] + Gf[0]
    G_FS[7] = gas_energies[0] + site_energies[0] + site_energies[2]
    dG_IS[7] = element_wise_addition([dGs[2] , dGs[0]])
    dG_FS[7] = element_wise_addition([[0]*8])
    G_TS[7] = Gf[11] + site_energies[2]
    dG_TS[7] = element_wise_addition([dGs[11]])
    
    G_IS[8] = Gf[7] + Gf[0]
    G_FS[8] = gas_energies[2] + site_energies[0] + site_energies[2]
    dG_IS[8] = element_wise_addition([dGs[7] , dGs[0]])
    dG_FS[8] = element_wise_addition([[0]*8])
    G_TS[8] = Gf[13] + site_energies[2]
    dG_TS[8] = element_wise_addition([dGs[13]])
    

    n_rxns = len(G_IS)
    for i in range(n_rxns):
        G_list = [G_IS[i],G_FS[i],G_TS[i]]
        G_TS_i = max(G_list)
        max_idx = G_list.index(G_TS_i)
        G_af[i] = G_TS_i - G_IS[i]
        G_ar[i] = G_TS_i - G_FS[i]
        kf = prefactor_list[i]*mpexp(-G_af[i]/(kB*T))
        kr = prefactor_list[i]*mpexp(-G_ar[i]/(kB*T))
        kfs.append(kf)
        krs.append(kr)
        if include_derivatives:
            dG_TS_i = [dG_IS[i],dG_FS[i],dG_TS[i]][max_idx]

            dG_TS_i = list(dG_TS_i[:n_adsorbates])
            dG_IS_i = list(dG_IS[i][:n_adsorbates])
            dG_FS_i = list(dG_FS[i][:n_adsorbates])

            dEf = [(-dG_TS_j + dG_IS_j) for dG_TS_j,dG_IS_j in zip(dG_TS_i,dG_IS_i)]
            dEr = [(-dG_TS_j + dG_FS_j) for dG_TS_j,dG_FS_j in zip(dG_TS_i,dG_FS_i)]
            dEfs.append(dEf)
            dErs.append(dEr)
        else:
            dEfs = dErs = None
    return kfs, krs, dEfs, dErs




def interaction_function(coverages,energies,epsilon,F,include_derivatives=True,include_integral=False):

    ##this dictionary is passed in via the "template" so that it can be compiled
    site_info_dict = {'s': [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1, {'slope': 1.0, 'cutoff': 0.66, 'smoothing': 0.05}], 's&h': [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1.0, {'slope': 1.0, 'cutoff': 2.0, 'smoothing': 0.05}], 'h': [[0], 1.0, {'slope': 1.0, 'cutoff': 2.0, 'smoothing': 0.05}]}

    N_ads = len(coverages)
    single_sites = [s for s in site_info_dict if '&' not in s]
    N_sites = len(single_sites)

    idx_lists = []
    f = [[0]*N_sites for x in range(N_sites)]
    df = [[0]*N_sites for x in range(N_sites)]
    d2f = [[0]*N_sites for x in range(N_sites)]

    ##form matrix of f, df, d2f for each site and cross-site
    for si,s in enumerate(single_sites):
        for qi,q in enumerate(single_sites):
            if s == q:
                idxs,max_cvg,F_params = site_info_dict[s]
                idx_lists.append(site_info_dict[s][0])
                theta_tot = 2*sum([coverages[i] for i in idxs])
            else:
                key1 = '&'.join([s,q])
                key2 = '&'.join([q,s])
                if key1 in site_info_dict:
                    idxs,max_cvg,F_params = site_info_dict[key1]
                elif key2 in site_info_dict:
                    idxs,max_cvg,F_params = site_info_dict[key2]
                else:
                    raise UserWarning(
                    ('No cross-site interactions'
                     ' specified for {s},{q}')
                      .format(**locals()))
                theta_tot = sum([coverages[i] for i in idxs])

            fs,dfs,d2fs = F(theta_tot,**F_params)
            f[si][qi] = f[qi][si] = fs
            df[si][qi] = df[qi][si] = dfs
            d2f[si][qi] = d2f[qi][si] = d2fs

    #initiate terms for first derivative
    term_1 = [0]*N_ads
    term_2 = [0]*N_ads
    
    #initate intermediate quantities
    eps_theta_theta = [[0]*N_sites for x in range(N_sites)]
    eps_theta = [[0]*N_sites for x in range(N_ads)]
    site_lookup = [0]*N_ads

    #form site_lookup and eps_theta_theta matrix.
    #site_lookup used to avoid empty sums over sites
    #eps_theta_theta used for term_2 and jacobian
    for s in range(N_sites):
        for i in idx_lists[s]:
            if i in idx_lists[s]:
                site_lookup[i] = s
            for q in range(N_sites):
                for j in idx_lists[q]:
                    ep_t_t = epsilon[i*N_ads+j]*coverages[i]*coverages[j]
                    eps_theta_theta[s][q] += ep_t_t

    #form term_1 and eps_theta matrix (eps_theta needed for jacobian)
    for k in range(N_ads):
        q_k = site_lookup[k]
        for s in range(N_sites):
            for i in idx_lists[s]:
                eps_theta[k][s] += epsilon[i*N_ads+k]*coverages[i]
            term_1[k] += ((f[s][q_k])**2)*eps_theta[k][s]

    #form term_2
    for k in range(N_ads):
        q_k = site_lookup[k]
        for s in range(N_sites):
            term_2[k] += 2*f[s][q_k]*df[s][q_k]*eps_theta_theta[s][q_k]

    #combine terms with constant energy to give E_diff
    E_diff = [a+b+c for a,b,c in zip(energies,term_1,term_2)]

    if include_derivatives:
        #compute the jacobian
        E_jacob = [[0]*N_ads for x in range(N_ads)]

        for k in range(N_ads):
            for l in range(N_ads):
                s_l = site_lookup[l]
                q_k = site_lookup[k]
                f_sq = f[s_l][q_k]
                df_sq = df[s_l][q_k]
                #obtained from tensor calculus, checked against numerical derivative
                E_jacob[l][k] += (     \
                         (f_sq**2)*epsilon[l*N_ads+k] +\
                         2*(f_sq*df_sq*(eps_theta[k][s_l]+eps_theta[l][q_k]) + \
                         eps_theta_theta[s_l][q_k]*((df_sq**2) + (f_sq*d2f[s_l][q_k]))))
                
                #term is only included if l and k are on the same site
                if l in idx_lists[q_k]:
                    for s in range(N_sites):
                        E_jacob[l][k] += 2*(     \
                             f[s][q_k]*df[s][q_k]*(eps_theta[k][s]+eps_theta[l][s]) + \
                             eps_theta_theta[s][q_k]*((df[s][q_k]**2) + \
                                                       f[s][q_k]*d2f[s][q_k]))
    else:
        E_jacob = None

    if include_integral:
        E = 0
        for i in range(N_ads):
            E += energies[i]*coverages[i]

        for s in range(N_sites):
            for q in range(N_sites):
                for i in idx_lists[s]:
                    for j in idx_lists[q]:
                        E += 0.5*(f[s][q]**2)*epsilon[i*N_ads+j]*coverages[i]*coverages[j]
    else:
        E = 0


    return E, E_diff, E_jacob




def interacting_mean_field_steady_state(rxn_parameters,theta,p,gas_energies,site_energies,T,F,mpf,matrix,mpexp,mpsqrt):

    
    def rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,mpf,matrix,mpexp,mpsqrt,include_derivatives=True):
        
        def interaction_function(coverages,energies,epsilon,F,include_derivatives=True,include_integral=False):
        
            ##this dictionary is passed in via the "template" so that it can be compiled
            site_info_dict = {'s': [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1, {'slope': 1.0, 'cutoff': 0.66, 'smoothing': 0.05}], 's&h': [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1.0, {'slope': 1.0, 'cutoff': 2.0, 'smoothing': 0.05}], 'h': [[0], 1.0, {'slope': 1.0, 'cutoff': 2.0, 'smoothing': 0.05}]}
        
            N_ads = len(coverages)
            single_sites = [s for s in site_info_dict if '&' not in s]
            N_sites = len(single_sites)
        
            idx_lists = []
            f = [[0]*N_sites for x in range(N_sites)]
            df = [[0]*N_sites for x in range(N_sites)]
            d2f = [[0]*N_sites for x in range(N_sites)]
        
            ##form matrix of f, df, d2f for each site and cross-site
            for si,s in enumerate(single_sites):
                for qi,q in enumerate(single_sites):
                    if s == q:
                        idxs,max_cvg,F_params = site_info_dict[s]
                        idx_lists.append(site_info_dict[s][0])
                        theta_tot = 2*sum([coverages[i] for i in idxs])
                    else:
                        key1 = '&'.join([s,q])
                        key2 = '&'.join([q,s])
                        if key1 in site_info_dict:
                            idxs,max_cvg,F_params = site_info_dict[key1]
                        elif key2 in site_info_dict:
                            idxs,max_cvg,F_params = site_info_dict[key2]
                        else:
                            raise UserWarning(
                            ('No cross-site interactions'
                             ' specified for {s},{q}')
                              .format(**locals()))
                        theta_tot = sum([coverages[i] for i in idxs])
        
                    fs,dfs,d2fs = F(theta_tot,**F_params)
                    f[si][qi] = f[qi][si] = fs
                    df[si][qi] = df[qi][si] = dfs
                    d2f[si][qi] = d2f[qi][si] = d2fs
        
            #initiate terms for first derivative
            term_1 = [0]*N_ads
            term_2 = [0]*N_ads
            
            #initate intermediate quantities
            eps_theta_theta = [[0]*N_sites for x in range(N_sites)]
            eps_theta = [[0]*N_sites for x in range(N_ads)]
            site_lookup = [0]*N_ads
        
            #form site_lookup and eps_theta_theta matrix.
            #site_lookup used to avoid empty sums over sites
            #eps_theta_theta used for term_2 and jacobian
            for s in range(N_sites):
                for i in idx_lists[s]:
                    if i in idx_lists[s]:
                        site_lookup[i] = s
                    for q in range(N_sites):
                        for j in idx_lists[q]:
                            ep_t_t = epsilon[i*N_ads+j]*coverages[i]*coverages[j]
                            eps_theta_theta[s][q] += ep_t_t
        
            #form term_1 and eps_theta matrix (eps_theta needed for jacobian)
            for k in range(N_ads):
                q_k = site_lookup[k]
                for s in range(N_sites):
                    for i in idx_lists[s]:
                        eps_theta[k][s] += epsilon[i*N_ads+k]*coverages[i]
                    term_1[k] += ((f[s][q_k])**2)*eps_theta[k][s]
        
            #form term_2
            for k in range(N_ads):
                q_k = site_lookup[k]
                for s in range(N_sites):
                    term_2[k] += 2*f[s][q_k]*df[s][q_k]*eps_theta_theta[s][q_k]
        
            #combine terms with constant energy to give E_diff
            E_diff = [a+b+c for a,b,c in zip(energies,term_1,term_2)]
        
            if include_derivatives:
                #compute the jacobian
                E_jacob = [[0]*N_ads for x in range(N_ads)]
        
                for k in range(N_ads):
                    for l in range(N_ads):
                        s_l = site_lookup[l]
                        q_k = site_lookup[k]
                        f_sq = f[s_l][q_k]
                        df_sq = df[s_l][q_k]
                        #obtained from tensor calculus, checked against numerical derivative
                        E_jacob[l][k] += (     \
                                 (f_sq**2)*epsilon[l*N_ads+k] +\
                                 2*(f_sq*df_sq*(eps_theta[k][s_l]+eps_theta[l][q_k]) + \
                                 eps_theta_theta[s_l][q_k]*((df_sq**2) + (f_sq*d2f[s_l][q_k]))))
                        
                        #term is only included if l and k are on the same site
                        if l in idx_lists[q_k]:
                            for s in range(N_sites):
                                E_jacob[l][k] += 2*(     \
                                     f[s][q_k]*df[s][q_k]*(eps_theta[k][s]+eps_theta[l][s]) + \
                                     eps_theta_theta[s][q_k]*((df[s][q_k]**2) + \
                                                               f[s][q_k]*d2f[s][q_k]))
            else:
                E_jacob = None
        
            if include_integral:
                E = 0
                for i in range(N_ads):
                    E += energies[i]*coverages[i]
        
                for s in range(N_sites):
                    for q in range(N_sites):
                        for i in idx_lists[s]:
                            for j in idx_lists[q]:
                                E += 0.5*(f[s][q]**2)*epsilon[i*N_ads+j]*coverages[i]*coverages[j]
            else:
                E = 0
        
        
            return E, E_diff, E_jacob
        
    
        kfs = []
        krs = []
        dEfs = []
        dErs = []
    
        n_adsorbates = 8
        n_transition_states = 7
        n_tot = n_adsorbates+n_transition_states
        # Account for the numbers solver where the extra coverage 
        # of the clean slab might be needed
        theta = theta[:n_adsorbates]
        if len(theta) == n_adsorbates:
            theta = list(theta) + [0]*n_transition_states #put transition-state coverages to 0
        elif len(theta) != n_adsorbates+n_transition_states:
            raise ValueError('Coverage vector was not correct length')
        energies = rxn_parameters[:n_tot]
        if len(rxn_parameters) == n_tot + n_tot**2:
            interaction_vector = rxn_parameters[-n_tot**2:]
        elif len(rxn_parameters) == n_tot:
            interaction_vector = [0]*n_tot**2
        else:
            raise ValueError('Length of reaction parameters is not correct. '+ str(rxn_parameters))
    
        G_int, Gf,dGs =  interaction_function(theta,energies,interaction_vector,F,include_derivatives=include_derivatives, include_integral=False)
    
        kB = mpf('0.00008617332478000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000029')
        h = mpf('0.000000000000004135667516000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000012')
        prefactor_list = [kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h]
    
        def element_wise_addition(lists):
            return [sum(L) for L in zip(*lists)]
    
        G_IS = [0]*9
        G_TS = [0]*9
        G_FS = [0]*9
        G_af = [0]*9
        G_ar = [0]*9
        G_IS[0] = gas_energies[3] + site_energies[2] + site_energies[2]
        G_FS[0] = Gf[0] + Gf[0]
        G_TS[0] = max([G_IS[0],G_FS[0]])
        
        G_IS[1] = site_energies[0] + gas_energies[1]
        G_FS[1] = Gf[6]
        G_TS[1] = max([G_IS[1],G_FS[1]])
        
        G_IS[2] = Gf[6] + Gf[0]
        G_FS[2] = Gf[4] + site_energies[2]
        G_TS[2] = Gf[12] + site_energies[2]
        
        G_IS[3] = Gf[4] + Gf[0]
        G_FS[3] = Gf[3] + site_energies[2]
        G_TS[3] = Gf[14] + site_energies[2]
        
        G_IS[4] = Gf[3] + site_energies[0]
        G_FS[4] = Gf[5] + Gf[7]
        G_TS[4] = Gf[9] + site_energies[0]
        
        G_IS[5] = Gf[5] + Gf[0]
        G_FS[5] = Gf[1] + site_energies[2]
        G_TS[5] = Gf[8] + site_energies[2]
        
        G_IS[6] = Gf[1] + Gf[0]
        G_FS[6] = Gf[2] + site_energies[2]
        G_TS[6] = Gf[10] + site_energies[2]
        
        G_IS[7] = Gf[2] + Gf[0]
        G_FS[7] = gas_energies[0] + site_energies[0] + site_energies[2]
        G_TS[7] = Gf[11] + site_energies[2]
        
        G_IS[8] = Gf[7] + Gf[0]
        G_FS[8] = gas_energies[2] + site_energies[0] + site_energies[2]
        G_TS[8] = Gf[13] + site_energies[2]
        
    
        n_rxns = len(G_IS)
        for i in range(n_rxns):
            G_list = [G_IS[i],G_FS[i],G_TS[i]]
            G_TS_i = max(G_list)
            max_idx = G_list.index(G_TS_i)
            G_af[i] = G_TS_i - G_IS[i]
            G_ar[i] = G_TS_i - G_FS[i]
            kf = prefactor_list[i]*mpexp(-G_af[i]/(kB*T))
            kr = prefactor_list[i]*mpexp(-G_ar[i]/(kB*T))
            kfs.append(kf)
            krs.append(kr)
            if include_derivatives:
                dG_TS_i = [dG_IS[i],dG_FS[i],dG_TS[i]][max_idx]
    
                dG_TS_i = list(dG_TS_i[:n_adsorbates])
                dG_IS_i = list(dG_IS[i][:n_adsorbates])
                dG_FS_i = list(dG_FS[i][:n_adsorbates])
    
                dEf = [(-dG_TS_j + dG_IS_j) for dG_TS_j,dG_IS_j in zip(dG_TS_i,dG_IS_i)]
                dEr = [(-dG_TS_j + dG_FS_j) for dG_TS_j,dG_FS_j in zip(dG_TS_i,dG_FS_i)]
                dEfs.append(dEf)
                dErs.append(dEr)
            else:
                dEfs = dErs = None
        return kfs, krs, dEfs, dErs
    
    
    kf, kr, junk, junk= rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,
            mpf,matrix,mpexp,mpsqrt,include_derivatives=False)
   
    r = [0]*len(kf)
    dtheta_dt = [0]*10
    
    s = [0]*3
    s[0] = (mpf('1.0') - theta[1] - theta[2] - theta[3] - theta[4] - theta[5] - theta[6] - theta[7])
    s[1] = (mpf('0.0'))
    s[2] = (mpf('1.0') - theta[0])
    r[0] = kf[0]*p[3]*theta[9]*theta[9] - kr[0]*theta[0]*theta[0]
    r[1] = kf[1]*p[1]*theta[8] - kr[1]*theta[6]
    r[2] = kf[2]*theta[6]*theta[0] - kr[2]*theta[4]*theta[9]
    r[3] = kf[3]*theta[4]*theta[0] - kr[3]*theta[3]*theta[9]
    r[4] = kf[4]*theta[3]*theta[8] - kr[4]*theta[5]*theta[7]
    r[5] = kf[5]*theta[5]*theta[0] - kr[5]*theta[1]*theta[9]
    r[6] = kf[6]*theta[1]*theta[0] - kr[6]*theta[2]*theta[9]
    r[7] = kf[7]*theta[2]*theta[0] - kr[7]*p[0]*theta[8]*theta[9]
    r[8] = kf[8]*theta[7]*theta[0] - kr[8]*p[2]*theta[8]*theta[9]
    dtheta_dt[0] =  + 2*r[0] + -1*r[2] + -1*r[3] + -1*r[5] + -1*r[6] + -1*r[7] + -1*r[8]
    dtheta_dt[1] =  + 1*r[5] + -1*r[6]
    dtheta_dt[2] =  + 1*r[6] + -1*r[7]
    dtheta_dt[3] =  + 1*r[3] + -1*r[4]
    dtheta_dt[4] =  + 1*r[2] + -1*r[3]
    dtheta_dt[5] =  + 1*r[4] + -1*r[5]
    dtheta_dt[6] =  + 1*r[1] + -1*r[2]
    dtheta_dt[7] =  + 1*r[4] + -1*r[8]
    dtheta_dt[8] =  + -1*r[1] + -1*r[4] + 1*r[7] + 1*r[8]
    dtheta_dt[9] =  + -2*r[0] + 1*r[2] + 1*r[3] + 1*r[5] + 1*r[6] + 1*r[7] + 1*r[8]
 
    r = matrix(r)
    dtheta_dt = matrix(dtheta_dt)

    if dtheta_dt.rows == len(theta)+1:
        dtheta_dt = dtheta_dt[:-1] 

    return dtheta_dt




def ideal_mean_field_steady_state(kf,kr,theta,p,mpf,matrix):

    r = [0]*len(kf)
    dtheta_dt = [0]*10
    s = [0]*3
    s[0] = (mpf('1.0') - theta[1] - theta[2] - theta[3] - theta[4] - theta[5] - theta[6] - theta[7])
    s[1] = (mpf('0.0'))
    s[2] = (mpf('1.0') - theta[0])
    r[0] = kf[0]*p[3]*theta[9]*theta[9] - kr[0]*theta[0]*theta[0]
    r[1] = kf[1]*p[1]*theta[8] - kr[1]*theta[6]
    r[2] = kf[2]*theta[6]*theta[0] - kr[2]*theta[4]*theta[9]
    r[3] = kf[3]*theta[4]*theta[0] - kr[3]*theta[3]*theta[9]
    r[4] = kf[4]*theta[3]*theta[8] - kr[4]*theta[5]*theta[7]
    r[5] = kf[5]*theta[5]*theta[0] - kr[5]*theta[1]*theta[9]
    r[6] = kf[6]*theta[1]*theta[0] - kr[6]*theta[2]*theta[9]
    r[7] = kf[7]*theta[2]*theta[0] - kr[7]*p[0]*theta[8]*theta[9]
    r[8] = kf[8]*theta[7]*theta[0] - kr[8]*p[2]*theta[8]*theta[9]
    dtheta_dt[0] =  + 2*r[0] + -1*r[2] + -1*r[3] + -1*r[5] + -1*r[6] + -1*r[7] + -1*r[8]
    dtheta_dt[1] =  + 1*r[5] + -1*r[6]
    dtheta_dt[2] =  + 1*r[6] + -1*r[7]
    dtheta_dt[3] =  + 1*r[3] + -1*r[4]
    dtheta_dt[4] =  + 1*r[2] + -1*r[3]
    dtheta_dt[5] =  + 1*r[4] + -1*r[5]
    dtheta_dt[6] =  + 1*r[1] + -1*r[2]
    dtheta_dt[7] =  + 1*r[4] + -1*r[8]
    dtheta_dt[8] =  + -1*r[1] + -1*r[4] + 1*r[7] + 1*r[8]
    dtheta_dt[9] =  + -2*r[0] + 1*r[2] + 1*r[3] + 1*r[5] + 1*r[6] + 1*r[7] + 1*r[8]

    r = matrix(r)
    dtheta_dt = matrix(dtheta_dt)

    if dtheta_dt.rows == len(theta)+1:
        dtheta_dt = dtheta_dt[:-1] 
    
    return dtheta_dt




def interacting_mean_field_jacobian(rxn_parameters,theta,p,gas_energies,site_energies,T,F,mpf,matrix,mpexp,mpsqrt):

#    print 'rxn_parameters = ', rxn_parameters
#    print 'theta = ', theta
#    print 'p = ', p
#    print 'gas_energies = ', gas_energies
#    print 'site_energies = ', site_energies
#    print 'T = ', T
    kB = mpf('0.00008617332478000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000029')
    n_adsorbates = 8
    kBT = kB*T


    J = [[0 for i in range(10)] for j in range(10)]

    
    def rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,mpf,matrix,mpexp,mpsqrt,include_derivatives=True):
        
        def interaction_function(coverages,energies,epsilon,F,include_derivatives=True,include_integral=False):
        
            ##this dictionary is passed in via the "template" so that it can be compiled
            site_info_dict = {'s': [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1, {'slope': 1.0, 'cutoff': 0.66, 'smoothing': 0.05}], 's&h': [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1.0, {'slope': 1.0, 'cutoff': 2.0, 'smoothing': 0.05}], 'h': [[0], 1.0, {'slope': 1.0, 'cutoff': 2.0, 'smoothing': 0.05}]}
        
            N_ads = len(coverages)
            single_sites = [s for s in site_info_dict if '&' not in s]
            N_sites = len(single_sites)
        
            idx_lists = []
            f = [[0]*N_sites for x in range(N_sites)]
            df = [[0]*N_sites for x in range(N_sites)]
            d2f = [[0]*N_sites for x in range(N_sites)]
        
            ##form matrix of f, df, d2f for each site and cross-site
            for si,s in enumerate(single_sites):
                for qi,q in enumerate(single_sites):
                    if s == q:
                        idxs,max_cvg,F_params = site_info_dict[s]
                        idx_lists.append(site_info_dict[s][0])
                        theta_tot = 2*sum([coverages[i] for i in idxs])
                    else:
                        key1 = '&'.join([s,q])
                        key2 = '&'.join([q,s])
                        if key1 in site_info_dict:
                            idxs,max_cvg,F_params = site_info_dict[key1]
                        elif key2 in site_info_dict:
                            idxs,max_cvg,F_params = site_info_dict[key2]
                        else:
                            raise UserWarning(
                            ('No cross-site interactions'
                             ' specified for {s},{q}')
                              .format(**locals()))
                        theta_tot = sum([coverages[i] for i in idxs])
        
                    fs,dfs,d2fs = F(theta_tot,**F_params)
                    f[si][qi] = f[qi][si] = fs
                    df[si][qi] = df[qi][si] = dfs
                    d2f[si][qi] = d2f[qi][si] = d2fs
        
            #initiate terms for first derivative
            term_1 = [0]*N_ads
            term_2 = [0]*N_ads
            
            #initate intermediate quantities
            eps_theta_theta = [[0]*N_sites for x in range(N_sites)]
            eps_theta = [[0]*N_sites for x in range(N_ads)]
            site_lookup = [0]*N_ads
        
            #form site_lookup and eps_theta_theta matrix.
            #site_lookup used to avoid empty sums over sites
            #eps_theta_theta used for term_2 and jacobian
            for s in range(N_sites):
                for i in idx_lists[s]:
                    if i in idx_lists[s]:
                        site_lookup[i] = s
                    for q in range(N_sites):
                        for j in idx_lists[q]:
                            ep_t_t = epsilon[i*N_ads+j]*coverages[i]*coverages[j]
                            eps_theta_theta[s][q] += ep_t_t
        
            #form term_1 and eps_theta matrix (eps_theta needed for jacobian)
            for k in range(N_ads):
                q_k = site_lookup[k]
                for s in range(N_sites):
                    for i in idx_lists[s]:
                        eps_theta[k][s] += epsilon[i*N_ads+k]*coverages[i]
                    term_1[k] += ((f[s][q_k])**2)*eps_theta[k][s]
        
            #form term_2
            for k in range(N_ads):
                q_k = site_lookup[k]
                for s in range(N_sites):
                    term_2[k] += 2*f[s][q_k]*df[s][q_k]*eps_theta_theta[s][q_k]
        
            #combine terms with constant energy to give E_diff
            E_diff = [a+b+c for a,b,c in zip(energies,term_1,term_2)]
        
            if include_derivatives:
                #compute the jacobian
                E_jacob = [[0]*N_ads for x in range(N_ads)]
        
                for k in range(N_ads):
                    for l in range(N_ads):
                        s_l = site_lookup[l]
                        q_k = site_lookup[k]
                        f_sq = f[s_l][q_k]
                        df_sq = df[s_l][q_k]
                        #obtained from tensor calculus, checked against numerical derivative
                        E_jacob[l][k] += (     \
                                 (f_sq**2)*epsilon[l*N_ads+k] +\
                                 2*(f_sq*df_sq*(eps_theta[k][s_l]+eps_theta[l][q_k]) + \
                                 eps_theta_theta[s_l][q_k]*((df_sq**2) + (f_sq*d2f[s_l][q_k]))))
                        
                        #term is only included if l and k are on the same site
                        if l in idx_lists[q_k]:
                            for s in range(N_sites):
                                E_jacob[l][k] += 2*(     \
                                     f[s][q_k]*df[s][q_k]*(eps_theta[k][s]+eps_theta[l][s]) + \
                                     eps_theta_theta[s][q_k]*((df[s][q_k]**2) + \
                                                               f[s][q_k]*d2f[s][q_k]))
            else:
                E_jacob = None
        
            if include_integral:
                E = 0
                for i in range(N_ads):
                    E += energies[i]*coverages[i]
        
                for s in range(N_sites):
                    for q in range(N_sites):
                        for i in idx_lists[s]:
                            for j in idx_lists[q]:
                                E += 0.5*(f[s][q]**2)*epsilon[i*N_ads+j]*coverages[i]*coverages[j]
            else:
                E = 0
        
        
            return E, E_diff, E_jacob
        
    
        kfs = []
        krs = []
        dEfs = []
        dErs = []
    
        n_adsorbates = 8
        n_transition_states = 7
        n_tot = n_adsorbates+n_transition_states
        # Account for the numbers solver where the extra coverage 
        # of the clean slab might be needed
        theta = theta[:n_adsorbates]
        if len(theta) == n_adsorbates:
            theta = list(theta) + [0]*n_transition_states #put transition-state coverages to 0
        elif len(theta) != n_adsorbates+n_transition_states:
            raise ValueError('Coverage vector was not correct length')
        energies = rxn_parameters[:n_tot]
        if len(rxn_parameters) == n_tot + n_tot**2:
            interaction_vector = rxn_parameters[-n_tot**2:]
        elif len(rxn_parameters) == n_tot:
            interaction_vector = [0]*n_tot**2
        else:
            raise ValueError('Length of reaction parameters is not correct. '+ str(rxn_parameters))
    
        G_int, Gf,dGs =  interaction_function(theta,energies,interaction_vector,F,include_derivatives=include_derivatives, include_integral=False)
    
        kB = mpf('0.00008617332478000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000029')
        h = mpf('0.000000000000004135667516000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000012')
        prefactor_list = [kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h]
    
        def element_wise_addition(lists):
            return [sum(L) for L in zip(*lists)]
    
        G_IS = [0]*9
        G_TS = [0]*9
        G_FS = [0]*9
        G_af = [0]*9
        G_ar = [0]*9
        dG_IS = [0]*9
        dG_TS = [0]*9
        dG_FS = [0]*9
        G_IS[0] = gas_energies[3] + site_energies[2] + site_energies[2]
        G_FS[0] = Gf[0] + Gf[0]
        dG_IS[0] = element_wise_addition([[0]*8])
        dG_FS[0] = element_wise_addition([dGs[0] , dGs[0]])
        G_TS[0] = max([G_IS[0],G_FS[0]])
        dG_TS[0] = None #determined later
        
        G_IS[1] = site_energies[0] + gas_energies[1]
        G_FS[1] = Gf[6]
        dG_IS[1] = element_wise_addition([[0]*8])
        dG_FS[1] = element_wise_addition([dGs[6]])
        G_TS[1] = max([G_IS[1],G_FS[1]])
        dG_TS[1] = None #determined later
        
        G_IS[2] = Gf[6] + Gf[0]
        G_FS[2] = Gf[4] + site_energies[2]
        dG_IS[2] = element_wise_addition([dGs[6] , dGs[0]])
        dG_FS[2] = element_wise_addition([dGs[4]])
        G_TS[2] = Gf[12] + site_energies[2]
        dG_TS[2] = element_wise_addition([dGs[12]])
        
        G_IS[3] = Gf[4] + Gf[0]
        G_FS[3] = Gf[3] + site_energies[2]
        dG_IS[3] = element_wise_addition([dGs[4] , dGs[0]])
        dG_FS[3] = element_wise_addition([dGs[3]])
        G_TS[3] = Gf[14] + site_energies[2]
        dG_TS[3] = element_wise_addition([dGs[14]])
        
        G_IS[4] = Gf[3] + site_energies[0]
        G_FS[4] = Gf[5] + Gf[7]
        dG_IS[4] = element_wise_addition([dGs[3]])
        dG_FS[4] = element_wise_addition([dGs[5] , dGs[7]])
        G_TS[4] = Gf[9] + site_energies[0]
        dG_TS[4] = element_wise_addition([dGs[9]])
        
        G_IS[5] = Gf[5] + Gf[0]
        G_FS[5] = Gf[1] + site_energies[2]
        dG_IS[5] = element_wise_addition([dGs[5] , dGs[0]])
        dG_FS[5] = element_wise_addition([dGs[1]])
        G_TS[5] = Gf[8] + site_energies[2]
        dG_TS[5] = element_wise_addition([dGs[8]])
        
        G_IS[6] = Gf[1] + Gf[0]
        G_FS[6] = Gf[2] + site_energies[2]
        dG_IS[6] = element_wise_addition([dGs[1] , dGs[0]])
        dG_FS[6] = element_wise_addition([dGs[2]])
        G_TS[6] = Gf[10] + site_energies[2]
        dG_TS[6] = element_wise_addition([dGs[10]])
        
        G_IS[7] = Gf[2] + Gf[0]
        G_FS[7] = gas_energies[0] + site_energies[0] + site_energies[2]
        dG_IS[7] = element_wise_addition([dGs[2] , dGs[0]])
        dG_FS[7] = element_wise_addition([[0]*8])
        G_TS[7] = Gf[11] + site_energies[2]
        dG_TS[7] = element_wise_addition([dGs[11]])
        
        G_IS[8] = Gf[7] + Gf[0]
        G_FS[8] = gas_energies[2] + site_energies[0] + site_energies[2]
        dG_IS[8] = element_wise_addition([dGs[7] , dGs[0]])
        dG_FS[8] = element_wise_addition([[0]*8])
        G_TS[8] = Gf[13] + site_energies[2]
        dG_TS[8] = element_wise_addition([dGs[13]])
        
    
        n_rxns = len(G_IS)
        for i in range(n_rxns):
            G_list = [G_IS[i],G_FS[i],G_TS[i]]
            G_TS_i = max(G_list)
            max_idx = G_list.index(G_TS_i)
            G_af[i] = G_TS_i - G_IS[i]
            G_ar[i] = G_TS_i - G_FS[i]
            kf = prefactor_list[i]*mpexp(-G_af[i]/(kB*T))
            kr = prefactor_list[i]*mpexp(-G_ar[i]/(kB*T))
            kfs.append(kf)
            krs.append(kr)
            if include_derivatives:
                dG_TS_i = [dG_IS[i],dG_FS[i],dG_TS[i]][max_idx]
    
                dG_TS_i = list(dG_TS_i[:n_adsorbates])
                dG_IS_i = list(dG_IS[i][:n_adsorbates])
                dG_FS_i = list(dG_FS[i][:n_adsorbates])
    
                dEf = [(-dG_TS_j + dG_IS_j) for dG_TS_j,dG_IS_j in zip(dG_TS_i,dG_IS_i)]
                dEr = [(-dG_TS_j + dG_FS_j) for dG_TS_j,dG_FS_j in zip(dG_TS_i,dG_FS_i)]
                dEfs.append(dEf)
                dErs.append(dEr)
            else:
                dEfs = dErs = None
        return kfs, krs, dEfs, dErs
    

    kf, kr, dEf, dEr = rate_constants(
                          rxn_parameters,theta,gas_energies,site_energies,T,F,
                          mpf,matrix,mpexp,mpsqrt,include_derivatives=True)
    s = [0]*3
    s[0] = (mpf('1.0') - theta[1] - theta[2] - theta[3] - theta[4] - theta[5] - theta[6] - theta[7])
    s[1] = (mpf('0.0'))
    s[2] = (mpf('1.0') - theta[0])
    kfkBT = [0]*9
    krkBT = [0]*9
    kfkBT[0] = kf[0]/kBT
    krkBT[0] = kr[0]/kBT
    kfkBT[1] = kf[1]/kBT
    krkBT[1] = kr[1]/kBT
    kfkBT[2] = kf[2]/kBT
    krkBT[2] = kr[2]/kBT
    kfkBT[3] = kf[3]/kBT
    krkBT[3] = kr[3]/kBT
    kfkBT[4] = kf[4]/kBT
    krkBT[4] = kr[4]/kBT
    kfkBT[5] = kf[5]/kBT
    krkBT[5] = kr[5]/kBT
    kfkBT[6] = kf[6]/kBT
    krkBT[6] = kr[6]/kBT
    kfkBT[7] = kf[7]/kBT
    krkBT[7] = kr[7]/kBT
    kfkBT[8] = kf[8]/kBT
    krkBT[8] = kr[8]/kBT
    J[0][0] = 0 + 2*(-2*kf[0]*p[3]*theta[9] + (kfkBT[0])*dEf[0][0]*p[3]*s[2]*s[2] - 2*kr[0]*theta[0] - (krkBT[0])*dEr[0][0]*theta[0]*theta[0]) + -1*(kf[2]*theta[6] + (kfkBT[2])*dEf[2][0]*theta[6]*theta[0] - -1*kr[2]*theta[4] - (krkBT[2])*dEr[2][0]*theta[4]*s[2]) + -1*(kf[3]*theta[4] + (kfkBT[3])*dEf[3][0]*theta[4]*theta[0] - -1*kr[3]*theta[3] - (krkBT[3])*dEr[3][0]*theta[3]*s[2]) + -1*(kf[5]*theta[5] + (kfkBT[5])*dEf[5][0]*theta[5]*theta[0] - -1*kr[5]*theta[1] - (krkBT[5])*dEr[5][0]*theta[1]*s[2]) + -1*(kf[6]*theta[1] + (kfkBT[6])*dEf[6][0]*theta[1]*theta[0] - -1*kr[6]*theta[2] - (krkBT[6])*dEr[6][0]*theta[2]*s[2]) + -1*(kf[7]*theta[2] + (kfkBT[7])*dEf[7][0]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[9] - (krkBT[7])*dEr[7][0]*p[0]*s[0]*s[2]) + -1*(kf[8]*theta[7] + (kfkBT[8])*dEf[8][0]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[9] - (krkBT[8])*dEr[8][0]*p[2]*s[0]*s[2])
    J[0][1] = 0 + 2*(0 + (kfkBT[0])*dEf[0][1]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][1]*theta[0]*theta[0]) + -1*(0 + (kfkBT[2])*dEf[2][1]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][1]*theta[4]*s[2]) + -1*(0 + (kfkBT[3])*dEf[3][1]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][1]*theta[3]*s[2]) + -1*(0 + (kfkBT[5])*dEf[5][1]*theta[5]*theta[0] - kr[5]*theta[8] - (krkBT[5])*dEr[5][1]*theta[1]*s[2]) + -1*(kf[6]*theta[0] + (kfkBT[6])*dEf[6][1]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][1]*theta[2]*s[2]) + -1*(0 + (kfkBT[7])*dEf[7][1]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][1]*p[0]*s[0]*s[2]) + -1*(0 + (kfkBT[8])*dEf[8][1]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][1]*p[2]*s[0]*s[2])
    J[0][2] = 0 + 2*(0 + (kfkBT[0])*dEf[0][2]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][2]*theta[0]*theta[0]) + -1*(0 + (kfkBT[2])*dEf[2][2]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][2]*theta[4]*s[2]) + -1*(0 + (kfkBT[3])*dEf[3][2]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][2]*theta[3]*s[2]) + -1*(0 + (kfkBT[5])*dEf[5][2]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][2]*theta[1]*s[2]) + -1*(0 + (kfkBT[6])*dEf[6][2]*theta[1]*theta[0] - kr[6]*theta[8] - (krkBT[6])*dEr[6][2]*theta[2]*s[2]) + -1*(kf[7]*theta[0] + (kfkBT[7])*dEf[7][2]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][2]*p[0]*s[0]*s[2]) + -1*(0 + (kfkBT[8])*dEf[8][2]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][2]*p[2]*s[0]*s[2])
    J[0][3] = 0 + 2*(0 + (kfkBT[0])*dEf[0][3]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][3]*theta[0]*theta[0]) + -1*(0 + (kfkBT[2])*dEf[2][3]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][3]*theta[4]*s[2]) + -1*(0 + (kfkBT[3])*dEf[3][3]*theta[4]*theta[0] - kr[3]*theta[8] - (krkBT[3])*dEr[3][3]*theta[3]*s[2]) + -1*(0 + (kfkBT[5])*dEf[5][3]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][3]*theta[1]*s[2]) + -1*(0 + (kfkBT[6])*dEf[6][3]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][3]*theta[2]*s[2]) + -1*(0 + (kfkBT[7])*dEf[7][3]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][3]*p[0]*s[0]*s[2]) + -1*(0 + (kfkBT[8])*dEf[8][3]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][3]*p[2]*s[0]*s[2])
    J[0][4] = 0 + 2*(0 + (kfkBT[0])*dEf[0][4]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][4]*theta[0]*theta[0]) + -1*(0 + (kfkBT[2])*dEf[2][4]*theta[6]*theta[0] - kr[2]*theta[8] - (krkBT[2])*dEr[2][4]*theta[4]*s[2]) + -1*(kf[3]*theta[0] + (kfkBT[3])*dEf[3][4]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][4]*theta[3]*s[2]) + -1*(0 + (kfkBT[5])*dEf[5][4]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][4]*theta[1]*s[2]) + -1*(0 + (kfkBT[6])*dEf[6][4]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][4]*theta[2]*s[2]) + -1*(0 + (kfkBT[7])*dEf[7][4]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][4]*p[0]*s[0]*s[2]) + -1*(0 + (kfkBT[8])*dEf[8][4]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][4]*p[2]*s[0]*s[2])
    J[0][5] = 0 + 2*(0 + (kfkBT[0])*dEf[0][5]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][5]*theta[0]*theta[0]) + -1*(0 + (kfkBT[2])*dEf[2][5]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][5]*theta[4]*s[2]) + -1*(0 + (kfkBT[3])*dEf[3][5]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][5]*theta[3]*s[2]) + -1*(kf[5]*theta[0] + (kfkBT[5])*dEf[5][5]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][5]*theta[1]*s[2]) + -1*(0 + (kfkBT[6])*dEf[6][5]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][5]*theta[2]*s[2]) + -1*(0 + (kfkBT[7])*dEf[7][5]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][5]*p[0]*s[0]*s[2]) + -1*(0 + (kfkBT[8])*dEf[8][5]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][5]*p[2]*s[0]*s[2])
    J[0][6] = 0 + 2*(0 + (kfkBT[0])*dEf[0][6]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][6]*theta[0]*theta[0]) + -1*(kf[2]*theta[0] + (kfkBT[2])*dEf[2][6]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][6]*theta[4]*s[2]) + -1*(0 + (kfkBT[3])*dEf[3][6]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][6]*theta[3]*s[2]) + -1*(0 + (kfkBT[5])*dEf[5][6]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][6]*theta[1]*s[2]) + -1*(0 + (kfkBT[6])*dEf[6][6]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][6]*theta[2]*s[2]) + -1*(0 + (kfkBT[7])*dEf[7][6]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][6]*p[0]*s[0]*s[2]) + -1*(0 + (kfkBT[8])*dEf[8][6]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][6]*p[2]*s[0]*s[2])
    J[0][7] = 0 + 2*(0 + (kfkBT[0])*dEf[0][7]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][7]*theta[0]*theta[0]) + -1*(0 + (kfkBT[2])*dEf[2][7]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][7]*theta[4]*s[2]) + -1*(0 + (kfkBT[3])*dEf[3][7]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][7]*theta[3]*s[2]) + -1*(0 + (kfkBT[5])*dEf[5][7]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][7]*theta[1]*s[2]) + -1*(0 + (kfkBT[6])*dEf[6][7]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][7]*theta[2]*s[2]) + -1*(0 + (kfkBT[7])*dEf[7][7]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][7]*p[0]*s[0]*s[2]) + -1*(kf[8]*theta[0] + (kfkBT[8])*dEf[8][7]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][7]*p[2]*s[0]*s[2])
    J[0][8] = 0 + 2*-1*kr[0]*-1*(theta[0]+theta[0]) + -1*(kf[2]*-1*(theta[0]+theta[6]) - kr[2]*-1*(1)*theta[8]) + -1*(kf[3]*-1*(theta[0]+theta[4]) - kr[3]*-1*(1)*theta[8]) + -1*(kf[5]*-1*(theta[0]+theta[5]) - kr[5]*-1*(1)*theta[8]) + -1*(kf[6]*-1*(theta[0]+theta[1]) - kr[6]*-1*(1)*theta[8]) + -1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[8]) + -1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[8])
    J[0][9] = 0 + 2*(2*kf[0]*p[3]*theta[9] - kr[0]*-1*(theta[0]+theta[0])) + -1*(kf[2]*-1*(theta[0]+theta[6]) - 1*kr[2]*theta[4]) + -1*(kf[3]*-1*(theta[0]+theta[4]) - 1*kr[3]*theta[3]) + -1*(kf[5]*-1*(theta[0]+theta[5]) - 1*kr[5]*theta[1]) + -1*(kf[6]*-1*(theta[0]+theta[1]) - 1*kr[6]*theta[2]) + -1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[9]) + -1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[9])
    J[1][0] = 0 + 1*(kf[5]*theta[5] + (kfkBT[5])*dEf[5][0]*theta[5]*theta[0] - -1*kr[5]*theta[1] - (krkBT[5])*dEr[5][0]*theta[1]*s[2]) + -1*(kf[6]*theta[1] + (kfkBT[6])*dEf[6][0]*theta[1]*theta[0] - -1*kr[6]*theta[2] - (krkBT[6])*dEr[6][0]*theta[2]*s[2])
    J[1][1] = 0 + 1*(0 + (kfkBT[5])*dEf[5][1]*theta[5]*theta[0] - kr[5]*theta[8] - (krkBT[5])*dEr[5][1]*theta[1]*s[2]) + -1*(kf[6]*theta[0] + (kfkBT[6])*dEf[6][1]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][1]*theta[2]*s[2])
    J[1][2] = 0 + 1*(0 + (kfkBT[5])*dEf[5][2]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][2]*theta[1]*s[2]) + -1*(0 + (kfkBT[6])*dEf[6][2]*theta[1]*theta[0] - kr[6]*theta[8] - (krkBT[6])*dEr[6][2]*theta[2]*s[2])
    J[1][3] = 0 + 1*(0 + (kfkBT[5])*dEf[5][3]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][3]*theta[1]*s[2]) + -1*(0 + (kfkBT[6])*dEf[6][3]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][3]*theta[2]*s[2])
    J[1][4] = 0 + 1*(0 + (kfkBT[5])*dEf[5][4]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][4]*theta[1]*s[2]) + -1*(0 + (kfkBT[6])*dEf[6][4]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][4]*theta[2]*s[2])
    J[1][5] = 0 + 1*(kf[5]*theta[0] + (kfkBT[5])*dEf[5][5]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][5]*theta[1]*s[2]) + -1*(0 + (kfkBT[6])*dEf[6][5]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][5]*theta[2]*s[2])
    J[1][6] = 0 + 1*(0 + (kfkBT[5])*dEf[5][6]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][6]*theta[1]*s[2]) + -1*(0 + (kfkBT[6])*dEf[6][6]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][6]*theta[2]*s[2])
    J[1][7] = 0 + 1*(0 + (kfkBT[5])*dEf[5][7]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][7]*theta[1]*s[2]) + -1*(0 + (kfkBT[6])*dEf[6][7]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][7]*theta[2]*s[2])
    J[1][8] = 0 + 1*(kf[5]*-1*(theta[0]+theta[5]) - kr[5]*-1*(1)*theta[8]) + -1*(kf[6]*-1*(theta[0]+theta[1]) - kr[6]*-1*(1)*theta[8])
    J[1][9] = 0 + 1*(kf[5]*-1*(theta[0]+theta[5]) - 1*kr[5]*theta[1]) + -1*(kf[6]*-1*(theta[0]+theta[1]) - 1*kr[6]*theta[2])
    J[2][0] = 0 + 1*(kf[6]*theta[1] + (kfkBT[6])*dEf[6][0]*theta[1]*theta[0] - -1*kr[6]*theta[2] - (krkBT[6])*dEr[6][0]*theta[2]*s[2]) + -1*(kf[7]*theta[2] + (kfkBT[7])*dEf[7][0]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[9] - (krkBT[7])*dEr[7][0]*p[0]*s[0]*s[2])
    J[2][1] = 0 + 1*(kf[6]*theta[0] + (kfkBT[6])*dEf[6][1]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][1]*theta[2]*s[2]) + -1*(0 + (kfkBT[7])*dEf[7][1]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][1]*p[0]*s[0]*s[2])
    J[2][2] = 0 + 1*(0 + (kfkBT[6])*dEf[6][2]*theta[1]*theta[0] - kr[6]*theta[8] - (krkBT[6])*dEr[6][2]*theta[2]*s[2]) + -1*(kf[7]*theta[0] + (kfkBT[7])*dEf[7][2]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][2]*p[0]*s[0]*s[2])
    J[2][3] = 0 + 1*(0 + (kfkBT[6])*dEf[6][3]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][3]*theta[2]*s[2]) + -1*(0 + (kfkBT[7])*dEf[7][3]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][3]*p[0]*s[0]*s[2])
    J[2][4] = 0 + 1*(0 + (kfkBT[6])*dEf[6][4]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][4]*theta[2]*s[2]) + -1*(0 + (kfkBT[7])*dEf[7][4]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][4]*p[0]*s[0]*s[2])
    J[2][5] = 0 + 1*(0 + (kfkBT[6])*dEf[6][5]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][5]*theta[2]*s[2]) + -1*(0 + (kfkBT[7])*dEf[7][5]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][5]*p[0]*s[0]*s[2])
    J[2][6] = 0 + 1*(0 + (kfkBT[6])*dEf[6][6]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][6]*theta[2]*s[2]) + -1*(0 + (kfkBT[7])*dEf[7][6]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][6]*p[0]*s[0]*s[2])
    J[2][7] = 0 + 1*(0 + (kfkBT[6])*dEf[6][7]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][7]*theta[2]*s[2]) + -1*(0 + (kfkBT[7])*dEf[7][7]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][7]*p[0]*s[0]*s[2])
    J[2][8] = 0 + 1*(kf[6]*-1*(theta[0]+theta[1]) - kr[6]*-1*(1)*theta[8]) + -1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[8])
    J[2][9] = 0 + 1*(kf[6]*-1*(theta[0]+theta[1]) - 1*kr[6]*theta[2]) + -1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[9])
    J[3][0] = 0 + 1*(kf[3]*theta[4] + (kfkBT[3])*dEf[3][0]*theta[4]*theta[0] - -1*kr[3]*theta[3] - (krkBT[3])*dEr[3][0]*theta[3]*s[2]) + -1*(0 + (kfkBT[4])*dEf[4][0]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][0]*theta[5]*theta[7])
    J[3][1] = 0 + 1*(0 + (kfkBT[3])*dEf[3][1]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][1]*theta[3]*s[2]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][1]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][1]*theta[5]*theta[7])
    J[3][2] = 0 + 1*(0 + (kfkBT[3])*dEf[3][2]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][2]*theta[3]*s[2]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][2]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][2]*theta[5]*theta[7])
    J[3][3] = 0 + 1*(0 + (kfkBT[3])*dEf[3][3]*theta[4]*theta[0] - kr[3]*theta[8] - (krkBT[3])*dEr[3][3]*theta[3]*s[2]) + -1*(kf[4]*(-1*theta[3] + 1**theta[8]) + (kfkBT[4])*dEf[4][3]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][3]*theta[5]*theta[7])
    J[3][4] = 0 + 1*(kf[3]*theta[0] + (kfkBT[3])*dEf[3][4]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][4]*theta[3]*s[2]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][4]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][4]*theta[5]*theta[7])
    J[3][5] = 0 + 1*(0 + (kfkBT[3])*dEf[3][5]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][5]*theta[3]*s[2]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][5]*theta[3]*s[0] - kr[4]*theta[7] - (krkBT[4])*dEr[4][5]*theta[5]*theta[7])
    J[3][6] = 0 + 1*(0 + (kfkBT[3])*dEf[3][6]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][6]*theta[3]*s[2]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][6]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][6]*theta[5]*theta[7])
    J[3][7] = 0 + 1*(0 + (kfkBT[3])*dEf[3][7]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][7]*theta[3]*s[2]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][7]*theta[3]*s[0] - kr[4]*theta[5] - (krkBT[4])*dEr[4][7]*theta[5]*theta[7])
    J[3][8] = 0 + 1*(kf[3]*-1*(theta[0]+theta[4]) - kr[3]*-1*(1)*theta[8]) + -1*(1*kf[4]*theta[3] - kr[4]*-1*(theta[7]+theta[5]))
    J[3][9] = 0 + 1*(kf[3]*-1*(theta[0]+theta[4]) - 1*kr[3]*theta[3]) + -1*(kf[4]*-1*(1)*theta[9] - kr[4]*-1*(theta[7]+theta[5]))
    J[4][0] = 0 + 1*(kf[2]*theta[6] + (kfkBT[2])*dEf[2][0]*theta[6]*theta[0] - -1*kr[2]*theta[4] - (krkBT[2])*dEr[2][0]*theta[4]*s[2]) + -1*(kf[3]*theta[4] + (kfkBT[3])*dEf[3][0]*theta[4]*theta[0] - -1*kr[3]*theta[3] - (krkBT[3])*dEr[3][0]*theta[3]*s[2])
    J[4][1] = 0 + 1*(0 + (kfkBT[2])*dEf[2][1]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][1]*theta[4]*s[2]) + -1*(0 + (kfkBT[3])*dEf[3][1]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][1]*theta[3]*s[2])
    J[4][2] = 0 + 1*(0 + (kfkBT[2])*dEf[2][2]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][2]*theta[4]*s[2]) + -1*(0 + (kfkBT[3])*dEf[3][2]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][2]*theta[3]*s[2])
    J[4][3] = 0 + 1*(0 + (kfkBT[2])*dEf[2][3]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][3]*theta[4]*s[2]) + -1*(0 + (kfkBT[3])*dEf[3][3]*theta[4]*theta[0] - kr[3]*theta[8] - (krkBT[3])*dEr[3][3]*theta[3]*s[2])
    J[4][4] = 0 + 1*(0 + (kfkBT[2])*dEf[2][4]*theta[6]*theta[0] - kr[2]*theta[8] - (krkBT[2])*dEr[2][4]*theta[4]*s[2]) + -1*(kf[3]*theta[0] + (kfkBT[3])*dEf[3][4]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][4]*theta[3]*s[2])
    J[4][5] = 0 + 1*(0 + (kfkBT[2])*dEf[2][5]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][5]*theta[4]*s[2]) + -1*(0 + (kfkBT[3])*dEf[3][5]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][5]*theta[3]*s[2])
    J[4][6] = 0 + 1*(kf[2]*theta[0] + (kfkBT[2])*dEf[2][6]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][6]*theta[4]*s[2]) + -1*(0 + (kfkBT[3])*dEf[3][6]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][6]*theta[3]*s[2])
    J[4][7] = 0 + 1*(0 + (kfkBT[2])*dEf[2][7]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][7]*theta[4]*s[2]) + -1*(0 + (kfkBT[3])*dEf[3][7]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][7]*theta[3]*s[2])
    J[4][8] = 0 + 1*(kf[2]*-1*(theta[0]+theta[6]) - kr[2]*-1*(1)*theta[8]) + -1*(kf[3]*-1*(theta[0]+theta[4]) - kr[3]*-1*(1)*theta[8])
    J[4][9] = 0 + 1*(kf[2]*-1*(theta[0]+theta[6]) - 1*kr[2]*theta[4]) + -1*(kf[3]*-1*(theta[0]+theta[4]) - 1*kr[3]*theta[3])
    J[5][0] = 0 + 1*(0 + (kfkBT[4])*dEf[4][0]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][0]*theta[5]*theta[7]) + -1*(kf[5]*theta[5] + (kfkBT[5])*dEf[5][0]*theta[5]*theta[0] - -1*kr[5]*theta[1] - (krkBT[5])*dEr[5][0]*theta[1]*s[2])
    J[5][1] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][1]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][1]*theta[5]*theta[7]) + -1*(0 + (kfkBT[5])*dEf[5][1]*theta[5]*theta[0] - kr[5]*theta[8] - (krkBT[5])*dEr[5][1]*theta[1]*s[2])
    J[5][2] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][2]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][2]*theta[5]*theta[7]) + -1*(0 + (kfkBT[5])*dEf[5][2]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][2]*theta[1]*s[2])
    J[5][3] = 0 + 1*(kf[4]*(-1*theta[3] + 1**theta[8]) + (kfkBT[4])*dEf[4][3]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][3]*theta[5]*theta[7]) + -1*(0 + (kfkBT[5])*dEf[5][3]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][3]*theta[1]*s[2])
    J[5][4] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][4]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][4]*theta[5]*theta[7]) + -1*(0 + (kfkBT[5])*dEf[5][4]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][4]*theta[1]*s[2])
    J[5][5] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][5]*theta[3]*s[0] - kr[4]*theta[7] - (krkBT[4])*dEr[4][5]*theta[5]*theta[7]) + -1*(kf[5]*theta[0] + (kfkBT[5])*dEf[5][5]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][5]*theta[1]*s[2])
    J[5][6] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][6]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][6]*theta[5]*theta[7]) + -1*(0 + (kfkBT[5])*dEf[5][6]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][6]*theta[1]*s[2])
    J[5][7] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][7]*theta[3]*s[0] - kr[4]*theta[5] - (krkBT[4])*dEr[4][7]*theta[5]*theta[7]) + -1*(0 + (kfkBT[5])*dEf[5][7]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][7]*theta[1]*s[2])
    J[5][8] = 0 + 1*(1*kf[4]*theta[3] - kr[4]*-1*(theta[7]+theta[5])) + -1*(kf[5]*-1*(theta[0]+theta[5]) - kr[5]*-1*(1)*theta[8])
    J[5][9] = 0 + 1*(kf[4]*-1*(1)*theta[9] - kr[4]*-1*(theta[7]+theta[5])) + -1*(kf[5]*-1*(theta[0]+theta[5]) - 1*kr[5]*theta[1])
    J[6][0] = 0 + 1*(0 + (kfkBT[1])*dEf[1][0]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][0]*theta[6]) + -1*(kf[2]*theta[6] + (kfkBT[2])*dEf[2][0]*theta[6]*theta[0] - -1*kr[2]*theta[4] - (krkBT[2])*dEr[2][0]*theta[4]*s[2])
    J[6][1] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][1]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][1]*theta[6]) + -1*(0 + (kfkBT[2])*dEf[2][1]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][1]*theta[4]*s[2])
    J[6][2] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][2]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][2]*theta[6]) + -1*(0 + (kfkBT[2])*dEf[2][2]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][2]*theta[4]*s[2])
    J[6][3] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][3]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][3]*theta[6]) + -1*(0 + (kfkBT[2])*dEf[2][3]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][3]*theta[4]*s[2])
    J[6][4] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][4]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][4]*theta[6]) + -1*(0 + (kfkBT[2])*dEf[2][4]*theta[6]*theta[0] - kr[2]*theta[8] - (krkBT[2])*dEr[2][4]*theta[4]*s[2])
    J[6][5] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][5]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][5]*theta[6]) + -1*(0 + (kfkBT[2])*dEf[2][5]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][5]*theta[4]*s[2])
    J[6][6] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][6]*p[1]*s[0] - kr[1] - (krkBT[1])*dEr[1][6]*theta[6]) + -1*(kf[2]*theta[0] + (kfkBT[2])*dEf[2][6]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][6]*theta[4]*s[2])
    J[6][7] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][7]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][7]*theta[6]) + -1*(0 + (kfkBT[2])*dEf[2][7]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][7]*theta[4]*s[2])
    J[6][8] = 0 + 1*(1*kf[1]*p[1] - kr[1]*-1*(1)) + -1*(kf[2]*-1*(theta[0]+theta[6]) - kr[2]*-1*(1)*theta[8])
    J[6][9] = 0 + 1*-1*kr[1]*-1*(1) + -1*(kf[2]*-1*(theta[0]+theta[6]) - 1*kr[2]*theta[4])
    J[7][0] = 0 + 1*(0 + (kfkBT[4])*dEf[4][0]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][0]*theta[5]*theta[7]) + -1*(kf[8]*theta[7] + (kfkBT[8])*dEf[8][0]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[9] - (krkBT[8])*dEr[8][0]*p[2]*s[0]*s[2])
    J[7][1] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][1]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][1]*theta[5]*theta[7]) + -1*(0 + (kfkBT[8])*dEf[8][1]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][1]*p[2]*s[0]*s[2])
    J[7][2] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][2]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][2]*theta[5]*theta[7]) + -1*(0 + (kfkBT[8])*dEf[8][2]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][2]*p[2]*s[0]*s[2])
    J[7][3] = 0 + 1*(kf[4]*(-1*theta[3] + 1**theta[8]) + (kfkBT[4])*dEf[4][3]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][3]*theta[5]*theta[7]) + -1*(0 + (kfkBT[8])*dEf[8][3]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][3]*p[2]*s[0]*s[2])
    J[7][4] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][4]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][4]*theta[5]*theta[7]) + -1*(0 + (kfkBT[8])*dEf[8][4]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][4]*p[2]*s[0]*s[2])
    J[7][5] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][5]*theta[3]*s[0] - kr[4]*theta[7] - (krkBT[4])*dEr[4][5]*theta[5]*theta[7]) + -1*(0 + (kfkBT[8])*dEf[8][5]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][5]*p[2]*s[0]*s[2])
    J[7][6] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][6]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][6]*theta[5]*theta[7]) + -1*(0 + (kfkBT[8])*dEf[8][6]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][6]*p[2]*s[0]*s[2])
    J[7][7] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][7]*theta[3]*s[0] - kr[4]*theta[5] - (krkBT[4])*dEr[4][7]*theta[5]*theta[7]) + -1*(kf[8]*theta[0] + (kfkBT[8])*dEf[8][7]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][7]*p[2]*s[0]*s[2])
    J[7][8] = 0 + 1*(1*kf[4]*theta[3] - kr[4]*-1*(theta[7]+theta[5])) + -1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[8])
    J[7][9] = 0 + 1*(kf[4]*-1*(1)*theta[9] - kr[4]*-1*(theta[7]+theta[5])) + -1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[9])
    J[8][0] = 0 + -1*(0 + (kfkBT[1])*dEf[1][0]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][0]*theta[6]) + -1*(0 + (kfkBT[4])*dEf[4][0]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][0]*theta[5]*theta[7]) + 1*(kf[7]*theta[2] + (kfkBT[7])*dEf[7][0]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[9] - (krkBT[7])*dEr[7][0]*p[0]*s[0]*s[2]) + 1*(kf[8]*theta[7] + (kfkBT[8])*dEf[8][0]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[9] - (krkBT[8])*dEr[8][0]*p[2]*s[0]*s[2])
    J[8][1] = 0 + -1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][1]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][1]*theta[6]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][1]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][1]*theta[5]*theta[7]) + 1*(0 + (kfkBT[7])*dEf[7][1]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][1]*p[0]*s[0]*s[2]) + 1*(0 + (kfkBT[8])*dEf[8][1]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][1]*p[2]*s[0]*s[2])
    J[8][2] = 0 + -1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][2]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][2]*theta[6]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][2]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][2]*theta[5]*theta[7]) + 1*(kf[7]*theta[0] + (kfkBT[7])*dEf[7][2]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][2]*p[0]*s[0]*s[2]) + 1*(0 + (kfkBT[8])*dEf[8][2]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][2]*p[2]*s[0]*s[2])
    J[8][3] = 0 + -1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][3]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][3]*theta[6]) + -1*(kf[4]*(-1*theta[3] + 1**theta[8]) + (kfkBT[4])*dEf[4][3]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][3]*theta[5]*theta[7]) + 1*(0 + (kfkBT[7])*dEf[7][3]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][3]*p[0]*s[0]*s[2]) + 1*(0 + (kfkBT[8])*dEf[8][3]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][3]*p[2]*s[0]*s[2])
    J[8][4] = 0 + -1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][4]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][4]*theta[6]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][4]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][4]*theta[5]*theta[7]) + 1*(0 + (kfkBT[7])*dEf[7][4]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][4]*p[0]*s[0]*s[2]) + 1*(0 + (kfkBT[8])*dEf[8][4]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][4]*p[2]*s[0]*s[2])
    J[8][5] = 0 + -1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][5]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][5]*theta[6]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][5]*theta[3]*s[0] - kr[4]*theta[7] - (krkBT[4])*dEr[4][5]*theta[5]*theta[7]) + 1*(0 + (kfkBT[7])*dEf[7][5]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][5]*p[0]*s[0]*s[2]) + 1*(0 + (kfkBT[8])*dEf[8][5]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][5]*p[2]*s[0]*s[2])
    J[8][6] = 0 + -1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][6]*p[1]*s[0] - kr[1] - (krkBT[1])*dEr[1][6]*theta[6]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][6]*theta[3]*s[0] - 0 - (krkBT[4])*dEr[4][6]*theta[5]*theta[7]) + 1*(0 + (kfkBT[7])*dEf[7][6]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][6]*p[0]*s[0]*s[2]) + 1*(0 + (kfkBT[8])*dEf[8][6]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][6]*p[2]*s[0]*s[2])
    J[8][7] = 0 + -1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][7]*p[1]*s[0] - 0 - (krkBT[1])*dEr[1][7]*theta[6]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][7]*theta[3]*s[0] - kr[4]*theta[5] - (krkBT[4])*dEr[4][7]*theta[5]*theta[7]) + 1*(0 + (kfkBT[7])*dEf[7][7]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][7]*p[0]*s[0]*s[2]) + 1*(kf[8]*theta[0] + (kfkBT[8])*dEf[8][7]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][7]*p[2]*s[0]*s[2])
    J[8][8] = 0 + -1*(1*kf[1]*p[1] - kr[1]*-1*(1)) + -1*(1*kf[4]*theta[3] - kr[4]*-1*(theta[7]+theta[5])) + 1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[8]) + 1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[8])
    J[8][9] = 0 + -1*-1*kr[1]*-1*(1) + -1*(kf[4]*-1*(1)*theta[9] - kr[4]*-1*(theta[7]+theta[5])) + 1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[9]) + 1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[9])
    J[9][0] = 0 + -2*(-2*kf[0]*p[3]*theta[9] + (kfkBT[0])*dEf[0][0]*p[3]*s[2]*s[2] - 2*kr[0]*theta[0] - (krkBT[0])*dEr[0][0]*theta[0]*theta[0]) + 1*(kf[2]*theta[6] + (kfkBT[2])*dEf[2][0]*theta[6]*theta[0] - -1*kr[2]*theta[4] - (krkBT[2])*dEr[2][0]*theta[4]*s[2]) + 1*(kf[3]*theta[4] + (kfkBT[3])*dEf[3][0]*theta[4]*theta[0] - -1*kr[3]*theta[3] - (krkBT[3])*dEr[3][0]*theta[3]*s[2]) + 1*(kf[5]*theta[5] + (kfkBT[5])*dEf[5][0]*theta[5]*theta[0] - -1*kr[5]*theta[1] - (krkBT[5])*dEr[5][0]*theta[1]*s[2]) + 1*(kf[6]*theta[1] + (kfkBT[6])*dEf[6][0]*theta[1]*theta[0] - -1*kr[6]*theta[2] - (krkBT[6])*dEr[6][0]*theta[2]*s[2]) + 1*(kf[7]*theta[2] + (kfkBT[7])*dEf[7][0]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[9] - (krkBT[7])*dEr[7][0]*p[0]*s[0]*s[2]) + 1*(kf[8]*theta[7] + (kfkBT[8])*dEf[8][0]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[9] - (krkBT[8])*dEr[8][0]*p[2]*s[0]*s[2])
    J[9][1] = 0 + -2*(0 + (kfkBT[0])*dEf[0][1]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][1]*theta[0]*theta[0]) + 1*(0 + (kfkBT[2])*dEf[2][1]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][1]*theta[4]*s[2]) + 1*(0 + (kfkBT[3])*dEf[3][1]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][1]*theta[3]*s[2]) + 1*(0 + (kfkBT[5])*dEf[5][1]*theta[5]*theta[0] - kr[5]*theta[8] - (krkBT[5])*dEr[5][1]*theta[1]*s[2]) + 1*(kf[6]*theta[0] + (kfkBT[6])*dEf[6][1]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][1]*theta[2]*s[2]) + 1*(0 + (kfkBT[7])*dEf[7][1]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][1]*p[0]*s[0]*s[2]) + 1*(0 + (kfkBT[8])*dEf[8][1]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][1]*p[2]*s[0]*s[2])
    J[9][2] = 0 + -2*(0 + (kfkBT[0])*dEf[0][2]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][2]*theta[0]*theta[0]) + 1*(0 + (kfkBT[2])*dEf[2][2]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][2]*theta[4]*s[2]) + 1*(0 + (kfkBT[3])*dEf[3][2]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][2]*theta[3]*s[2]) + 1*(0 + (kfkBT[5])*dEf[5][2]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][2]*theta[1]*s[2]) + 1*(0 + (kfkBT[6])*dEf[6][2]*theta[1]*theta[0] - kr[6]*theta[8] - (krkBT[6])*dEr[6][2]*theta[2]*s[2]) + 1*(kf[7]*theta[0] + (kfkBT[7])*dEf[7][2]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][2]*p[0]*s[0]*s[2]) + 1*(0 + (kfkBT[8])*dEf[8][2]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][2]*p[2]*s[0]*s[2])
    J[9][3] = 0 + -2*(0 + (kfkBT[0])*dEf[0][3]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][3]*theta[0]*theta[0]) + 1*(0 + (kfkBT[2])*dEf[2][3]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][3]*theta[4]*s[2]) + 1*(0 + (kfkBT[3])*dEf[3][3]*theta[4]*theta[0] - kr[3]*theta[8] - (krkBT[3])*dEr[3][3]*theta[3]*s[2]) + 1*(0 + (kfkBT[5])*dEf[5][3]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][3]*theta[1]*s[2]) + 1*(0 + (kfkBT[6])*dEf[6][3]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][3]*theta[2]*s[2]) + 1*(0 + (kfkBT[7])*dEf[7][3]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][3]*p[0]*s[0]*s[2]) + 1*(0 + (kfkBT[8])*dEf[8][3]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][3]*p[2]*s[0]*s[2])
    J[9][4] = 0 + -2*(0 + (kfkBT[0])*dEf[0][4]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][4]*theta[0]*theta[0]) + 1*(0 + (kfkBT[2])*dEf[2][4]*theta[6]*theta[0] - kr[2]*theta[8] - (krkBT[2])*dEr[2][4]*theta[4]*s[2]) + 1*(kf[3]*theta[0] + (kfkBT[3])*dEf[3][4]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][4]*theta[3]*s[2]) + 1*(0 + (kfkBT[5])*dEf[5][4]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][4]*theta[1]*s[2]) + 1*(0 + (kfkBT[6])*dEf[6][4]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][4]*theta[2]*s[2]) + 1*(0 + (kfkBT[7])*dEf[7][4]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][4]*p[0]*s[0]*s[2]) + 1*(0 + (kfkBT[8])*dEf[8][4]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][4]*p[2]*s[0]*s[2])
    J[9][5] = 0 + -2*(0 + (kfkBT[0])*dEf[0][5]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][5]*theta[0]*theta[0]) + 1*(0 + (kfkBT[2])*dEf[2][5]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][5]*theta[4]*s[2]) + 1*(0 + (kfkBT[3])*dEf[3][5]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][5]*theta[3]*s[2]) + 1*(kf[5]*theta[0] + (kfkBT[5])*dEf[5][5]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][5]*theta[1]*s[2]) + 1*(0 + (kfkBT[6])*dEf[6][5]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][5]*theta[2]*s[2]) + 1*(0 + (kfkBT[7])*dEf[7][5]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][5]*p[0]*s[0]*s[2]) + 1*(0 + (kfkBT[8])*dEf[8][5]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][5]*p[2]*s[0]*s[2])
    J[9][6] = 0 + -2*(0 + (kfkBT[0])*dEf[0][6]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][6]*theta[0]*theta[0]) + 1*(kf[2]*theta[0] + (kfkBT[2])*dEf[2][6]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][6]*theta[4]*s[2]) + 1*(0 + (kfkBT[3])*dEf[3][6]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][6]*theta[3]*s[2]) + 1*(0 + (kfkBT[5])*dEf[5][6]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][6]*theta[1]*s[2]) + 1*(0 + (kfkBT[6])*dEf[6][6]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][6]*theta[2]*s[2]) + 1*(0 + (kfkBT[7])*dEf[7][6]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][6]*p[0]*s[0]*s[2]) + 1*(0 + (kfkBT[8])*dEf[8][6]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][6]*p[2]*s[0]*s[2])
    J[9][7] = 0 + -2*(0 + (kfkBT[0])*dEf[0][7]*p[3]*s[2]*s[2] - 0 - (krkBT[0])*dEr[0][7]*theta[0]*theta[0]) + 1*(0 + (kfkBT[2])*dEf[2][7]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][7]*theta[4]*s[2]) + 1*(0 + (kfkBT[3])*dEf[3][7]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][7]*theta[3]*s[2]) + 1*(0 + (kfkBT[5])*dEf[5][7]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][7]*theta[1]*s[2]) + 1*(0 + (kfkBT[6])*dEf[6][7]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][7]*theta[2]*s[2]) + 1*(0 + (kfkBT[7])*dEf[7][7]*theta[2]*theta[0] - -1*kr[7]*p[0]*theta[8] - (krkBT[7])*dEr[7][7]*p[0]*s[0]*s[2]) + 1*(kf[8]*theta[0] + (kfkBT[8])*dEf[8][7]*theta[7]*theta[0] - -1*kr[8]*p[2]*theta[8] - (krkBT[8])*dEr[8][7]*p[2]*s[0]*s[2])
    J[9][8] = 0 + -2*-1*kr[0]*-1*(theta[0]+theta[0]) + 1*(kf[2]*-1*(theta[0]+theta[6]) - kr[2]*-1*(1)*theta[8]) + 1*(kf[3]*-1*(theta[0]+theta[4]) - kr[3]*-1*(1)*theta[8]) + 1*(kf[5]*-1*(theta[0]+theta[5]) - kr[5]*-1*(1)*theta[8]) + 1*(kf[6]*-1*(theta[0]+theta[1]) - kr[6]*-1*(1)*theta[8]) + 1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[8]) + 1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[8])
    J[9][9] = 0 + -2*(2*kf[0]*p[3]*theta[9] - kr[0]*-1*(theta[0]+theta[0])) + 1*(kf[2]*-1*(theta[0]+theta[6]) - 1*kr[2]*theta[4]) + 1*(kf[3]*-1*(theta[0]+theta[4]) - 1*kr[3]*theta[3]) + 1*(kf[5]*-1*(theta[0]+theta[5]) - 1*kr[5]*theta[1]) + 1*(kf[6]*-1*(theta[0]+theta[1]) - 1*kr[6]*theta[2]) + 1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[9]) + 1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[9])
    
    J = matrix(J)
    return J




def ideal_mean_field_jacobian(kf,kr,theta,p,mpf,matrix):
    n_adsorbates = 8
    J = [[0 for i in range(10)] for j in range(10)]

    s = [0]*3
    s[0] = (mpf('1.0') - theta[1] - theta[2] - theta[3] - theta[4] - theta[5] - theta[6] - theta[7])
    s[1] = (mpf('0.0'))
    s[2] = (mpf('1.0') - theta[0])
    J[0][0] = 0 + 2*(-2*kf[0]*p[3]*theta[9] - 2*kr[0]*theta[0]) + -1*(kf[2]*theta[6] - -1*kr[2]*theta[4]) + -1*(kf[3]*theta[4] - -1*kr[3]*theta[3]) + -1*(kf[5]*theta[5] - -1*kr[5]*theta[1]) + -1*(kf[6]*theta[1] - -1*kr[6]*theta[2]) + -1*(kf[7]*theta[2] - -1*kr[7]*p[0]*theta[9]) + -1*(kf[8]*theta[7] - -1*kr[8]*p[2]*theta[9])
    J[0][1] = 0 + -1*-1*kr[5]*theta[8] + -1*kf[6]*theta[0] + -1*-1*-1*kr[7]*p[0]*theta[8] + -1*-1*-1*kr[8]*p[2]*theta[8]
    J[0][2] = 0 + -1*-1*kr[6]*theta[8] + -1*(kf[7]*theta[0] - -1*kr[7]*p[0]*theta[8]) + -1*-1*-1*kr[8]*p[2]*theta[8]
    J[0][3] = 0 + -1*-1*kr[3]*theta[8] + -1*-1*-1*kr[7]*p[0]*theta[8] + -1*-1*-1*kr[8]*p[2]*theta[8]
    J[0][4] = 0 + -1*-1*kr[2]*theta[8] + -1*kf[3]*theta[0] + -1*-1*-1*kr[7]*p[0]*theta[8] + -1*-1*-1*kr[8]*p[2]*theta[8]
    J[0][5] = 0 + -1*kf[5]*theta[0] + -1*-1*-1*kr[7]*p[0]*theta[8] + -1*-1*-1*kr[8]*p[2]*theta[8]
    J[0][6] = 0 + -1*kf[2]*theta[0] + -1*-1*-1*kr[7]*p[0]*theta[8] + -1*-1*-1*kr[8]*p[2]*theta[8]
    J[0][7] = 0 + -1*-1*-1*kr[7]*p[0]*theta[8] + -1*(kf[8]*theta[0] - -1*kr[8]*p[2]*theta[8])
    J[0][8] = 0 + 2*-1*kr[0]*-1*(theta[0]+theta[0]) + -1*(kf[2]*-1*(theta[0]+theta[6]) - kr[2]*-1*(1)*theta[8]) + -1*(kf[3]*-1*(theta[0]+theta[4]) - kr[3]*-1*(1)*theta[8]) + -1*(kf[5]*-1*(theta[0]+theta[5]) - kr[5]*-1*(1)*theta[8]) + -1*(kf[6]*-1*(theta[0]+theta[1]) - kr[6]*-1*(1)*theta[8]) + -1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[8]) + -1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[8])
    J[0][9] = 0 + 2*(2*kf[0]*p[3]*theta[9] - kr[0]*-1*(theta[0]+theta[0])) + -1*(kf[2]*-1*(theta[0]+theta[6]) - 1*kr[2]*theta[4]) + -1*(kf[3]*-1*(theta[0]+theta[4]) - 1*kr[3]*theta[3]) + -1*(kf[5]*-1*(theta[0]+theta[5]) - 1*kr[5]*theta[1]) + -1*(kf[6]*-1*(theta[0]+theta[1]) - 1*kr[6]*theta[2]) + -1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[9]) + -1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[9])
    J[1][0] = 0 + 1*(kf[5]*theta[5] - -1*kr[5]*theta[1]) + -1*(kf[6]*theta[1] - -1*kr[6]*theta[2])
    J[1][1] = 0 + 1*-1*kr[5]*theta[8] + -1*kf[6]*theta[0]
    J[1][2] = 0 + -1*-1*kr[6]*theta[8]
    J[1][3] = 0
    J[1][4] = 0
    J[1][5] = 0 + 1*kf[5]*theta[0]
    J[1][6] = 0
    J[1][7] = 0
    J[1][8] = 0 + 1*(kf[5]*-1*(theta[0]+theta[5]) - kr[5]*-1*(1)*theta[8]) + -1*(kf[6]*-1*(theta[0]+theta[1]) - kr[6]*-1*(1)*theta[8])
    J[1][9] = 0 + 1*(kf[5]*-1*(theta[0]+theta[5]) - 1*kr[5]*theta[1]) + -1*(kf[6]*-1*(theta[0]+theta[1]) - 1*kr[6]*theta[2])
    J[2][0] = 0 + 1*(kf[6]*theta[1] - -1*kr[6]*theta[2]) + -1*(kf[7]*theta[2] - -1*kr[7]*p[0]*theta[9])
    J[2][1] = 0 + 1*kf[6]*theta[0] + -1*-1*-1*kr[7]*p[0]*theta[8]
    J[2][2] = 0 + 1*-1*kr[6]*theta[8] + -1*(kf[7]*theta[0] - -1*kr[7]*p[0]*theta[8])
    J[2][3] = 0 + -1*-1*-1*kr[7]*p[0]*theta[8]
    J[2][4] = 0 + -1*-1*-1*kr[7]*p[0]*theta[8]
    J[2][5] = 0 + -1*-1*-1*kr[7]*p[0]*theta[8]
    J[2][6] = 0 + -1*-1*-1*kr[7]*p[0]*theta[8]
    J[2][7] = 0 + -1*-1*-1*kr[7]*p[0]*theta[8]
    J[2][8] = 0 + 1*(kf[6]*-1*(theta[0]+theta[1]) - kr[6]*-1*(1)*theta[8]) + -1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[8])
    J[2][9] = 0 + 1*(kf[6]*-1*(theta[0]+theta[1]) - 1*kr[6]*theta[2]) + -1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[9])
    J[3][0] = 0 + 1*(kf[3]*theta[4] - -1*kr[3]*theta[3])
    J[3][1] = 0 + -1*-1*kf[4]*theta[3]
    J[3][2] = 0 + -1*-1*kf[4]*theta[3]
    J[3][3] = 0 + 1*-1*kr[3]*theta[8] + -1*kf[4]*(-1*theta[3] + 1**theta[8])
    J[3][4] = 0 + 1*kf[3]*theta[0] + -1*-1*kf[4]*theta[3]
    J[3][5] = 0 + -1*(-1*kf[4]*theta[3] - kr[4]*theta[7])
    J[3][6] = 0 + -1*-1*kf[4]*theta[3]
    J[3][7] = 0 + -1*(-1*kf[4]*theta[3] - kr[4]*theta[5])
    J[3][8] = 0 + 1*(kf[3]*-1*(theta[0]+theta[4]) - kr[3]*-1*(1)*theta[8]) + -1*(1*kf[4]*theta[3] - kr[4]*-1*(theta[7]+theta[5]))
    J[3][9] = 0 + 1*(kf[3]*-1*(theta[0]+theta[4]) - 1*kr[3]*theta[3]) + -1*(kf[4]*-1*(1)*theta[9] - kr[4]*-1*(theta[7]+theta[5]))
    J[4][0] = 0 + 1*(kf[2]*theta[6] - -1*kr[2]*theta[4]) + -1*(kf[3]*theta[4] - -1*kr[3]*theta[3])
    J[4][1] = 0
    J[4][2] = 0
    J[4][3] = 0 + -1*-1*kr[3]*theta[8]
    J[4][4] = 0 + 1*-1*kr[2]*theta[8] + -1*kf[3]*theta[0]
    J[4][5] = 0
    J[4][6] = 0 + 1*kf[2]*theta[0]
    J[4][7] = 0
    J[4][8] = 0 + 1*(kf[2]*-1*(theta[0]+theta[6]) - kr[2]*-1*(1)*theta[8]) + -1*(kf[3]*-1*(theta[0]+theta[4]) - kr[3]*-1*(1)*theta[8])
    J[4][9] = 0 + 1*(kf[2]*-1*(theta[0]+theta[6]) - 1*kr[2]*theta[4]) + -1*(kf[3]*-1*(theta[0]+theta[4]) - 1*kr[3]*theta[3])
    J[5][0] = 0 + -1*(kf[5]*theta[5] - -1*kr[5]*theta[1])
    J[5][1] = 0 + 1*-1*kf[4]*theta[3] + -1*-1*kr[5]*theta[8]
    J[5][2] = 0 + 1*-1*kf[4]*theta[3]
    J[5][3] = 0 + 1*kf[4]*(-1*theta[3] + 1**theta[8])
    J[5][4] = 0 + 1*-1*kf[4]*theta[3]
    J[5][5] = 0 + 1*(-1*kf[4]*theta[3] - kr[4]*theta[7]) + -1*kf[5]*theta[0]
    J[5][6] = 0 + 1*-1*kf[4]*theta[3]
    J[5][7] = 0 + 1*(-1*kf[4]*theta[3] - kr[4]*theta[5])
    J[5][8] = 0 + 1*(1*kf[4]*theta[3] - kr[4]*-1*(theta[7]+theta[5])) + -1*(kf[5]*-1*(theta[0]+theta[5]) - kr[5]*-1*(1)*theta[8])
    J[5][9] = 0 + 1*(kf[4]*-1*(1)*theta[9] - kr[4]*-1*(theta[7]+theta[5])) + -1*(kf[5]*-1*(theta[0]+theta[5]) - 1*kr[5]*theta[1])
    J[6][0] = 0 + -1*(kf[2]*theta[6] - -1*kr[2]*theta[4])
    J[6][1] = 0 + 1*-1*kf[1]*p[1]
    J[6][2] = 0 + 1*-1*kf[1]*p[1]
    J[6][3] = 0 + 1*-1*kf[1]*p[1]
    J[6][4] = 0 + 1*-1*kf[1]*p[1] + -1*-1*kr[2]*theta[8]
    J[6][5] = 0 + 1*-1*kf[1]*p[1]
    J[6][6] = 0 + 1*(-1*kf[1]*p[1] - kr[1]) + -1*kf[2]*theta[0]
    J[6][7] = 0 + 1*-1*kf[1]*p[1]
    J[6][8] = 0 + 1*(1*kf[1]*p[1] - kr[1]*-1*(1)) + -1*(kf[2]*-1*(theta[0]+theta[6]) - kr[2]*-1*(1)*theta[8])
    J[6][9] = 0 + 1*-1*kr[1]*-1*(1) + -1*(kf[2]*-1*(theta[0]+theta[6]) - 1*kr[2]*theta[4])
    J[7][0] = 0 + -1*(kf[8]*theta[7] - -1*kr[8]*p[2]*theta[9])
    J[7][1] = 0 + 1*-1*kf[4]*theta[3] + -1*-1*-1*kr[8]*p[2]*theta[8]
    J[7][2] = 0 + 1*-1*kf[4]*theta[3] + -1*-1*-1*kr[8]*p[2]*theta[8]
    J[7][3] = 0 + 1*kf[4]*(-1*theta[3] + 1**theta[8]) + -1*-1*-1*kr[8]*p[2]*theta[8]
    J[7][4] = 0 + 1*-1*kf[4]*theta[3] + -1*-1*-1*kr[8]*p[2]*theta[8]
    J[7][5] = 0 + 1*(-1*kf[4]*theta[3] - kr[4]*theta[7]) + -1*-1*-1*kr[8]*p[2]*theta[8]
    J[7][6] = 0 + 1*-1*kf[4]*theta[3] + -1*-1*-1*kr[8]*p[2]*theta[8]
    J[7][7] = 0 + 1*(-1*kf[4]*theta[3] - kr[4]*theta[5]) + -1*(kf[8]*theta[0] - -1*kr[8]*p[2]*theta[8])
    J[7][8] = 0 + 1*(1*kf[4]*theta[3] - kr[4]*-1*(theta[7]+theta[5])) + -1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[8])
    J[7][9] = 0 + 1*(kf[4]*-1*(1)*theta[9] - kr[4]*-1*(theta[7]+theta[5])) + -1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[9])
    J[8][0] = 0 + 1*(kf[7]*theta[2] - -1*kr[7]*p[0]*theta[9]) + 1*(kf[8]*theta[7] - -1*kr[8]*p[2]*theta[9])
    J[8][1] = 0 + -1*-1*kf[1]*p[1] + -1*-1*kf[4]*theta[3] + 1*-1*-1*kr[7]*p[0]*theta[8] + 1*-1*-1*kr[8]*p[2]*theta[8]
    J[8][2] = 0 + -1*-1*kf[1]*p[1] + -1*-1*kf[4]*theta[3] + 1*(kf[7]*theta[0] - -1*kr[7]*p[0]*theta[8]) + 1*-1*-1*kr[8]*p[2]*theta[8]
    J[8][3] = 0 + -1*-1*kf[1]*p[1] + -1*kf[4]*(-1*theta[3] + 1**theta[8]) + 1*-1*-1*kr[7]*p[0]*theta[8] + 1*-1*-1*kr[8]*p[2]*theta[8]
    J[8][4] = 0 + -1*-1*kf[1]*p[1] + -1*-1*kf[4]*theta[3] + 1*-1*-1*kr[7]*p[0]*theta[8] + 1*-1*-1*kr[8]*p[2]*theta[8]
    J[8][5] = 0 + -1*-1*kf[1]*p[1] + -1*(-1*kf[4]*theta[3] - kr[4]*theta[7]) + 1*-1*-1*kr[7]*p[0]*theta[8] + 1*-1*-1*kr[8]*p[2]*theta[8]
    J[8][6] = 0 + -1*(-1*kf[1]*p[1] - kr[1]) + -1*-1*kf[4]*theta[3] + 1*-1*-1*kr[7]*p[0]*theta[8] + 1*-1*-1*kr[8]*p[2]*theta[8]
    J[8][7] = 0 + -1*-1*kf[1]*p[1] + -1*(-1*kf[4]*theta[3] - kr[4]*theta[5]) + 1*-1*-1*kr[7]*p[0]*theta[8] + 1*(kf[8]*theta[0] - -1*kr[8]*p[2]*theta[8])
    J[8][8] = 0 + -1*(1*kf[1]*p[1] - kr[1]*-1*(1)) + -1*(1*kf[4]*theta[3] - kr[4]*-1*(theta[7]+theta[5])) + 1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[8]) + 1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[8])
    J[8][9] = 0 + -1*-1*kr[1]*-1*(1) + -1*(kf[4]*-1*(1)*theta[9] - kr[4]*-1*(theta[7]+theta[5])) + 1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[9]) + 1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[9])
    J[9][0] = 0 + -2*(-2*kf[0]*p[3]*theta[9] - 2*kr[0]*theta[0]) + 1*(kf[2]*theta[6] - -1*kr[2]*theta[4]) + 1*(kf[3]*theta[4] - -1*kr[3]*theta[3]) + 1*(kf[5]*theta[5] - -1*kr[5]*theta[1]) + 1*(kf[6]*theta[1] - -1*kr[6]*theta[2]) + 1*(kf[7]*theta[2] - -1*kr[7]*p[0]*theta[9]) + 1*(kf[8]*theta[7] - -1*kr[8]*p[2]*theta[9])
    J[9][1] = 0 + 1*-1*kr[5]*theta[8] + 1*kf[6]*theta[0] + 1*-1*-1*kr[7]*p[0]*theta[8] + 1*-1*-1*kr[8]*p[2]*theta[8]
    J[9][2] = 0 + 1*-1*kr[6]*theta[8] + 1*(kf[7]*theta[0] - -1*kr[7]*p[0]*theta[8]) + 1*-1*-1*kr[8]*p[2]*theta[8]
    J[9][3] = 0 + 1*-1*kr[3]*theta[8] + 1*-1*-1*kr[7]*p[0]*theta[8] + 1*-1*-1*kr[8]*p[2]*theta[8]
    J[9][4] = 0 + 1*-1*kr[2]*theta[8] + 1*kf[3]*theta[0] + 1*-1*-1*kr[7]*p[0]*theta[8] + 1*-1*-1*kr[8]*p[2]*theta[8]
    J[9][5] = 0 + 1*kf[5]*theta[0] + 1*-1*-1*kr[7]*p[0]*theta[8] + 1*-1*-1*kr[8]*p[2]*theta[8]
    J[9][6] = 0 + 1*kf[2]*theta[0] + 1*-1*-1*kr[7]*p[0]*theta[8] + 1*-1*-1*kr[8]*p[2]*theta[8]
    J[9][7] = 0 + 1*-1*-1*kr[7]*p[0]*theta[8] + 1*(kf[8]*theta[0] - -1*kr[8]*p[2]*theta[8])
    J[9][8] = 0 + -2*-1*kr[0]*-1*(theta[0]+theta[0]) + 1*(kf[2]*-1*(theta[0]+theta[6]) - kr[2]*-1*(1)*theta[8]) + 1*(kf[3]*-1*(theta[0]+theta[4]) - kr[3]*-1*(1)*theta[8]) + 1*(kf[5]*-1*(theta[0]+theta[5]) - kr[5]*-1*(1)*theta[8]) + 1*(kf[6]*-1*(theta[0]+theta[1]) - kr[6]*-1*(1)*theta[8]) + 1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[8]) + 1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[8])
    J[9][9] = 0 + -2*(2*kf[0]*p[3]*theta[9] - kr[0]*-1*(theta[0]+theta[0])) + 1*(kf[2]*-1*(theta[0]+theta[6]) - 1*kr[2]*theta[4]) + 1*(kf[3]*-1*(theta[0]+theta[4]) - 1*kr[3]*theta[3]) + 1*(kf[5]*-1*(theta[0]+theta[5]) - 1*kr[5]*theta[1]) + 1*(kf[6]*-1*(theta[0]+theta[1]) - 1*kr[6]*theta[2]) + 1*(kf[7]*-1*(theta[0]+theta[2]) - 1*kr[7]*p[0]*theta[9]) + 1*(kf[8]*-1*(theta[0]+theta[7]) - 1*kr[8]*p[2]*theta[9])

    J = matrix(J)
    return J




def constrain_coverage_function(cvgs,mpf,c_min):
    cvgs = [max(ci,c_min) for ci in cvgs]

    n_adsorbates = 8
    site_info_dict = {'s': [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1, {'slope': 1.0, 'cutoff': 0.66, 'smoothing': 0.05}], 's&h': [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1.0, {'slope': 1.0, 'cutoff': 2.0, 'smoothing': 0.05}], 'h': [[0], 1.0, {'slope': 1.0, 'cutoff': 2.0, 'smoothing': 0.05}]}
    max_coverage_list = [1.0, 1, 1, 1, 1, 1, 1, 1]

    #enforce explicit maxima
    cvgs = [min(ci,maxi) for ci,maxi in zip(cvgs,max_coverage_list)]

    #enforce site conservation
    single_sites = [s for s in site_info_dict if '&' not in s]
    for s in single_sites:
        idxs = [idx for idx in site_info_dict[s][0] if idx < n_adsorbates]
        tot = site_info_dict[s][1]
        sum_a = sum([cvgs[id] for id in idxs])
        if sum_a > tot:
            for id in idxs:
                cvgs[id] = cvgs[id]/sum_a
                cvgs[id] = cvgs[id]*tot
    return cvgs




def elementary_rates(rate_constants,theta,p,mpf,matrix):

    # using // for integer division
    # ensuring integer result under Py2 and Py3
    kf = rate_constants[0:len(rate_constants)//2]
    kr = rate_constants[len(rate_constants)//2:]

    r = matrix([0]*len(kf))
    dtheta_dt = matrix([0]*10)
    
    s = [0]*3
    s[0] = (mpf('1.0') - theta[1] - theta[2] - theta[3] - theta[4] - theta[5] - theta[6] - theta[7])
    s[1] = (mpf('0.0'))
    s[2] = (mpf('1.0') - theta[0])
    r[0] = kf[0]*p[3]*theta[9]*theta[9] - kr[0]*theta[0]*theta[0]
    r[1] = kf[1]*p[1]*theta[8] - kr[1]*theta[6]
    r[2] = kf[2]*theta[6]*theta[0] - kr[2]*theta[4]*theta[9]
    r[3] = kf[3]*theta[4]*theta[0] - kr[3]*theta[3]*theta[9]
    r[4] = kf[4]*theta[3]*theta[8] - kr[4]*theta[5]*theta[7]
    r[5] = kf[5]*theta[5]*theta[0] - kr[5]*theta[1]*theta[9]
    r[6] = kf[6]*theta[1]*theta[0] - kr[6]*theta[2]*theta[9]
    r[7] = kf[7]*theta[2]*theta[0] - kr[7]*p[0]*theta[8]*theta[9]
    r[8] = kf[8]*theta[7]*theta[0] - kr[8]*p[2]*theta[8]*theta[9]
    dtheta_dt[0] =  + 2*r[0] + -1*r[2] + -1*r[3] + -1*r[5] + -1*r[6] + -1*r[7] + -1*r[8]
    dtheta_dt[1] =  + 1*r[5] + -1*r[6]
    dtheta_dt[2] =  + 1*r[6] + -1*r[7]
    dtheta_dt[3] =  + 1*r[3] + -1*r[4]
    dtheta_dt[4] =  + 1*r[2] + -1*r[3]
    dtheta_dt[5] =  + 1*r[4] + -1*r[5]
    dtheta_dt[6] =  + 1*r[1] + -1*r[2]
    dtheta_dt[7] =  + 1*r[4] + -1*r[8]
    dtheta_dt[8] =  + -1*r[1] + -1*r[4] + 1*r[7] + 1*r[8]
    dtheta_dt[9] =  + -2*r[0] + 1*r[2] + 1*r[3] + 1*r[5] + 1*r[6] + 1*r[7] + 1*r[8]
    
    return r



