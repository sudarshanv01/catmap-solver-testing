
def rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,mpf,matrix,mpexp,mpsqrt,include_derivatives=True):
    
    def interaction_function(coverages,energies,interaction_vector,F,include_derivatives=True,include_integral=False): 
    
    #    Function for evaluating coverage-dependent intearction energies. 
    #
    #    -coverages is a vector of n_coverages coverages (thetas).
    #    -energies is a vector of energies (length n_coverages).
    #    -interaction_vector is a raveled matrix of self and cross interaction parameters
    #    (length n_coverages^2).
    #    -include_derivatives: True or False
    #
    #    Function evaluates:
    #
    #    E_i = E_0 + F(theta_tot)*sum_j(epsilon_ij*theta_j)
    #
    #    dE_i/dtheta_n = sum_j(dF/dtheta_tot*dtheta_tot/dtheta_n*epsilon_ij*theta_j + F*epsilon_ij*delta_jn)
    #
    #    where sum_j is summation over index j, epsilon_ij is the interaction parameter between 
    #    adsorbate i and adsorbat j, theta_tot is a sum on theta, delta_nj is the kronecker delta,
    #    theta_j is the coverage of adsorbate j, and F is the "response function".
    
        Es = []
        dEs = []
        n_cvg = len(coverages)
    
    #    Dictionary with site names as keys and values a list of:
    #    [indices_corresponding_to_site, max_coverage_of_site, piecewise_linearity_threshold]
    
        site_info_dict = {'s': [[0, 1, 2, 3], 1, {'cutoff': 0.25, 'smoothing': 0.01}]}
    
        f_vector = [0]*len(coverages)
        df_vector = [0]*len(coverages)
        sites = [0]*len(coverages)
        c_tots = [0]*len(coverages)
        for s in site_info_dict:
            idxs, max_cvg, F_params = site_info_dict[s]
            cvgs = [coverages[j] for j in idxs]
            c_tot_i = sum(cvgs)
            f, df = F(c_tot_i,**F_params)
            for idx in idxs:
                f_vector[idx] = f
                df_vector[idx] = df
                sites[idx] = s
                c_tots[idx] = c_tot_i
        Es = []
        dEs = []
    
        for i,theta_i in enumerate(coverages):
            E_0 = energies[i]
            epsilon_vector = interaction_vector[i*n_cvg:(i+1)*n_cvg]
            f_epsilon_dot_theta = 0
            for ep,nc,f in zip(epsilon_vector,coverages,f_vector):
                f_epsilon_dot_theta += f*ep*nc
    
            E_i = E_0 + f_epsilon_dot_theta
            Es.append(E_i)
    
            if include_derivatives:
                diff_vector = []
                for n,theta_n in enumerate(coverages):
                    df_sum = 0
                    f_sum = 0
                    for j,theta_j in enumerate(coverages):
                        df= df_vector[j]
                        f = f_vector[j]
                        ep = epsilon_vector[j]
                        if j==n:
                            delta = 1
                        else:
                            delta = 0
                        if sites[j] == sites[n]:
                            dsum = 1
                        else:
                            dsum = 0
                        df_sum += df*dsum*ep*theta_j
                        f_sum += f*ep*delta
                    diff_vector.append(df_sum+f_sum)
                dEs.append(diff_vector)
            else:
                dEs = None
        E_int = None
        if include_integral:
            raise UserWarning('First-order interactions are not '
                                'compatible with integral energies')
        return E_int, Es, dEs
        

    kfs = []
    krs = []
    dEfs = []
    dErs = []

    n_adsorbates = 2
    n_transition_states = 2
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

    kB = mpf('0.000086173324779999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999992')
    h = mpf('0.0000000000000041356675159999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999968')
    prefactor_list = [kB*T/h, kB*T/h, kB*T/h]

    def element_wise_addition(lists):
        return [sum(L) for L in zip(*lists)]

    G_IS = [0]*3
    G_TS = [0]*3
    G_FS = [0]*3
    G_af = [0]*3
    G_ar = [0]*3
    dG_IS = [0]*3
    dG_TS = [0]*3
    dG_FS = [0]*3
    G_IS[0] = site_energies[1] + gas_energies[1]
    G_FS[0] = Gf[0]
    dG_IS[0] = element_wise_addition([[0]*2])
    dG_FS[0] = element_wise_addition([dGs[0]])
    G_TS[0] = max([G_IS[0],G_FS[0]])
    dG_TS[0] = None #determined later
    
    G_IS[1] = site_energies[1] + site_energies[1] + gas_energies[2]
    G_FS[1] = Gf[1] + Gf[1]
    dG_IS[1] = element_wise_addition([[0]*2])
    dG_FS[1] = element_wise_addition([dGs[1] , dGs[1]])
    G_TS[1] = Gf[3] + site_energies[1]
    dG_TS[1] = element_wise_addition([dGs[3]])
    
    G_IS[2] = Gf[0] + Gf[1]
    G_FS[2] = gas_energies[0] + site_energies[1] + site_energies[1]
    dG_IS[2] = element_wise_addition([dGs[0] , dGs[1]])
    dG_FS[2] = element_wise_addition([[0]*2])
    G_TS[2] = Gf[2] + site_energies[1]
    dG_TS[2] = element_wise_addition([dGs[2]])
    

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




def interaction_function(coverages,energies,interaction_vector,F,include_derivatives=True,include_integral=False): 

#    Function for evaluating coverage-dependent intearction energies. 
#
#    -coverages is a vector of n_coverages coverages (thetas).
#    -energies is a vector of energies (length n_coverages).
#    -interaction_vector is a raveled matrix of self and cross interaction parameters
#    (length n_coverages^2).
#    -include_derivatives: True or False
#
#    Function evaluates:
#
#    E_i = E_0 + F(theta_tot)*sum_j(epsilon_ij*theta_j)
#
#    dE_i/dtheta_n = sum_j(dF/dtheta_tot*dtheta_tot/dtheta_n*epsilon_ij*theta_j + F*epsilon_ij*delta_jn)
#
#    where sum_j is summation over index j, epsilon_ij is the interaction parameter between 
#    adsorbate i and adsorbat j, theta_tot is a sum on theta, delta_nj is the kronecker delta,
#    theta_j is the coverage of adsorbate j, and F is the "response function".

    Es = []
    dEs = []
    n_cvg = len(coverages)

#    Dictionary with site names as keys and values a list of:
#    [indices_corresponding_to_site, max_coverage_of_site, piecewise_linearity_threshold]

    site_info_dict = {'s': [[0, 1, 2, 3], 1, {'cutoff': 0.25, 'smoothing': 0.01}]}

    f_vector = [0]*len(coverages)
    df_vector = [0]*len(coverages)
    sites = [0]*len(coverages)
    c_tots = [0]*len(coverages)
    for s in site_info_dict:
        idxs, max_cvg, F_params = site_info_dict[s]
        cvgs = [coverages[j] for j in idxs]
        c_tot_i = sum(cvgs)
        f, df = F(c_tot_i,**F_params)
        for idx in idxs:
            f_vector[idx] = f
            df_vector[idx] = df
            sites[idx] = s
            c_tots[idx] = c_tot_i
    Es = []
    dEs = []

    for i,theta_i in enumerate(coverages):
        E_0 = energies[i]
        epsilon_vector = interaction_vector[i*n_cvg:(i+1)*n_cvg]
        f_epsilon_dot_theta = 0
        for ep,nc,f in zip(epsilon_vector,coverages,f_vector):
            f_epsilon_dot_theta += f*ep*nc

        E_i = E_0 + f_epsilon_dot_theta
        Es.append(E_i)

        if include_derivatives:
            diff_vector = []
            for n,theta_n in enumerate(coverages):
                df_sum = 0
                f_sum = 0
                for j,theta_j in enumerate(coverages):
                    df= df_vector[j]
                    f = f_vector[j]
                    ep = epsilon_vector[j]
                    if j==n:
                        delta = 1
                    else:
                        delta = 0
                    if sites[j] == sites[n]:
                        dsum = 1
                    else:
                        dsum = 0
                    df_sum += df*dsum*ep*theta_j
                    f_sum += f*ep*delta
                diff_vector.append(df_sum+f_sum)
            dEs.append(diff_vector)
        else:
            dEs = None
    E_int = None
    if include_integral:
        raise UserWarning('First-order interactions are not '
                            'compatible with integral energies')
    return E_int, Es, dEs
    



def interacting_mean_field_steady_state(rxn_parameters,theta,p,gas_energies,site_energies,T,F,mpf,matrix,mpexp,mpsqrt):

    
    def rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,mpf,matrix,mpexp,mpsqrt,include_derivatives=True):
        
        def interaction_function(coverages,energies,interaction_vector,F,include_derivatives=True,include_integral=False): 
        
        #    Function for evaluating coverage-dependent intearction energies. 
        #
        #    -coverages is a vector of n_coverages coverages (thetas).
        #    -energies is a vector of energies (length n_coverages).
        #    -interaction_vector is a raveled matrix of self and cross interaction parameters
        #    (length n_coverages^2).
        #    -include_derivatives: True or False
        #
        #    Function evaluates:
        #
        #    E_i = E_0 + F(theta_tot)*sum_j(epsilon_ij*theta_j)
        #
        #    dE_i/dtheta_n = sum_j(dF/dtheta_tot*dtheta_tot/dtheta_n*epsilon_ij*theta_j + F*epsilon_ij*delta_jn)
        #
        #    where sum_j is summation over index j, epsilon_ij is the interaction parameter between 
        #    adsorbate i and adsorbat j, theta_tot is a sum on theta, delta_nj is the kronecker delta,
        #    theta_j is the coverage of adsorbate j, and F is the "response function".
        
            Es = []
            dEs = []
            n_cvg = len(coverages)
        
        #    Dictionary with site names as keys and values a list of:
        #    [indices_corresponding_to_site, max_coverage_of_site, piecewise_linearity_threshold]
        
            site_info_dict = {'s': [[0, 1, 2, 3], 1, {'cutoff': 0.25, 'smoothing': 0.01}]}
        
            f_vector = [0]*len(coverages)
            df_vector = [0]*len(coverages)
            sites = [0]*len(coverages)
            c_tots = [0]*len(coverages)
            for s in site_info_dict:
                idxs, max_cvg, F_params = site_info_dict[s]
                cvgs = [coverages[j] for j in idxs]
                c_tot_i = sum(cvgs)
                f, df = F(c_tot_i,**F_params)
                for idx in idxs:
                    f_vector[idx] = f
                    df_vector[idx] = df
                    sites[idx] = s
                    c_tots[idx] = c_tot_i
            Es = []
            dEs = []
        
            for i,theta_i in enumerate(coverages):
                E_0 = energies[i]
                epsilon_vector = interaction_vector[i*n_cvg:(i+1)*n_cvg]
                f_epsilon_dot_theta = 0
                for ep,nc,f in zip(epsilon_vector,coverages,f_vector):
                    f_epsilon_dot_theta += f*ep*nc
        
                E_i = E_0 + f_epsilon_dot_theta
                Es.append(E_i)
        
                if include_derivatives:
                    diff_vector = []
                    for n,theta_n in enumerate(coverages):
                        df_sum = 0
                        f_sum = 0
                        for j,theta_j in enumerate(coverages):
                            df= df_vector[j]
                            f = f_vector[j]
                            ep = epsilon_vector[j]
                            if j==n:
                                delta = 1
                            else:
                                delta = 0
                            if sites[j] == sites[n]:
                                dsum = 1
                            else:
                                dsum = 0
                            df_sum += df*dsum*ep*theta_j
                            f_sum += f*ep*delta
                        diff_vector.append(df_sum+f_sum)
                    dEs.append(diff_vector)
                else:
                    dEs = None
            E_int = None
            if include_integral:
                raise UserWarning('First-order interactions are not '
                                    'compatible with integral energies')
            return E_int, Es, dEs
            
    
        kfs = []
        krs = []
        dEfs = []
        dErs = []
    
        n_adsorbates = 2
        n_transition_states = 2
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
    
        kB = mpf('0.000086173324779999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999992')
        h = mpf('0.0000000000000041356675159999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999968')
        prefactor_list = [kB*T/h, kB*T/h, kB*T/h]
    
        def element_wise_addition(lists):
            return [sum(L) for L in zip(*lists)]
    
        G_IS = [0]*3
        G_TS = [0]*3
        G_FS = [0]*3
        G_af = [0]*3
        G_ar = [0]*3
        G_IS[0] = site_energies[1] + gas_energies[1]
        G_FS[0] = Gf[0]
        G_TS[0] = max([G_IS[0],G_FS[0]])
        
        G_IS[1] = site_energies[1] + site_energies[1] + gas_energies[2]
        G_FS[1] = Gf[1] + Gf[1]
        G_TS[1] = Gf[3] + site_energies[1]
        
        G_IS[2] = Gf[0] + Gf[1]
        G_FS[2] = gas_energies[0] + site_energies[1] + site_energies[1]
        G_TS[2] = Gf[2] + site_energies[1]
        
    
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
    dtheta_dt = [0]*3
    
    s = [0]*2
    s[0] = (mpf('0.0'))
    s[1] = (mpf('1.0') - theta[0] - theta[1])
    r[0] = kf[0]*p[1]*s[1] - kr[0]*theta[0]
    r[1] = kf[1]*p[2]*s[1]*s[1] - kr[1]*theta[1]*theta[1]
    r[2] = kf[2]*theta[0]*theta[1] - kr[2]*p[0]*s[1]*s[1]
    dtheta_dt[0] =  + 1*r[0] + -1*r[2]
    dtheta_dt[1] =  + 2*r[1] + -1*r[2]
    dtheta_dt[2] =  + -1*r[0] + -2*r[1] + 2*r[2]
 
    r = matrix(r)
    dtheta_dt = matrix(dtheta_dt)

    if dtheta_dt.rows == len(theta)+1:
        dtheta_dt = dtheta_dt[:-1] 

    return dtheta_dt




def ideal_mean_field_steady_state(kf,kr,theta,p,mpf,matrix):

    r = [0]*len(kf)
    dtheta_dt = [0]*3
    s = [0]*2
    s[0] = (mpf('0.0'))
    s[1] = (mpf('1.0') - theta[0] - theta[1])
    r[0] = kf[0]*p[1]*s[1] - kr[0]*theta[0]
    r[1] = kf[1]*p[2]*s[1]*s[1] - kr[1]*theta[1]*theta[1]
    r[2] = kf[2]*theta[0]*theta[1] - kr[2]*p[0]*s[1]*s[1]
    dtheta_dt[0] =  + 1*r[0] + -1*r[2]
    dtheta_dt[1] =  + 2*r[1] + -1*r[2]
    dtheta_dt[2] =  + -1*r[0] + -2*r[1] + 2*r[2]

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
    kB = mpf('0.000086173324779999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999992')
    n_adsorbates = 2
    kBT = kB*T


    J = [[0 for i in range(3)] for j in range(3)]

    
    def rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,mpf,matrix,mpexp,mpsqrt,include_derivatives=True):
        
        def interaction_function(coverages,energies,interaction_vector,F,include_derivatives=True,include_integral=False): 
        
        #    Function for evaluating coverage-dependent intearction energies. 
        #
        #    -coverages is a vector of n_coverages coverages (thetas).
        #    -energies is a vector of energies (length n_coverages).
        #    -interaction_vector is a raveled matrix of self and cross interaction parameters
        #    (length n_coverages^2).
        #    -include_derivatives: True or False
        #
        #    Function evaluates:
        #
        #    E_i = E_0 + F(theta_tot)*sum_j(epsilon_ij*theta_j)
        #
        #    dE_i/dtheta_n = sum_j(dF/dtheta_tot*dtheta_tot/dtheta_n*epsilon_ij*theta_j + F*epsilon_ij*delta_jn)
        #
        #    where sum_j is summation over index j, epsilon_ij is the interaction parameter between 
        #    adsorbate i and adsorbat j, theta_tot is a sum on theta, delta_nj is the kronecker delta,
        #    theta_j is the coverage of adsorbate j, and F is the "response function".
        
            Es = []
            dEs = []
            n_cvg = len(coverages)
        
        #    Dictionary with site names as keys and values a list of:
        #    [indices_corresponding_to_site, max_coverage_of_site, piecewise_linearity_threshold]
        
            site_info_dict = {'s': [[0, 1, 2, 3], 1, {'cutoff': 0.25, 'smoothing': 0.01}]}
        
            f_vector = [0]*len(coverages)
            df_vector = [0]*len(coverages)
            sites = [0]*len(coverages)
            c_tots = [0]*len(coverages)
            for s in site_info_dict:
                idxs, max_cvg, F_params = site_info_dict[s]
                cvgs = [coverages[j] for j in idxs]
                c_tot_i = sum(cvgs)
                f, df = F(c_tot_i,**F_params)
                for idx in idxs:
                    f_vector[idx] = f
                    df_vector[idx] = df
                    sites[idx] = s
                    c_tots[idx] = c_tot_i
            Es = []
            dEs = []
        
            for i,theta_i in enumerate(coverages):
                E_0 = energies[i]
                epsilon_vector = interaction_vector[i*n_cvg:(i+1)*n_cvg]
                f_epsilon_dot_theta = 0
                for ep,nc,f in zip(epsilon_vector,coverages,f_vector):
                    f_epsilon_dot_theta += f*ep*nc
        
                E_i = E_0 + f_epsilon_dot_theta
                Es.append(E_i)
        
                if include_derivatives:
                    diff_vector = []
                    for n,theta_n in enumerate(coverages):
                        df_sum = 0
                        f_sum = 0
                        for j,theta_j in enumerate(coverages):
                            df= df_vector[j]
                            f = f_vector[j]
                            ep = epsilon_vector[j]
                            if j==n:
                                delta = 1
                            else:
                                delta = 0
                            if sites[j] == sites[n]:
                                dsum = 1
                            else:
                                dsum = 0
                            df_sum += df*dsum*ep*theta_j
                            f_sum += f*ep*delta
                        diff_vector.append(df_sum+f_sum)
                    dEs.append(diff_vector)
                else:
                    dEs = None
            E_int = None
            if include_integral:
                raise UserWarning('First-order interactions are not '
                                    'compatible with integral energies')
            return E_int, Es, dEs
            
    
        kfs = []
        krs = []
        dEfs = []
        dErs = []
    
        n_adsorbates = 2
        n_transition_states = 2
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
    
        kB = mpf('0.000086173324779999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999992')
        h = mpf('0.0000000000000041356675159999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999968')
        prefactor_list = [kB*T/h, kB*T/h, kB*T/h]
    
        def element_wise_addition(lists):
            return [sum(L) for L in zip(*lists)]
    
        G_IS = [0]*3
        G_TS = [0]*3
        G_FS = [0]*3
        G_af = [0]*3
        G_ar = [0]*3
        dG_IS = [0]*3
        dG_TS = [0]*3
        dG_FS = [0]*3
        G_IS[0] = site_energies[1] + gas_energies[1]
        G_FS[0] = Gf[0]
        dG_IS[0] = element_wise_addition([[0]*2])
        dG_FS[0] = element_wise_addition([dGs[0]])
        G_TS[0] = max([G_IS[0],G_FS[0]])
        dG_TS[0] = None #determined later
        
        G_IS[1] = site_energies[1] + site_energies[1] + gas_energies[2]
        G_FS[1] = Gf[1] + Gf[1]
        dG_IS[1] = element_wise_addition([[0]*2])
        dG_FS[1] = element_wise_addition([dGs[1] , dGs[1]])
        G_TS[1] = Gf[3] + site_energies[1]
        dG_TS[1] = element_wise_addition([dGs[3]])
        
        G_IS[2] = Gf[0] + Gf[1]
        G_FS[2] = gas_energies[0] + site_energies[1] + site_energies[1]
        dG_IS[2] = element_wise_addition([dGs[0] , dGs[1]])
        dG_FS[2] = element_wise_addition([[0]*2])
        G_TS[2] = Gf[2] + site_energies[1]
        dG_TS[2] = element_wise_addition([dGs[2]])
        
    
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
    s = [0]*2
    s[0] = (mpf('0.0'))
    s[1] = (mpf('1.0') - theta[0] - theta[1])
    kfkBT = [0]*3
    krkBT = [0]*3
    kfkBT[0] = kf[0]/kBT
    krkBT[0] = kr[0]/kBT
    kfkBT[1] = kf[1]/kBT
    krkBT[1] = kr[1]/kBT
    kfkBT[2] = kf[2]/kBT
    krkBT[2] = kr[2]/kBT
    J[0][0] = 0 + 1*(-1*kf[0]*p[1] + (kfkBT[0])*dEf[0][0]*p[1]*s[1] - kr[0] - (krkBT[0])*dEr[0][0]*theta[0]) + -1*(kf[2]*theta[1] + (kfkBT[2])*dEf[2][0]*theta[0]*theta[1] - -2*kr[2]*p[0]*s[1] - (krkBT[2])*dEr[2][0]*p[0]*s[1]*s[1])
    J[0][1] = 0 + 1*(-1*kf[0]*p[1] + (kfkBT[0])*dEf[0][1]*p[1]*s[1] - 0 - (krkBT[0])*dEr[0][1]*theta[0]) + -1*(kf[2]*theta[0] + (kfkBT[2])*dEf[2][1]*theta[0]*theta[1] - -2*kr[2]*p[0]*s[1] - (krkBT[2])*dEr[2][1]*p[0]*s[1]*s[1])
    J[0][2] = 0 + 1*(1*kf[0]*p[1] - kr[0]*-1*(1)) + -1*(kf[2]*-1*(theta[1]+theta[0]) - 2*kr[2]*p[0]*s[1])
    J[1][0] = 0 + 2*(-2*kf[1]*p[2]*s[1] + (kfkBT[1])*dEf[1][0]*p[2]*s[1]*s[1] - 0 - (krkBT[1])*dEr[1][0]*theta[1]*theta[1]) + -1*(kf[2]*theta[1] + (kfkBT[2])*dEf[2][0]*theta[0]*theta[1] - -2*kr[2]*p[0]*s[1] - (krkBT[2])*dEr[2][0]*p[0]*s[1]*s[1])
    J[1][1] = 0 + 2*(-2*kf[1]*p[2]*s[1] + (kfkBT[1])*dEf[1][1]*p[2]*s[1]*s[1] - 2*kr[1]*theta[1] - (krkBT[1])*dEr[1][1]*theta[1]*theta[1]) + -1*(kf[2]*theta[0] + (kfkBT[2])*dEf[2][1]*theta[0]*theta[1] - -2*kr[2]*p[0]*s[1] - (krkBT[2])*dEr[2][1]*p[0]*s[1]*s[1])
    J[1][2] = 0 + 2*(2*kf[1]*p[2]*s[1] - kr[1]*-1*(theta[1]+theta[1])) + -1*(kf[2]*-1*(theta[1]+theta[0]) - 2*kr[2]*p[0]*s[1])
    J[2][0] = 0 + -1*(-1*kf[0]*p[1] + (kfkBT[0])*dEf[0][0]*p[1]*s[1] - kr[0] - (krkBT[0])*dEr[0][0]*theta[0]) + -2*(-2*kf[1]*p[2]*s[1] + (kfkBT[1])*dEf[1][0]*p[2]*s[1]*s[1] - 0 - (krkBT[1])*dEr[1][0]*theta[1]*theta[1]) + 2*(kf[2]*theta[1] + (kfkBT[2])*dEf[2][0]*theta[0]*theta[1] - -2*kr[2]*p[0]*s[1] - (krkBT[2])*dEr[2][0]*p[0]*s[1]*s[1])
    J[2][1] = 0 + -1*(-1*kf[0]*p[1] + (kfkBT[0])*dEf[0][1]*p[1]*s[1] - 0 - (krkBT[0])*dEr[0][1]*theta[0]) + -2*(-2*kf[1]*p[2]*s[1] + (kfkBT[1])*dEf[1][1]*p[2]*s[1]*s[1] - 2*kr[1]*theta[1] - (krkBT[1])*dEr[1][1]*theta[1]*theta[1]) + 2*(kf[2]*theta[0] + (kfkBT[2])*dEf[2][1]*theta[0]*theta[1] - -2*kr[2]*p[0]*s[1] - (krkBT[2])*dEr[2][1]*p[0]*s[1]*s[1])
    J[2][2] = 0 + -1*(1*kf[0]*p[1] - kr[0]*-1*(1)) + -2*(2*kf[1]*p[2]*s[1] - kr[1]*-1*(theta[1]+theta[1])) + 2*(kf[2]*-1*(theta[1]+theta[0]) - 2*kr[2]*p[0]*s[1])
    
    J = matrix(J)
    return J




def ideal_mean_field_jacobian(kf,kr,theta,p,mpf,matrix):
    n_adsorbates = 2
    J = [[0 for i in range(3)] for j in range(3)]

    s = [0]*2
    s[0] = (mpf('0.0'))
    s[1] = (mpf('1.0') - theta[0] - theta[1])
    J[0][0] = 0 + 1*(-1*kf[0]*p[1] - kr[0]) + -1*(kf[2]*theta[1] - -2*kr[2]*p[0]*s[1])
    J[0][1] = 0 + 1*-1*kf[0]*p[1] + -1*(kf[2]*theta[0] - -2*kr[2]*p[0]*s[1])
    J[0][2] = 0 + 1*(1*kf[0]*p[1] - kr[0]*-1*(1)) + -1*(kf[2]*-1*(theta[1]+theta[0]) - 2*kr[2]*p[0]*s[1])
    J[1][0] = 0 + 2*-2*kf[1]*p[2]*s[1] + -1*(kf[2]*theta[1] - -2*kr[2]*p[0]*s[1])
    J[1][1] = 0 + 2*(-2*kf[1]*p[2]*s[1] - 2*kr[1]*theta[1]) + -1*(kf[2]*theta[0] - -2*kr[2]*p[0]*s[1])
    J[1][2] = 0 + 2*(2*kf[1]*p[2]*s[1] - kr[1]*-1*(theta[1]+theta[1])) + -1*(kf[2]*-1*(theta[1]+theta[0]) - 2*kr[2]*p[0]*s[1])
    J[2][0] = 0 + -1*(-1*kf[0]*p[1] - kr[0]) + -2*-2*kf[1]*p[2]*s[1] + 2*(kf[2]*theta[1] - -2*kr[2]*p[0]*s[1])
    J[2][1] = 0 + -1*-1*kf[0]*p[1] + -2*(-2*kf[1]*p[2]*s[1] - 2*kr[1]*theta[1]) + 2*(kf[2]*theta[0] - -2*kr[2]*p[0]*s[1])
    J[2][2] = 0 + -1*(1*kf[0]*p[1] - kr[0]*-1*(1)) + -2*(2*kf[1]*p[2]*s[1] - kr[1]*-1*(theta[1]+theta[1])) + 2*(kf[2]*-1*(theta[1]+theta[0]) - 2*kr[2]*p[0]*s[1])

    J = matrix(J)
    return J




def constrain_coverage_function(cvgs,mpf,c_min):
    cvgs = [max(ci,c_min) for ci in cvgs]

    n_adsorbates = 2
    site_info_dict = {'s': [[0, 1, 2, 3], 1, {'cutoff': 0.25, 'smoothing': 0.01}]}
    max_coverage_list = [1, 1]

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
    dtheta_dt = matrix([0]*3)
    
    s = [0]*2
    s[0] = (mpf('0.0'))
    s[1] = (mpf('1.0') - theta[0] - theta[1])
    r[0] = kf[0]*p[1]*s[1] - kr[0]*theta[0]
    r[1] = kf[1]*p[2]*s[1]*s[1] - kr[1]*theta[1]*theta[1]
    r[2] = kf[2]*theta[0]*theta[1] - kr[2]*p[0]*s[1]*s[1]
    dtheta_dt[0] =  + 1*r[0] + -1*r[2]
    dtheta_dt[1] =  + 2*r[1] + -1*r[2]
    dtheta_dt[2] =  + -1*r[0] + -2*r[1] + 2*r[2]
    
    return r



