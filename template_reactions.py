templates = {}

templates[
    "mkm_job"
] = r"""

    from catmap import ReactionModel
    import os
    import json

    mkm_file = 'input.mkm'
    model = ReactionModel(setup_file=mkm_file)

    # Decide on the solver based on the current
    # working directory.
    folder_name = os.getcwd().split('/')[-1]
    if 'coverage' in folder_name:
        # This is default solver.
        model.use_numbers_solver = False
        model.fix_x_star = True 
    elif 'numbers_fix_xstar' in folder_name:
        # This is the numbers solver with fixed x*=1
        model.use_numbers_solver = True
        model.fix_x_star = True
    elif 'numbers_free_xstar' in folder_name:
        # This is the numbers solver with free x*
        model.use_numbers_solver = True
        model.fix_x_star = False
    else:
        raise ValueError('Unknown folder name. '
                        'Use coverage or numbers_fix_xstar or numbers_free_xstar.')

    # Import the json file with the solver specifics
    with open('solver_specifics.json') as f:
        solver_specifics = json.load(f)

    # Use the same solver specifications for all test cases.
    model.decimal_precision = solver_specifics['decimal_precision'] 
    model.tolerance = solver_specifics['tolerance']
    model.max_rootfinding_iterations = solver_specifics['max_rootfinding_iterations']
    model.max_bisections = solver_specifics['max_bisections']
    model.max_damping_iterations = solver_specifics['max_damping_iterations']

    # Run in DEBUG mode
    model.DEBUG = True

    model.run()

"""

templates[
    "mkm_CO_oxidation"
] = r"""

    rxn_expressions = [
                '*_s + CO_g -> CO*', 
                '2*_s + O2_g <-> O-O* + *_s -> 2O*',
                'CO* +  O* <-> O-CO* + * -> CO2_g + 2*',
                    ]

    surface_names = ['Pt', 'Ag', 'Cu','Rh','Pd','Au','Ru','Ni'] 

    descriptor_names= ['O_s','CO_s'] #descriptor names

    descriptor_ranges = [[${desc1}, ${desc1}],[${desc2}, ${desc2}]]

    temperature = 500 #Temperature of the reaction

    species_definitions = {}
    species_definitions['CO_g'] = {'pressure':1.} #define the gas pressures
    species_definitions['O2_g'] = {'pressure':1./3.}
    species_definitions['CO2_g'] = {'pressure':0}
    species_definitions['s'] = {'site_names': ['111'], 'total':1} #define the sites

    data_file = 'data.pkl'

    # Parser parameters

    input_file = 'energies.txt'

    # Scaler parameters
    gas_thermo_mode = "shomate_gas"
    adsorbate_thermo_mode = "frozen_adsorbate"
    scaling_constraint_dict = {
                            'O_s':['+',0,None],
                            'CO_s':[0,'+',None],
                            'O-CO_s':'initial_state',
                            'O-O_s':'final_state',
                            }
    """
