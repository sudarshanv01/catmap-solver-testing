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
folder_name = os.getcwd().split('/')[-2]
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

resolution = 1

species_definitions = {}
species_definitions['CO_g'] = {'pressure':1.} #define the gas pressures
species_definitions['O2_g'] = {'pressure':1./3.}
species_definitions['CO2_g'] = {'pressure':0}
species_definitions['s'] = {'site_names': ['111'], 'total':1} #define the sites

max_initial_guesses = 1

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

templates[
    "mkm_CO_oxidation_ads_ads"
] = r"""

rxn_expressions = [

               '*_s + CO_g -> CO*', 
               '2*_s + O2_g <-> O-O* + *_s -> 2O*',
               'CO* +  O* <-> O-CO* + * -> CO2_g + 2*',

                   ]


surface_names = ['Pt', 'Ag', 'Cu','Rh','Pd','Au','Ru','Ni'] #surfaces to include in scaling (need to have descriptors defined for each)

descriptor_names= ['O_s','CO_s'] #descriptor names

descriptor_ranges = [[${desc1}, ${desc1}],[${desc2}, ${desc2}]]

temperature = 500 #Temperature of the reaction

resolution = 1

species_definitions = {}
species_definitions['CO_g'] = {'pressure':1.} #define the gas pressures
species_definitions['O2_g'] = {'pressure':1./3.}
species_definitions['CO2_g'] = {'pressure':0}

max_initial_guesses = 1

species_definitions['s'] = {'site_names': ['111'], 'total':1} #define the sites

data_file = 'data.pkl'

input_file = 'energies.txt' #input data

gas_thermo_mode = "shomate_gas"

adsorbate_thermo_mode = "frozen_adsorbate"

scaling_constraint_dict = {
                           'O_s':['+',0,None],
                           'CO_s':[0,'+',None],
                           'O-CO_s':'initial_state',
                           'O-O_s':'final_state',
                           }

adsorbate_interaction_model = 'first_order' #use "first order" interaction model
interaction_response_function = 'smooth_piecewise_linear' #use "smooth piecewise linear" interactions
species_definitions['s']['interaction_response_parameters'] = {'cutoff':0.25,'smoothing':0.01}
max_self_interaction = 'Pd' #self interaction parameters cannot be higher than the parameter for Pd
transition_state_cross_interaction_mode = 'transition_state_scaling' #use TS scaling for TS interaction
cross_interaction_mode = 'geometric_mean' #use geometric mean for cross parameters
species_definitions['CO_s'] = {'self_interaction_parameter':[3.248, 0.965, 3.289, 3.209, 3.68, None, None, None]} #3.1
species_definitions['O_s'] = {'self_interaction_parameter':[3.405, 5.252, 6.396, 2.708, 3.87, None, None, None]} #3.1
"""


templates[
    "mkm_CO_hydrogenation"
] = r"""
scaler = 'GeneralizedLinearScaler'

rxn_expressions = [

                   'H2_g + 2*_h -> 2H_h',                                #1
                   '*_s + CO_g -> CO*',                                  #2           
                   'CO* + H_h <-> H-CO* + *_h -> CHO* + *_h',            #3
                   'CHO* + H_h <-> HCO-H* + *_h -> CHOH* + *_h',         #4             
                   'CHOH* + * <-> CH-OH* + * -> CH* + OH*',              #5

                   'CH* + H_h <-> CH-H* + *_h -> CH2* + *_h',            #6
                   'CH2* + H_h <-> CH2-H* + *_h ->  CH3* + *_h',          #7       
                   'CH3* + H_h <-> CH3-H* + *_h ->  CH4_g + * + *_h',     #8      

                   'OH* + H_h <->  H-OH* + *_h ->  H2O_g + * + *_h ',      #9

                   'CH*  +  CO* <-> CH-CO* + * -> CHCO* + *',             #10

                   'CHCO* + H_h <-> H-CHCO* + *_h -> CH2CO* + *_h ',     #11
                   'CH2CO* + H_h <-> H-CH2CO* + *_h -> CH3CO* + *_h',    #12

                   'CH3CO* + H_h  <-> H-CH3CO* + *_h -> CH3CHO* + *_h',  #13
                   'CH3CHO* -> CH3CHO_g + *',                            #14

                   'CH3CHO* + H_h  <-> CH3CHO-H* + *_h -> CH3CHOH* + *_h',  #15
                   'CH3CHOH* + H_h <-> CH3CHOH-H* + *_h -> CH3CH2OH* + *_h', #16
                   'CH3CH2OH* -> CH3CH2OH_g + *',

                   # MeOH production from JSY
                   'CHO_s + H_h <-> H-CHO_s + *_h <-> CH2O_s + *_h', #18
                   'CH2O_s + H_h <-> H-CH2O_s + *_h <-> CH3O_s + *_h', #19
                   'CH3O_s + H_h <-> CH3O-H_s + *_h <-> CH3OH_g + *_s + *_h', #20
                   'CHOH_s + H_h <-> H-CHOH_s + *_h <-> CH2OH_s + *_h', #21
                   'CH2OH_s + H_h <-> H-CH2OH_s + *_h <-> CH3OH_g + *_s + *_h', #22

                   # additional Pathway C-O bond breaking and C2+oxy formation
                   'CH3O_s + H_h + *_s <-> CH3-OH_s + *_h + *_s -> CH3_s + OH_s + *_h',
                   'CH3_s + CO_s <-> CH3-CO_s + *_s -> CH3CO_s + *_s',
                   ]


surface_names = ['Rh','Ag','Cu','Pd','Pt']

descriptor_names = ['CO_s','H-CO_s']

descriptor_ranges = [[${desc1}, ${desc1}],[${desc2}, ${desc2}]]

resolution = 1

temperature = 523
P = 20

species_definitions = {}
delta = 1e-20
species_definitions['CO_g'] = {'pressure':P*(1./3. - 4*delta)} 
species_definitions['H2_g'] = {'pressure':P*(2./3. - 4*delta)}
species_definitions['CH4_g'] = {'pressure':P*delta}
species_definitions['H2O_g'] = {'pressure':P*delta}
species_definitions['CH3CH2OH_g'] = {'pressure':P*delta}
species_definitions['CH3CHO_g'] = {'pressure':P*delta}
species_definitions['CH3OH_g'] = {'pressure':P*delta}

max_initial_guesses = 1

interaction_fitting_mode = None
numerical_delta_theta = 0.11
adsorbate_interaction_model = 'second_order'
default_self_interaction_parameter = 0

interaction_response_function = 'smooth_piecewise_linear'

interaction_response_parameters = {'slope':1,
                                   'cutoff':0.66,
                                   'smoothing':0.05,
                                   }
transition_state_cross_interaction_mode = 'transition_state_scaling' #use TS sc

species_definitions['s'] = {'site_names': ['111'], 'total':1} #upper step site
species_definitions['h'] = {'site_names': ['111'], 
                            'total':1,
                            'cross_interaction_response_parameters':{'s':{
                                                               'slope':1.0,
                                                               'cutoff':0.66/2.,
                                                               'smoothing':0.05,
                                                               }},
                            } #hydrogen reservoir

CO_cross_int = {\
                'CH2O_s':[-2.05,None,None,None,None],
                'CH3O_s':[-3.22,None,None,None,None],
                'CH2OH_s':[-2.13,None,None,None,None],
                'CH3O-H_s':[-3.88,None,None,None,None],
                'H-CH2OH_s':[-6.35,None,None,None,None],
                'H-CHO_s':[1.578,None,None,None,None],
                'H-CH2O_s':[-1.11,None,None,None,None],
                'H-CHOH_s':[-0.46,None,None,None,None],
                'CH3CO_s': [-1.6945625,None,None,None,None],
                'CHO_s': [-1.3309375,None,None,None,None],
                'CH3CHOH-H_s': [-4.5309375,None,None,None,None],
                'CH3CHO-H_s': [-5.91275,None,None,None,None],
                'CH3CHOH_s': [-4.749125,None,None,None,None],
                'CHOH_s': [-3.949125,None,None,None,None],
                'CH-OH_s': [-1.1854375,None,None,None,None],
                'H-CH3CO_s': [-0.31275,None,None,None,None],
                'CH2-H_s': [1.8690625,None,None,None,None],
                'HCO-H_s': [-0.676375,None,None,None,None],
                'H_h': [0.87275,None,None,None,None],
                'CH3CH2OH_s': [-5.6945625,None,None,None,None],
                'CH3-H_s': [-0.0945625,None,None,None,None],
                'H-OH_s': [1.0690625,None,None,None,None],
                'CH2_s': [-3.149125,None,None,None,None],
                'CH_s': [-2.64,None,None,None,None],
                'CH3CHO_s': [-5.2581875,None,None,None,None],
                'CH3_s': [1.0690625,None,None,None,None],
                'CH-CO_s': [2.23275,None,None,None,None],
                'CHCO_s': [0.48725,None,None,None,None],
                'OH_s': [-2.7855,None,None,None,None],
                'CH2CO_s': [-1.3309375,None,None,None,None],
                'H-CH2CO_s': [0.63275,None,None,None,None],
                'H-CHCO_s': [2.08725,None,None,None,None],
                'CH-H_s': [0.56,None,None,None,None],
                'H-CO_s': [0.2690625,None,None,None,None],
                'CH3-OH_s':[-1.26,None,None,None,None],
                'CH3-CO_s':[2.09,None,None,None,None],
        }


species_definitions['CO_s'] = { 'self_interaction_parameter':[5.6,None,None,None,None],
                                'cross_interaction_parameters':CO_cross_int}
species_definitions['CHOH_s'] = { 'self_interaction_parameter':[8.6,None,None,None,None]}
species_definitions['CH_s'] = {'self_interaction_parameter':[2.8,None,None,None,None]}
species_definitions['OH_s'] = {'self_interaction_parameter':[1.0,None,None,None,None]}
species_definitions['CHO_s'] = {'self_interaction_parameter':[8.4,None,None,None,None]}

scaling_constraint_dict = {
                          'C_s':['+',0,None],
                          'O_s':[None,None,None],
                          'CO_s':['+','0',None],
                          'H_h':['+','0',None],
                          'CHO_s':['+','+',None],
                          'CHOH_s':['+','0',None],
                          'OH_s':[None,None,None],
                          'CH_s':['+',0,None],
                          'CH2_s':['+',0,None],
                          'CH3_s':['+',0,None],
                          'CHCO_s':['+',0,None],
                          'CH2CO_s':['+',0,None],
                          'CH3CO_s':['+',0,None],
                          'CH3CHO_s':['+','+',None],
                          'CH3CHOH_s':['+',0,None],
                          'CH3CH2OH_s':[None,None,None],
                          'CH2OH_s':['+',0,None],
                          'CH3O_s': [None,None,None],
                          'CH2O_s': ['+','+',None],
                          'HCO-H_s': 'final_state',
                          'CH-OH_s': 'final_state',
                          'CH-H_s': 'BEP',
                          'CH2-H_s':'BEP',
                          'CH3-H_s': 'BEP',
                          'H-OH_s': 'BEP',
                          'CH-CO_s': 'BEP',
                          'H-CHCO_s': 'initial_state',
                          'H-CH2CO_s': 'final_state',
                          'H-CH3CO_s': 'BEP',
                          'CH3CHO-H_s':'BEP',
                          'CH3CHOH-H_s': 'BEP',
                          'H-CO_s': 'BEP',
                          'H-CHOH_s':'initial_state',
                          'H-CH2O_s':'BEP',
                          'H-CH2OH_s':'BEP',
                          'CH3O-H_s':'BEP',
                          'H-CHO_s':'BEP',
                          
                        'CH3-OH_s':'BEP',
                        'CH3-CO_s':'BEP',
        }

interaction_scaling_constraint_dict = {
                           'C_s':['+',0,None],
                           'O_s':[0,'+',None],
                           'CO_s':['+',0,None],
                           'OH_s':[0,'+',None],
                           }

data_file = 'data.pkl'
input_file = 'energies.txt' #input data

gas_thermo_mode = "shomate_gas"
adsorbate_thermo_mode = "harmonic_adsorbate"
max_entropy_per_mode = 3*8.617e-5

"""


templates[
    "mkm_methanol_synthesis"
] = r"""

rxn_expressions = [

    # methanation
    'H2_g + 2*_h <-> H-H_h + *_h -> 2H_h',  # 0
    '*_s + CO_g -> CO_s',  # 1
    'CO_s + *_f + H_h <-> C-OH_s + *_f + *_h -> C_f + OH_s + *_h',  # 2
    'C_f + H_h <-> H-C_f + *_h -> CH_f + *_h',  # 3
    'CH_f + H_h + *_t <-> H-CH_f + *_t + *_h -> CH2_t + *_h + *_f',  # 4
    'CH2_t + H_h <-> H-CH2_t + *_h -> CH3_t + *_h',  # 5
    'CH3_t + H_h <-> H-CH3_t + *_h -> CH4_g + *_t + *_h',  # 6
    'O_s + H_h <-> O-H_s + *_h -> OH_s + *_h ',  # 7
    'OH_s + H_h <-> H-OH_s + *_h -> H2O_g + *_s+ *_h ',  # 8
    '2OH_s <-> O-H2O_s + *_s-> H2O_g + O_s + *_s',  # 211only                #9

    # MeOH
    'CO_s + H_h <-> H-CO_s + *_h -> HCO_s + *_h',  # 10
    'HCO_s + H_h <-> H-HCO_s + *_h -> CH2O_s + *_h',  # 11
    'CH2O_s + H_h <-> H-CH2O_s + *_h -> CH3O_s + *_h',  # 12
    'CH3O_s + H_h <-> CH3O-H_s + *_h -> CH3OH_g + *_h + *_s',  # 13

    'HCO_s + *_f <-> O-CH_s + *_f -> CH_f + O_s',  # 211only          #14
    'CH2O_s + *_s<-> O-CH2_s + *_s-> CH2_s + O_s',  # 211only                #15
    'CH3O_s + *_s<-> O-CH3_s + *_s-> CH3_s + O_s',  # 211only                #16
    'CH3_s + *_t -> CH3_t + *_s',  # 17
    'CH2_s + *_t -> CH2_t + *_s',  # 18

    # C-CO coupling
    'C_f + CO_s  <-> C-CO_f + *_s-> CCO_f + *_s',        # 19
    'CH_f  + CO_s <-> CH-CO_f + *_s-> CHCO_f + *_s',      # 20
    'CH2_s + CO_s  <-> CH2-CO_s + *_s-> CH2CO_s + *_s',  # 21
    'CH3_s + CO_s  <-> CH3-CO_s + *_s-> CH3CO_s + *_s',  # 22

    # CCO hydrogenation
    'CCO_f + H_h <-> H-CCO_f + *_h -> CHCO_f + *_h',
    'CCO_f + H_h <-> H-CCO_f + *_h -> CCHO_f + *_h',
    'CHCO_f + H_h + *_s<-> H-CHCO_f + *_h + *_s-> CH2CO_s + *_h + *_f',
    'CHCO_f + H_h + *_s<-> CHCO-H_f + *_h + *_s-> CHCOH_s + *_h + *_f',
    'CHCO_f + H_h + *_s<-> CHCO-H_f + *_h + *_s-> CHCHO_s + *_h + *_f',

    'CH2CO_s + H_h <-> CH2CO-H_s + *_h -> CH2COH_s + *_h',
    'CHCOH_s + H_h <-> H-CHCOH_s + *_h -> CH2COH_s + *_h',
    'CHCOH_s + H_h <-> CHCOH-H_s + *_h -> CHCHOH_s + *_h',
    'CH2COH_s + H_h <-> CH2COH-H_s + *_h -> CH2CHOH_s + *_h',
    'CHCHOH_s + H_h <-> H-CHCHOH_s + *_h -> CH2CHOH_s + *_h',
    'CHCHOH_s + H_h <-> CHCHOH-H_s + *_h -> CHCH2OH_s + *_h',
    'CH2CHOH_s + H_h <-> H-CH2CHOH_s + *_h -> CH3CHOH_s + *_h',
    'CH2CHOH_s + H_h <-> CH2CHOH-H_s + *_h -> CH2CH2OH_s + *_h',
    'CHCH2OH_s + H_h <-> H-CHCH2OH_s + *_h -> CH2CH2OH_s + *_h',
    'CH2CH2OH_s + H_h <-> H-CH2CH2OH_s + *_h -> CH3CH2OH_g + *_s+ *_h',
    'CH3CHOH_s + H_h <-> H-CH3CHOH_s + *_h -> CH3CH2OH_g + *_s+ *_h',

    'CCO_f + H_h <-> CCO-H_f + *_h -> CCOH_f + *_h',
    'CCHO_f + H_h + *_s<-> H-CCHO_f + *_h + *_s -> CHCHO_s + *_h + *_f',
    'CCHO_f + H_h + *_s<-> CCHO-H_f + *_h + *_s -> CCHOH_s + *_h + *_f',
    'CCOH_f + H_h + *_s<-> H-CCOH_f + *_h + *_s -> CCHOH_s + *_h + *_f',
    'CCOH_f + H_h + *_s<-> CCOH-H_f + *_h + *_s -> CHCOH_s + *_h + *_f',
    'CCHOH_s + H_h <-> H-CCHOH_s + *_h -> CHCHOH_s + *_h',
    'CHCHO_s + H_h <-> CHCHO-H_s + *_h -> CHCHOH_s + *_h',
    'CH2COH_s + H_h <-> CH2COH-H_s + *_h -> CH2CHOH_s + *_h',

    'CH3CO_s + H_h  <-> CH3CO-H_s + *_h -> CH3COH_s + *_h',
    'CH3CO_s + H_h  <-> H-CH3CO_s + *_h -> CH3CHO_s + *_h',
    'CH3CHO_s + H_h  <-> CH3CHO-H_s + *_h -> CH3CHOH_s + *_h',
    'CH3COH_s + H_h  <-> H-CH3COH_s + *_h -> CH3CHOH_s + *_h',
]


transition_state_scaling_parameters = {'C2': [0.856923076923077, 1.0213846153846156],
                                       'C1': [0.856923076923077, 1.0213846153846156],
                                       'HA': [0.856923076923077, 1.0213846153846156], 
                                       'O': [0.856923076923077, 1.0213846153846156],
                                       'HC': [0.856923076923077, 1.0213846153846156],
                                       'HCHF': [0.9040078556953679, 1.0850540423932544],
                                       'HCF': [0.92666528478110022, 1.1472740172658014],
                                       'HHH': [0.93434150378549374, 0.8706685500086766],
                                       'CH3OHS': [0.79769953720375431, 1.3288816719091958],
                                       'CH2COS': [0.66153883256667578, 2.2480723262738671],
                                       'HCH2OS': [0.65661108228923626, 1.1682969482947245],
                                       'CCOF': [0.76566126306385485, 1.8603218796633014],
                                       'CH3COS': [0.94412985940579097, 0.95213477042121075],
                                       'HCH2T': [1.0117818769374605, 0.62184801171077664],
                                       'HOHS':[0.75553852148731815, 0.87353763486229896],
                                       'OCH3S': [0.63348882023953523, 2.0188057837259343],
                                       'OCHS':[0.81393143870201701, 2.3124332180026932],
                                       'CHCOF':[0.86418687117997095, 1.4098531921772275],
                                       'OHS':[0.86707852180806666, 0.84084205423986202],
                                       'HCH3T':[0.95662009478240573, 0.73742912372397385],
                                       'HHCOS':[0.80681391439157879, 0.98031578916602335],
                                       'HCOS':[0.77006109667770817, 1.1203407392851128],
                                       'COHS':[1.0693928920752096, 0.63985949945071552],
                                       'OH2OS':[0.88672406942095738, 1.0944069371952425],
                                       'OCH2S': [0.78928878444791162, 2.3594008814522867],
                                       }

surface_names = ['Ni', 'Cu', 'Pd', 'Pt', 'Ag', 'Rh', 'Ru', 'CuCo']

descriptor_names = ['CO_s', 'OH_s']

descriptor_ranges = [[${desc1}, ${desc1}],[${desc2}, ${desc2}]]

resolution = 1

scaling_constraint_dict = {

    # Methanation intermediates
    'C_f': ['+', 0, None],
    'CO_s': ['+', 0, None],
    'H_h': ['+', 0, None],
    'CH_f': ['+', 0, None],
    'CH2_t': ['+', 0, None],
    'CH2_s': ['+', 0, None],
    'CH3_t': ['+', 0, None],
    'CH3_s': ['+', 0, None],
    'OH_s': [0, '+', None],
    'O_s': [0, '+', None],

    # Methanol intermediates (from J. Cat. MeOH paper)
    'HCO_s': ['+', 0, None],
    'CH2O_s': ['+', '+', None],
    'CH3O_s': [0, '+', None],

    # CC intermediates
    'CCO_f': ['+', 0, None],
    'CHCO_f': ['+', '+', None],
    'CCHO_f': ['+', '+', None],
    'CCOH_f': ['+', 0, None],
    'CH2CO_s': ['+', '+', None],
    'CHCHO_s': ['+', '+', None],
    'CHCOH_s': ['+', 0, None],
    'CCHOH_s': ['+', 0, None],
    'CH2COH_s': ['+', 0, None],
    'CHCHOH_s': ['+', 0, None],
    'CH2CHOH_s': ['+', 0, None],
    'CHCH2OH_s': ['+', 0, None],
    'CH2CH2OH_s': ['+', 0, None],
    'CH3CHOH_s': ['+', 0, None],
    'CH3CO_s': ['+', '+', None],
    'CH3CHO_s': ['+', 0, None],
    'CH3COH_s': ['+', 0, None],

    # formic acid and CO pathway
    'HCOO_s': [0, '+', None],
    'COOH_s': ['+', 0, None],
    'HCOOH_s': [0, 0, 2.01],  # Studt approximation, Nat. Chem.
    'H-COO_s': 'final_state',
    'COO-H_s': 'initial_state',
    'HCOO-H_s': 'initial_state',
    'H-COOH_s': 'final_state',
    'OC-OH_s': 'final_state',

    # acetaldehyde pathway
    'H2COOH_s': [0, '+', None],
    'H-HCOOH_s': 'final_state',
    'H2CO-OH_s': 'final_state',

    'C-CO_f': 'TS(C_f+CO):CCOF',
    'H-CH3COH_s': 'TS(CH3COH+H_h):HC',
    'CH2CHOH-H_s': 'TS(CH2CHOH+H_h):HA',
    'H-CHCH2OH_s': 'TS(CHCH2OH+H_h):HC',
    'H-CH3CHOH_s': 'TS(CH3CHOH+H_h):HC',
    'H-CCHO_f': 'TS(CCHO_f+H_h):HC',
    'CHCHOH-H_s': 'TS(CHCHOH+H_h):HA',
    'CH2-CO_s': 'TS(CH2+CO):CH2COS',
    'H-CH2CHOH_s': 'TS(CH2CHOH+H_h):HC',
    'CH3CHO-H_s': 'TS(CH3CHO+H_h):HA',
    'H-CHCHO_s': 'TS(CHCHO+H_h):HC',
    'CH3-CO_s': 'TS(CH3+CO):CH3COS',
    'H-CH2_t': 'TS(CH2_t+H_h):HCH2T',
    'CH3O-H_s': 'TS(CH3O+H_h):CH3OHS',
    'H-CCHOH_s': 'TS(CCHOH+H_h):HC',
    'H-CHCHOH_s': 'TS(CHCHOH+H_h):HC',
    'H-CCOH_f': 'TS(CCOH_f+H_h):HC',
    'H-CH2O_s': 'TS(CH2O+H_h):HCH2OS',
    'H-OH_s': 'TS(OH+H_h):HOHS',
    'O-CH3_s': 'TS(O+CH3):OCH3S',
    'H-CH_f': 'TS(CH_f+H_h):HCHF',
    'H-C_f': 'TS(C_f+H_h):HCF',
    'O-CH_s': 'TS(O+CH_f):OCHS',
    'CH-CO_f': 'TS(CH_f+CO):CHCOF',
    'CHCOH-H_s': 'TS(CHCOH+H_h):HA',
    'O-H_s': 'TS(O+H_h):OHS',
    'H-CCO_f': 'TS(CCO_f+H_h):HA',
    'H-CH3_t': 'TS(CH3_t+H_h):HCH3T',
    'CCOH-H_f': 'TS(CCOH_f+H_h):HA',
    'H-CHCOH_s': 'TS(CHCOH+H_h):HC',
    'H-H_h': 'TS(H_h+H_h):HHH',
    'CH2COH-H_s': 'TS(CH2COH+H_h):HA',
    'H-HCO_s': 'TS(HCO+H_h):HHCOS',
    'CCO-H_f': 'TS(CCO_f+H_h):HC',
    'H-CO_s': 'TS(CO+H_h):HCOS',
    'CH3CO-H_s': 'TS(CH3COH+H_h):HA',
    'CHCO-H_f': 'TS(CHCO_f+H_h):HC',
    'C-OH_s': 'TS(OH+C_f):COHS',
    'H-CHCO_f': 'TS(CHCO_f+H_h):HC',
    'H-CH2CH2OH_s': 'TS(CH2CH2OH+H_h):HC',
    'O-H2O_s': 'TS(OH+OH):OH2OS',
    'CHCHO-H_s': 'TS(CHCHO+H_h):HA',
    'CCHO-H_f': 'TS(CCHO_f+H_h):HA',
    'O-CH2_s': 'TS(O+CH2):OCH2S',
    'H-CH3CO_s': 'TS(CH3CO+H_h):HC',
    'CH2CO-H_s': 'TS(CH2CO+H_h):HA',
}

data_file = 'data.pkl'

total_C = 6.7
species_definitions = {}
species_definitions['CO_g'] = {'pressure': total_C}
species_definitions['H2_g'] = {'pressure': 13.3}
species_definitions['CH4_g'] = {'pressure': 1e-6}
species_definitions['CH3OH_g'] = {'pressure': 1e-9}
species_definitions['CH3CH2OH_g'] = {'pressure': 1e-9}
species_definitions['H2O_g'] = {'pressure': 1e-2}

species_definitions['s'] = {'site_names': ['211'], 'total': 1}
species_definitions['t'] = {'site_names': ['111'], 'total': 1}
species_definitions['f'] = {'site_names': ['211'], 'total': 1}
species_definitions['h'] = {'site_names': ['111'], 'total': 1}

max_initial_guesses = 1

temperature = 523.15
adsorbate_thermo_mode = 'harmonic_adsorbate'
gas_thermo_mode = 'shomate_gas'

input_file = 'energies.txt'

"""
