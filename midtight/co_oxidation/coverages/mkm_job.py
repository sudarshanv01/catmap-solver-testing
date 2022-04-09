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

# Use the same resolution
model.resolution = solver_specifics['resolution']

# Store the production rate
model.output_variables += ['production_rate']

model.run()