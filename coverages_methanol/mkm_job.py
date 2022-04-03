from catmap import ReactionModel
import numpy as np

mkm_file = 'EtOH.mkm'
model = ReactionModel(setup_file=mkm_file)

# overbinding
overbinding = 0.25
CO_energies = model.species_definitions['CO_s']['formation_energy']
CO_energies = [E+overbinding for E in CO_energies if E != None]
model.species_definitions['CO_s']['formation_energy'] = CO_energies

model.run()
