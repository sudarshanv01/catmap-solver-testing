from catmap import ReactionModel
from catmap import analyze
import numpy as np

mkm_file = 'CH3OH.mkm'
model = ReactionModel(setup_file=mkm_file)

model.run()