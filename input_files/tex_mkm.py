from catmap import ReactionModel
import os
import json

mkm_file = 'input.mkm'
model = ReactionModel(setup_file=mkm_file)
model.model_summary()
