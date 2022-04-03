from catmap import ReactionModel

mkm_file = 'EtOH.mkm'
model = ReactionModel(setup_file=mkm_file)
model.run()

# from catmap import analyze
# vm = analyze.VectorMap(model)
# vm.plot_variable = 'rate' #tell the model which output to plot
# vm.log_scale = True #rates should be plotted on a log-scale
# vm.min = 1e-25 #minimum rate to plot
# vm.max = 1e3 #maximum rate to plot
# vm.plot(save='rate.pdf') #draw the plot and save it as "rate.pdf"

# vm.unique_only = False
# vm.plot(save='all_rates.pdf')
# vm.unique_only = True

# model.output_variables += ['production_rate']
# model.run()
# vm.production_rate_map = model.production_rate_map #attach map
# vm.threshold = 1e-30 #do not plot rates below this
# vm.plot_variable = 'production_rate'
# vm.plot(save='production_rate.pdf')

# vm.plot_variable = 'coverage'
# vm.log_scale = False
# vm.min = 0
# vm.max = 1
# vm.plot(save='coverage.pdf')

# vm.include_labels = ['CO_s']
# vm.plot(save='CO_coverage.pdf')

# sa = analyze.ScalingAnalysis(model)
# sa.plot(save='scaling.pdf')
