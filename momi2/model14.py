import momi
import numpy as np
import matplotlib
import logging
from autograd.numpy import log

'''
tree structure with bidirectional southeast/central and central/east migration, 
and unidirectional west into central and east into southeast migration events 
with southeast population size change
'''

logging.basicConfig(level=logging.INFO,
                    filename="model14.log")

model14 = momi.DemographicModel(N_e=1e5,gen_time=3,muts_per_gen=2.43e-9)

# read data
sfs = momi.Sfs.load("outfile_sfs")
model14.set_data(sfs,length=411642)

# population size parameters
model14.add_size_param("n_east")
model14.add_size_param("n_west")
model14.add_size_param("n_central")
model14.add_size_param("n_southeast")
model14.add_size_param("n_southeast_bn")

# split time parameters
# constrain split times so that central-west split is most recent
model14.add_time_param("tdiv_central_west")
model14.add_time_param("tdiv_east_west",lower_constraints=["tdiv_central_west"])
# constrain split times so that east-southeast is earliest split in the tree
model14.add_time_param("tdiv_west_southeast",lower_constraints=["tdiv_central_west","tdiv_east_west"])

# migration parameters
model14.add_time_param("tm_southeast_central", upper_constraints=["tdiv_central_west"])
model14.add_time_param("tm_central_southeast", upper_constraints=["tdiv_central_west"])
model14.add_time_param("tm_east_southeast", upper_constraints=["tdiv_central_west"])
model14.add_time_param("tm_east_central", upper_constraints=["tdiv_central_west"])
model14.add_time_param("tm_central_east", upper_constraints=["tdiv_central_west"])
model14.add_time_param("tm_west_central", upper_constraints=["tdiv_central_west"])
model14.add_time_param("tb_southeast", upper_constraints=["tdiv_west_southeast"])

model14.add_pulse_param("p_southeast_central")
model14.add_pulse_param("p_central_southeast")
model14.add_pulse_param("p_east_southeast")
model14.add_pulse_param("p_east_central")
model14.add_pulse_param("p_central_east")
model14.add_pulse_param("p_west_central")

# add tips
model14.add_leaf("Southeast", N="n_southeast", t=0)
model14.set_size("Southeast", N="n_southeast_bn", t="tb_southeast")
model14.add_leaf("Central", N="n_central", t=0)
model14.add_leaf("East", N="n_east", t=0)
model14.add_leaf("West", N="n_west", t=0)

# migration events
model14.move_lineages("Southeast", "Central", t="tm_southeast_central", p="p_southeast_central")
model14.move_lineages("Central", "Southeast", t="tm_central_southeast", p="p_central_southeast")
model14.move_lineages("East", "Southeast", t="tm_east_southeast", p="p_east_southeast")
model14.move_lineages("East", "Central", t="tm_east_central", p="p_east_central")
model14.move_lineages("Central", "East", t="tm_central_east", p="p_central_east")
model14.move_lineages("West", "Central", t="tm_west_central", p="p_west_central")

# split events
# Central joins onto West
model14.move_lineages("Central", "West", t="tdiv_central_west")
# East joins onto West
model14.move_lineages("East", "West",  t="tdiv_east_west")
# West joins onto Southeast
model14.move_lineages("West", "Southeast", t="tdiv_west_southeast")

# figure
yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5, 2e6]

fig = momi.DemographyPlot(
    model14, ["West", "Central", "East", "Southeast"],
    figsize=(6,8),
    major_yticks=yticks,
    linthreshy=1e5, pulse_color_bounds=(0,1.0), draw=False)

# # plot bootstraps onto the canvas in transparency
# for params in bootstrap_results:
#     fig.add_bootstrap(
#         params,
#         # alpha=0: totally transparent. alpha=1: totally opaque
#         alpha=1/10)

# now draw the inferred demography on top of the bootstraps
fig.draw()
fig.draw_N_legend(loc="upper left")

from matplotlib import pyplot as plt
plt.savefig("model14_fig1.png")

#set up for n runs
results = []
n_runs = 5

#do the  runs
for i in range(n_runs):
    print(f"Starting run {i+1} out of {n_runs}...")
    model14.set_params(randomize=True)
    results.append(model14.optimize(method="L-BFGS-B", options={"maxiter":100}))
# sort results according to log likelihood, pick the best (largest) one
best_result = sorted(results, key=lambda r: r.log_likelihood)[-1]
model14.set_params(best_result.parameters)
best_result
print("Best result: {}".format(best_result))

#calculate the AIC

for model in [model14]:
    lik = model.log_likelihood()
    nparams = len(model.get_params())
    aic = 2*nparams - 2*lik
    print("AIC {}".format(aic))
