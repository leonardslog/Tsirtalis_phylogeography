import momi
import numpy as np
import matplotlib
import logging
from autograd.numpy import log

# null model, tree structure with no migration

logging.basicConfig(level=logging.INFO,
                    filename="model0.log")

model0 = momi.DemographicModel(N_e=1e5,gen_time=3,muts_per_gen=2.43e-9)

# read data
sfs = momi.Sfs.load("outfile_sfs")
model0.set_data(sfs,length=411642)

# population size parameters
model0.add_size_param("n_east")
model0.add_size_param("n_west")
model0.add_size_param("n_central")
model0.add_size_param("n_southeast")

# split time parameters
# constrain split times so that central-west split is most recent
model0.add_time_param("tdiv_central_west")
model0.add_time_param("tdiv_east_west",lower_constraints=["tdiv_central_west"])
# constrain split times so that east-southeast is earliest split in the tree
model0.add_time_param("tdiv_west_southeast",lower_constraints=["tdiv_central_west","tdiv_east_west"])


# add tips
model0.add_leaf("Southeast", N="n_southeast", t=0)
model0.add_leaf("Central", N="n_central", t=0)
model0.add_leaf("East", N="n_east", t=0)
model0.add_leaf("West", N="n_west", t=0)

# split events
# Central joins onto West
model0.move_lineages("Central", "West", t="tdiv_central_west")
# East joins onto West
model0.move_lineages("East", "West",  t="tdiv_east_west")
# West joins onto Southeast
model0.move_lineages("West", "Southeast", t="tdiv_west_southeast")

# figure
yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5, 2e6]

fig = momi.DemographyPlot(
    model0, ["West", "Central", "East", "Southeast"],
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
plt.savefig("model0_fig1.png")

#set up for n runs
results = []
n_runs = 5

#do the  runs
for i in range(n_runs):
    print(f"Starting run {i+1} out of {n_runs}...")
    model0.set_params(randomize=True)
    results.append(model0.optimize(method="L-BFGS-B", options={"maxiter":100}))
# sort results according to log likelihood, pick the best (largest) one
best_result = sorted(results, key=lambda r: r.log_likelihood)[-1]
model0.set_params(best_result.parameters)
best_result
print("Best result: {}".format(best_result))

#calculate the AIC

for model in [model0]:
    lik = model.log_likelihood()
    nparams = len(model.get_params())
    aic = 2*nparams - 2*lik
    print("AIC {}".format(aic))
