import os
import pandas as pd

if not os.path.exists('../results/cactus/'):
    os.makedirs('../results/cactus/')

f = open( '../results/cactus/' + snakemake.config["alignment_name"] + ".txt", "a")

f.writelines('"' + snakemake.config["speciesTree"] + '"\n\n')

for x in snakemake.params["genomes"]:
  f.writelines(x + ' data/' + x + '.fa.gz\n')

f.close()

if not os.path.exists("../results/checkpoints/"):
    os.makedirs("../results/checkpoints/")

g = open( "../results/checkpoints/jobstore_setup.txt", "a")
g.writelines('setup complete')
g.close()