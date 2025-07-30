import yaml
import snakemake.io

with open('config.yaml', 'r') as file:
    config = yaml.safe_load(file)

from parse_config import *

cfg = PrepConfig(config)

bett,fun=cfg.get_possible_output('True')



print(bett['method'])

wildcards={'dataset':'Biase'}
adat = lambda wildcards: cfg.get_from_dataset(wildcards.dataset,'file')




params=cfg.get_clustering_params('leiden', 'res')
print(params)

b= 'False'

print("True"==str(b))