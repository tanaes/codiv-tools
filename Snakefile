import pandas as pd
from os.path import join

configfile: "config.yaml"

clades_fp = config['clades']

clades_df = pd.read_csv(clades_fp, index_col=0)

clades_df.index = clades_df.index.astype(str)

print(clades_df)

include: "resources/snakefiles/phylo.smk"
include: "resources/snakefiles/align.smk"

rule all:
    input:
        expand("output/phylo/raxml/raxml_{clade}/RAxML_bestTree.{clade}_out",
               clade=clades_df.index)
