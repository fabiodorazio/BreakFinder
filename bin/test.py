import pandas as pd
import numpy as np

g = pd.read_csv('/Users/fabiodorazio/Github/BreakFinder/results/merged_breaks/Sample14.breakends_merged.bed',
                sep='\t', header=None, 
                               usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'counts'])
print(g)
coord_col="start"
count_col="counts"
step=1
g = g[[coord_col, count_col]].groupby(coord_col, as_index=True).sum()
start = int(g.index.min())
stop  = int(g.index.max())
# full coordinate index
full = np.arange(start, stop + 1, step)
out = g.reindex(full, fill_value=0).rename_axis(coord_col).reset_index()

print(out.sort_values('counts'))
print(out[out['counts'] >= 2])