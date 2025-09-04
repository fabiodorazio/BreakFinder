#!/usr/bin/env python3

import pandas as pd
import sys
import glob
from pathlib import Path
import utils

def main():
    """
    Reads columns (chr | start   | end) of Sample beds and reference AsiSI
    Iterates over sample bed rows (small dfs) and checks for intersections with reference
    Updates intersections count
    Retains indexes of overlapping regions -> slices df for overlapping and non-overlapping regions
    Calculates intersection rate as n intersections / n rows reference 
    """
    # ensure all arguments are passed
    if len(sys.argv) != 5:
        print("Usage: python Count_intersections.py <merged_bed> <reference_bed> <output_prefix>")
        sys.exit(1)
    
    root_path = Path(__file__).parent.parent
    # load arguments
    merged_bed = sys.argv[1]
    reference_bed = sys.argv[2]
    sample_prefix = sys.argv[3]
    binsize = int(sys.argv[4])

    # extract original read counts equivalent to sample
    s_basename = merged_bed.split('/')[-1].split('.')[0]
    non_intersected_path = f'{root_path}/results/intersections/{s_basename}.*nonintersected*.bed'
    non_intersected = glob.glob(non_intersected_path)

    try:
        # Read files
        merged_df = pd.read_csv(merged_bed, sep='\t', header=None, 
                               usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'counts'])
        ref_df = pd.read_csv(reference_bed, sep='\t', header=None,
                           usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
        non_intersected_df = pd.read_csv(non_intersected[0], sep='\t', header=None,
                           usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'counts'])

        # add flags
        merged_df['Flag'] = 'Asisi Overlap'
        non_intersected_df['Flag'] = 'Asisi No Overlap'
        # recreate original file for chr size estimation
        original_df = pd.concat([merged_df, non_intersected_df])

        # 1) Count intersections
        intersections = len(merged_df)
        print(merged_bed)
        
        # 2) estimate chr size
        #genome_size = 3_200_000_000
        # in this case chromosome size is approximated
        if len(non_intersected_df['chrom'].unique()) >= 23:
            # full genome size can be expected
            region_size = 3_200_000_000
        else:
            region_size = 0
            # expect less chromosomes
            for chr in non_intersected_df['chrom'].unique():
                chr_subset = non_intersected_df[non_intersected_df['chrom'] == chr]
                chr_size = max(chr_subset['end']) - min(chr_subset['start'])

                region_size += chr_size


        # Other metrics:
        if merged_df.shape[0] > 0:
            # 3) alculate Asisi intersection rate
            Asisi_intersection_rate = intersections / len(ref_df)
            # 4) calculate Over total intersection rate
            Tot_intersection_rate = intersections / len(original_df)
        else:
            Asisi_intersection_rate = 0
            Tot_intersection_rate = 0
        # 5) all breaks rate
        all_break_rate = len(original_df) / region_size
        
        # create output data frame
        summary_df = pd.DataFrame.from_dict({
            'Metrics':['Reference_sites', 'Asisi_Intersections_N', 'Asis_Intersection_rate',
                       'Tot_intersection_rate', 'All_break_rate'],
            'Values': [len(ref_df), intersections, Asisi_intersection_rate, 
                       Tot_intersection_rate, all_break_rate]
        })

        summary_df.to_csv(f"{sample_prefix}.summary.txt", sep='\t', index=False)
        
        # save
        with open(f"{sample_prefix}.count.txt", 'w') as f:
            f.write(str(intersections))
        
        # Try to plot
        try:
            import matplotlib.pyplot
            import seaborn
            import visuals
            visuals.breaks_location(sample_prefix, original_df, ref_df, root_path, bin_size=binsize)
        except ImportError:
            print("Matplotlib not available - skipping plots. Try running nextflow with `-profile conda`")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()