import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import utils
import numpy as np


def breaks_location(sample_name, breaks_all, asisi_df, output, bin_size = 5000):
    '''
    The aim is to plot read counts in bins over each chromosome in bed files
    Asisi sites are mapped to the input data and highlighted
    '''

    # use only the start of the coordinate (1bp reads) 
    breaks_all['location'] = breaks_all['start'].astype(int)

    # extract the asisi sites from input file
    asisi_sites = breaks_all[breaks_all['Flag'] == 'Asisi Overlap']
    # subset and rename col start for clarity
    asi = asisi_sites[['chrom','start']].rename(columns={'start':'pos'})
    # standardize types
    asi['chrom'] = asi['chrom'].astype(str)
    asi['pos'] = asi['pos'].astype('int64')
    # bin positions
    asi['bin_start'] = (asi['pos'] // bin_size) * bin_size
    # create asisi df
    asisi_bins = (asi[['chrom','bin_start']]
                .drop_duplicates()
                .assign(has_asisi=True))

    asisi_positions = set(asisi_df['start']).union(set(asisi_df['end']))

    # Create bins for breaks
    breaks_all['bin_start'] = (breaks_all['location'] // bin_size) * bin_size
    # ensure no overlaps with bins
    breaks_all['bin_end'] = breaks_all['bin_start'] + bin_size - 1
    # initialise column to store booleans
    breaks_all['has_asisi'] = False
    

    # Group by chromosome AND bin, aggregate counts and merge with asisi bins
    binned_df = breaks_all.groupby(['chrom', 'bin_start']).agg({
        'counts': 'sum',
    }).reset_index().merge(asisi_bins, on=['chrom','bin_start'], how='left').assign(has_asisi=lambda d: d['has_asisi'].fillna(False),
                     Flag=lambda d: np.where(d['has_asisi'],
                                             'Asisi Overlap', 'Asisi No Overlap'))


    # run through chromosomes
    for chr, chr_data in binned_df.groupby('chrom'):
        palette = {'Asisi Overlap': 'red', 'Asisi No Overlap': 'blue'}
        colors = chr_data['Flag'].map(palette).to_list()
        
        # Create the plot
        ax = sns.barplot(data=chr_data, x='bin_start', y='counts', 
                        hue='Flag',
                        palette={'Asisi Overlap': 'red', 'Asisi No Overlap': 'blue'},
                        dodge=False,
                        alpha=0.7,
                        linewidth=0.5)
        
        # avoid x axis crowding
        for i, (patch, flag) in enumerate(zip(ax.patches, chr_data['Flag'])):
            if flag == 'Asisi Overlap':
                # Thicker line for overlaps
                patch.set_linewidth(3)
                patch.set_edgecolor('darkred')
        plt.xlabel(f'Genomic Bin Start Position (chr {chr})')
        plt.ylabel('Break Counts')
        plt.title(f'DNA Breaks with Asisi Overlap - Chromosome {chr}')
        plt.legend(title='Asisi Overlap')
        plt.xticks([])
        plt.tight_layout()
        output_path = utils.check_output_dir(f"{output}/results/Plots/")
        plt.savefig(f"{output_path}{sample_name}_{chr}_Breaks_location.png")
        plt.close()


def intersection_rate_per_sample(df, output):
    '''
    Simple bar chart of intersection rates by sample
    '''
    # subset for metrics
    df = df[df['Metrics'] == 'Asis_Intersection_rate']

    # sort by sample number
    df['sample_num'] = df['Sample'].str.split('Sample').str[-1].astype(int)
    df_sorted = df.sort_values('sample_num')

    plt.figure(figsize=(10, 6))
    plt.bar(df_sorted['sample_num'], df_sorted['Values'])
    plt.xlabel('Samples')
    plt.ylabel('AsiSI Intersection Rate')
    plt.title('AsiSI Intersection Rate Across Samples')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('intersection_rate_comparison.png', dpi=300)
    output_path = utils.check_output_dir(f"{output}/results/Plots/")
    plt.savefig(f"{output_path}Intersection_Rate.png")