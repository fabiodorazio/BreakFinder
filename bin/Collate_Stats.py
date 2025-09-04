import pandas as pd
import sys
from pathlib import Path

def main():
    '''
    Loops through file list, adds Sample name column and bind all files into 1 table
    Add visuals
    '''
    root_path = Path(__file__).parent.parent

    if len(sys.argv) < 3:
        print("Usage: python concatenate_stats.py <output_file> <input_file1> [<input_file2> ...]")
        sys.exit(1)
    
    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    
    if not input_files:
        raise ValueError("No statistics files found")

    # read and combine all files
    dfs = []
    for file in input_files:
        try:
            # input file basenames
            file_basename = file.split('.')[0]
            df = pd.read_csv(file, sep='\\t')
            # add sample name
            df['Sample'] = file_basename
            dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not read {file}: {e}")

    if not dfs:
        raise ValueError("No valid data found in any files")

    # concatenate 
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # save
    combined_df.to_csv(output_file, index=False)

    # Try to plot
    try:
        import matplotlib.pyplot
        import seaborn
        import visuals
        visuals.intersection_rate_per_sample(combined_df, root_path)
    except ImportError:
        print("Matplotlib not available - skipping plots")

    

if __name__ == "__main__":
    main()