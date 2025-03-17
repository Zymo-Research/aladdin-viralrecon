import argparse
import pandas as pd
import os

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(
        description="Merge multiple featureCounts TSV files into a single CSV for DESeq2 input."
    )

    # Arguments for input files, output, and customization
    parser.add_argument(
        "-f", "--featureCounts",
        default=None,
        nargs='+',  
        help="Path to the featureCounts files."
    )
    parser.add_argument(
        "-o", "--output",
        default="merged_featurecounts.csv",
        help="Output CSV filename (default: 'merged_featurecounts.csv')."
    )
    parser.add_argument(
        "-d", "--delimiter",
        default="\t",  
        help="Delimiter used in the input files (default: tab)."
    )
    parser.add_argument(
        "--gene_col",
        default="Geneid",  
        help="Name of the gene identifier column (default: 'Geneid')."
    )
    parser.add_argument(
        "--gene_name_col",
        default="gene_name",  
        help="Name of the gene name column (default: 'gene_name')."
    )
    parser.add_argument(
        "--count_col",
        type=int,
        default=-1,
        help=("Zero-indexed column for the counts. Use -1 for the last column "
              "(default: -1).")
    )

    # Parsing the arguments
    args = parser.parse_args()

    files = args.featureCounts
    if not files:
        print("No files found:", args.featureCounts)
        return

    dfs = []  
    gene_name_mapping = {}  # Dictionary to store gene_name mapping

    for f in files:
        try:
            print(f"Reading file: {f}")
            df = pd.read_csv(f, sep=args.delimiter, comment="#", engine='python')

            print(f"Columns in {f}: {df.columns.tolist()}")
            print(df.head())

            if args.gene_col not in df.columns:
                print(f"Warning: '{args.gene_col}' column not found in {f}")
                continue

            if args.count_col == -1:
                count_col_name = df.columns[-1]
            else:
                count_col_name = df.columns[args.count_col]

            if count_col_name not in df.columns:
                print(f"Warning: Count column {count_col_name} not found in {f}")
                continue

            sample_name = os.path.splitext(os.path.basename(f))[0].replace(".featureCounts", "")

            # Extract gene_name column if present
            if args.gene_name_col in df.columns:
                df_gene_name = df[[args.gene_col, args.gene_name_col]]
                gene_name_mapping.update(df_gene_name.set_index(args.gene_col).to_dict()[args.gene_name_col])

            df = df[[args.gene_col, count_col_name]].rename(columns={count_col_name: sample_name})
            dfs.append(df)

        except Exception as e:
            print(f"Error reading {f}: {e}")

    if not dfs:
        print("No valid featureCounts files processed.")
        return

    # Merge all DataFrames on the gene identifier
    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on=args.gene_col, how="outer")

    # Fill missing counts with 0
    merged_df.fillna(0, inplace=True)

    # Add gene_name column
    merged_df = pd.merge(merged_df, pd.DataFrame(gene_name_mapping.items(), columns=[args.gene_col, args.gene_name_col]), on=args.gene_col, how="left")
    # merged_df.rename(columns={'gene': 'gene_name'}, inplace=True)
    # merged_df.fillna('Unknown', inplace=True)
    merged_df.dropna(inplace=True)

    # Save to CSV
    merged_df.to_csv(args.output, index=False)
    print(f"Merged counts file saved as {args.output}")

if __name__ == "__main__":
    main()