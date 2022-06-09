import os
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(
    description = """Convert TCRrep data from 10X, Adaptive, or AIRR format to IMGT.
    \n**Adaptive/AIRR**: Can only take single-chain data and won't output all columns by default. Specify extra columns with -e/--extra
    """,
    usage='%(prog)s [-h HELP] -i INPUT -c CHAINS -o OUTPUT [-f INPUT FORMAT] [-e EXTRA COLUMNS] [-s SPECIES]',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

parser.add_argument('-i', '--input', 
                    help='Input TCR filepath', metavar='', type=str, required=True)
parser.add_argument('-c', '--chains', 
                    choices=['alpha-beta', 'alpha', 'beta', 'gamma-delta', 'gamma', 'delta'], 
                    help='Input TCR chain(s). Options are: "alpha-beta", "alpha", "beta", "gamma-delta", "gamma", "delta". Default is "beta"', 
                    metavar='', type=str, required=True)
parser.add_argument('-o', '--output', 
                    help='Output TCR filepath. Extension should be .tsv', metavar='', type=str, required=True)
parser.add_argument('-f', '--format',
                    choices=['10x', 'adaptive', 'airr'],
                    default='10x',
                    help='Format of input TCR data. Options are: "10x", "adaptive", "airr". Default is "10x"', metavar='', type=str, required=False)
parser.add_argument('-e', '--extra',
                    default=[''],
                    help='List of extra input Adaptive/AIRR columns to keep in addition to the CDR3/V/J ones', metavar='', type=list, required=False)
parser.add_argument('-s', '--species', 
                    choices=['human', 'mouse'], 
                    default='human', 
                    help='Options are: "human", "mouse". Default is "human"', metavar='', type=str, required=False)

args = parser.parse_args()


def main():
    input = args.input
    chains = args.chains
    type = args.format
    outpath = args.output
    extra_columns = args.extra
    organism = args.species

    df = load_data(input)

    if type == '10x':
        out_df = tenx_to_imgt(df, chains)
    elif type == 'adaptive':
        out_df = adaptive_to_imgt(df, extra_columns, chains, organism)
    elif type == 'airr':
        out_df = airr_to_imgt(df, extra_columns, chains)
    
    out_df.to_csv(outpath, sep='\t', index=False)


def load_data(in_file):
    '''Load into a dataframe'''
    loading_method = {'.xlsx': pd.read_excel,
                      '.csv': pd.read_csv,
                      '.tsv': pd.read_table}
    extension = str(os.path.splitext(in_file)[1])
    raw_df = loading_method[extension](in_file)
    return raw_df


def tenx_to_imgt(df, chains):
    '''Convert from 10X to IMGT by adding alleles'''
    # TODO: Check if allele already exists (check for * in column)
    if chains == 'alpha-beta':
            df['v_a_gene'] = df['v_a_gene'].apply(lambda x: f"{x}*01")
            df['v_b_gene'] = df['v_b_gene'].apply(lambda x: f"{x}*01")
            df['j_a_gene'] = df['j_a_gene'].apply(lambda x: f"{x}*01")
            df['j_b_gene'] = df['j_b_gene'].apply(lambda x: f"{x}*01")
    elif chains == 'gamma-delta':
            df['v_g_gene'] = df['v_g_gene'].apply(lambda x: f"{x}*01")
            df['v_d_gene'] = df['v_d_gene'].apply(lambda x: f"{x}*01")
            df['j_g_gene'] = df['j_g_gene'].apply(lambda x: f"{x}*01")
            df['j_d_gene'] = df['j_d_gene'].apply(lambda x: f"{x}*01")
    else:
        letter = chains[0]
        v_gene = f'v_{letter}_gene'
        j_gene = f'j_{letter}_gene'
        df[v_gene] = df[v_gene].apply(lambda x: f"{x}*01")
        df[j_gene] = df[j_gene].apply(lambda x: f"{x}*01")
    return df


#############################
#   Modifed from tcrdist3   #
#############################

header_dict = {'alpha':["cdr3_b_aa", "v_a_gene", "j_a_gene", "cdr3_b_nucseq"], 
               'beta' :["cdr3_b_aa", "v_b_gene", "j_b_gene", "cdr3_b_nucseq"],
               'gamma':["cdr3_g_aa", "v_g_gene", "j_g_gene", "cdr3_g_nucseq"],
               'delta':["cdr3_d_aa", "v_d_gene", "j_d_gene", "cdr3_d_nucseq"]}


def valid_cdr3(cdr3):
    '''Return True if all amino acids are part of standard amino acid list'''
    if not isinstance(cdr3, str):
        return False
    else:
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 
                        'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        valid = np.all([aa in amino_acids for aa in cdr3])
    return valid


def adaptive_to_imgt(df, extra_columns, chains, organism):
    '''Convert from Adaptive to IMGT'''
    item_names = header_dict[chains]
    # Parse bio-identity
    ns= {0:"cdr3_aa", 1:"v_gene", 2:"j_gene"}
    # Expand bio_idenitty column to 3-column cdr3,v,j
    new_df = df.copy()
    cdr_v_j = new_df['bio_identity'].str.split("+", expand = True).\
        rename(columns = lambda x: ns[x])
    new_df[[item_names[0], 'v_gene', 'j_gene']] = cdr_v_j
    # Convert Names from Adapative to IMGT
    new_df[item_names[1]] = new_df['v_gene'].\
      apply(lambda x : adaptive_to_imgt[organism].get(x))
    new_df[item_names[2]] = new_df['j_gene'].\
      apply(lambda x : adaptive_to_imgt[organism].get(x))
    # Define columns
    new_df['valid_cdr3'] = new_df[item_names[0]].apply(lambda cdr3: valid_cdr3(cdr3)) 
    new_df = new_df[new_df['valid_cdr3']]
    new_df['productive_frequency'] = pd.to_numeric(new_df['productive_frequency'], errors='coerce')
    new_df['count'] = pd.to_numeric(new_df['templates'], errors='coerce')
    new_df[item_names[3]] = new_df['rearrangement']
    new_df = new_df[[ item_names[0], item_names[1], item_names[2], item_names[3],
                      'productive_frequency', 'count'] + extra_columns]
    return new_df


def airr_to_imgt(df, extra_columns, chains):
    '''Convert from AIRR to IMGT'''
    item_names = header_dict[chains]
    df = df.rename(columns = {'junction_aa':item_names[0],
                              'v_call':item_names[1],
                              'j_call':item_names[2],
                              'junction': item_names[3]})
    df[item_names[1]] = df[item_names[1]].\
      apply(lambda x: x.replace("*00", "*01"))
    df[item_names[2]] = df[item_names[2]].\
      apply(lambda x: x.replace("*00", "*01"))
    df['valid_cdr3'] = df[item_names[0]].apply(lambda cdr3: valid_cdr3(cdr3)) 
    df = df[df['valid_cdr3']]
    # Subset to only necessary columns 
    item_names = item_names + extra_columns
    df = df[item_names]
    return df


if __name__ == "__main__":
    main()
