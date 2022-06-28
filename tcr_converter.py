import os
import argparse
import pandas as pd
import numpy as np
import sys
from tcrdist.swap_gene_name import adaptive_to_imgt

parser = argparse.ArgumentParser(
    description = """Convert TCRrep data from Adaptive v4 or AIRR format to IMGT.
    \n**Adaptive/AIRR**: Can only take single-chain data and won't output all columns 
    by default. Specify extra columns with -e/--extra
    """,
    usage='%(prog)s [-h HELP] -i INPUT [-c CHAIN] [-o OUTPUT] [-e EXTRA COLUMNS] [-s SPECIES]',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

parser.add_argument('-i', '--input', 
                    help='Input TCR filepath', metavar='', type=str, required=True)
parser.add_argument('-c', '--chain', 
                    choices=['alpha', 'beta', 'gamma', 'delta'],
                    default=['beta'],
                    help='Input TCR chain. Options are: "alpha", "beta", "gamma", "delta". Default is "beta"', 
                    metavar='', type=str, required=False)
parser.add_argument('-o', '--output', 
                    default='./tcr_converter_out.csv',
                    help='Output TCR filepath. Extension should be .csv', metavar='', type=str, required=False)
parser.add_argument('-e', '--extra',
                    default=[],
                    help='List of extra input Adaptive/AIRR columns to keep in addition to the CDR3/V/J ones', metavar='', type=list, required=False)
parser.add_argument('-s', '--species', 
                    choices=['human', 'mouse'], 
                    default='human', 
                    help='Options are: "human", "mouse". Default is "human"', metavar='', type=str, required=False)

args = parser.parse_args()


def main():
    input = args.input
    chains = args.chain
    outpath = args.output
    extra_columns = args.extra
    organism = args.species

    format = get_input_format(in_file=input)

    if format == 'adaptive_v4':
        df = load_data(input)
        out_df = adapt_to_imgt(df, extra_columns, chains, organism)
    elif format == 'airr':
        df = load_data(input)
        out_df = airr_to_imgt(df, extra_columns, chains)
    else:
        print('Only AIRR and Adaptive v4 supported at this time.')
        sys.exit(0)
    
    print(f'Writing output to: "{outpath}"')
    out_df.to_csv(outpath, index=False)

      
# Columns unique to these formats
headers = {
    # https://clients.adaptivebiotech.com/assets/downloads/
    # immunoSEQ_AnalyzerManual.pdf 
    'adaptive_v2': {'nucleotide', 'aminoAcid', 'count (templates/reads)',
                'frequencyCount (%)', 'cdr3Length', 'vMaxResolved', 
                'vFamilyName', 'vGeneName', 'vGeneAllele', 'vFamilyTies', 
                'vGeneNameTies', 'vGeneAlleleTies', 'dMaxResolved', 
                'dFamilyName', 'dGeneName', 'dGeneAllele', 'dFamilyTies', 
                'dGeneNameTies', 'dGeneAlleleTies', 'jMaxResolved', 
                'jFamilyName', 'jGeneName', 'jGeneAllele', 'jFamilyTies', 
                'jGeneNameTies', 'jGeneAlleleTies', 'vDeletion', 'n1Insertion', 
                'd5Deletion', 'd3Deletion', 'n2Insertion', 'jDeletion',
                'vIndex', 'n1Index', 'dIndex', 'n2Index', 'jIndex',
                'estimatedNumberGenomes', 'sequenceStatus', 'cloneResolved', 
                'vOrphon', 'dOrphon', 'jOrphon', 'vFunction', 'dFunction', 
                'jFunction', 'fractionNucleated', 'vAlignLength', 
                'vAlignSubstitutionCount', 'vAlignSubstitutionIndexes', 
                'vAlignSubstitutionGeneThreePrimeIndexes', 'vSeqWithMutations'},

    # https://www.adaptivebiotech.com/wp-content/uploads/2019/07/
    # MRK-00342_immunoSEQ_TechNote_DataExport_WEB_REV.pdf
    'adaptive_v4': {'rearrangement', 'extended_rearrangement', 'bio_identity', 
                'amino_acid','templates', 'frame_type', 'rearrangement_type', 
                'productive_frequency', 'cdr1_start_index', 
                'cdr1_rearrangement_length', 'cdr2_start_index', 
                'cdr2_rearrangement_length', 'cdr3_start_index', 'cdr3_length',
                'v_index', 'n1_index', 'd_index', 'n2_index', 'j_index', 
                'v_deletions', 'n2_insertions', 'd3_deletions', 'd5_deletions', 
                'n1_insertions', 'j_deletions', 'chosen_j_allele', 
                'chosen_j_family', 'chosen_j_gene', 'chosen_v_allele', 
                'chosen_v_family', 'chosen_v_gene', 'd_allele', 
                'd_allele_ties', 'd_family', 'd_family_ties', 'd_gene_ties',
                'd_resolved', 'j_allele', 'j_allele_ties', 'j_family', 
                'j_family_ties', 'j_gene_ties', 'j_resolved', 'v_allele', 
                'v_allele_ties', 'v_family', 'v_family_ties', 'v_gene_ties', 
                'v_resolved', 'Frequency', 'Reads', 'Rearrangement', 
                'v_shm_count', 'v_shm_indexes', 'Antibody'},

    # https://support.10xgenomics.com/single-cell-vdj/software/pipelines/
    # latest/output/annotation#contig-annotation
    '10x': {'barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 
            'chain', 'c_gene', 'full_length', 'fwr1_nt', 'cdr1_nt', 'fwr2_nt', 
            'cdr2_nt', 'fwr3_nt', 'cdr3_nt', 'fwr4_nt', 'reads', 'umis', 
            'raw_clonotype_id', 'raw_consensus_id', 'exact_subclonotype_id'},
    
    # From a real file
    'airr': {'sequence_id', 'sequence', 'sequence_aa', 'rev_comp', 'vj_in_frame', 
            'stop_codon', 'complete_vdj', 'locus', 'v_call', 'd_call', 'd2_call', 
            'j_call', 'c_call', 'sequence_alignment', 'sequence_alignment_aa', 
            'germline_alignment', 'germline_alignment_aa', 'junction', 'junction_aa', 
            'np1', 'np1_aa', 'np2', 'np2_aa', 'np3', 'np3_aa', 'cdr1_aa', 'cdr2_aa', 
            'cdr3_aa', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa', 'v_score', 
            'v_identity', 'v_support', 'v_cigar', 'd_score', 'd_identity',
            'd_support', 'd_cigar', 'd2_score', 'd2_identity', 'd2_support',
            'd2_cigar', 'j_score', 'j_identity', 'j_support', 'j_cigar', 'c_score',
            'c_identity', 'c_support', 'c_cigar', 'v_sequence_start',
            'v_sequence_end', 'v_germline_start', 'v_germline_end',
            'v_alignment_start', 'v_alignment_end', 'd_sequence_start',
            'd_sequence_end', 'd_germline_start', 'd_germline_end',
            'd_alignment_start', 'd_alignment_end', 'd2_sequence_start',
            'd2_sequence_end', 'd2_germline_start', 'd2_germline_end',
            'd2_alignment_start', 'd2_alignment_end', 'j_sequence_start',
            'j_sequence_end', 'j_germline_start', 'j_germline_end',
            'j_alignment_start', 'j_alignment_end', 'cdr1_start', 'cdr1_end',
            'cdr2_start', 'cdr2_end', 'cdr3_start', 'cdr3_end', 'fwr1_start',
            'fwr1_end', 'fwr2_start', 'fwr2_end', 'fwr3_start', 'fwr3_end',
            'fwr4_start', 'fwr4_end', 'v_sequence_alignment',
            'v_sequence_alignment_aa', 'd_sequence_alignment',
            'd_sequence_alignment_aa', 'd2_sequence_alignment',
            'd2_sequence_alignment_aa', 'j_sequence_alignment',
            'j_sequence_alignment_aa', 'c_sequence_alignment',
            'c_sequence_alignment_aa', 'v_germline_alignment',
            'v_germline_alignment_aa', 'd_germline_alignment',
            'd_germline_alignment_aa', 'd2_germline_alignment',
            'd2_germline_alignment_aa', 'j_germline_alignment',
            'j_germline_alignment_aa', 'c_germline_alignment',
            'c_germline_alignment_aa', 'junction_length', 'junction_aa_length',
            'np1_length', 'np2_length', 'np3_length', 'n1_length', 'n2_length',
            'n3_length', 'p3v_length', 'p5d_length', 'p3d_length', 'p5d2_length',
            'p3d2_length', 'p5j_length', 'consensus_count', 'duplicate_count',
            'cell_id', 'clone_id', 'rearrangement_id', 'repertoire_id',
            'sample_processing_id', 'data_processing_id', 'germline_database',
            'rearrangement_set_id'}
}


def load_header(in_file):
    '''Returns set of column headers'''
    load = {'.xlsx': pd.read_excel,
            '.csv': pd.read_csv,
            '.tsv': pd.read_table}
    ext = str(os.path.splitext(in_file)[1])
    in_cols = load[ext](in_file, index_col=0, nrows=0).columns
    in_cols = set(in_cols)
    return in_cols


def get_input_format(in_file):
    '''Detects standard format based on column headers'''
    in_headers = load_header(in_file=in_file)
    format = []

    # Check for intersections (that aren't empty sets) with known headers
    if not in_headers.isdisjoint(headers['10x']):
        format.append('10x')
    if not in_headers.isdisjoint(headers['adaptive_v2']):
        format.append('adaptive_v2')
    if not in_headers.isdisjoint(headers['adaptive_v4']):
        format.append('adaptive_v4')
    if not in_headers.isdisjoint(headers['airr']):
        format.append('airr')

    # If doesn't intersect with any
    if len(format) == 0:
        print('Input file does not have any standard headers, please convert by '\
            'hand. Exiting now.')
        sys.exit(0)
    # If intersects with multiple
    elif len(format) > 1:
        print('Found column names from multiple formats:', format, '\nInput file'\
            'is not a standard format, please convert by hand. Exiting now.')
        sys.exit(0)
    else:
        print('Detected format:', format[0])

    return format[0]


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
