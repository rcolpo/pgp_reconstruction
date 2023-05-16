import pandas as pd
import os
from subprocess import Popen, PIPE, STDOUT

def parse_diamond_output(filepath, remove_extra_columns=True):
    """
    Parses a diamond output file and loads it into a pandas DataFrame.

    Parameters:
    filepath (str): Path to the input file
    remove_extra_columns (bool): If True, unnecessary columns are dropped from the DataFrame (default: True)

    Returns:
    pandas.DataFrame: DataFrame representing the parsed data
    """
    col_names = ['source_gene', 'target_gene', 'pident', 'alignment_length', 'mismatches', 'gap_opens',
                 'query_start', 'query_end', 'sstart', 'send', 'e_value', 'score']

    df = pd.read_csv(filepath, sep='\t', names=col_names)

    if remove_extra_columns:
        df = df[['source_gene', 'target_gene', 'score']]

    return df


def execute_diamond_blast(input_seq_file, seq_type, output_file, db_path, extra_args=None, show_output=True):
    """
    Executes a Diamond BLAST alignment.

    Parameters:
    input_seq_file (str): Path to the input sequence file in FASTA format
    seq_type (str): Type of the sequence ('protein' or 'dna')
    output_file (str): Path to the output file
    db_path (str): Path to the diamond protein database file
    extra_args (str): Additional command-line arguments for diamond (optional)
    show_output (bool): If True, diamond's output is printed to stdout (default: True)

    Returns:
    int: The exit code from diamond

    Notes:
    Default arguments are: --min-score 50 -k0
    """

    assert seq_type in ['protein', 'dna'], "Sequence type must be 'protein' or 'dna'"

    command = ['diamond']

    if seq_type == 'protein':
        command.append('blastp')
    elif seq_type == 'dna':
        command.append('blastx')

    command.extend(['-d', db_path, '-q', input_seq_file, '-o', output_file, '-p', '4'])

    if not extra_args:
        extra_args = "--min-score 50 -k0"

    command.extend(extra_args.split())

    if show_output:
        print(' '.join(command))

    with open(os.devnull, 'w') as null_file:
        try:
            proc = Popen(command, stdout=null_file if not show_output else PIPE, stderr=STDOUT)
            exit_status = proc.wait()
        except OSError:
            exit_status = None

    return exit_status