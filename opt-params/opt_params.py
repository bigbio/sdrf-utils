import logging
import random
import re
import uuid
from collections import Counter
from typing import Tuple, List, Dict
import json
import os

import numpy as np
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
from pandas import DataFrame
import click

logging.basicConfig(level=logging.INFO)

folder_plots_uui = ""
folder_sperator = os.path.sep

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

def run_sage(fragment_tolerance: int = 0, precursor_tolerance: int = 0, fragment_type: str = "ppm",
             mzml_files: list = [], fasta_path: str = "", sage_config_file: str = None, use_file_values: bool = True ) -> DataFrame:

    if sage_config_file is None:
        raise ValueError("The sage config file is required.")

    with open(sage_config_file) as f:
        data = json.load(f)

    if fragment_tolerance != 0 and precursor_tolerance != 0 or not use_file_values:
        data["precursor_tol"]["ppm"] = [int(-1 * (precursor_tolerance)), int(precursor_tolerance)]
        if fragment_type == "ppm":
            data["fragment_tol"][fragment_type] = [int(-1 * (fragment_tolerance)), int(fragment_tolerance)]
        else:
            data["fragment_tol"][fragment_type] = [-1 * (fragment_tolerance), fragment_tolerance]
    else:
        logging.info("Using the values from the file.")
        if "ppm" in data["precursor_tol"]:
            precursor_tolerance = data["precursor_tol"]["ppm"][1]
        else:
            precursor_tolerance = data["precursor_tol"]["da"][1]

        if "ppm" in data["fragment_tol"]:
            fragment_tolerance = data["fragment_tol"]["ppm"][1]
        else:
            fragment_tolerance = data["fragment_tol"]["da"][1]

    data["mzml_paths"] = mzml_files

    if os.path.exists(fasta_path):
        data["database"]["fasta"] = fasta_path
    else:
        logging.error(f"File {fasta_path} does not exist.")
        raise FileNotFoundError(f"File {fasta_path} does not exist.")

    temp_sage_file = str(uuid.uuid4()) + ".json"
    with open(temp_sage_file, "w") as f:
        json.dump(data, f, indent=4)

    logging.info("Running SAGE with fragment tolerance: {} and precursor tolerance: {}".format(fragment_tolerance,precursor_tolerance))

    result = subprocess.run(["sage", temp_sage_file, "--write-pin"], capture_output=True, text=True)
    os.remove(temp_sage_file)

    if result.returncode != 0:
        logging.error("Error running SAGE.")
        logging.error(result.stderr)
        raise ValueError("Error running SAGE.")

    sage_table = pd.read_csv("results.sage.tsv", sep="\t")
    sage_table = compute_entrapment_qvalues(sage_table)
    return sage_table


def extract_ptms(sage_table_target):

    def extract_modifications(peptide):
        return re.findall(r'([A-Z]\[\+[0-9.]+\])', peptide)

    # Apply the function to the peptides column and flatten the list of lists
    modifications = [mod for peptide in sage_table_target['peptide'] for mod in extract_modifications(peptide)]

    # Count the occurrences of each modification
    modification_counts = Counter(modifications)

    # convert modification counts to a dictionary key is the mod name, value is the count
    modification_counts_dict = dict(modification_counts)

    return modification_counts_dict


def get_stats_from_sage(sage_table: pd.DataFrame, fragment_tolerance: int = 0, precursor_tolerance: int = 0,
                        number_psms: int = 0) -> dict:
    sage_table = sage_table[sage_table['spectrum_q'] <= 0.01]
    sage_table = sage_table[sage_table['entrapment_qvalue'] <= 0.01]
    decoy_filter = sage_table['proteins'].str.contains("DECOY_")
    sage_table_target = sage_table[~decoy_filter]
    entrap_peptides = sage_table_target[sage_table_target['proteins'].str.contains("ENTRAP")]
    sage_table_decoy = sage_table[decoy_filter]
    precursor_std_error = sage_table_target['precursor_ppm'].std() * 4
    fragment_std_error = sage_table_target['fragment_ppm'].std() * 4

    plot_distribution(sage_table_target, sage_table_decoy, entrap_peptides, fragment_tolerance, precursor_tolerance)

    ptms_dist = extract_ptms(sage_table_target)

    stats = {
        "fragment_tolerance": fragment_tolerance,
        "precursor_tolerance": precursor_tolerance,
        "total_peptides": len(sage_table['peptide'].unique()),
        "total_entrap_peptides": len(entrap_peptides['peptide'].unique()),
        "total_decoys": len(sage_table_decoy['peptide'].unique()),
        "total_target": len(sage_table_target['peptide'].unique()),
        "precursor_std_error_plus4": precursor_std_error,
        "fragment_std_error_plus4": fragment_std_error,
        "number_psms": number_psms,
    }

    for ky, val in ptms_dist.items():
        stats[f"{ky}"] = val
        stats[f"{ky}_percentage"] = (val / len(sage_table_target['peptide'])) * 100

    return stats


def plot_distribution(sage_table_target, sage_table_decoy, entrap_peptides, fragment_tolerance, precursor_tolerance):
    plt.hist(sage_table_target['spectrum_q'], bins=20, alpha=0.5, label='Target')
    plt.hist(sage_table_decoy['spectrum_q'], bins=20, alpha=0.5, label='Decoy')
    plt.hist(entrap_peptides['spectrum_q'], bins=20, alpha=0.5, label='Entrap')
    number_psms = len(sage_table_target['peptide'])
    number_peptides = len(sage_table_target['peptide'].unique())
    plt.title("Q value distribution - Number of peptides: {} - Number of PSMs: {}".format(number_peptides, number_psms))
    plt.legend(loc='upper right')
    plt.xlabel(
        "Q value - fragment tolerance: {} - precursor tolerance: {}".format(fragment_tolerance, precursor_tolerance))
    file_name = folder_plots_uui + folder_sperator + "q_value_distribution-{}-{}.png".format(fragment_tolerance, precursor_tolerance)
    plt.savefig(file_name)
    logging.info("Showing Q value distribution plot.")
    plt.close()
    # plt.show()


def compute_entrapment_qvalues(sage_data_frame: pd.DataFrame) -> pd.DataFrame:
    """ compute the q values using the entrapments peptides instead of decoys.
    :param sage_data_frame: pandas data frame with the sage results
    """

    # Use inplace=True for sorting
    sage_data_frame.sort_values(by='sage_discriminant_score', ascending=False, inplace=True)

    # Use a more descriptive name for columns
    sage_data_frame['is_entrap'] = sage_data_frame['proteins'].apply(
        lambda x: 'entrap' if x.startswith('ENTRAP_') else 'target')
    sage_data_frame['is_decoy'] = sage_data_frame['proteins'].apply(
        lambda x: 'decoy' if x.startswith('DECOY_') else 'target')

    # Compute the q values using entrapment peptides
    sage_data_frame['spectrum_entrapment_q'] = 0
    sage_data_frame['cumulative_target'] = (sage_data_frame['is_entrap'] == 'target').cumsum()
    sage_data_frame['cumulative_entrap'] = (sage_data_frame['is_entrap'] == 'entrap').cumsum()
    sage_data_frame['FDR_ENTRAP'] = sage_data_frame['cumulative_entrap'] / sage_data_frame['cumulative_target']

    # Initialize the q-values with a large number
    sage_data_frame['entrapment_qvalue'] = 0

    # Use vectorized operation for calculating q-values
    sage_data_frame['entrapment_qvalue'] = sage_data_frame['FDR_ENTRAP'].expanding().min()

    return sage_data_frame


def compute_best_combination(sage_table: pd.DataFrame) -> int:
    """
    Compute the difference between the decoys and entrapmet qvalues.
    Sum all the differences and return the value.
    :param sage_data_frame:
    :return:
    """
    # sage_data_frame['q_diff'] = (sage_data_frame['spectrum_q'] - sage_data_frame['spectrum_entrapment_q']).abs()
    # make copy of the data frame
    current_data_frame = sage_table.copy()
    current_data_frame = current_data_frame[current_data_frame['spectrum_q'] <= 0.01]
    current_data_frame = current_data_frame[current_data_frame['entrapment_qvalue'] <= 0.01]
    filter = current_data_frame['proteins'].str.contains("DECOY_")
    sage_table_target = current_data_frame[~filter]
    return len(sage_table_target['peptide'])


def combined_search(results: list = [], start_fragment_tolerance: int = 0, start_precursor_tolerance: int = 0,
                    min_fragment_tolerance: int = 1, max_fragment_tolerance: int = 50,
                    min_precursor_tolerance: int = 10,
                    max_precursor_tolerance: int = 100, num_psms: int = 0, fragment_type: str = "ppm",
                    mzml_files: list = [],
                    grid_frag_steps: int=10,
                    grid_prec_steps: int=5, max_iterations=10, search_radius=1, fasta_file="", initial_temp=100,
                    cooling_rate=0.95, sage_config_file=None) -> Tuple[int, int, List[Dict]]:
    def acceptance_probability(old_value, new_value, temperature):
        if new_value > old_value:
            return 1.0
        return np.exp((new_value - old_value) / temperature)

    best_fragment_tolerance = int(start_fragment_tolerance)
    best_precursor_tolerance = int(start_precursor_tolerance)
    best_value = int(num_psms)

    # Coarse Grid Search
    fragment_tolerances = np.linspace(min_fragment_tolerance, max_fragment_tolerance, grid_frag_steps).astype(int)
    precursor_tolerances = np.linspace(min_precursor_tolerance, max_precursor_tolerance, grid_prec_steps).astype(int)

    for ft in fragment_tolerances:
        for pt in precursor_tolerances:
            sage_table = run_sage(fragment_tolerance=ft, precursor_tolerance=pt,
                                  fragment_type=fragment_type, mzml_files=mzml_files, fasta_path=fasta_file, sage_config_file=sage_config_file, use_file_values=False)
            new_value = compute_best_combination(sage_table)
            results.append(get_stats_from_sage(sage_table, ft, pt, new_value))

            if new_value > best_value:
                best_fragment_tolerance = ft
                best_precursor_tolerance = pt
                best_value = new_value

    # Simulated Annealing
    temperature = initial_temp
    current_fragment_tolerance = best_fragment_tolerance
    current_precursor_tolerance = best_precursor_tolerance

    for _ in range(max_iterations):
        fragment_tolerance = np.random.randint(max(min_fragment_tolerance, current_fragment_tolerance - search_radius),
                                               min(max_fragment_tolerance, current_fragment_tolerance + search_radius))
        precursor_tolerance = np.random.randint(
            max(min_precursor_tolerance, current_precursor_tolerance - search_radius),
            min(max_precursor_tolerance, current_precursor_tolerance + search_radius))

        sage_table = run_sage(fragment_tolerance=fragment_tolerance, precursor_tolerance=precursor_tolerance,
                              fragment_type=fragment_type, mzml_files=mzml_files, fasta_path=fasta_file, use_file_values=False)
        new_value = compute_best_combination(sage_table)
        results.append(get_stats_from_sage(sage_table, fragment_tolerance, precursor_tolerance, new_value))

        if acceptance_probability(best_value, new_value, temperature) > random.random():
            current_fragment_tolerance = fragment_tolerance
            current_precursor_tolerance = precursor_tolerance
            best_fragment_tolerance = fragment_tolerance
            best_precursor_tolerance = precursor_tolerance
            best_value = new_value

        temperature *= cooling_rate

    return best_fragment_tolerance, best_precursor_tolerance, results


@click.command("tolerances", help="Optimize the fragment and precursor tolerances using the SAGE algorithm.")
@click.option("--fragment-type", default="ppm", help="The type of fragment tolerance to use (ppm or da)")
@click.option("--mzml-path", default=".", help="The path to the mzML files to use for the SAGE analysis.")
@click.option("--initial-fragment-tolerance", help="The initial fragment tolerance to use for the optimization.", type=int)
@click.option("--initial-precursor-tolerance", help="The initial precursor tolerance to use for the optimization.", type=int)
@click.option("--min-fragment-tolerance", default=1, help="The minimum fragment tolerance to consider.", type=int)
@click.option("--max-fragment-tolerance", default=50, help="The minimum precursor tolerance to consider.", type=int)
@click.option("--min-precursor-tolerance", default=10, help="The minimum precursor tolerance to consider.", type=int)
@click.option("--max-precursor-tolerance", default=50, help="The maximum precursor tolerance to consider.", type=int)
@click.option("--fasta-file", default="Homo-sapiens-uniprot-reviewed-contaminants-entrap-decoy-20240615.fasta", help="The path to the fasta file to use for the SAGE analysis.")
@click.option("--sage-config-file", default="sage-general.json", help="The path to the Sage config file to use for the SAGE analysis.")
@click.option("--max-iterations", default=10, help="The maximum number of iterations to run the optimization.", type=int)
def tolerances(fragment_type:str, mzml_path: str, initial_fragment_tolerance: int, initial_precursor_tolerance: int,
         min_fragment_tolerance: int, max_fragment_tolerance: int, min_precursor_tolerance: int, max_precursor_tolerance: int,
         fasta_file, sage_config_file, max_iterations: int):
    results = []

    # detect absolute all the mzML files in the mzml-path
    mzml_files = [os.path.join(mzml_path, f) for f in os.listdir(mzml_path) if f.endswith(".mzML")]

    # generate uui unique identifier for the folder with the plots
    global folder_plots_uui
    folder_plots_uui = str(uuid.uuid4())
    os.makedirs(folder_plots_uui)

    num_psms = 0
    best_precursor_tolerance = 0
    best_fragment_tolerance = 0

    if sage_config_file is None:
        sage_config_file = "general-sage.json"


    if initial_fragment_tolerance is not None and initial_precursor_tolerance is not None:
        sage_table = run_sage(int(initial_fragment_tolerance), int(initial_precursor_tolerance), fragment_type,
                              mzml_files, fasta_path = fasta_file, sage_config_file=sage_config_file, use_file_values=False)
        num_psms = compute_best_combination(sage_table)
        results.append(
            get_stats_from_sage(sage_table, initial_fragment_tolerance, initial_precursor_tolerance, num_psms))
        best_precursor_tolerance = initial_precursor_tolerance
        best_fragment_tolerance = initial_fragment_tolerance

    best_fragment_tolerance, best_precursor_tolerance, results = combined_search(results=results,
                                                                               start_fragment_tolerance=best_fragment_tolerance,
                                                                               start_precursor_tolerance=best_precursor_tolerance,
                                                                               min_fragment_tolerance=min_fragment_tolerance,
                                                                               max_fragment_tolerance=max_fragment_tolerance,
                                                                               min_precursor_tolerance=min_precursor_tolerance,
                                                                               max_precursor_tolerance=max_precursor_tolerance,
                                                                               num_psms=num_psms,
                                                                               fragment_type=fragment_type,
                                                                               mzml_files=mzml_files,
                                                                               max_iterations=max_iterations, fasta_file=fasta_file,
                                                                                 sage_config_file=sage_config_file)

    print("Best tolerances found: Fragment tolerance: {} - Precursor tolerance: {}".format(best_fragment_tolerance,
                                                                                           best_precursor_tolerance))

    # Write the results to a file
    results_df = pd.DataFrame(results)

    # Aggregate the data to handle duplicate entries by averaging number_psms
    df_agg = results_df.groupby(['fragment_tolerance', 'precursor_tolerance']).agg({'number_psms': 'mean'}).reset_index()

    # Pivot the aggregated dataframe to prepare for plotting
    df_pivot = df_agg.pivot(index='fragment_tolerance', columns='precursor_tolerance', values='number_psms')

    # Plot the results
    ax = df_pivot.plot(kind='line', marker='o', figsize=(10, 6))

    # Set plot title and labels
    ax.set_title('Number of PSMs vs. Fragment Tolerances')
    ax.set_xlabel('Fragment Tolerance')
    ax.set_ylabel('Number of PSMs')

    # Set legend title
    ax.legend(title='Precursor Tolerances')

    # Save the plot to a file
    output_file = f"{folder_plots_uui}/final_results_tolerances.png"
    plt.savefig(output_file)
    plt.close()

    print(f"Plot saved to {output_file}")

    results_df.to_csv("sage_results_tolerances.tsv", sep="\t", index=False)


@click.command("ptms", help="Extract PTMs from SAGE results.")
@click.option("--mzml-path", default=".", help="The path to the mzML files to use for the SAGE analysis.")
@click.option("--fasta-file", default="Homo-sapiens-uniprot-reviewed-contaminants-entrap-decoy-20240615.fasta", help="The path to the fasta file to use for the SAGE analysis.")
@click.option("--sage-config-file", default="general-sage-ptms.json", help="The path to the Sage config file to use for the SAGE analysis.", required = True)
@click.option("--open-search", help="Open search analysis", is_flag=True)
def ptms(mzml_path: str, fasta_file: str, sage_config_file: str, open_search: bool):
    # detect absolute all the mzML files in the mzml-path

    if not os.path.exists(mzml_path):
        logging.error(f"Folder {mzml_path} does not exist.")
        raise FileNotFoundError(f"Folder {mzml_path} does not exist.")

    if not os.path.exists(fasta_file):
        logging.error(f"File {fasta_file} does not exist.")
        raise FileNotFoundError(f"File {fasta_file} does not exist.")

    if not os.path.exists(sage_config_file):
        logging.error(f"File {sage_config_file} does not exist.")
        raise FileNotFoundError(f"File {sage_config_file} does not exist.")
    else:
        sage_params = json.load(open(sage_config_file))
        if "ppm" in sage_params["precursor_tol"]:
            precursor_tolerance = sage_params["precursor_tol"]["ppm"][1]
        else:
            precursor_tolerance = sage_params["precursor_tol"]["da"][1]

        if "ppm" in sage_params["fragment_tol"]:
            fragment_tolerance = sage_params["fragment_tol"]["ppm"][1]
        else:
            fragment_tolerance = sage_params["fragment_tol"]["da"][1]


    mzml_files = [os.path.join(mzml_path, f) for f in os.listdir(mzml_path) if f.endswith(".mzML")]

    global folder_plots_uui
    folder_plots_uui = str(uuid.uuid4())
    os.makedirs(folder_plots_uui)

    sage_table = run_sage(fragment_type="ppm", mzml_files=mzml_files, fasta_path=fasta_file, sage_config_file=sage_config_file, use_file_values=True)
    num_psms = compute_best_combination(sage_table)
    stats = get_stats_from_sage(sage_table, fragment_tolerance=fragment_tolerance, precursor_tolerance=precursor_tolerance, number_psms=num_psms)

    # check verbose the PTMs that percentace bigger than 1% of the psms
    logging.info("PTMs with percentage bigger than 1% of the psms:")
    for key, value in stats.items():
        if "_percentage" in key and value > 1:
            logging.info(f"{key} - {value}")

    # print dictionary stats to a file as key value.
    with open("sage_stats_ptms.tsv", "w") as f:
        for key, value in stats.items():
            f.write(f"{key}\t{value}\n")


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    Main function to run the CLI
    """
    pass

cli.add_command(tolerances)
cli.add_command(ptms)

if __name__ == "__main__":
    cli()
