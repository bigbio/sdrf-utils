import logging
import math
import random
import re
import uuid
from collections import Counter
from typing import Tuple, List, Dict
import json
import os
from sagepy.core import EnzymeBuilder, SageSearchConfiguration, validate_mods, validate_var_mods, \
    Scorer, RawSpectrum, SpectrumProcessor, Precursor, Tolerance, IonType
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import DataFrame
import click
from pyteomics import mzml
import concurrent.futures as cf

from tqdm import tqdm

logging.basicConfig(level=logging.INFO)

folder_plots_uui = ""
folder_sperator = os.path.sep

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

indexed_db = None


def generate_indexed_db(fasta_path, sage_config_file):
    if sage_config_file is None:
        raise ValueError("The sage config file is required.")

    with open(sage_config_file) as f:
        param_data = json.load(f)
        cleave_at = param_data["database"]["enzyme"]["cleave_at"]
        missed_cleavages = param_data["database"]["enzyme"]["missed_cleavages"]
        min_len = param_data["database"]["enzyme"]["min_len"]
        max_len = param_data["database"]["enzyme"]["max_len"]
        restrict = param_data["database"]["enzyme"]["restrict"]
        static_mods = param_data["database"]["static_mods"]
        variable_mods = param_data["database"]["variable_mods"]
        fragment_min_mz = param_data["database"]["fragment_min_mz"]
        fragment_max_mz = param_data["database"]["fragment_max_mz"]
        peptide_min_mass = param_data["database"]["peptide_min_mass"]
        peptide_max_mass = param_data["database"]["peptide_max_mass"]
        ion_kinds = param_data["database"]["ion_kinds"]
        min_ion_index = param_data["database"]["min_ion_index"]
        max_variable_mods = param_data["database"]["max_variable_mods"]
        generate_decoys = param_data["database"]["generate_decoys"]
        decoy_tag = param_data["database"]["decoy_tag"]

    # configure a trypsin-like digestor of fasta files
    enzyme_builder = EnzymeBuilder(
        missed_cleavages=missed_cleavages,
        min_len=min_len,
        max_len=max_len,
        cleave_at=cleave_at,
        restrict=restrict,
        c_terminal=True
    )

    # # generate static cysteine modification
    # static_mods = {k: v for k, v in [SAGE_KNOWN_MODS.cysteine_static()]}
    #
    # # generate variable methionine modification
    # variable_mods = {k: v for k, v in [SAGE_KNOWN_MODS.methionine_variable()]}

    # generate SAGE compatible mod representations
    static = validate_mods(static_mods)
    variab = validate_var_mods(variable_mods)

    with open(fasta_path, 'r') as infile:
        fasta = infile.read()

    # generate IonType class
    ion_kinds = [IonType(ion).get_py_ptr() for ion in ion_kinds]

    # set-up a config for a sage-database
    sage_config = SageSearchConfiguration(
        fasta=fasta,
        static_mods=static,
        variable_mods=variab,
        enzyme_builder=enzyme_builder,
        fragment_min_mz=fragment_min_mz,
        fragment_max_mz=fragment_max_mz,
        peptide_min_mass=peptide_min_mass,
        peptide_max_mass=peptide_max_mass,
        ion_kinds=ion_kinds,
        generate_decoys=generate_decoys,
        decoy_tag=decoy_tag,
        min_ion_index=min_ion_index,
        max_variable_mods=max_variable_mods,
        bucket_size=int(np.power(2, 14))
    )

    # generate the database for searching against
    global indexed_db
    indexed_db = sage_config.generate_indexed_database()


def mass_to_mod(mass: float) -> str:
    """ Convert a mass to a UNIMOD modification annotation.

    Args:
        mass: a mass in Da

    Returns:
        a UNIMOD modification annotation
    """
    maybe_key = int(np.round(mass))
    mod_dict = {
        1: '[UNIMOD:7]',
        42: '[UNIMOD:1]',
        57: '[UNIMOD:4]',
        80: '[UNIMOD:21]',
        16: '[UNIMOD:35]',
        119: '[UNIMOD:312]',
    }
    # try to translate to UNIMOD annotation
    try:
        return mod_dict[maybe_key]
    except KeyError:
        raise KeyError(f"Rounded mass not in dict: {maybe_key}")


# Adapted from https://github.com/theGreatHerrLebert/sagepy/blob/v0.2.20-alpha/sagepy/sagepy/core/peptide.py#L125C1
# Will update in next release https://github.com/theGreatHerrLebert/sagepy/pull/12/files
def to_unimod_sequence(modifications, sequence) -> str:
    """ Get Peptide sequence with UNIMOD modification annotations.

    Returns:
        str: Peptide sequence with UNIMOD modification annotations.
    """

    seq = ''

    for i, (s, m) in enumerate(zip(sequence, modifications)):
        if m != 0:
            if i == 0:
                if mass_to_mod(m) == '[UNIMOD:1]':
                    seq += f'{mass_to_mod(m)}{s}'
                else:
                    seq += f'{s}{mass_to_mod(m)}'
            else:
                seq += f'{s}{mass_to_mod(m)}'
        else:
            seq += s

    return seq


def run_one_raw(spec_file, param_data):
    global indexed_db

    raw_spectrums = mzml.MzML(spec_file)
    res = []

    # loading parameters
    fragment_min_mz = param_data["database"]["fragment_min_mz"]
    fragment_max_mz = param_data["database"]["fragment_max_mz"]
    deisotope_bool = param_data["deisotope"]
    min_matched_peaks = param_data["min_matched_peaks"]
    min_isotope_err = param_data["min_isotope_err"]
    max_isotope_err = param_data["max_isotope_err"]
    min_precursor_charge = param_data["min_precursor_charge"]
    max_precursor_charge = param_data["max_precursor_charge"]
    chimera = param_data["chimera"]
    report_psms = param_data["report_psms"]
    annotate_matches = param_data["annotate_matches"]
    max_fragment_charge = param_data["max_fragment_charge"]

    if "ppm" in param_data["precursor_tol"]:
        precursor_tolerance = Tolerance(
            ppm=(param_data["precursor_tol"]["ppm"][0], param_data["precursor_tol"]["ppm"][1]))
    else:
        precursor_tolerance = Tolerance(
            ppm=(param_data["precursor_tol"]["da"][0], param_data["precursor_tol"]["da"][1]))

    if "ppm" in param_data["fragment_tol"]:
        fragment_tolerance = Tolerance(ppm=(param_data["fragment_tol"]["ppm"][0], param_data["fragment_tol"]["ppm"][1]))
    else:
        fragment_tolerance = Tolerance(ppm=(param_data["fragment_tol"]["da"][0], param_data["fragment_tol"]["da"][1]))

    # Begin search and score
    for spectrum in tqdm(raw_spectrums):
        if spectrum["ms level"] == 1:
            continue
        selected_ion = spectrum["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0][
            "selected ion m/z"]
        charge = spectrum["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["charge state"]
        mz = spectrum["m/z array"].astype(np.float32)
        intensity = spectrum["intensity array"].astype(np.float32)
        precursor = Precursor(
            charge=charge,
            mz=selected_ion,
        )
        raw_spectrum = RawSpectrum(
            file_id=1,
            spec_id=spectrum["id"],
            total_ion_current=spectrum["total ion current"],
            precursors=[precursor],
            mz=mz,
            intensity=intensity
        )

        spec_processor = SpectrumProcessor(take_top_n=75, min_fragment_mz=fragment_min_mz,
                                           max_fragment_mz=fragment_max_mz,
                                           deisotope=deisotope_bool)
        query = spec_processor.process(raw_spectrum)
        scorer = Scorer(report_psms=report_psms, min_matched_peaks=min_matched_peaks,
                        precursor_tolerance=precursor_tolerance,
                        fragment_tolerance=fragment_tolerance, min_isotope_err=min_isotope_err,
                        max_isotope_err=max_isotope_err,
                        min_precursor_charge=min_precursor_charge, max_precursor_charge=max_precursor_charge,
                        min_fragment_mass=fragment_min_mz, max_fragment_mass=fragment_max_mz, chimera=chimera,
                        annotate_matches=annotate_matches, max_fragment_charge=max_fragment_charge)

        results = scorer.score(db=indexed_db, spectrum=query)
        for feature in results:
            peptide = indexed_db[feature.peptide_idx]
            sequence = peptide.sequence
            modifications = peptide.modifications
            unimod_sequence = to_unimod_sequence(modifications, sequence)
            precursor_ppm = abs(feature.expmass - feature.calcmass - feature.isotope_error) * 2e6 / (feature.expmass + feature.calcmass- feature.isotope_error)
            res.append({"file_id": feature.file_id, "spec_id": feature.spec_id, "peptide": unimod_sequence,
                        "proteins": ";".join(peptide.proteins), "rank": feature.rank,
                        "label": feature.label, "exp.mass": feature.expmass, "cal.mass": feature.calcmass,
                        "charge": feature.charge, "retention time": feature.rt,
                        "precursor_ppm": precursor_ppm, "fragment_ppm": feature.average_ppm,
                        "delta mass": feature.delta_mass, "isotope error": feature.isotope_error,
                        "hyperscore": feature.hyperscore, "delta_next": feature.delta_next,
                        "delta_best": feature.delta_best,
                        "matched peaks": feature.matched_peaks, "missed cleavages": feature.missed_cleavages,
                        "scored candidates": feature.scored_candidates,
                        "sage_discriminant_score": feature.discriminant_score,
                        "posterior error": feature.posterior_error,
                        "spectrum_q": feature.spectrum_q, "peptide q": feature.peptide_q,
                        "protein q": feature.protein_q,
                        "ms2 intensity": feature.ms2_intensity})

    res = pd.DataFrame(res, index=None)
    return res


def run_sage(fragment_tolerance: int = 0, precursor_tolerance: int = 0, fragment_type: str = "ppm",
             mzml_files: list = [], sage_config_file: str = None,
             use_file_values: bool = True, n_process=4) -> DataFrame:
    if sage_config_file is None:
        raise ValueError("The sage config file is required.")

    with open(sage_config_file) as f:
        param_data = json.load(f)

        # loading and set default parameters
        param_data["database"]["fragment_min_mz"] = param_data["database"]["fragment_min_mz"] if "fragment_min_mz" in \
                                                                                                 param_data[
                                                                                                     "database"] else 150
        param_data["database"]["fragment_max_mz"] = param_data["database"]["fragment_max_mz"] if "fragment_max_mz" in \
                                                                                                 param_data[
                                                                                                     "database"] else 2000
        param_data["deisotope"] = param_data["deisotope"] if "deisotope" in param_data else True
        param_data["min_matched_peaks"] = param_data["min_matched_peaks"] if "min_matched_peaks" in param_data else 6
        param_data["min_isotope_err"] = param_data["min_isotope_err"] if "min_isotope_err" in param_data else -1
        param_data["max_isotope_err"] = param_data["max_isotope_err"] if "max_isotope_err" in param_data else 3
        param_data["min_precursor_charge"] = param_data[
            "min_precursor_charge"] if "min_precursor_charge" in param_data else 2
        param_data["max_precursor_charge"] = param_data[
            "max_precursor_charge"] if "max_precursor_charge" in param_data else 4
        param_data["chimera"] = param_data["chimera"] if "chimera" in param_data else False
        param_data["report_psms"] = param_data["report_psms"] if "report_psms" in param_data else 1
        param_data["annotate_matches"] = param_data["annotate_matches"] if "annotate_matches" in param_data else False
        param_data["max_fragment_charge"] = param_data[
            "max_fragment_charge"] if "max_fragment_charge" in param_data else 1

    if fragment_tolerance != 0 and precursor_tolerance != 0 or not use_file_values:
        param_data["precursor_tol"]["ppm"] = [int(-1 * precursor_tolerance), int(precursor_tolerance)]
        precursor_tolerance = Tolerance(ppm=(int(-1 * precursor_tolerance), int(precursor_tolerance)))
        if fragment_type == "ppm":
            param_data["fragment_tol"][fragment_type] = [int(-1 * fragment_tolerance), int(fragment_tolerance)]
            fragment_tolerance = Tolerance(ppm=(int(-1 * fragment_tolerance), int(fragment_tolerance)))
        else:
            param_data["fragment_tol"][fragment_type] = [-1 * fragment_tolerance, fragment_tolerance]
            fragment_tolerance = Tolerance(da=(-1 * fragment_tolerance, fragment_tolerance))

    else:
        logging.info("Using the values from the file.")
        if "ppm" in param_data["precursor_tol"]:
            precursor_tolerance = param_data["precursor_tol"]["ppm"][1]
        else:
            precursor_tolerance = param_data["precursor_tol"]["da"][1]

        if "ppm" in param_data["fragment_tol"]:
            fragment_tolerance = param_data["fragment_tol"]["ppm"][1]
        else:
            fragment_tolerance = param_data["fragment_tol"]["da"][1]

    param_data["mzml_paths"] = mzml_files

    temp_sage_file = str(uuid.uuid4()) + ".json"
    with open(temp_sage_file, "w") as f:
        json.dump(param_data, f, indent=4)

    logging.info("Running SAGE with fragment tolerance: {} and precursor tolerance: {}".format(fragment_tolerance,
                                                                                               precursor_tolerance))

    sage_table = pd.DataFrame()
    tasks = []

    with cf.ProcessPoolExecutor(n_process) as pp:
        for m in mzml_files:
            tasks.append(pp.submit(run_one_raw, m, param_data))

        for future in cf.as_completed(tasks):
            sage_table = pd.concat([sage_table, future.result()], ignore_index=True)

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
    file_name = folder_plots_uui + folder_sperator + "q_value_distribution-{}-{}.png".format(fragment_tolerance,
                                                                                             precursor_tolerance)
    plt.savefig(file_name)
    logging.info("Showing Q value distribution plot.")
    plt.close()
    # plt.show()


def compute_entrapment_qvalues(sage_data_frame: pd.DataFrame) -> pd.DataFrame:
    """ compute the q values using the entrapment peptides instead of decoys.
    :param sage_data_frame: panda data frame with the sage results
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
    :param sage_table: pandas data frame with the sage results
    """
    # sage_data_frame['q_diff'] = (sage_data_frame['spectrum_q'] - sage_data_frame['spectrum_entrapment_q']).abs()
    # make copy of the data frame
    current_data_frame = sage_table.copy()
    current_data_frame = current_data_frame[current_data_frame['spectrum_q'] <= 0.01]
    current_data_frame = current_data_frame[current_data_frame['entrapment_qvalue'] <= 0.01]
    filter_decoy = current_data_frame['proteins'].str.contains("DECOY_")
    sage_table_target = current_data_frame[~filter_decoy]
    return len(sage_table_target['peptide'])


def combined_search(results: list = [], start_fragment_tolerance: int = 0, start_precursor_tolerance: int = 0,
                    min_fragment_tolerance: int = 1, max_fragment_tolerance: int = 50,
                    min_precursor_tolerance: int = 10,
                    max_precursor_tolerance: int = 100, num_psms: int = 0, fragment_type: str = "ppm",
                    mzml_files: list = [],
                    grid_frag_steps: int = 10,
                    grid_prec_steps: int = 5, max_iterations=10, search_radius=1, fasta_file="", initial_temp=100,
                    cooling_rate=0.95, sage_config_file=None) -> Tuple[int, int, List[Dict]]:
    def acceptance_probability(old_value, new_value, temperature):
        if new_value > old_value:
            return 1.0
        return np.exp((new_value - old_value) / temperature)

    best_fragment_tolerance = int(start_fragment_tolerance)
    best_precursor_tolerance = int(start_precursor_tolerance)
    best_value = int(num_psms)

    # Coarse Grid Search
    precursor_tolerances = np.linspace(min_precursor_tolerance, max_precursor_tolerance, grid_prec_steps).astype(int)

    grid_best_value = 1
    for ft in precursor_tolerances:
        sage_table = run_sage(fragment_tolerance=start_fragment_tolerance, precursor_tolerance=ft,
                              fragment_type=fragment_type, mzml_files=mzml_files,
                              sage_config_file=sage_config_file, use_file_values=False, n_process=4)
        new_value = compute_best_combination(sage_table)
        results.append(get_stats_from_sage(sage_table, start_fragment_tolerance, ft, new_value))

        if (new_value - grid_best_value) / grid_best_value > 0.01:  # 1% improvement
            best_fragment_tolerance = start_fragment_tolerance
            best_precursor_tolerance = ft
            grid_best_value = new_value
            logging.info(
                "New Best value for precursor tolerance {}: {}".format(best_precursor_tolerance, grid_best_value))
        else:
            logging.info("Current value for precursor tolerance worst than previous {}: {}".format(ft, grid_best_value))
            break

    if best_precursor_tolerance == max_precursor_tolerance:
        logging.info(
            "Best precursor tolerance is the maximum value {}. You have to consider locking for for unknown modifications.".format(
                best_precursor_tolerance))

    fragment_tolerances = np.linspace(start_fragment_tolerance, max_fragment_tolerance, grid_frag_steps).astype(int)

    for ft in fragment_tolerances:
        for pt in precursor_tolerances:
            sage_table = run_sage(fragment_tolerance=ft, precursor_tolerance=pt,
                                  fragment_type=fragment_type, mzml_files=mzml_files,
                                  sage_config_file=sage_config_file, use_file_values=False)
            new_value = compute_best_combination(sage_table)
            results.append(get_stats_from_sage(sage_table, ft, pt, new_value))

            if (new_value - grid_best_value) / grid_best_value > 0.01:  # 1% improvement
                best_fragment_tolerance = ft
                best_precursor_tolerance = pt
                best_value = new_value
                logging.info("New Best value for fragment tolerance {}: {}".format(best_fragment_tolerance, best_value))
            else:
                logging.info(
                    "Current value for fragment tolerance worst than previous {}: {}".format(ft, grid_best_value))
                break

    # Simulated Annealing
    if grid_best_value < best_value:
        best_value = int(num_psms)
        logging.info("Best value from grid search is better than the initial value. Using it.")
        best_fragment_tolerance = start_fragment_tolerance
        best_precursor_tolerance = start_precursor_tolerance

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
                              fragment_type=fragment_type, mzml_files=mzml_files,
                              sage_config_file=sage_config_file,
                              use_file_values=False)
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
@click.option("--initial-fragment-tolerance", help="The initial fragment tolerance to use for the optimization.",
              type=int, default=20)
@click.option("--initial-precursor-tolerance", help="The initial precursor tolerance to use for the optimization.",
              type=int, default=20)
@click.option("--min-fragment-tolerance", default=1, help="The minimum fragment tolerance to consider.", type=int)
@click.option("--max-fragment-tolerance", default=50, help="The minimum precursor tolerance to consider.", type=int)
@click.option("--min-precursor-tolerance", default=10, help="The minimum precursor tolerance to consider.", type=int)
@click.option("--max-precursor-tolerance", default=50, help="The maximum precursor tolerance to consider.", type=int)
@click.option("--fasta-file", default="Homo-sapiens-uniprot-reviewed-contaminants-entrap-decoy-20240615.fasta",
              help="The path to the fasta file to use for the SAGE analysis.")
@click.option("--sage-config-file", default="general-sage.json",
              help="The path to the Sage config file to use for the SAGE analysis.")
@click.option("--max-iterations", default=10, help="The maximum number of iterations to run the optimization.",
              type=int)
def tolerances(fragment_type: str, mzml_path: str, initial_fragment_tolerance: int, initial_precursor_tolerance: int,
               min_fragment_tolerance: int, max_fragment_tolerance: int, min_precursor_tolerance: int,
               max_precursor_tolerance: int,
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

        if os.path.exists(fasta_file):
            logging.info(f"Generating indexed database from {fasta_file}.")
            generate_indexed_db(fasta_file, sage_config_file)
        else:
            logging.error(f"File {fasta_file} does not exist.")
            raise FileNotFoundError(f"File {fasta_file} does not exist.")

        sage_table = run_sage(int(initial_fragment_tolerance), int(initial_precursor_tolerance), fragment_type,
                              mzml_files, sage_config_file=sage_config_file,
                              use_file_values=False, n_process=4)
        sage_table.to_csv("sage_table.csv", index=False)
        num_psms = compute_best_combination(sage_table)
        stats = get_stats_from_sage(sage_table, initial_fragment_tolerance, initial_precursor_tolerance, num_psms)
        results.append(stats)
        initial_fragment_tolerance = math.ceil(stats["fragment_std_error_plus4"])
        best_precursor_tolerance = initial_precursor_tolerance
        best_fragment_tolerance = initial_fragment_tolerance

    best_fragment_tolerance, best_precursor_tolerance, results = combined_search(results=results,
                                                                                 start_fragment_tolerance=best_fragment_tolerance,
                                                                                 start_precursor_tolerance=best_precursor_tolerance,
                                                                                 min_fragment_tolerance=initial_fragment_tolerance,
                                                                                 max_fragment_tolerance=max_fragment_tolerance,
                                                                                 min_precursor_tolerance=min_precursor_tolerance,
                                                                                 max_precursor_tolerance=max_precursor_tolerance,
                                                                                 num_psms=num_psms,
                                                                                 fragment_type=fragment_type,
                                                                                 mzml_files=mzml_files,
                                                                                 max_iterations=max_iterations,
                                                                                 fasta_file=fasta_file,
                                                                                 sage_config_file=sage_config_file)

    print("Best tolerances found: Fragment tolerance: {} - Precursor tolerance: {}".format(best_fragment_tolerance,
                                                                                           best_precursor_tolerance))

    # Write the results to a file
    results_df = pd.DataFrame(results)

    # Aggregate the data to handle duplicate entries by averaging number_psms
    df_agg = results_df.groupby(['fragment_tolerance', 'precursor_tolerance']).agg(
        {'number_psms': 'mean'}).reset_index()

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
@click.option("--fasta-file", default="Homo-sapiens-uniprot-reviewed-contaminants-entrap-decoy-20240615.fasta",
              help="The path to the fasta file to use for the SAGE analysis.")
@click.option("--sage-config-file", default="general-sage-ptms.json",
              help="The path to the Sage config file to use for the SAGE analysis.", required=True)
def ptms(mzml_path: str, fasta_file: str, sage_config_file: str):
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

    sage_table = run_sage(fragment_type="ppm", mzml_files=mzml_files, fasta_path=fasta_file,
                          sage_config_file=sage_config_file, use_file_values=True)

    num_psms = compute_best_combination(sage_table)
    stats = get_stats_from_sage(sage_table, fragment_tolerance=fragment_tolerance,
                                precursor_tolerance=precursor_tolerance, number_psms=num_psms)

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
