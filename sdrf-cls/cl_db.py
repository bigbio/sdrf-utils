import re
from typing import Union
import click
import pandas
import pandas as pd
import glob
import numpy as np
import spacy
from sklearn.metrics.pairwise import cosine_similarity


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

nlp = spacy.load("en_core_web_md")  # Load the spacy model

def string_if_not_empty(param: list) -> Union[None, str]:
    """
    Return a string if the list is not empty
    :param param: List
    :return: None if the list is empty, the string otherwise
    """
    if param == "None":
        param = []
    if param and len(param) > 0:
        l = [
            x
            for x in param
            if isinstance(x, float)
            and ~np.isnan(x)
            or not isinstance(x, float)
            and x != None
        ]
        return "; ".join(l)
    return "no available"

def get_cell_line_code(sdrf_file):
    """
    Extracts a list of unique cell line codes from an SDRF file.

    Reads the specified SDRF file and attempt to retrieve unique values
    from the 'characteristics[cell line]' column. If the column is not
    found, an error message is printed.

    Parameters:
        sdrf_file (str): The path to the SDRF file to be processed.

    Returns:
        list: A list of unique cell line codes, or None if the column is not found.
    """
    sdrf = pd.read_csv(sdrf_file, sep="\t")
    try:
        cl_list = sdrf["characteristics[cell line]"].unique().tolist()
        return cl_list
    except KeyError:
        print("The SDRF file does not have a column named 'characteristics[cell line]' -- {}".format(sdrf_file))
    return None

def modo_dict_to_context(obo_list: list) -> str:
    context = ""
    for entry in obo_list:
        if "obsolete" not in entry["name"]:
            context += f"{entry['id']}: {entry['name']}\n"
    return context


def calculate_similarity(cellosaurus_text, synonyms):
    query_vec = nlp(" ".join(cellosaurus_text)).vector
    synonyms_vec = [nlp(" ".join(synonym)).vector for synonym in synonyms]
    similarities = cosine_similarity([query_vec], synonyms_vec)
    return max(similarities[0])


def map_celllines(cellosaurus_text: str, context: list):
    # Use the LLM to find the correct MONDO term
    max_similarity = 0
    closest_match = None
    for entry in context:
        synonyms = entry.split(";")
        similarity = calculate_similarity(cellosaurus_text, synonyms)
        if similarity > max_similarity:
            closest_match = entry
            max_similarity = similarity
    return closest_match


def read_cell_line_database(database) -> dict:
    """
    The database is a tab-delimited with the following structure. The information for each cell lines is:

    cell line: Selected cell line name (e.g., A549, HELA, etc.)
    cellosaurus name: Cellosaurus name
    cellosaurus accession: Cellosaurus accession
    bto cell line: BTO cell line
    organism: Organism of the cell line
    organism part: Organism part of the cell line
    sampling site:	Sampling site of the cell line
    age: Age of the cell line
    developmental stage: Developmental stage of the cell line
    sex: Sex of the cell line
    ancestry category: Ancestry category of the cell line
    disease: Disease associated with the cell line
    cell type: Cell type of the cell line
    Material type: Material used to grow the cell line
    synonyms: Synonyms for the cell line
    curated: The cell line has been curated or not. Possible values (curated, not curated, ai curated)

    If multiple values are present for a give field; they are separated by ;

    :param database: Database file path
    :return: List of dictionaries with the database content
    """

    database_df = pd.read_csv(database, sep="\t", comment="#", header=0, dtype=str)

    # keep only curated cell lines
    database_df = database_df[database_df["curated"] == "curated"]

    # Convert the dataframe to a list of dictionaries
    database_list = database_df.to_dict(orient="records")

    # convert disease and sampling site to list divided by ;
    for entry in database_list:
        entry["disease"] = entry["disease"].split(";")
        entry["disease"] = [disease.strip() for disease in entry["disease"]]
        entry["sampling site"] = entry["sampling site"].split(";")
        entry["sampling site"] = [site.strip() for site in entry["sampling site"]]
        entry["synonyms"] = entry["synonyms"].split(";")
        # remove spaces in the synonyms
        entry["synonyms"] = [synonym.strip() for synonym in entry["synonyms"]]
    database_list = {entry["cell line"]: entry for entry in database_list}

    return database_list


def find_cell_line(old_cl: str, current_cl_database: dict) -> Union[dict, None]:
    """
    Find a given cell line annotated in an SDRF in a standarized cell line database
    :param old_cl: Code (e.g., HELA, A549, etc.)
    :param current_cl_database: Database of all cell lines for the multiomics configuration
    :return:
    """

    # Normalize the cell line name to lower case and remove spaces
    old_cl = old_cl.lower().strip()

    for key, entry in current_cl_database.items():
        if "cell line" in entry:
            if entry["cell line"].lower().strip() == old_cl:
                return entry
        if "cellosaurus name" in entry:
            if entry["cellosaurus name"].lower().strip() == old_cl:
                return entry
        if "bto cell line" in entry:
            if entry["bto cell line"].lower().strip() == old_cl:
                return entry
        for synonym in entry["synonyms"]:
            if synonym.lower().strip() == old_cl:
                return entry
    return None


def is_in_synonyms(old_cl, cellosaurus) -> bool:
    """
    Check if a cell line is in the synonyms list
    :param old_cl: Old cell line code
    :param sysnonyms: Synonyms list
    :return:
    """
    if "synonyms" not in cellosaurus:
        return False
    for synonym in cellosaurus["synonyms"]:
        if synonym.lower().strip() == old_cl.lower().strip():
            return True
    return False


def get_cell_line_bto(bto_code: str, bto_list: list):
    for entry in bto_list:
        if entry["id"] == bto_code:
            return entry
    return None

def validate_ages_as_sdrf(age_string: str) -> bool:
    """
    Validate the age string from the SDRF. The age should be in multiple format:
     - Year format: 1Y, 10Y, 100Y, etc.
     - Year and Month: 40Y5M, 10Y10M, etc.
     - Year, Month, Day: 10Y10M10D, 100Y1M3D, etc.
     - Weeks: 8W, etc
    All the ages could also include intervals like 10Y-20Y, 10Y-20Y5M, etc.
    @param age_string: Age string
    @return: True if the age is valid, False otherwise
    """
    # Regular expression to match the age format
    pattern = r"(\d+Y)?(\d+M)?(\d+D)?(\d+W)?(-(\d+Y)?(\d+M)?(\d+D)?(\d+W)?)?"
    match = re.match(pattern, age_string)
    if match:
        return True
    print(f"Age {age_string} is not valid")

    return False


def get_age_consensus(cell_passport_entry, cellosaurus_entry, ae_entry):
    """
    The Age in SDRF could be in multiple formats, we will use the following rules to get the age:
    Year format: 1Y, 10Y, 100Y, etc.
    Year and Month: 40Y5M, 10Y10M, etc.
    Year, Month, Day: 10Y10M10D, 100Y1M3D, 0Y9M etc.
    Weeks: 8W, etc
    All the ages could also include intervals like 10Y-20Y, 10Y-20Y5M, etc.

    """
    if (
        cell_passport_entry is not None
        and "age" in cell_passport_entry
        and cell_passport_entry["age"] != "no available"
    ) and int(cell_passport_entry["age"]) > 0:
        return str(cell_passport_entry["age"]) + "Y"
    if (
        cellosaurus_entry is not None
        and "age" in cellosaurus_entry
        and cellosaurus_entry["age"] != "no available"
    ):
        return cellosaurus_entry["age"]
    if (
        ae_entry is not None
        and "age" in ae_entry
        and ae_entry["age"] != "no available"
        and ae_entry["age"] != "nan"
        and "available" not in ae_entry["age"]
    ):
        age = ae_entry["age"].upper().replace("YEAR", "").strip()
        return str(age) + "Y"
    return "no available"


def estimate_developmental_stage(age_string: str) -> str:
    """
    Estimate the developmental stage from the age string
    """
    # remove Y from age and check if is integer
    age = age_string.replace("Y", "")
    if age.isdigit():
        age = int(age)
        if 1 <= age <= 2:
            return "Infant"
        elif 3 <= age < 12:
            return "Children"
        elif 12 <= age < 18:
            return "Juvenile"
        elif 18 <= age < 65:
            return "Adult"
        elif age >= 65:
            return "Elderly"
    return "no available"


def create_new_entry(
    cellosaurus_entry, cell_passport_entry, ae_entry
) -> Union[dict, None]:
    """
    The entry is a dictionary with the following fields:
    cell line
    cellosaurus name
    cellosaurus accession
    bto cell line
    organism
    organism part
    sampling site
    age
    developmental stage
    sex
    ancestry category
    disease
    cell type
    Material type
    synonyms
    curated
    """
    # Create a new entry
    entry = {
        "cell line": "no available",
        "cellosaurus name": "no available",
        "cellosaurus accession": "no available",
        "bto cell line": "no available",
        "organism": "no available",
        "organism part": "no available",
        "sampling site": ["no available", "no available"],
        "age": "no available",
        "developmental stage": "no available",
        "sex": "no available",
        "ancestry category": "no available",
        "disease": ["no available", "no available"],
        "cell type": "no available",
        "Material type": "cell",
        "synonyms": [],
        "curated": "not curated",
    }
    original = entry.copy()
    # The cell passport is the reference for the database, we will use it for the cell line name.
    if cell_passport_entry is not None:
        entry["cell line"] = cell_passport_entry["cell line"]
    elif cellosaurus_entry is not None:
        entry["cell line"] = cellosaurus_entry["cellosaurus name"]

    if cellosaurus_entry is not None:
        entry["cellosaurus name"] = cellosaurus_entry["cellosaurus name"]
        entry["cellosaurus accession"] = cellosaurus_entry["cellosaurus accession"]
        if "bto cell line" in cellosaurus_entry:
            entry["bto cell line"] = cellosaurus_entry["bto cell line"]

    # Set the organism using the cellosaurus entry and cell passport entry
    if cellosaurus_entry is not None and cell_passport_entry is not None:
        if (
            cellosaurus_entry["organism"].lower()
            != cell_passport_entry["organism"].lower()
            and cellosaurus_entry["organism"] != "no available"
            and cell_passport_entry["organism"] != "no available"
        ):
            raise ValueError(
                f"Organism mismatch: {cellosaurus_entry['organism']} vs {cell_passport_entry['organism']}"
            )
        else:
            entry["organism"] = cell_passport_entry["organism"].capitalize()
    elif (
        cell_passport_entry is not None
        and cell_passport_entry["organism"].lower() != "no available"
    ):
        entry["organism"] = cell_passport_entry["organism"].capitalize()
    elif (
        cellosaurus_entry is not None
        and cellosaurus_entry["organism"].lower() != "no available"
    ):
        entry["organism"] = cellosaurus_entry["organism"].capitalize()
    else:
        entry["organism"] = "no available"

    # Set the sampling site using the cell passport entry, cell

    if (
        cell_passport_entry is not None
        and cell_passport_entry["sampling site"].lower() != "no available"
        and cell_passport_entry["sampling site"].lower() != "unknown"
    ):
        entry["sampling site"][0] = cell_passport_entry["sampling site"].strip().capitalize()
    if (
        cellosaurus_entry is not None
        and cellosaurus_entry["sampling site"].lower() != "no available"
        and cellosaurus_entry["sampling site"].lower() != "unknown"
    ):
        entry["sampling site"][1] = cellosaurus_entry["sampling site"].strip().capitalize()

    if (
        cell_passport_entry is not None
        and cell_passport_entry["disease"].lower() != "no available"
    ):
        entry["disease"][0] = cell_passport_entry["disease"].strip().capitalize()
    if (
        cellosaurus_entry is not None
        and cellosaurus_entry["disease"].lower() != "no available"
    ):
        entry["disease"][1] = cellosaurus_entry["disease"].strip().capitalize()

    # Set organism part using the cell passport entry
    if (
        cell_passport_entry is not None
        and cell_passport_entry["organism part"].lower() != "no available"
    ):
        entry["organism part"] = cell_passport_entry["organism part"]
    elif ae_entry is not None and ae_entry["organism part"].lower() != "no available":
        entry["organism part"] = ae_entry["organism part"].strip().capitalize()

    # Set the age using cell passports, cellosaurus, and ae entries
    entry["age"] = get_age_consensus(cell_passport_entry, cellosaurus_entry, ae_entry)

    # Set sex using the cell passport entry
    if (
        cell_passport_entry is not None
        and cell_passport_entry["sex"].lower() != "no available"
    ):
        entry["sex"] = cell_passport_entry["sex"].capitalize()
    elif (
        cellosaurus_entry is not None
        and cellosaurus_entry["sex"].lower() != "no available"
    ):
        entry["sex"] = cellosaurus_entry["sex"].capitalize()
    elif ae_entry is not None and "available" not in ae_entry["sex"].lower():
        entry["sex"] = ae_entry["sex"].capitalize()

    # Set ancestry category using the cell passport entry
    if (
        cellosaurus_entry is not None
        and cellosaurus_entry["ancestry category"].lower() != "no available"
    ):
        entry["ancestry category"] = cellosaurus_entry["ancestry category"]
    elif (
        cell_passport_entry is not None
        and cell_passport_entry["ancestry category"].lower() != "no available"
    ):
        entry["ancestry category"] = cell_passport_entry["ancestry category"]

    # Set cell type using the cellosaurus entry
    if (
        cellosaurus_entry is not None
        and "available" not in cellosaurus_entry["cell type"].lower()
    ):
        entry["cell type"] = cellosaurus_entry["cell type"]

    # Synonyms are the union of the cell passport and cellosaurus synonyms and each one of them
    # should be unique
    if (
        cell_passport_entry is not None
        and "available" not in cell_passport_entry["synonyms"]
    ):
        entry["synonyms"] += [
            cell_line.upper().strip()
            for cell_line in cell_passport_entry["synonyms"].split(";")
        ]
    if cellosaurus_entry is not None and "available" not in cellosaurus_entry:
        entry["synonyms"] += [
            cell_line.upper().strip()
            for cell_line in cellosaurus_entry["synonyms"].split(";")
        ]
    if ae_entry is not None and "available" not in ae_entry["synonyms"]:
        entry["synonyms"] += [
            cell_line.upper().strip() for cell_line in ae_entry["synonyms"].split(";")
        ]
    # Remove duplicates in the synonym list
    entry["synonyms"] = list(set(entry["synonyms"]))

    # development stage
    if (
        cellosaurus_entry is not None
        and "available" not in cellosaurus_entry["developmental stage"].lower()
    ):
        entry["developmental stage"] = cellosaurus_entry["developmental stage"].capitalize()
    elif entry["age"] != "no available":
        entry["developmental stage"] = estimate_developmental_stage(entry["age"])

    if entry == original or entry["organism"] == "no available":
        return None

    entry["curated"] = "not curated"
    return entry

def write_database(current_cl_database: list, database: str) -> None:
    """
    Write the database objects to the database file
    :param current_cl_database: current cell line database list
    :param database: database file path
    :return:
    """
    def get_string_available(list_values: list)-> str:

        # split some of the words in the list by , and add the list to the values
        list_values = [value.split(",") for value in list_values]
        list_values = [item for sublist in list_values for item in sublist]
        # remove duplicates
        list_values = list(set(list_values))
        # remove the no available values from list
        list_values = [value.capitalize() for value in list_values if value != "no available"]
        if not list_values:
            return "no available"
        return "; ".join(list_values)

    with open(database, "w") as file:
        headers = [
            "cell line",
            "cellosaurus name",
            "cellosaurus accession",
            "bto cell line",
            "organism",
            "organism part",
            "sampling site",
            "age",
            "developmental stage",
            "sex",
            "ancestry category",
            "disease",
            "cell type",
            "Material type",
            "synonyms",
            "curated",
        ]
        # Write the header row
        file.write("\t".join(headers) + "\n")

        for key, entry in current_cl_database.items():
            row = [
                entry.get("cell line", "no available"),
                entry.get("cellosaurus name", "no available"),
                entry.get("cellosaurus accession", "no available"),
                entry.get("bto cell line", "no available"),
                entry.get("organism", "no available"),
                entry.get("organism part", "no available"),
                get_string_available(entry.get("sampling site", ["no available", "no available"])),
                entry.get("age", "no available"),
                entry.get("developmental stage", "no available"),
                entry.get("sex", "no available"),
                entry.get("ancestry category", "no available"),
                get_string_available(entry.get("disease", ["no available", "no available"])),
                entry.get("cell type", "no available"),
                entry.get("Material type", "no available"),
                string_if_not_empty(entry.get("synonyms", [])),
                entry.get("curated", "not curated"),
            ]

            row = ["no available" if item is None else str(item) for item in row]
            file.write("\t".join(row) + "\n")

@click.command(
    "cl-database",
    short_help="Create a cell lines metadata database for annotating cell lines SDRFs",
)
@click.option(
    "--database",
    help="Current database file with cell lines",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "--cellosaurus-database",
    help="Cellosaurus database file",
    required=True,
    type=click.Path(exists=True),
    default="cellosaurus-db.tsv",
)
@click.option(
    "--ea-database",
    help="EA Atlas database file",
    required=True,
    type=click.Path(exists=True),
    default="ea-cls-db.tsv",
)
@click.option(
    "--cell-passports-database",
    help="Cell passports database file",
    required=True,
    type=click.Path(exists=True),
    default="cell-passports-db.tsv",
)
@click.option(
    "--sdrf-path",
    help="SDRF folder with all existing SDRF files",
    required=False,
    type=click.Path(exists=True),
)
@click.option(
    "--include-all-cellpassports",
    help="Include all cell passports cell lines",
    is_flag=True,
)
@click.option(
    "--ai-synonyms",
    help="AI synonyms file",
    required=True,
    type=click.Path(exists=True),
    default="ai-synonyms.tsv",
)
@click.option(
    "--unknown", help="Output for unknown cell lines in cellosaurus", required=True
)
def cl_database(
    database: str,
    cellosaurus_database: str,
    ea_database: str,
    cell_passports_database: str,
    sdrf_path: str,
    include_all_cellpassports: bool,
    ai_synonyms: str,
    unknown: str,
) -> None:
    """
    Creates and updates a cell lines metadata database for annotating cell lines in SDRFs.

    This command-line tool processes various database files and SDRF files to compile
    a comprehensive cell line database. It checks for existing cell lines, adds new ones,
    and writes the updated database to a file. It also logs any cell lines that could not
    be found in the provided databases.

    Parameters:
        database (str): Path to the current cell line database file.
        cellosaurus_database (str): Path to the CelloSaurus database file.
        ea_database (str): Path to the EA Atlas database file.
        cell_passports_database (str): Path to the cell passports database file.
        sdrf_path (str): Path to the folder containing SDRF files.
        include_all_cellpassports (bool): Flag to include all cell passports cell lines.
        ai_synonyms (str): Path to the AI synonyms file.
        unknown (str): Path to the output file for unknown cell lines.

    Returns:
        None
    """

    cls = []  # List of cell lines
    if sdrf_path is None and not include_all_cellpassports:
        raise ValueError(
            "The cell lines that wants to be added search from existing SDRF must be provided"
        )
    if sdrf_path is not None:
        sdrf_files = glob.glob(sdrf_path + "/**/*.tsv", recursive=True)
        for sdrf_file in sdrf_files:
            cls += get_cell_line_code(sdrf_file)
        print("Number of cell lines in the SDRF files: ", len(cls))

    # Read the current cell line database
    current_cl_database = read_cell_line_database(database)

    # Parse the CelloSaurus file and transform the list to dictionary of cellosaurus where key is cellosaurus name
    cellosaurus = pandas.read_csv(cellosaurus_database, sep="\t", header=0, dtype=str)
    cellosaurus = cellosaurus.to_dict(orient="records")
    cellosaurus = [{k: str(v) for k, v in record.items()} for record in cellosaurus]
    cellosaurus = {entry["cellosaurus name"]: entry for entry in cellosaurus}

    # Parse the EA Atlas file and transform list to dictionary of ea atlas where key is cell line
    ea_atlas = pandas.read_csv(ea_database, sep="\t", header=0, dtype=str)
    ea_atlas = ea_atlas.to_dict(orient="records")
    ea_atlas = [{k: str(v) for k, v in record.items()} for record in ea_atlas]
    ea_atlas = {entry["cell line"]: entry for entry in ea_atlas}

    # Parse the cell passports file and transform list to dictionary of cell passports where key is cell line
    cell_passports = pandas.read_csv(
        cell_passports_database, sep="\t", header=0, dtype=str
    )
    cell_passports = cell_passports.to_dict(orient="records")
    cell_passports = [
        {k: str(v) for k, v in record.items()} for record in cell_passports
    ]
    cell_passports = {entry["cell line"]: entry for entry in cell_passports}

    if include_all_cellpassports:
        # get cell lines names from cell pass
        cls += [value["cell line"] for key, value in cell_passports.items()]

    cls = list(set(cls))
    print("Final number of cell lines to annotated -- {}".format(str(len(cls))))

    # Add the cell lines that are not in the current cell line database
    non_found_cl = []

    ai_synonyms_dic = None
    if ai_synonyms is not None:
        ai_synonyms_dic = pandas.read_csv(ai_synonyms, sep="\t", header=0, dtype=str)
        ai_synonyms_dic = ai_synonyms_dic.to_dict(orient="records")
        ai_synonyms_dic = [{k: str(v) for k, v in record.items()} for record in ai_synonyms_dic]
        ai_synonyms_dic = {entry["cell line"]: entry for entry in ai_synonyms_dic}

    def find_cell_line_cellosaurus(cl: str, cellosaurus: dict) -> Union[dict, None]:
        for key, cellosaurus_entry in cellosaurus.items():
            if cellosaurus_entry["cellosaurus name"].lower() == cl.lower():
                return cellosaurus_entry
            if "synonyms" in cellosaurus_entry:
                for synonym in cellosaurus_entry["synonyms"]:
                    if cl.lower() in synonym.lower():
                        return cellosaurus_entry
            if "cellosaurus accession" in cellosaurus_entry:
                if cellosaurus_entry["cellosaurus accession"].lower() == cl.lower():
                    return cellosaurus_entry
        return None

    def find_cell_line_cell_passports(
        cl: str, cell_passports: dict
    ) -> Union[dict, None]:
        for key, cell_passports_entry in cell_passports.items():
            if cell_passports_entry["cell line"].lower() == cl.lower():
                return cell_passports_entry
            if "synonyms" in cell_passports_entry:
                for synonym in cell_passports_entry["synonyms"]:
                    if cl.lower() in synonym.lower():
                        return cell_passports_entry
            if "cellosaurus accession" in cell_passports_entry:
                if cell_passports_entry["cellosaurus accession"].lower() == cl.lower():
                    return cell_passports_entry
        return None

    def find_cell_line_ea_atlas(cl: str, ea_atlas: dict) -> Union[dict, None]:
        """
        Searches for a cell line in the EA Atlas database.

        Iterates through the EA Atlas entries to find a match for the given cell line
        name or its synonyms. Returns the corresponding entry if found, otherwise
        returns None.

        Parameters:
            cl (str): The cell line name to search for.
            ea_atlas (dict): The EA Atlas database represented as a dictionary.

        Returns:
            dict or None: The EA Atlas entry for the cell line if found, otherwise None.
        """
        for key, ea_atlas_entry in ea_atlas.items():
            if ea_atlas_entry["cell line"].lower() == cl.lower():
                return ea_atlas_entry
            if "synonyms" in ea_atlas_entry:
                for synonym in ea_atlas_entry["synonyms"]:
                    if cl.lower() in synonym.lower():
                        return ea_atlas_entry
        return None

    def find_in_synonyms_table(cl: str, ai_synonyms: dict) -> Union[str, None]:
        """
        Searches for a cell line in the ai_synonyms dictionary that matches the given
        cell line identifier `cl`. If an exact match is found, it returns the cell line.
        If not, it checks the synonyms for a match and returns the corresponding cell
        line if found. If no match is found, it returns the original `cl`.

        Parameters:
            cl (str): The cell line identifier to search for.
            ai_synonyms (dict): A dictionary containing cell line entries with possible synonyms.

        Returns:
            str: The matched cell line from the dictionary or the original `cl` if no match is found.
        """
        for key, ai_synonyms_entry in ai_synonyms.items():
            if ai_synonyms_entry["cell line"].lower() == cl.lower():
                return ai_synonyms_entry["cell line"]
            if "synonyms" in ai_synonyms_entry:
                for synonym in ai_synonyms_entry["synonyms"].split(";"):
                    if cl.lower() in synonym.lower():
                        return ai_synonyms_entry["cell line"]
        return cl

    for cl in cls:

        if find_cell_line(cl, current_cl_database) is None:

            if ai_synonyms_dic is not None:
                cl = find_in_synonyms_table(cl, ai_synonyms_dic)

            cellosaurus_entry = find_cell_line_cellosaurus(cl, cellosaurus)
            cell_passports_entry = find_cell_line_cell_passports(cl, cell_passports)
            ea_atlas_entry = find_cell_line_ea_atlas(cl, ea_atlas)

            if cell_passports_entry is not None:
                for key, value in cellosaurus.items():
                    # override the cellosaurus entry with the cell passports entry link to cellosaurus
                    if (
                        value["cellosaurus accession"].lower()
                        == cell_passports_entry["cellosaurus accession"].lower()
                    ):
                        cellosaurus_entry = value
                        break

            if cellosaurus_entry is not None:
                for key, value in cell_passports.items():
                    # if the cell passports are not found using the cell line name try to find it througth cellosaurus
                    if (
                        value["cellosaurus accession"].lower()
                        == cellosaurus_entry["cellosaurus accession"].lower()
                    ):
                        cell_passports_entry = value
                        break

            if (
                cellosaurus_entry is None
                and cell_passports_entry is None
                and ea_atlas_entry is None
            ):
                non_found_cl.append(cl)
            else:
                new_cl_entry = create_new_entry(
                    cellosaurus_entry, cell_passports_entry, ea_atlas_entry
                )
                if new_cl_entry is not None:
                    if new_cl_entry["cell line"] not in current_cl_database:
                        current_cl_database[new_cl_entry["cell line"]] = new_cl_entry
                    else:
                        print(f"Cell line {cl} already in the database")
                else:
                    non_found_cl.append(cl)
        else:
            print(f"Cell line {cl} already in the database")

    write_database(current_cl_database, database)
    with open(unknown, "w") as file:
        for cl in non_found_cl:
            file.write(cl + "\n")

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    Main function to run the CLI
    """
    pass

cli.add_command(cl_database)

if __name__ == "__main__":
    cli()
