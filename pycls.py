import gzip
import os
from typing import Union

import click
import pandas
import pandas as pd
import glob

import spacy
import numpy as np
from huggingface_hub import login
from owlready2 import *
from sklearn.metrics.pairwise import cosine_similarity
from transformers import pipeline

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

nlp = spacy.load("en_core_web_md")  # Load the spacy model


def get_cell_line_code(sdrf_file):
    sdrf = pd.read_csv(sdrf_file, sep="\t")
    try:
        cl_list = sdrf["characteristics[cell line]"].unique().tolist()
    except KeyError:
        print("The SDRF file does not have a column named 'characteristics[cell line]' -- {}".format(sdrf_file))
    return cl_list


def get_sampling_site(cellosaurus_comment: str) -> Union[None, str]:
    """
    Extract the sampling site from the cellosaurus comment field
    :param cellosaurus_comment: Cellosaurus comment field
    :return: Sampling site
    """
    # The regular expression pattern to match "In situ", "Colon", and the UBERON accession
    pattern = r"Derived from site:\s*(.+?);\s*(.+?);\s*UBERON=(UBERON_\d{7})"

    # Search for the pattern in the text
    match = re.search(pattern, cellosaurus_comment)

    if match:
        tissue = match.group(2)
        return tissue
    return None


def get_cell_type(cellosaurus_comment: str) -> Union[None, str]:
    """
    Extract the cell type from the cellosaurus comment field
    :param cellosaurus_comment: Cellosaurus comment field
    :return: Cell type
    """
    # The regular expression pattern to match the cell type and the CL accession
    pattern = r"Cell type:\s*(.+?);\s*CL=(CL_\d{7})"

    # Search for the pattern in the text
    match = re.search(pattern, cellosaurus_comment)

    # Check if a match is found and print the results
    if match:
        cell_type = match.group(2)
        cell_type = cell_type.replace("_", ":").strip()
        return cell_type
    return None


def parse_cellosaurus_file(file_path, bto: dict, cl_type: dict):
    """
    Parse the CelloSaurus file and return a list of dictionaries with the parsed content
    :param file_path: CelloSaurus file path
    :param bto: BTO ontology list
    :return: List of CelloSaurus dictionaries
    """

    def parse_entry(entry, bto: dict, cl_type: dict):
        # init the data dictionary also add None to all the fields
        data = {
            "cellosaurus name": "no available",
            "cellosaurus accession": "no available",
            "bto cell line": "no available",
            "efo": "no available",
            "organism": "no available",
            "age": "no available",
            "developmental stage": "no available",
            "sex": "no available",
            "ancestry category": "no available",
            "disease": "no available",
            "cell type": "no available",
            "sampling site": "no available",
            "synonyms": [],
        }

        lines = entry.strip().split("\n")
        for line in lines:
            if line.startswith("ID"):
                data["cellosaurus name"] = line.split("ID ")[1].strip()
            elif line.startswith("AC"):
                data["cellosaurus accession"] = line.split("AC ")[1].strip()
            elif line.startswith("SY"):
                data["synonyms"] = line.split("SY ")[1].strip().split("; ")
            elif line.startswith("DR   BTO"):
                bto_accession = line.split("; ")[1]
                if bto_accession in bto:
                    data["bto cell line"] = bto[bto_accession]["name"]
                    if "synonyms" in bto[bto_accession]:
                        if "synonyms" in data:
                            data["synonyms"] += bto[bto_accession]["synonyms"]
                        else:
                            data["synonyms"] = bto[bto_accession]["synonyms"]
            elif line.startswith("DR   EFO"):
                data["efo"] = line.split("; ")[1]
            elif line.startswith("OX"):
                data["organism"] = line.split("OX ")[1].strip()
                scientific_name, tax = parse_cellosaurus_taxonomy(data["organism"])
                data["organism"] = scientific_name
            elif line.startswith("SX"):
                data["sex"] = line.split()[1]
            elif line.startswith("AG"):
                data["age"] = line.split("AG ")[1].strip()
            elif line.startswith("CC") and "Population" in line:
                data["ancestry category"] = line.split(": ")[1].strip().replace(".", "")
            elif line.startswith("DI"):
                if "NCIt" in line:
                    # Regular expression to match all disease annotations and capture the disease name
                    pattern = r"NCIt;\s*C\d+;\s*([^;]+)"
                    match = re.search(pattern, line)
                    if match:
                        data["disease"] = match.group(1).strip()
            elif line.startswith("CC") and "Derived from site" in line:
                data["sampling site"] = get_sampling_site(line)
            elif line.startswith("CC") and "Cell type" in line:
                code_cl = get_cell_type(line)
                if code_cl is not None and code_cl in cl_type:
                    data["cell type"] = cl_type[code_cl]["name"]

        return data

    # Read the file and split into entries, the file is gzipped
    with gzip.open(file_path, "r") as file:
        content = file.read().decode("utf-8")

    # Split the content by entries
    entries = content.split("//\n")

    # Parse each entry
    parsed_data = [
        parse_entry(entry, bto, cl_type) for entry in entries if entry.strip()
    ]
    # remove empty entries
    parsed_data = [entry for entry in parsed_data if entry]
    return parsed_data


def read_obo_file(file_path) -> dict:
    """
    Read an obo file and return a list of dictionaries with the parsed content
    :param file_path: OBO file path
    :return: List of OBO dictionaries
    """
    with open(file_path, "r") as file:
        content = file.read()

    # Split the content by entries
    entries = content.split("\n\n")

    def parse_obo_term(entry: str) -> dict:
        obo_dict = {}
        lines = entry.strip().split("\n")
        for line in lines:
            if line.startswith("id:"):
                # extract the id from pattern id: BTO:0000001 using regex
                obo_dict["id"] = line.split("id: ")[1].strip()
            elif line.startswith("name:"):
                obo_dict["name"] = line.split("name: ")[1].strip()
            elif line.startswith("def:"):
                obo_dict["def"] = line.split("def: ")[1].strip()
            elif line.startswith("is_a:"):
                if "is_a" not in obo_dict:
                    obo_dict["is_a"] = []
                obo_dict["is_a"].append(line.split("is_a: ")[1].strip())
            elif line.startswith("synonym:"):
                if "synonyms" not in obo_dict:
                    obo_dict["synonyms"] = []
                obo_dict["synonyms"].append(line.split("synonym: ")[1].strip())
                obo_dict["synonyms"] = [
                    synonym.replace("RELATED []", "")
                    .replace("RELATED MS []", "")
                    .strip()
                    .strip('"')
                    for synonym in obo_dict["synonyms"]
                ]
        if "synonyms" not in obo_dict:
            obo_dict["synonyms"] = []
        return obo_dict

    # Create dictionary of OBO terms, the key is the id of the obo term
    obo_dict = {
        parse_obo_term(entry)["id"]: parse_obo_term(entry)
        for entry in entries
        if entry.strip() and "id:" in entry
    }

    return obo_dict


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


def parse_cellosaurus_taxonomy(organism_text: str):
    """
    Parse the organism text from the cellosaurus database
    :param organism_text: organism text from the cellosaurus database
    :return: return the species name and the taxonomy id
    """
    pattern = r"NCBI_TaxID=(\d+); ! ([\w\s]+) \(([\w\s]+)\)"
    match = re.search(pattern, organism_text)

    if match:
        ncbi_tax_id = match.group(1)
        species_name = match.group(2)
        return species_name, ncbi_tax_id
    return None, None


def is_age_in_text(age_text: str) -> bool:
    """
    Check if the age field contains a number is an Age if not is a developmental stage.
    :param age_text: String text
    :return: True if the age is a number, False otherwise
    """
    return any(char.isdigit() for char in age_text)


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


def create_new_entry_from_cellosaurus(cellosaurus):
    """
    Create a new entry for a cell line not found in the database
    :param old_cl:
    :param bto:
    :param cellosaurus_list:
    :param modo_context:
    :return:
    """
    entry = {}
    entry["cellosaurus name"] = cellosaurus["cellosaurus name"]
    if "bto cell line" in cellosaurus:
        entry["bto cell line"] = cellosaurus["bto cell line"]
    else:
        entry["bto cell line"] = None
    if "cellosaurus accession" in cellosaurus:
        entry["cellosaurus accession"] = cellosaurus["cellosaurus accession"]
    if "organism" in cellosaurus:
        entry["organism"] = cellosaurus["organism"]
    entry["organism part"] = None
    if "age" in cellosaurus:
        if is_age_in_text(cellosaurus["age"]):
            entry["age"] = cellosaurus["age"]
            entry["developmental stage"] = None
        else:
            entry["developmental stage"] = cellosaurus["age"]
            entry["age"] = None
    if "sex" in cellosaurus:
        entry["sex"] = cellosaurus["sex"]
    else:
        entry["sex"] = None
    if "ancestry category" in cellosaurus:
        entry["ancestry category"] = cellosaurus["ancestry category"]
    else:
        entry["ancestry category"] = None
    if "disease" in cellosaurus:
        entry["disease"] = cellosaurus["disease"]
    else:
        entry["disease"] = None
    if "synonyms" in cellosaurus:
        entry["synonyms"] = cellosaurus["synonyms"]
    else:
        entry["synonyms"] = []
    if "cell type" in cellosaurus:
        entry["cell type"] = cellosaurus["cell type"]
    if "sampling site" in cellosaurus:
        entry["sampling site"] = cellosaurus["sampling site"]
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


def write_database_cellosaurus(current_cl_database: list, database: str) -> None:
    """
    Write the database objects to the database file
    :param current_cl_database: current cell line database list
    :param database: database file path
    :return:
    """

    with open(database, "w") as file:
        headers = [
            "cellosaurus name",
            "cellosaurus accession",
            "bto cell line",
            "organism",
            "age",
            "developmental stage",
            "sex",
            "ancestry category",
            "disease",
            "cell type",
            "sampling site",
            "synonyms",
        ]
        # Write the header row
        file.write("\t".join(headers) + "\n")

        for entry in current_cl_database:
            row = [
                entry.get("cellosaurus name", "no available"),
                entry.get("cellosaurus accession", "no available"),
                entry.get("bto cell line", "no available"),
                entry.get("organism", "no available"),
                entry.get("age", "no available"),
                entry.get("developmental stage", "no available"),
                entry.get("sex", "no available"),
                entry.get("ancestry category", "no available"),
                entry.get("disease", "no available"),
                entry.get("cell type", "no available"),
                entry.get("sampling site", "no available"),
                string_if_not_empty(entry.get("synonyms", [])),
            ]
            row = ["no available" if item is None else str(item) for item in row]
            file.write("\t".join(row) + "\n")


@click.command(
    "cellosaurus-database",
    short_help="Create the cellosaurus database from the cellosaurus file",
)
@click.option(
    "--cellosaurus",
    help="CelloSaurus database file, the file is gzipped",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "--output", help="Output file with the cellosaurus database", required=True
)
@click.option(
    "--bto", help="BTO ontology file", required=True, type=click.Path(exists=True)
)
@click.option(
    "--cl", help="Cell type ontology file", required=True, type=click.Path(exists=True)
)
@click.option(
    "--filter-species", help="Include only the following species", required=False
)
def cellosaurus_db(
    cellosaurus: str, output: str, bto: str, cl: str, filter_species
) -> None:
    """
    The following function creates a celloSaurus database from the cellosaurus file, it does some mapping to bto
    also parse organisms, diseases, etc.
    :param cellosaurus: CelloSaurus database file
    :param output: Output file with the cellosaurus database
    :param bto: BTO ontology file
    :param cl: Cell type ontology file
    :return:
    """
    # Read the BTO files
    bto = read_obo_file(bto)
    cl_type = read_obo_file(cl)

    # Parse the CelloSaurus file
    cellosaurus_list = parse_cellosaurus_file(cellosaurus, bto, cl_type)
    if filter_species:
        filter_species = filter_species.split(",")
        cellosaurus_list = [
            entry for entry in cellosaurus_list if entry["organism"] in filter_species
        ]

    current_cl_database = []
    for cellosaurus_cl in cellosaurus_list:
        cl_new_entry = create_new_entry_from_cellosaurus(cellosaurus_cl)
        current_cl_database.append(cl_new_entry)

    write_database_cellosaurus(current_cl_database, output)


@click.command(
    "mistral-recommendation",
    short_help="Looks for cell lines in cellosaurus using Mixtral LLM because the name do not match the cellosaurus name",
)
@click.option("--unknown", help="unknown cell lines", required=True)
@click.option("--output", help="File with the recommendations", required=True)
def mistral_recommendation(unknown: str, output: str) -> None:
    """
    The following function creates a vector database using LLMs using the CelloSaurus database, BTO and EFO ontologies
    :param unknown: File with the unknown cell lines
    :param output: Output file with a vector database constructed using LLMs
    :return:
    """

    login(token=os.environ.get("HUGGINGFACE_TOKEN"))

    with open(unknown, "r") as file:
        content = file.read()

    unknown_cl = content.split("\n")

    model_name = "mistralai/Mistral-7B-Instruct-v0.3"  # Replace with the actual path to the model
    generator = pipeline("text-generation", model=model_name, truncation=True)

    def get_cell_line_name(cell_code):
        prompt = f"Provide the Cellosaurus cell line name for the code {cell_code}:"
        response = generator(prompt, max_length=50, num_return_sequences=1)
        return response[0]["generated_text"].strip()

    cell_line_names = {}

    for cell_line in unknown_cl:
        cell_line_names[cell_line] = get_cell_line_name(cell_line)
        print(f"{cell_line}: {cell_line_names[cell_line]}")

    # Print the results
    with open(output, "w") as file:
        for code, name in cell_line_names.items():
            print(f"{code}: {name}")
            file.write(f"{code} - {name}\n")


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


@click.command(
    "cell-passports-database",
    short_help="Create a database from cell passports files",
)
@click.option("--cell-passports", help="Cell passports file", required=True)
@click.option("--output", help="Output file with the database", required=True)
def cell_passports_to_database(cell_passports: str, output: str) -> None:
    """
    The following function creates a database of celllines from cell passports files with the following information:
    each cell line will contain the following information:
    model name -> cell line
    synonyms
    tissue -> organism part
    cancer_type_detail -> disease second
    sample_site -> sampling site
    RRID -> cellosaurus accession
    species -> organism
    cancer_type -> disease
    gender -> sex
    ethnicity -> ancestry category
    age
    model_id
    sample_id
    patient_id

    :param cell_passports: path to the folder containing the cell passport files
    :param output: path to the output file
    :return:
    """
    cell_passports = pd.read_csv(cell_passports, sep=",", header=0)

    # Filter by model_type = Cell Line
    cell_passports = cell_passports[cell_passports["model_type"] == "Cell Line"]
    print(
        "The number of cell lines in the cell passports file is: ", len(cell_passports)
    )
    columns = [
        "model_name",
        "synonyms",
        "tissue",
        "cancer_type",
        "sample_site",
        "cancer_type_detail",
        "RRID",
        "species",
        "gender",
        "ethnicity",
        "age_at_sampling",
        "model_id",
        "sample_id",
        "patient_id",
    ]
    # sublect columns
    cell_passports = cell_passports[columns]
    cell_passports = cell_passports.fillna("no available")
    # convert age_at_sampling to no decimal places
    cell_passports["age_at_sampling"] = cell_passports["age_at_sampling"].apply(
        lambda x: int(x) if x != "no available" else x
    )
    # write pandas dataframe to file

    # rename some columns to match the database
    cell_passports = cell_passports.rename(
        columns={
            "model_name": "cell line",
            "tissue": "organism part",
            "cancer_type": "disease",
            "sample_site": "sampling site",
            "gender": "sex",
            "cancer_type_detail": "cancer type detail",
            "species": "organism",
            "age_at_sampling": "age",
            "ethnicity": "ancestry category",
            "RRID": "cellosaurus accession",
        }
    )
    cell_passports.to_csv(output, sep="\t", index=False)


@click.command(
    "ea-database", short_help="Create a database from big expression atlas files"
)
@click.option("--ea-folder", help="Expression Atlas folder", required=True)
@click.option(
    "--ea-cl-catalog", help="Expression Atlas cell line catalog", required=True
)
@click.option(
    "--output",
    help="Output file with the database",
    required=True,
    type=click.Path(exists=False),
    default="ea-cls-db.tsv",
)
def ea_create_database(ea_folder: str, ea_cl_catalog: str, output: str) -> None:
    """
    The following function creates a database of celllines file from expression atlas experiments with the following information:
    each cell line will contain the following information:
    - cell line name
    - organism
    - organism part
    - age
    - developmental stage
    - sex
    - ancestry category
    - disease

    :param ea_folder: Expression Atlas folder
    :param ea_cl_catalog: Expression Atlas cell line catalog, this is a list of cell lines curated by expression atlas.
    :param output: Output file with the database
    :return:
    """

    ea_files = glob.glob(ea_folder + "/**/*.tsv", recursive=True)

    cell_lines_dict = {}
    for file in ea_files:
        # read tab-delimited file
        data = pd.read_csv(file, sep="\t")

        # remove duplicates
        data = data.drop_duplicates(
            subset=[
                "Sample Characteristic[organism]",
                "Sample Characteristic[organism part]",
                "Sample Characteristic[cell line]",
                "Sample Characteristic[disease]",
            ]
        )
        columns_data = list(data.columns)

        # add to dictionary with cell line as key
        for i, row in data.iterrows():
            cell_line = row["Sample Characteristic[cell line]"]
            if cell_line not in cell_lines_dict:
                cell_lines_dict[cell_line] = {}

                cell_lines_dict[cell_line]["organism"] = []
                cell_lines_dict[cell_line]["organism"].append(
                    row["Sample Characteristic[organism]"]
                )
                cell_lines_dict[cell_line]["organism part"] = []
                cell_lines_dict[cell_line]["organism part"].append(
                    row["Sample Characteristic[organism part]"]
                )
                cell_lines_dict[cell_line]["disease"] = []
                cell_lines_dict[cell_line]["disease"].append(
                    row["Sample Characteristic[disease]"]
                )

                # check if the other fields are present
                cell_lines_dict[cell_line]["age"] = []
                if "Sample Characteristic[age]" in columns_data:
                    cell_lines_dict[cell_line]["age"].append(
                        row["Sample Characteristic[age]"]
                    )

                cell_lines_dict[cell_line]["developmental stage"] = []
                if "Sample Characteristic[developmental stage]" in columns_data:
                    cell_lines_dict[cell_line]["developmental stage"].append(
                        row["Sample Characteristic[developmental stage]"]
                    )

                cell_lines_dict[cell_line]["sex"] = []
                if "Sample Characteristic[sex]" in columns_data:
                    cell_lines_dict[cell_line]["sex"].append(
                        row["Sample Characteristic[sex]"]
                    )

                cell_lines_dict[cell_line]["ancestry category"] = []
                if "Sample Characteristic[ancestry category]" in columns_data:
                    cell_lines_dict[cell_line]["ancestry category"].append(
                        row["Sample Characteristic[ancestry category]"]
                    )
            else:
                # check that all the fields are the same, if not raise error:
                if (
                    cell_lines_dict[cell_line]["organism"]
                    != row["Sample Characteristic[organism]"]
                ):
                    print(f"Organism is different for cell line {cell_line}")
                if (
                    cell_lines_dict[cell_line]["organism part"]
                    != row["Sample Characteristic[organism part]"]
                ):
                    print(f"Organism part is different for cell line {cell_line}")
                if (
                    row["Sample Characteristic[disease]"]
                    not in cell_lines_dict[cell_line]["disease"]
                ):
                    cell_lines_dict[cell_line]["disease"].append(
                        row["Sample Characteristic[disease]"]
                    )
                    print(
                        f"Disease is different for cell line {cell_line} - values are {cell_lines_dict[cell_line]['disease']} and {row['Sample Characteristic[disease]']}"
                    )

                if (
                    "Sample Characteristic[age]" in columns_data
                    and row["Sample Characteristic[age]"]
                    not in cell_lines_dict[cell_line]["age"]
                ):
                    cell_lines_dict[cell_line]["age"].append(
                        row["Sample Characteristic[age]"]
                    )
                    print(f"Age is different for cell line {cell_line}")

                if (
                    "Sample Characteristic[developmental stage]" in columns_data
                    and row["Sample Characteristic[developmental stage]"]
                    not in cell_lines_dict[cell_line]["developmental stage"]
                ):
                    cell_lines_dict[cell_line]["developmental stage"].append(
                        row["Sample Characteristic[developmental stage]"]
                    )
                    print(f"Developmental stage is different for cell line {cell_line}")

                if (
                    "Sample Characteristic[sex]" in columns_data
                    and row["Sample Characteristic[sex]"]
                    not in cell_lines_dict[cell_line]["sex"]
                ):
                    cell_lines_dict[cell_line]["sex"].append(
                        row["Sample Characteristic[sex]"]
                    )

                if (
                    "Sample Characteristic[ancestry category]" in columns_data
                    and row["Sample Characteristic[ancestry category]"]
                    not in cell_lines_dict[cell_line]["ancestry category"]
                ):
                    cell_lines_dict[cell_line]["ancestry category"].append(
                        row["Sample Characteristic[ancestry category]"]
                    )

                print(f"Cell line {cell_line} already in database")

    # read the cell line catalog
    ae_cl_catalog = pd.read_csv(ea_cl_catalog, sep=",", header=0)

    # check if the cell lines in the catalog are in the database
    for i, row in ae_cl_catalog.iterrows():
        if row["cell line"] in cell_lines_dict:
            print(f"Cell line {row['cell line']} found in the database")
            if row["organism"] not in cell_lines_dict[row["cell line"]]["organism"]:
                cell_lines_dict[row["cell line"]]["organism"].append(row["organism"])
            if (
                row["organism part"]
                not in cell_lines_dict[row["cell line"]]["organism part"]
            ):
                cell_lines_dict[row["cell line"]]["organism part"].append(
                    row["organism part"]
                )
            if row["disease"] not in cell_lines_dict[row["cell line"]]["disease"]:
                cell_lines_dict[row["cell line"]]["disease"].append(row["disease"])
            if row["age"] not in cell_lines_dict[row["cell line"]]["age"]:
                cell_lines_dict[row["cell line"]]["age"].append(row["age"])
            if (
                row["developmental stage"]
                not in cell_lines_dict[row["cell line"]]["developmental stage"]
            ):
                cell_lines_dict[row["cell line"]]["developmental stage"].append(
                    row["developmental stage"]
                )
            if row["sex"] not in cell_lines_dict[row["cell line"]]["sex"]:
                cell_lines_dict[row["cell line"]]["sex"].append(row["sex"])
            cell_lines_dict[row["cell line"]]["synonyms"] = [row["synonyms"]]
        else:
            print(f"Cell line {row['cell line']} not found in the database")
            cell_lines_dict[row["cell line"]] = {}
            cell_lines_dict[row["cell line"]]["organism"] = [row["organism"]]
            cell_lines_dict[row["cell line"]]["organism part"] = [row["organism part"]]
            cell_lines_dict[row["cell line"]]["disease"] = [row["disease"]]
            cell_lines_dict[row["cell line"]]["age"] = [row["age"]]
            cell_lines_dict[row["cell line"]]["developmental stage"] = [
                row["developmental stage"]
            ]
            cell_lines_dict[row["cell line"]]["sex"] = [row["sex"]]
            cell_lines_dict[row["cell line"]]["synonyms"] = [row["synonyms"]]

    # write the ea atlas database to file as a comma separated file.
    with open(output, "w", newline="") as file:
        # Define the CSV headers
        headers = [
            "cell line",
            "organism",
            "organism part",
            "disease",
            "age",
            "developmental stage",
            "sex",
            "ancestry category",
            "synonyms",
        ]

        # Write the header row
        file.write("\t".join(headers) + "\n")

        for cell_line, data in cell_lines_dict.items():
            # Construct the row
            row = [
                cell_line,
                string_if_not_empty(data.get("organism")),
                string_if_not_empty(data.get("organism part")),
                string_if_not_empty(data.get("disease")),
                string_if_not_empty(data.get("age")),
                string_if_not_empty(data.get("developmental stage")),
                string_if_not_empty(data.get("sex")),
                string_if_not_empty(data.get("ancestry category", [])),
                string_if_not_empty(data.get("synonyms", [])),
            ]
            # Write the row
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
    help="CelloSaurus database file",
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
    The following function creates a vector database using LLMs using the CelloSaurus database, BTO and EFO ontologies
    :param database: Current database file with cell lines
    :param cellosaurus_database: CelloSaurus database file
    :param ea_database: EA Atlas database file
    :param cell_passports_database: Cell passports database file
    :param sdrf_path: SDRF folder with all existing SDRF files
    :param include_all_cellpassports: Include all cell passports cell lines
    :param ai_synonyms: AI synonyms file
    :param unknown: Output for unknown cell lines in cellosaurus
    :return:
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
        for key, ea_atlas_entry in ea_atlas.items():
            if ea_atlas_entry["cell line"].lower() == cl.lower():
                return ea_atlas_entry
            if "synonyms" in ea_atlas_entry:
                for synonym in ea_atlas_entry["synonyms"]:
                    if cl.lower() in synonym.lower():
                        return ea_atlas_entry
        return None

    def find_in_synonyms_table(cl: str, ai_synonyms: dict) -> Union[str, None]:
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
cli.add_command(mistral_recommendation)
cli.add_command(ea_create_database)
cli.add_command(cellosaurus_db)
cli.add_command(cell_passports_to_database)


if __name__ == "__main__":
    cli()
