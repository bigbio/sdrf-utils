import gzip
import re
from typing import Union, List, Dict, Tuple, Optional
import click
import numpy as np

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

def read_obo_file(file_path: str) -> Dict[str, dict]:
    """
    Reads an OBO file and returns a dictionary with the parsed content.

    Parameters:
        file_path (str): Path to the OBO file.

    Returns:
        Dict[str, dict]: A dictionary where keys are OBO term IDs and values are dictionaries of term attributes.
    """
    with open(file_path, "r") as file:
        content = file.read()

    entries = content.split("\n\n")

    def parse_obo_term(entry: str) -> dict:
        obo_dict = {}
        lines = entry.strip().split("\n")
        for line in lines:
            if line.startswith("id:"):
                obo_dict["id"] = line.split("id: ")[1].strip()
            elif line.startswith("name:"):
                obo_dict["name"] = line.split("name: ")[1].strip()
            elif line.startswith("def:"):
                obo_dict["def"] = line.split("def: ")[1].strip()
            elif line.startswith("is_a:"):
                obo_dict.setdefault("is_a", []).append(line.split("is_a: ")[1].strip())
            elif line.startswith("synonym:"):
                obo_dict.setdefault("synonyms", []).append(line.split("synonym: ")[1].strip())
                obo_dict["synonyms"] = [
                    synonym.replace("RELATED []", "")
                    .replace("RELATED MS []", "")
                    .strip()
                    .strip('"')
                    for synonym in obo_dict["synonyms"]
                ]
        obo_dict.setdefault("synonyms", [])
        return obo_dict

    return {
        parse_obo_term(entry)["id"]: parse_obo_term(entry)
        for entry in entries
        if entry.strip() and "id:" in entry
    }

def string_if_not_empty(param: List[Union[str, float]]) -> str:
    """
    Returns a string if the list is not empty, otherwise returns 'no available'.

    Parameters:
        param (List[Union[str, float]]): List of strings or floats.

    Returns:
        str: Concatenated string of list elements or 'no available'.
    """
    if param == "None":
        param = []
    if param:
        l = [
            x for x in param
            if (isinstance(x, float) and not np.isnan(x)) or (not isinstance(x, float) and x is not None)
        ]
        return "; ".join(l)
    return "no available"

def write_database_cellosaurus(current_cl_database: List[dict], database: str) -> None:
    """
    Writes the database objects to the specified database file.

    Parameters:
        current_cl_database (List[dict]): Current cell line database list.
        database (str): Path to the database file.
    """
    headers = [
        "cellosaurus name", "cellosaurus accession", "bto cell line", "organism",
        "age", "developmental stage", "sex", "ancestry category", "disease",
        "cell type", "sampling site", "synonyms"
    ]

    with open(database, "w") as file:
        file.write("\t".join(headers) + "\n")
        for entry in current_cl_database:
            row = [
                entry.get(header, "no available") for header in headers[:-1]
            ] + [string_if_not_empty(entry.get("synonyms", []))]
            file.write("\t".join(row) + "\n")

def is_age_in_text(age_text: str) -> bool:
    """
    Checks if the age field contains a number.

    Parameters:
        age_text (str): Age text to check.

    Returns:
        bool: True if the age is a number, False otherwise.
    """
    return any(char.isdigit() for char in age_text)

def get_sampling_site(cellosaurus_comment: str) -> Optional[str]:
    """
    Extracts the sampling site from a Cellosaurus comment string.

    Parameters:
        cellosaurus_comment (str): The comment string from which to extract the sampling site information.

    Returns:
        Optional[str]: The tissue sampling site if found, otherwise None.
    """
    pattern = r"Derived from site:\s*(.+?);\s*(.+?);\s*UBERON=(UBERON_\d{7})"
    match = re.search(pattern, cellosaurus_comment)
    return match.group(2) if match else None

def get_cell_type(cellosaurus_comment: str) -> Optional[str]:
    """
    Extracts the cell type from the Cellosaurus comment field.

    Parameters:
        cellosaurus_comment (str): Cellosaurus comment field.

    Returns:
        Optional[str]: Cell type if found, otherwise None.
    """
    pattern = r"Cell type:\s*(.+?);\s*CL=(CL_\d{7})"
    match = re.search(pattern, cellosaurus_comment)
    return match.group(2).replace("_", ":").strip() if match else None

def parse_cellosaurus_taxonomy(organism_text: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parses the organism text from the Cellosaurus database.

    Parameters:
        organism_text (str): Organism text from the Cellosaurus database.

    Returns:
        Tuple[Optional[str], Optional[str]]: Species name and taxonomy ID if found, otherwise (None, None).
    """
    pattern = r"NCBI_TaxID=(\d+); ! ([\w\s]+) \(([\w\s]+)\)"
    match = re.search(pattern, organism_text)
    return (match.group(2), match.group(1)) if match else (None, None)

def parse_cellosaurus_file(file_path: str, bto: Dict[str, dict], cl_type: Dict[str, dict]) -> List[dict]:
    """
    Parses the CelloSaurus file and returns a list of dictionaries with the parsed content.

    Parameters:
        file_path (str): CelloSaurus file path.
        bto (Dict[str, dict]): BTO ontology dictionary.
        cl_type (Dict[str, dict]): Cell type ontology dictionary.

    Returns:
        List[dict]: List of CelloSaurus dictionaries.
    """
    def parse_entry(entry: str, bto: Dict[str, dict], cl_type: Dict[str, dict]) -> dict:
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
                    data["synonyms"].extend(bto[bto_accession].get("synonyms", []))
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
                    pattern = r"NCIt;\s*C\d+;\s*([^;]+)"
                    match = re.search(pattern, line)
                    if match:
                        data["disease"] = match.group(1).strip()
            elif line.startswith("CC") and "Derived from site" in line:
                data["sampling site"] = get_sampling_site(line)
            elif line.startswith("CC") and "Cell type" in line:
                code_cl = get_cell_type(line)
                if code_cl and code_cl in cl_type:
                    data["cell type"] = cl_type[code_cl]["name"]

        return data

    with gzip.open(file_path, "r") as file:
        content = file.read().decode("utf-8")

    entries = content.split("//\n")
    parsed_data = [parse_entry(entry, bto, cl_type) for entry in entries if entry.strip()]
    return [entry for entry in parsed_data if entry]

def create_new_entry_from_cellosaurus(cellosaurus: dict) -> dict:
    """
    Creates a new entry for a cell line not found in the database.

    Parameters:
        cellosaurus (dict): Dictionary containing cellosaurus data.

    Returns:
        dict: New entry dictionary.
    """
    entry = {
        "cellosaurus name": cellosaurus.get("cellosaurus name", "no available"),
        "bto cell line": cellosaurus.get("bto cell line"),
        "cellosaurus accession": cellosaurus.get("cellosaurus accession"),
        "organism": cellosaurus.get("organism"),
        "organism part": None,
        "age": None,
        "developmental stage": None,
        "sex": cellosaurus.get("sex"),
        "ancestry category": cellosaurus.get("ancestry category"),
        "disease": cellosaurus.get("disease"),
        "synonyms": cellosaurus.get("synonyms", []),
        "cell type": cellosaurus.get("cell type"),
        "sampling site": cellosaurus.get("sampling site"),
    }

    if "age" in cellosaurus:
        if is_age_in_text(cellosaurus["age"]):
            entry["age"] = cellosaurus["age"]
        else:
            entry["developmental stage"] = cellosaurus["age"]

    return entry

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
    cellosaurus: str, output: str, bto: str, cl: str, filter_species: Optional[str]
) -> None:
    """
    Creates a CelloSaurus database from the cellosaurus file, mapping to BTO and parsing organisms, diseases, etc.

    Parameters:
        cellosaurus (str): CelloSaurus database file.
        output (str): Output file with the cellosaurus database.
        bto (str): BTO ontology file.
        cl (str): Cell type ontology file.
        filter_species (Optional[str]): Comma-separated list of species to include.
    """
    bto_data = read_obo_file(bto)
    cl_type_data = read_obo_file(cl)

    cellosaurus_list = parse_cellosaurus_file(cellosaurus, bto_data, cl_type_data)
    if filter_species:
        filter_species_list = filter_species.split(",")
        cellosaurus_list = [
            entry for entry in cellosaurus_list if entry["organism"] in filter_species_list
        ]

    current_cl_database = [
        create_new_entry_from_cellosaurus(cellosaurus_cl) for cellosaurus_cl in cellosaurus_list
    ]

    write_database_cellosaurus(current_cl_database, output)

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    Main function to run the CLI.
    """
    pass

cli.add_command(cellosaurus_db)

if __name__ == "__main__":
    cli()