import glob
from typing import Union

import click
import numpy as np
import pandas as pd

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

def string_if_not_empty(param: list) -> Union[None, str]:
    """
    Return a string if the list is not empty
    :param param: List
    :return: None if the list is empty, the string otherwise
    """
    if param == "None":
        param = []
    if param and len(param) > 0:
        filtered_elements = [
            x
            for x in param
            if isinstance(x, float)
            and not np.isnan(x)
            or not isinstance(x, float)
            and x is not None
        ]
        return "; ".join(filtered_elements)
    return "no available"


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

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    Main function to run the CLI
    """
    pass

cli.add_command(ea_create_database)

if __name__ == "__main__":
    cli()