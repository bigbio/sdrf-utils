import click
import pandas as pd

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

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


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    Main function to run the CLI
    """
    pass

cli.add_command(cell_passports_to_database)

if __name__ == "__main__":
    cli()
