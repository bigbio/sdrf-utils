import os
import click
from transformers import pipeline, login
import os

import click
from huggingface_hub import login
from transformers import pipeline

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

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
    # Log in using the Hugging Face token from environment variables
    login(token=os.environ.get("HUGGINGFACE_TOKEN"))

    # Read unknown cell lines from the input file
    with open(unknown, "r") as file:
        unknown_cl = file.read().splitlines()

    # Initialize the text generation model
    model_name = "mistralai/Mistral-7B-Instruct-v0.3"
    generator = pipeline("text-generation", model=model_name, truncation=True)

    # Function to get the Cellosaurus cell line name
    def get_cell_line_name(cell_code):
        prompt = f"Provide the Cellosaurus cell line name for the code {cell_code}:"
        response = generator(prompt, max_length=50, num_return_sequences=1)
        return response[0]["generated_text"].strip()

    # Generate cell line names and store them in a dictionary
    cell_line_names = {cell_line: get_cell_line_name(cell_line) for cell_line in unknown_cl}

    # Write the results to the output file
    with open(output, "w") as file:
        for code, name in cell_line_names.items():
            print(f"{code}: {name}")
            file.write(f"{code} - {name}\n")


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    Main function to run the CLI
    """
    pass


if __name__ == "__main__":
    cli()
