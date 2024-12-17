# Cell Line Metadata Database

This repository contains scripts for creating and managing a cell line metadata database to enable the annotation of SDRFs for cell lines datasets. The main use case is the annotation of SDRF datasets for the [quantms.org resource](quantms.org). This repo uses multiple ontologies and natural language processing (NLP) to annotate cell lines in SDRF files. 

## Motivation

Cell lines are a fundamental part of biological research, and they are used in a wide range of experiments. However, cell line metadata can be inconsistent and difficult to manage. Here we are creating a database that can be used to annotate/validate proteomics SDRF for cell lines studies. These are the major sources of cell line metadata:

- [Cellosaurus](https://web.expasy.org/cellosaurus/): Cellosaurus is the **main source** of metadata in our database. The source of the metadata can be downloaded from [cellosaurus.txt](https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt). We converted the file to a shorter version with only the fields that we are interested and the taxonomy. We use the script `[cellosaurus_db.py](cellosaurus_db.py)` to create the database.
- [Cell model passports](https://cog.sanger.ac.uk/cmp/download/model_list_20240110.csv): The cell model passports are a collection of cell lines from multiple sources. We use the file [model_list_20240110.csv](model_list_20240110.csv) to create a database extracting only the cell lines information `[cellpassports_db.py](cellpassports_db.py)`.
- [EA](https://https://www.ebi.ac.uk/gxa): Expression Atlas has been curating for more than 10 years the metadata of multiple RNA experiments. We collect multiple cell lines experiments from EA in folder [ea](ea); and try to create a catalog of cell lines metadata as an extra source. We use the following script `[ea_db.py](ea_db.py)` to create the database.

## Ontologies

- [MONDO](https://bioportal.bioontology.org/ontologies/MONDO): The Monarch Disease Ontology (MONDO) is used to annotate the disease of the cell line.
- [BTO](https://bioportal.bioontology.org/ontologies/BTO): The BRENDA Tissue Ontology (BTO) is used to annotate an extra reference for the cell line ID. 

> **Note**: Additionally, we use other resources such as [Coriell cell line Catalog](https://www.coriell.org/), [cell bank riken](https://cell.brc.riken.jp/en/) and [atcc](https://www.atcc.org/) for manual annotation of cell lines in the database. 

## Features for every cell line

The database is created using SQLite and contains the following fields:

- **cell line**: The cell line name as defined by the curation team (ai or manual).
- **cellosaurus name**: The cell line name as annotated in Cellosaurus `ID` 
- **cellosaurus accession**: The cell line accession as annotated in Cellsaurus `AC`
- **bto cell line**: The cell line name as annotated in BTO
- **organism**: The organism of the cell line as annotated in Cellosaurus
- **organism part**: This information is not available in Cellosaurus, we use other sources to _annotate_ this field.
- **sampling site**: The sampling site of the cell line as annotated in Cellosaurus. If the information is not available, we use other sources to _annotate_ this field.
- **age**: The age of the cell line as annotated in Cellosaurus. If the age is not available (empty), we annotated the age from other sources such as [atcc](https://www.atcc.org/) or [Coriell cell line Catalog](https://www.coriell.org/)
- **developmental stage**: The developmental stage of the cell line as annotated in Cellosaurus; if the information is not available is inferred from the age of the cell line. 
- **sex**: Sex as provided by Cellosaurus
- **ancestry category**: The ancestry category of the cell line as annotated in Cellosaurus. If not available we use other sources. 
- **disease**: The disease is _"agreed"_ among sources.  
- **cell type**: The cell type is _"agreed"_ among sources.
- **Material type**: The material is _"agreed"_ among sources.
- **synonyms**: This field is built using all the accessions and synonyms from all sources.
- **curated**: This field is used to annotate if the cell line has been curated by the team, the classes are _not curated_, _ai curated_, _manual curated_.

> **Note**: The database is a tab-delimited file that can be easily read and search using pandas or GitHub table rendering. 

## Requirements

- Python 3.x
- Libraries: `pandas`, `spacy`, `click`, `owlready2`, `scikit-learn`
- Spacy model: `en_core_web_md`

## Installation

Install the required Python packages using pip: