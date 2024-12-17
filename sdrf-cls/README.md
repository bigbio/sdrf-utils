# SDRF Cell Line Metadata Database

This repository provides tools to create and manage a **cell line metadata database** for annotating SDRFs (Sample and Data Relationship Format) in proteomics studies. The primary use case is enhancing annotation consistency for [quantms.org](https://quantms.org) datasets. The scripts integrate multiple ontologies and **natural language processing (NLP)** methods to standardize cell line metadata.

---

## Table of Contents

1. [Motivation](#motivation)  
2. [Metadata Sources](#metadata-sources)  
3. [Ontologies](#ontologies)  
4. [Database Structure](#database-structure)  
5. [Features](#features)  
6. [Requirements](#requirements)  
7. [Installation](#installation)  
8. [Usage](#usage)  
9. [Contribution](#contribution)  
10. [License](#license)

---

## Motivation

Cell lines are essential in biological research but often lack standardized metadata, leading to inconsistencies. This repository aims to:

- **Create a centralized database** for cell line metadata.
- **Facilitate annotation and validation** of cell line SDRFs, particularly in proteomics datasets.

---

## Metadata Sources

We integrate metadata from **three main sources** and additional curation efforts:

1. **[Cellosaurus](https://web.expasy.org/cellosaurus/):**  
   The primary metadata source.  
   - Download: [cellosaurus.txt](https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt)  
   - Script: [`cellosaurus_db.py`](cellosaurus_db.py) extracts relevant fields.  

2. **[Cell Model Passports](https://cog.sanger.ac.uk/cmp):**  
   A collection of cell lines from various sources.  
   - Input file: [model_list_20240110.csv](model_list_20240110.csv)  
   - Script: [`cellpassports_db.py`](cellpassports_db.py) processes this data.  

3. **[Expression Atlas (EA)](https://www.ebi.ac.uk/gxa):**  
   Metadata curated over RNA experiments for over 10 years.  
   - Collected data: Stored in the `ea/` folder.  
   - Script: [`ea_db.py`](ea_db.py) processes this source.  

> **Additional Curation**: Manual annotation is performed using data from:  
> - [Coriell Cell Line Catalog](https://www.coriell.org/)  
> - [Cell Bank RIKEN](https://cell.brc.riken.jp/en/)  
> - [ATCC](https://www.atcc.org/)

---

## Ontologies

The following ontologies are used for annotation:

1. **[MONDO](https://bioportal.bioontology.org/ontologies/MONDO):**  
   Used to annotate the disease associated with a cell line.

2. **[BTO](https://bioportal.bioontology.org/ontologies/BTO):**  
   Provides additional references for cell line IDs.

---

## Database Structure

The database is implemented using **SQLite** and contains the following key fields:

| Field Name              | Description                                                                 |
|-------------------------|-----------------------------------------------------------------------------|
| **cell line**           | Cell line name (curated: AI or manual).                                     |
| **cellosaurus name**    | Name as annotated in Cellosaurus `ID`.                                      |
| **cellosaurus accession** | Accession ID from Cellosaurus `AC`.                                       |
| **bto cell line**       | Name as annotated in BTO.                                                   |
| **organism**            | Organism species (from Cellosaurus).                                        |
| **organism part**       | Annotated from supplementary sources.                                       |
| **sampling site**       | Sampling site of the cell line.                                             |
| **age**                 | Age of the cell line (from Cellosaurus or additional sources).              |
| **developmental stage** | Developmental stage (inferred from age if missing).                        |
| **sex**                 | Sex information (from Cellosaurus).                                         |
| **ancestry category**   | Ancestry classification (from Cellosaurus or supplementary sources).        |
| **disease**             | Agreed-upon disease annotation across sources.                             |
| **cell type**           | Agreed-upon cell type annotation across sources.                           |
| **material type**       | Agreed-upon material classification.                                       |
| **synonyms**            | Consolidated synonyms and accessions from all sources.                     |
| **curated**             | Curation status: `_not curated_`, `_AI curated_`, or `_manual curated_`.    |

> **Note**: The final database is provided as a **tab-delimited file** for easy integration. It can be loaded into tools like **Pandas** or viewed directly via GitHub's table renderer.

---

## Features

- Standardizes metadata from **multiple sources**.
- Uses **ontologies** to annotate diseases and tissue information.
- Supports **AI-based curation** and manual validation for accuracy.
- Provides **easy-to-query** tab-delimited outputs.

---

## Requirements

To use the scripts, ensure the following are installed:

- **Python 3.x**
- Required libraries:  
