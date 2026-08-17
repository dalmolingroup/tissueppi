# tissuePPI-app

## Name
Tissue-specific interactome

## Description
Get a protein-protein interaction (PPI) network on the context of your target tissue. 
Chose a gene symbol and a tissue, and get the PPI Network for the protein and their interactions on the tissue selected.

The results are available as an interative network or a downloadable HTML/CSV file.

## Application diagram

```mermaid
flowchart TD
    A[User Shiny Interface] -->|Selects gene, tissue| B[Front-end Shiny]
    B -->|HTTP request| C[Plumber API]
    C -->|Query to SQLite DB| D[SQLite: Protein-protein interactions]
    D -->|Returns filtered data| C
    C -->|Processes data and generates outputs| E[Output]
    E -->|Interactive network HTML| B
    E -->|Interactive image HTML| F[Temporary storage]
    E -->|Edgelist file CSV| F
    B -->|Displays network| A
    A -->|User downloads files| B
    B -->|GET requests to download files| F
```

## Installation

To run the application locally:
 - First create a directory
 - Download the data files to the created directory
 - Change to created directory
 - Clone de application with command: git clone 
 - Create the database with the command: sqlite3 ./dados/interacoes.sqlite < create-database.sql.
 
 Then, make sure you have docker and docker compose installed, and you only need to run
 
 `docker compose up -d`
 
And open http://localhost:3838/ in your browser of preference.

## Authors and acknowledgment
This app is authored by Julia Apolonio and Fábio Silva, with supervision of dr. Patrick Terremate and dr. Rodrigo Dalmolin
