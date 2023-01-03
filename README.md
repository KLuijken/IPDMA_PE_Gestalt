# IPDMA_PE_Gestalt

IN PROGRESS

This repository accomplishes the manuscript 'Accuracy of the Physicians’ Intuitive Risk Estimation in the Diagnostic Management of Pulmonary Embolism: An Individual Patient Data Meta-Analysis' by Rosanne van Maanen, Emily S.L. Martens, Toshihiko Takada, Pierre-Marie Roy, Kerstin de Wit, Sameer Parpia, Noémie Kraaijpoel, Menno V. Huisman, Philip S. Wells, Gregoire le Gal, Marc Righini, Yonathan Freund, Javier Galipienzo, Jeanet W. Blom, Karel G.M. Moons, Frans H. Rutten, Maarten van Smeden, Frederikus A. Klok, Geert-Jan Geersing, Kim Luijken

## Purpose of scripts
The scripts in the current repository can be used to replicate the results in the main text of the above mentioned study, as well as the additional analyses presented in the supplementary materials.

## Run analysis
An explainer of the code can be found [here]()
The script in the current repository can be run using a simple simulated dataset that does not contain meaningful information. 
The file ./rcode/exe runs the analysis and generates figures at once. This takes around XX time.  

## Generating tables and figures
The R code in ./rcode/visualisation creates figures depicting the results of the analysis. The script produces .pdf and .txt files. 

Attached packages:  

   
 

## Project organization

```
.
├── .gitignore
├── CITATION.md
├── LICENSE.md
├── README.md
├── data                      - contains a simulated dataset
├── docs
│   ├── code_explainer        - html file explaining all code
├── results
│   ├── figures               - figures shown in manuscript
│   ├── tables                - tables shown in manuscript
└── rcode                     - source code for this project
    ├── analysis              - helper scripts and script for main analysis
    ├── exe                   - main script from which analysis is executed
    └── packages              - required dependencies

```

## Citation


