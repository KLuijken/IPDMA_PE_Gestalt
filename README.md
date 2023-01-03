# IPDMA_PE_Gestalt

IN PROGRESS

This repository accomplishes the manuscript 'Accuracy of the Physicians’ Intuitive Risk Estimation in the Diagnostic Management of Pulmonary Embolism: An Individual Patient Data Meta-Analysis' by Rosanne van Maanen, Emily S.L. Martens, Toshihiko Takada, Pierre-Marie Roy, Kerstin de Wit, Sameer Parpia, Noémie Kraaijpoel, Menno V. Huisman, Philip S. Wells, Gregoire le Gal, Marc Righini, Yonathan Freund, Javier Galipienzo, Jeanet W. Blom, Karel G.M. Moons, Frans H. Rutten, Maarten van Smeden, Frederikus A. Klok, Geert-Jan Geersing, Kim Luijken

## Purpose of scripts
The scripts in the current repository can be used to replicate the results in the main text of the above mentioned study, as well as the additional analyses presented in the supplementary materials.

## Run analysis
An explainer of the code can be found [here](https://kluijken.github.io/IPDMA_PE_Gestalt/)
The script in the current repository can be run using a simple simulated dataset that does not contain meaningful information. 
The file ./1_exe.R runs the analysis and generates figures at once. This takes around 15 seconds.  

## Generating tables and figures
The R code in ./5_analysis_main.R and ./6_analysis_descriptive.R create figures and results presenting the results of the analysis. The script produces .pdf and .docx files. 

Attached packages:  
gt_0.8.0  
rstanarm_2.21.3  
lme4_1.1-31  
dplyr_1.0.10  
metafor_3.8-1  
mice_3.15.0  
MASS_7.3-58.1  

## Project organization

```
.
├── .gitignore
├── CITATION.md
├── LICENSE.md
├── README.md
├── IPDMA-PE.Proj
├── data                      - contains a simulated dataset
├── docs
│   ├── code_explainer        - html file explaining all code
├── figures                   - outputted figures
├── tables                    - outputted tables
├── 1_exe.R                   - main script from which analysis is executed
├── 2_packages.R
├── 3_simulate_data.R
├── 4_helpers.R               - user defined functions for analysis
├── 5_analysis_main.R
└── 6_analysis_descriptive.R

```

## Citation


