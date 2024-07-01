# A case study of applying a pretrained ML model to predict types of acute leukemia for a German hospital cohort

Artificial intelligence (AI) and machine learning (ML) algorithms have shown great promise in clinical medicine, offering potential improvements in diagnostic accuracy and patient outcomes. Despite the increasing number of published algorithms, most remain unvalidated in real-world clinical settings. This study aims to simulate the practical implementation challenges of a recently developed ML algorithm, AI-PAL, designed for the diagnosis of acute leukemia and report on its performance.  

We conducted a detailed simulation of the AI-PAL algorithm's implementation at the University Hospital Essen. The study involved building a cohort from our clinical research database, identifying all initially diagnosed patients with acute leukemia. The algorithm's performance was assessed by reproducing the original study's results and recalibrating diagnosis thresholds to fit the local patient population.


## AI-PAL machine learning model
The pretrained AI-PAL model used in this worrk was developed and published by **Alcazer et al., 2024**. It is a free and open-source software package built in R.

**Evaluation of a machine-learning model based on laboratory parameters for the prediction of acute leukaemia subtypes: a multicentre model development and validation study in France**<br>
Alcazer, Vincent et al.<br>
The Lancet Digital Health, Volume 6, Issue 5, e323 - e333<br>
DOI:https://doi.org/10.1016/S2589-7500(24)00044-X<br>
https://github.com/VincentAlcazer/AIPAL, 01.07.2024

## Getting started
The cohort of patients with acute leukemia was extratced from the FHIR server of the University Medicine Essen.
The Jupter notebook for extracting the cohort is provided [here](cohort-extraction-fhir/cohort_ume_aipal.ipynb).<br>
The R scripts for calculating the AUROC, confidence tables and PPV/NPV optimizations for predictions made by the pretrained AI-PAL model based for cohort files are provided [here](AIPAL_predictions_R).


## Requirements
- Python 3.6 or larger
- R 4.3 or larger
