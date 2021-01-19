## Sco-GEM: The consensus genome-scale metabolic model of _Streptomyces coelicolor_

[![memote tested](https://img.shields.io/badge/memote-tested-blue.svg?style=plastic)](https://sysbiochalmers.github.io/Sco-GEM/)
[![Build Status](https://travis-ci.org/SysBioChalmers/Sco-GEM.svg?branch=master)](https://travis-ci.org/SysBioChalmers/Sco-GEM)
[![Join the chat at https://gitter.im/SysBioChalmers/Sco-GEM](https://badges.gitter.im/SysBioChalmers/Sco-GEM.svg)](https://gitter.im/SysBioChalmers/Sco-GEM?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
![Version release (latest by date)](https://img.shields.io/github/v/release/SysBioChalmers/Sco-GEM?style=plastic)

### Description

This repository contains the consensus genome-scale metabolic model **Sco-GEM** for the antibiotic producer _Streptomyces coelicolor_ A3(2), a representative species of soil-dwelling, filamentous and gram-positive actinobacterium harbouring enriched secondary metabolite biosynthesis gene clusters. As a well-known pharmaceutical and bioactive compound producer, _S. coelicolor_ has been exploited for antibiotic and secondary metabolite production.

### Citation

> S. Sulheim, T. Kumelj, D. van Dissel, A. Salehzadeh-Yazdi, C. Du, G.P. van Wezel, K. Nieselt, E. Almaas, A. Wentzel, E.J. Kerkhoven (2020). Enzyme-constrained models and omics analysis of _Streptomyces coelicolor_ reveal metabolic changes that enhance heterologous production. iScience, 23(9), 101525 [doi:10.1016/j.isci.2020.101525](https://doi.org/10.1016/j.isci.2020.101525), [pmid:32942174](https://pubmed.ncbi.nlm.nih.gov/32942174/)

The Sco-GEM model distributed on this GitHub repository is continuously updated. To get access to the models, data and code associated to the Sulheim _et al_. (2020) publication, use [Sco-GEM release 1.2.1](https://github.com/SysBioChalmers/Sco-GEM/releases/tag/v1.2.1).

### Keywords

**Utilisation:** Predictive simulation; Multi-omics integrative analysis  
**Field:** Metabolic-network reconstruction  
**Type of model:** Curated reconstruction  
**Model source:** [iKS1317](http://doi.org/10.1002/biot.201800180)
**Omic source:** [Genomics](http://dx.doi.org/10.1038/417141a); [Transcriptomics](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132487); [Proteomics](http://dx.doi.org/10.6019/PXD013178) 
**Taxonomic name:** _Streptomyces coelicolor_
**Taxonomy ID:** [taxonomy:100226](https://identifiers.org/taxonomy:100226)  
**Genome ID:** [insdc.gca:GCA_000203835.1](https://identifiers.org/insdc.gca:GCA_000203835.1)
**Metabolic system:** General metabolism (primary and secondary)
**Strain:** A3(2)  
**Condition:** Complex medium  

### Model Overview

| Taxonomy | Template Model | Reactions | Metabolites| Genes | Memote score |
| ------------- |:-------------:|:-------------:|:-------------:|-----:|-----:|
| _Streptomyces coelicolor_ A3(2) | iKS1317 | 2612 | 2073 | 1777 | 77%|

### Installation & Usage

Sco-GEM is primarily distributed in SBML L3V1 FBCv1 format, and therefore  works well with RAVEN Toolbox, cobrapy, COBRA Toolbox, etc.

Development of the model is preferably done via cobrapy, using the [export.py](.code/export.py) function distributed in /code to save new model version, to reduce the number of changes in the model and aid comparison of model versions.

### Contributing

Contributions are always welcome! Please read the [contributing guideline](.github/CONTRIBUTING.md) to get started.

### Contributors
* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Sweden
* [Snorre Sulheim](https://www.sintef.no/en/all-employees/employee/?empId=5675) ([@sulheim](https://github.com/sulheim)), SINTEF Industry, Norway / Norwegian University of Science and Technology, Norway
* [Tjasa Kumelj](https://www.ntnu.edu/employees/tjasa.kumelj) ([@tjasakumelj](https://github.com/tjasakumelj)), Norwegian University of Science and Technology, Norway
* [Ali Salehzadeh-Yazdi](https://www.sbi.uni-rostock.de/team/detail/ali-salehzadeh-yazdi) ([@alisalehzadeh-yazdi](https://github.com/alisalehzadeh-yazdi)), University of Rostock, Germany