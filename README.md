
[![memote tested](https://img.shields.io/badge/memote-tested-blue.svg?style=plastic)](https://sysbiochalmers.github.io/Sco-GEM/)
[![Build Status](https://travis-ci.org/SysBioChalmers/Sco-GEM.svg?branch=master)](https://travis-ci.org/SysBioChalmers/Sco-GEM)
[![Join the chat at https://gitter.im/SysBioChalmers/Sco-GEM](https://badges.gitter.im/SysBioChalmers/Sco-GEM.svg)](https://gitter.im/SysBioChalmers/Sco-GEM?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/SysBioChalmers/Sco-GEM?style=plastic)
# Sco-GEM: The consensus genome-scale metabolic model of _Streptomyces coelicolor_

### Brief Repository Description
This repository contains the consensus genome-scale metabolic model **Sco-GEM** for the antibiotic producer _Streptomyces coelicolor_ A3(2), a representative species of _Actinomycetales_.

### Abstract
_Streptomyces coelicolor_ is a representative species of soil-dwelling, filamentous and gram-positive actinobacterium harbouring enriched secondary metabolite biosynthesis gene clusters. As a well-known pharmaceutical and bioactive compound producer, _S. coelicolor_ has been exploited for antibiotic and secondary metabolite production.

### Model KeyWords
**GEM Category:** Species; **Utilisation:** Predictive simulation; **Field:** Metabolic-network reconstruction; **Type of Model:** Reconstruction and refinement; **Model Source:** [iKS1317](http://dx.doi.org/); **Omic Source:** [Genomics](http://dx.doi.org/10.1038/417141a); **Taxonomy:** _Streptomyces coelicolor_; **Metabolic System:** General Metabolism; **Strain:** A3(2); **Condition:** Complex medium;

### Reference
S. Sulheim, T. Kumelj, D. van Dissel, A. Salehzadeh-Yazdi, C. Du, G. P. van Wezel, K. Nieselt, E. Almaas, A. Wentzel, E. J. Kerkhoven (2020). Enzyme-constrained models and omics analysis of Streptomyces coelicolor reveal metabolic changes that enhance heterologous production. Iscience, 23(9), 101525


- Pubmed ID: https://pubmed.ncbi.nlm.nih.gov/32942174/
- doi: https://doi.org/10.1016/j.isci.2020.101525
- Last update: 2020-07-18

- The model contains:

| Taxonomy | Template Model | Reactions | Metabolites| Genes | Memote score |
| ------------- |:-------------:|:-------------:|:-------------:|-----:|-----:|
| _Streptomyces coelicolor_ A3(2) | iKS1317 | 2612 | 2073 | 1777 | 77%|


This repository is administered by Eduard Kerkhoven ([@SysBioChalmers](https://github.com/SysBioChalmers)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology

## Installation

### Recommended Software:
* A functional Matlab installation (MATLAB 7.3 or higher).
* Python 3.6 or 3.7 including standard packages and additional packages described in requirements.txt
* [RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN) for MATLAB (required for contributing to development). 
* libSBML MATLAB API ([version 5.16.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/)  is recommended).
* [Gurobi Optimizer for MATLAB](http://www.gurobi.com/registration/download-reg).
* For contributing to development: a [git wrapper](https://github.com/manur/MATLAB-git) added to the search path.

### Installation Instructions
* Clone the [master](https://github.com/SysBioChalmers/sco-GEM) branch from [SysBioChalmers GitHub](https://github.com/SysBioChalmers).
* Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com).

### Contributors
* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Sweden
* [Snorre Sulheim](https://www.sintef.no/en/all-employees/employee/?empId=5675) ([@sulheim](https://github.com/sulheim)), SINTEF Industry, Norway / Norwegian University of Science and Technology, Norway
* [Tjasa Kumelj](https://www.ntnu.edu/employees/tjasa.kumelj) ([@tjasakumelj](https://github.com/tjasakumelj)), Norwegian University of Science and Technology, Norway
* [Ali Salehzadeh-Yazdi](https://www.sbi.uni-rostock.de/team/detail/ali-salehzadeh-yazdi) ([@alisalehzadeh-yazdi](https://github.com/alisalehzadeh-yazdi)), University of Rostock, Germany

