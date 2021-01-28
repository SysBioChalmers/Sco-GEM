## Sco-GEM: The consensus genome-scale metabolic model of _Streptomyces coelicolor_

[![memote tested](https://img.shields.io/badge/memote-tested-blue.svg?style=plastic)](https://sysbiochalmers.github.io/Sco-GEM/)
[![Build Status](https://travis-ci.org/SysBioChalmers/Sco-GEM.svg?branch=master)](https://travis-ci.org/SysBioChalmers/Sco-GEM)
![Version release (latest by date)](https://img.shields.io/github/v/release/SysBioChalmers/Sco-GEM?style=plastic)
[![DOI](https://zenodo.org/badge/145685631.svg)](https://zenodo.org/badge/latestdoi/145685631)

### Description

This repository contains the consensus genome-scale metabolic model **Sco-GEM** for the antibiotic producer _Streptomyces coelicolor_ A3(2), a representative species of soil-dwelling, filamentous and gram-positive actinobacterium harbouring enriched secondary metabolite biosynthesis gene clusters. As a well-known pharmaceutical and bioactive compound producer, _S. coelicolor_ has been exploited for antibiotic and secondary metabolite production.

### Citation

> S. Sulheim, T. Kumelj, D. van Dissel, A. Salehzadeh-Yazdi, C. Du, G.P. van Wezel, K. Nieselt, E. Almaas, A. Wentzel, E.J. Kerkhoven (2020). Enzyme-constrained models and omics analysis of _Streptomyces coelicolor_ reveal metabolic changes that enhance heterologous production. iScience, 23(9), 101525 [doi:10.1016/j.isci.2020.101525](https://doi.org/10.1016/j.isci.2020.101525), [pmid:32942174](https://pubmed.ncbi.nlm.nih.gov/32942174/)

The Sco-GEM model distributed on this GitHub repository is continuously updated, with the latest releases available [here](https://github.com/SysBioChalmers/Sco-GEM/releases). To get access to the models, data and code associated to the Sulheim _et al_. (2020) publication, use [Sco-GEM release 1.2.1](https://github.com/SysBioChalmers/Sco-GEM/releases/tag/v1.2.1).

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

## Installation & Usage

### **User:**

To obtain Sco-GEM, clone it from [`master`](https://github.com/sysbiochalmers/Sco-GEM) in the GitHub repository, or just download the [latest release](https://github.com/sysbiochalmers/Sco-GEM/releases).

Sco-GEM is distributed in SBML L3V1 FBCv1 format (`model/Sco-GEM.xml`), and therefore works well with any appropriate constraint-based modelling package, such as [cobrapy](https://github.com/opencobra/cobrapy), [RAVEN Toolbox](https://github.com/sysbiochalmers/raven/) and [COBRA Toolbox](https://github.com/opencobra/cobratoolbox). Installation instructions for each package are provided on their website, after which you can use their default functions for loading and exporting of the models:

***cobrapy***
```python
import cobra
model = cobra.io.read_sbml_model('Sco-GEM.xml')
cobra.io.write_sbml_model(model, 'Sco-GEM.xml')
```

***RAVEN Toolbox*** \* 
```matlab
model = importModel('Sco-GEM.xml')
exportModel(model, 'Sco-GEM.xml')
```

***COBRA Toolbox*** \*
```matlab
model = readCbModel('Sco-GEM.xml')
writeCbModel(model, 'Sco-GEM.xml')
```

\* note that some annotation might be lost when exporting the model from RAVEN and COBRA Toolboxes

### **Contributor:**

Development of the model is done via cobrapy, to ensure that model content is retained as much as possible (I/O through other software might result in undesired loss of annotation). Therefore, make sure you have Python 3.6+ installed.

[Fork](https://github.com/sysbiochalmers/sco-gem/fork) the Sco-GEM repository to your own GitHub account, and create a new branch from `devel`.

In Python, create an environment with all requirements:

```python
pip install -r requirements.txt  # installs all dependencies
touch .env                       # creates a .env file for locating the root, required
```

Load the model using the default code specified [above](#user). Export the model with the `export` function provided in the repository:
```python
cd ./code
import export
export.export(model, formats = ['xml', 'yml'])
```

More information on contributing to Sco-GEM can be found in the [contributing guidelines](.github/CONTRIBUTING.md), read these to get started. Contributions are always welcome!

### Contributors
* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Sweden
* [Snorre Sulheim](https://www.sintef.no/en/all-employees/employee/?empId=5675) ([@sulheim](https://github.com/sulheim)), SINTEF Industry, Norway / Norwegian University of Science and Technology, Norway
* [Tjasa Kumelj](https://www.ntnu.edu/employees/tjasa.kumelj) ([@tjasakumelj](https://github.com/tjasakumelj)), Norwegian University of Science and Technology, Norway
* [Ali Salehzadeh-Yazdi](https://www.sbi.uni-rostock.de/team/detail/ali-salehzadeh-yazdi) ([@alisalehzadeh-yazdi](https://github.com/alisalehzadeh-yazdi)), University of Rostock, Germany