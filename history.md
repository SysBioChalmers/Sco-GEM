# History

### Sco-GEM v1.3.0:

- fix:
  - correct biocyc annotation ZCAROTDH2 (closes #33)
  - UniProt ID, PFAM, Panther, GO term and refseq annotations as gathered from Uniprot (closes #44)
  - correct invalid KEGG metabolite IDs identified in Sco4 paper (closes #62)
  - gene names are gathered from the Sanger genome annotation, strepDB, UniProt and additional literature (closes #64)
  - directionality of ICDHyr (closes #87)
  - remove reaction with incorrect co-factor (G6PDH1b) (closes #88)
  - remove unnecessary anomers of carbohydrates, only keeping the generic forms (closes #89)
  - directionality of ferredoxin-NAPD reductase (closes #90)
  - remove incorrect in GPR for CAT reaction (closes #100)
  - curated duplicated reactions, involved GABTA, ABPYRATA, CA2abc1 and AMEt (closes #104)
  - correct and deprecate multiple reaction and metabolite annotations (closes #105)
  - remove reaction GND2, wrong co-factor (closes #106)
  - remove multiple annotations when only one is correct (KEGG, ec-code) (closes #111)
  - acetyl-CoA rewiring from iKS1317 (closes #127)
  - remove unused genes, not annotated to any reaction
  - remove inconsistent and unnecessary `<notes>` entries, with information already covered elsewhere
  - `export.py` recognizes `sbo` terms, adheres to [standard-GEM](https://github.com/MetabolicAtlas/standard-GEM), include subroutine to load previous releases of the model (`export.get_earlier_model_unversioned`) 
  - remove two unnecessary boundary metabolites (closes #136)

- feat:
  - protein sequence, length and mass are gathered from UniProt (closes #44)
  - `increaseVersion.py` function that should be run on master branch to prepare a new version release (closes #133)
  - template curation script `code/curation/vx_x_x.py`

- refactor:
  - renamed folders and reformatted README.md to adhere to [standard-GEM](https://github.com/MetabolicAtlas/standard-GEM)
  - all model files (`xml` and `yml` for now) are located in the `model/` folder (closes #110)
  - moved scripts, data and ec-models specific for Sulheim _et al._ (2020) to dedicated folders in /code and /data (closes #110)
  - separate requirements.txt for `code/sulheim2020/`, reducing the packages in `/requirements.txt`
  - I/O by latest cobrapy (0.20) adds zero charge for metabolites whose charge was previously not specified (hence, metabolite charges should be curated, see #79)

- doc:
  - add Zenodo batch (closes #15, closes #65)
  - point to GitHub Discussions instead of Gitter (closes #124)
  - update contributing guidelines, move to `.github` folder (closes #133)

- chore:
  - add model.id field (closes #128)

### Sco-GEM v1.2.1:

This is the version of the Sco-GEM model(s) and repository used in the iScience publication [Enzyme-Constrained Models and Omics Analysis of Streptomyces coelicolor Reveal Metabolic Changes that Enhance Heterologous Production](https://doi.org/10.1016/j.isci.2020.101525). 

This release contains mostly minor edits related to annotations and additional plotting functions used to make figures for the aforementioned paper.

### Sco-GEM v1.2.0:

The curation and addition of transport reactions are the main improvement from Sco-GEM v1.1 to Sco-GEM v1.2.0
Additionally, branches for plotting and analysing data have been incorporated. The basis for this work was done by Tjasa Kumelj, and only the scripting adding the changes to the model reconstruction pipeline was performed by Snorre Sulheim.


| Taxonomy | Template Model | Reactions | Metabolites| Genes | Memote score |
| ------------- |:-------------:|:-------------:|:-------------:|-----:|-----:|
| _Streptomyces coelicolor_ A3(2) | iKS1317 | 2612 | 2073 | 1777 | 77%|

The manuscript describing the model reconstruction and strain analysis of _Streptomyces coelicolor_ M145 and M1152 is now available on BioRxiv: https://doi.org/10.1101/796722

### Sco-GEM v1.1.0:

A minor release of Sco-GEM with updated Ec-ScoGEM models and updated random sampling results.

Basic stats
- 1613 genes
- 2540 reactions
- 2023 metabolites

Sco-GEM update
- Subsystem and pathway annotation of reactions
- Deletion of 10 reactions that didn't have any gene annotation
- Updated reaction bounds of 40 reactions, 38 of them by ATP-driven reactions which in general are assumed to be irreversible
- Added exchange reaction for Coelimycin P1

EcSco-GEM update
- Updated from Sco-GEM 1.1
- New random sampling data with 5000 iterations
- Updated GECKO version

### Sco-GEM v1.0.0:

First release of scoGEM, the community-developed genome-scale reconstruction of _Streptomyces coelicolor_.
The reconstruction is created by using [iKS1317](https://onlinelibrary.wiley.com/doi/full/10.1002/biot.201800180) as a template, additional genes, reactions and metabolites have been added from [Sco4](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006541) and [iAA1259](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4905-5)

Basic stats
- 1613 genes
- 2544 reactions
- 2017 metabolites
- Memote score of 73

New features
- Increased coverage compared to previous models
- Divided the biomass-reaction into pseudo-metabolites:
  - lipid
  - rna
  - dna
  - cell wall
  - protein
  - carbohydrate
  - misc
- Added SBO annotations to all genes, reactions and metabolites
- Curated biomass reaction with respect to prosthetic groups
- Improved coverage and curated annotation of 
  - genes: GO, UniProt, AA-sequence, mass, EC-number, PANTHER (all in separate csv-file)
  - reactions: KEGG, BioCyc, MetaNetX
  - metabolites: KEGG, BioCyc, CHEBI, MetaNetX
- Reviewed and curated reaction bounds by using calculations of the change in Gibbs free energy from [eQuilibrator](http://equilibrator.weizmann.ac.il/)
- Created an enzyme-constrained model from scoGEM 1.0 by using [GECKO](http://msb.embopress.org/content/13/8/935.long)
  - Generated generic ec-models with protein pools for both strain M145 and M1152
  - Generated condition-specific proteome-constrained ec-models for both M145 and M1152
