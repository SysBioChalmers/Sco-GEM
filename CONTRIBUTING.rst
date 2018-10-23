.. highlight:: shell

==============================================================
Contributing
==============================================================
Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

**[TODO: mention Gitter room when implemented]**

You can contribute in many ways, but please follow the guidelines given below.

--------------------------------------------------------------

Types of Contributions
--------------------------------------------------------------

Report problems or suggestions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Report problems with the metabolic model at https://github.com/edkerk/scoGEM/issues.

If you are reporting a problem, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting, such as COBRA/RAVEN/cobrapy version.
* Detailed steps to reproduce the problem.

Modify the model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Besides reporting issues, you are of course welcome to directly contribute to development of the model by making modifications yourself. This can be related to your own issues, but you are also welcome to contribute to any of the other issues that are raised by others at the [`GitHub issues <https://github.com/edkerk/scoGEM/issues>`_] page. Please make sure you read the Contribution Guidelines below.

--------------------------------------------------------------

Contribution Guidelines
--------------------------------------------------------------

Getting started
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Regular developers can work directly on the main repo and should clone the ``scoGEM`` repo locally using your favourite Git software::

    $ git clone git@github.com:edkerk/scoGEM.git

2. If you've not been added as regular contributor, you can still fork the ``scoGEM`` repo on GitHub and clone your fork locally::

    $ git clone git@github.com:your_name_here/scoGEM.git

3. Regardless, develop the model using the guidelines below, make pull requests to the main repo using the guidelines below.

Software
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To ensure consistency between model versions, allowing for easy diffing and merging of changes, all modifications of ``scoGEM`` should be for the moment be performed using `cobrapy <https://opencobra.github.io/cobrapy/>`_. If you're not fluent in using `cobrapy` but still want to directly submit changes to the model, please contact the repo maintainers on Gitter **[LINK]**.

Branches
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The ``master`` branch is protected, which means that direct commits are not allowed. Generally all pull requests should come from the ``devel`` branch and should be carefully reviewed before merged. **[Mention memote results]**. An update in the ``master`` branch triggers an update in model version.
* The ``devel`` branch is protected, which means that direct commits are not allowed. All pull requests should come from semantically named feature branches.
* When starting development, either contribute to an existing (relevant) feature branch, or make a new branch that is connected to the latest commit on the ``devel`` branch.
* Try to avoid merging feature branches. Rather, make a pull request to ``devel`` for the first branch and once accepted then merge ``devel`` with the second branch. This way, pull requests only contain changes from one branch and should be more focused. If needed, temporary `test` branches can be used to evaluate the effect of merging different feature branches.

Semantic naming of feature branches, commits and pull requests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Every feature branch name, commit title and pull request title should follow semanting naming containing an "action" and "keyword/phrase". An "action" is one of the following:

========  ==============================  ============================================================
Action    Type of changes                 Examples
========  ==============================  ============================================================
feat      new feature added               new annotations, new reactions, new dataset
fix       fixing of bugs, errors          remove duplicate metabolites, fix reaction reversibility
refactor  code refactoring                changes in scripts that don't change output
style     style/formatting of files       
docs      documentation                   update readme, this file
chore     maintenance chores              export model with new version of cobrapy
========  ==============================  ============================================================

* Format of feature branch name: ``action/keyword``, where keyword is descriptive, but short and has no spaces

  * Examples: ``feat/newPathway`` or ``fix/geneAssocation`` or ``docs/updateReadme``

* Format of commit and pull request titles: ``action: keyphrase``, where keyphrase is a short descriptive phrase

  * Examples: ``feat: new antibiotic pathway`` or ``fix: correct genes in reaction FBP`` or ``docs: attribute new developers``
  * Besides the title, commits and pull requests have descriptions that are more descriptive of the changes, and should use the templates that are automatically provided when making a new commit or pull request.

Reproducible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To strive for reproducible science, ``scoGEM`` aims to provide scripts describing all changes to the model. Such scripts are stored in the `ComplementaryScripts` folder (detailed below), in a relevant subfolder. Please pay attention to existing scripts to get an idea of how these scripts should look like. Any data files required for this, such as tables with new annotations, are stored in the relevant subfolder of `ComplementaryData`.

Model consistency and conventions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- For metabolite and reaction IDs, ``scoGEM`` aims to use `BIGG <http://bigg.ucsd.edu>`_ IDs were possible, as they are human readable, relatively well defined and compatible with e.g. SBML, short and memorable. If no BIGG ID exists, ``scoGEM`` uses BIGG-ish IDs.
- ... (other conventions being followed?)

Directory structure, file formats and names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The repo contains of a number of folders and subfolders, as specified here:

* ComplementaryData: contains data files required for model curation and simulation, all in CSV or tab-delimited format
  
  - curation: files required for model curation
  - essentiality: data for essentiality testing by memote
  - growth: experimentally measured growth data
  - media: details on media compositions for testing by memote
  - models: earlier *S. coelicolor* models that were used to construct the consensus ``scoGEM``

* ComplementaryScripts: contains scripts required for model curation and simulation, in either python or 
  
  - consensusModel: scripts to generate the consensus model ``scoGEM`` version 1.0
  - gecko: scripts to generate the enzyme-constrained version ``ecScoGEM``, using the `Gecko toolbox <https://github.com/SysBioChalmers/GECKO>`_

* ModelFiles

  - txt: ``scoGEM`` in text format, to facilitate diffing changes between models, automatically generated by cobrapy
  - xml: ``scoGEM`` in SBML L3V1 FBCv2 format as stored by cobrapy, ready for use in simulations by any other SBML-compatible software package, or further curation
  - yml: ``scoGEM`` in YAML format, to facilitate diffing changes between models, automatically generated by cobrapy

--------------------------------------------------------------

Contributors
--------------------------------------------------------------
-
-
-

In addition, ``scoGEM`` leverages the hard labour that has previously been performed in the development of genome-scale models of *Streptomyces coelicolor*, as published in the following papers:

- Borodina I, Krabben P, Nielsen J. Genome Res. 2005;15: 820–9. `doi <http://doi.org/10.1101/gr.3364705>`_
- Alam MT, Merlo ME, Hodgson DA, Wellington EMH, Takano E, Breitling R. BMC Genomics. 2010;11: 202. `doi <http://doi.org/10.1186/1471-2164-11-202>`_
- Kim M, Sang Yi J, Kim J, Kim J-N, Kim MW, Kim B-G. Biotechnol J. 2014;9: 1185–94. `doi <http://doi.org/doi:10.1002/biot.201300539>`_
- Amara A, Takano E, Breitling R. BMC Genomics. 2018;19: 519. `doi <https://doi.org/10.1186/s12864-018-4905-5>`_
- Wang H, Marcišauskas S, Sánchez BJ, Domenzain I, Hermansson D, Agren R, Nielsen J, Kerkhoven EJ. bioRxiv. 2018; 321067. `doi <http://doi.org/10.1101/321067>`_ **[change to PLOS Comp Biol]**
- Sulheim S, Kumelj T, Wentzel A, Almaas E. **[further details]**
