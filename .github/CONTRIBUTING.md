## Contributor guidelines

First of all, thank you for contributing to Sco-GEM!! Anybody is welcome to contribute, but please abide by the following guidelines.

You can contribute in 2 main ways: by creating issues, and by sending pull requests (PRs) with additions, deletions, corrections, etc. to the model.

### Reporting issues in the model

Report an issue at https://github.com/SysBioChalmers/Sco-GEM/issues if you note any of the following:

* Incorrect annotation for any model component.
* Missing feature or field you would like the model to have.
* Bug/weird simulation results.
* Lacking documentation.
* Any type of feedback.

If you are unsure about the issue, consider asking first in the [Discussions](https://github.com/SysBioChalmers/Sco-GEM/discussions) forum.

When creating the issue, please make sure:

* You tested your code (if any) with all requirements for running the model.
* You did your analysis in the `master` branch of the repository.
* You provide any necessary files/links needed for understanding the issue.
* You checked that a similar issue does not exist already

Feel free to also comment on any of the [open issues](https://github.com/SysBioChalmers/Sco-GEM/issues).

Finally, if you like Sco-GEM please remember to 'star' our Github page (click on the star at top right corner), that way we also have an idea of who is using Sco-GEM!

### Contributing to the model

Want to contribute to the model with some additions or improvements? Consider starting by raising an issue and assign it to yourself to describe what you want to achieve. This way, we reduce the risk of duplicated efforts and you may also get suggestions on how to best proceed, e.g. there may be half-finished work in some branch that you could start with. Also, feel free to browse our [open issues](https://github.com/SysBioChalmers/Sco-GEM/issues): anything tagged with "help wanted" is open to whoever wants to implement it!

Here's how to set up Sco-GEM for local development to contribute smaller features or changes that you can implement yourself:

1. Firstly, make sure that you have all [requirements](https://github.com/SysBioChalmers/Sco-GEM#contributor) for contributing to Sco-GEM. Note that model development is done with cobrapy, to ensure that model content is retained as much as possible (I/O through other software might result in undesired loss of annotation).

2. [Fork](https://github.com/SysBioChalmers/Sco-GEM/fork) the Sco-GEM repository on GitHub.

3. Clone your fork locally:
    ```
    git clone https://github.com/<your Github name>/Sco-GEM.git
    ```

4. Check out the branch that you want to contribute to. Most likely that will be ``devel``:
    ```
    git checkout devel
    ```

5. Create a branch for local development based on the previously checked out branch ([see below](#branching-model) for details on the branching model and how to name your branch):
    ```
    git checkout -b name-of-your-branch
    ```

6. Now you can make your changes locally!
    * Never directly edit the model files, but make changes through cobrapy and document them in scripts.
    * A **template script to use** is available from **`code/curation/vx.x.x_template.py`**. Modify this script, define subroutines for your curations, or refer to other scripts in the repository if required. Particularly if larger curations or changes are made to the model it is probably more suitable to define and call dedicate scripts. Rename the template script to the intended new model version, this can be changed later if another version number is decided.
    * As example, see [`code/curation/v1.3.0.py`](https://github.com/SysBioChalmers/Sco-GEM/blob/master/code/curation/v1.3.0.py).
    * Scripts are placed in `code/` and data (as `.tsv` or `.csv`) are placed in `data/`, using relevant subfolders. Note that binary data such as `.mat` structures or `.xls` tables cannot be stored in the repo (as they cannot be version-controlled, and they increment too much the size of the repo).
    * For new metabolites and/or reactions, use the BiGG nomenclature for identifiers whenever possible ([metabolites](http://bigg.ucsd.edu/universal/metabolites), [reactions](http://bigg.ucsd.edu/universal/reactions), [download](http://bigg.ucsd.edu/data_access)). If no BiGG identifier is specified for the new metabolite/reaction, then clearly indicate this in the script. 
    * When you are done making changes, review locally your changes with `git diff` or any git client, to make sure you are modifying the model as you intended.

7. Commit your changes and push your branch to GitHub.
    ```
    git add .
    git commit -m "Title of your commit"
    git push origin name-of-your-branch
    ```
    [See below](#semantic-commits) for recommendations on how to name your commits. In case of larger updates, you can of course make several commits on a single contribution. However, if you need to do too many commits, consider if your contribution could be instead split into separate branches (making it easier for reviewing later).

8. Submit a pull request through the GitHub website (https://help.github.com/articles/creating-a-pull-request-from-a-fork/) to the `devel` branch of the original SysBioChalmers repo (not your fork). We recommend ticking the box "Allow edits from maintainers" if you wish for us to be able to contribute directly to your branch (speeding-up the reviewing process).

Note that steps 3, 4, 5 and 7 can be done, if you prefer, with any git client, such as [Github Desktop](https://desktop.github.com/), [Fork](https://git-fork.com/) or [GitKraken](https://www.gitkraken.com/git-client).

Finally, and for larger features that you want to work on collaboratively, you may consider to first request to join our development team to get write access to the repository so that you can create a branch directly in the main repository (or simply ask the administrator to create a branch for you). Once you have a new branch, you can push your changes directly to the main repository and when finished, submit a pull request from that branch to ``devel``. [See below](#development-team-guidelines) for more details.

Thank you very much for contributing to Sco-GEM!

#### Branching model

* `devel`: Is the branch all pull-requests should be based on.

* `master`: Is only touched by the administrator and is the branch with the tested & reviewed model that is released or ready for the next release.

* `gh-pages`: Is only touched by the administrator and for maintaining the [landing page](http://sysbiochalmers.github.io/Sco-GEM/) of Sco-GEM.

* `{chore, doc, feat, fix, refactor, style}/descriptive-name`: Any other branch created in the model. If you work on a fix, start the branch name with `fix/`, if you work on a feature, start the branch name with `feat/`. Examples: `fix/format_reactions` or `feat/new_algorithms`. [See below](#semantic-commits) for more details on the possible actions you can use.

#### Semantic commits

Please use concise descriptive commit messages. Ideally, use semantic commit messages to make it easier to show what you are aiming to do:

`action: brief description`

`action` refers to what exactly are you doing in the commit, following a [standard definition](http://karma-runner.github.io/2.0/dev/git-commit-msg.html) in software development:
* `chore`: updating toolbox, data files, etc.
* `doc`: updating documentation or explanatory comments in functions.
* `feat`: new feature added, e.g. new reaction / metabolite / function / etc.
* `fix`: something that was incorrect in the model and now has been corrected.
* `refactor`: see [code refactoring](https://en.wikipedia.org/wiki/Code_refactoring).
* `style`: minor format changes of model, functions or data (spaces, semi-colons, etc., no code change).

Some examples:

|commit|commit message|
|:---:|:---:|
|Add new rxns|`feat: add methanol pathway`|
|Remove a metabolite|`fix: remove duplicated citrate`|
|Correct metabolite formula|`fix: correct carbohydrate formulas`|
|Fix rxn stoichiometry|`fix: complex V stoich coeffs`|
|Update gene IDs|`fix: update gene IDs from SWISSPROT`|
|Format name of compartment|`style: remove uppercases compartment name`|
|Split a rxn in 2|`refactor: split isomerase in 2 steps`|
|Add some data|`feat: metabolomics data`|
|Update documentation of function|`doc: improved comments v1_3_0.py`|
|Update toolbox|`chore: update cobrapy version`|

More examples [here](https://github.com/SysBioChalmers/Sco-GEM/commits/master). A more detailed explanation or comments is encouraged to be left in the commit description.

## Development team guidelines

This section is meant for the development team of Sco-GEM. As a member of the development team, you should comply with all [previous contributor guidelines](#contributor-guidelines) as well. Besides that, please also consider the following guidelines.

### Creating pull requests

Changes to the model _cannot_ be directly committed to the `master` or `devel` branches (in fact they are protected). Commits are made to side-branches, after which pull requests are made for merging with `devel`. For this, follow the [instructions](#contributing-to-the-model) for contributors, but consider that as members of the development team have write access to the repository, you can create a branch directly in the main repository without needing to fork, for your convenience. This means that you can:

* Skip step 2 of the contribution process.
* In step 3 of the contribution process, clone directly the original repo:
  ```
  git clone https://github.com/SysBioChalmers/Sco-GEM.git
  ```

Follow all other steps in the same way. Also, when creating your pull request (or after you have created it):
* Choose 1 or more members of the team (ideally with knowledge on the pull request) as reviewers. Note that the person making the pull request and the reviewer _cannot_ be the same person.
* Assign appropriate [labels](https://github.com/SysBioChalmers/Sco-GEM/issues/labels).


### Reviewing pull requests

Every pull request must be approved by at least one reviewer before it can be merged. When reviewing someone else's pull request, keep in mind the following aspects:
* **Compatibility:** First of all, make sure that the model is still compatible with the loading/saving wrappers (`loadYeastModel.m` & `saveYeastModel.m`) and that no errors appear. Check also that [`requirements.txt`](https://github.com/SysBioChalmers/Sco-GEM/blob/master/requirements.txt) does not change in any unexpected ways (e.g. an "unknown" toolbox version). Finally, ensure that the SBML fields `model metaid`, `model id` and `model name` never change, as if they change it would create a conflict in the next release.
* **Documentation:** Every change should be justified with a reference/link/argument. This can be provided as data in `/data`, or directly as a comment in the pull request.
* **Reproducibility:** Run the script that can update the previous model release to the next intended model version, by checking out the branch from where the pull request is made. If the script runs without problems, then check with `git diff` or similar whether any additional changes were introduced.
* **Style:** Ensure that the changes to the model are compliant with the model's rxn/met/gene naming conventions (when unsure, take a look at a similar field in the model). Also, make sure that scripts have a compliant style, and datasets are straight-forward to understand.
* Avoid vague comments and try to be as explicit as possible (e.g.: _"Please include X here"_ instead of _"X could be included here"_).

### Releasing a new version

* A merge of `devel` with `master` invokes a new release.
* A new release should be made as soon as there is substantial new work in `devel` (as rule of thumb, after around 3 pull request merges).

Sco-GEM follows [semantic versioning](https://semver.org/), adapted to GEMs:
* A `major` release is seldom used and only meant for a new publication. Backwards compatibility should be, ideally, always preserved.
* A `minor` release involves a substantial change in the model (several new reactions/metabolites/genes), such as:
  * Addition of genes/reactions/metabolites from a whole genome annotation.
  * Addition of several annotation fields.
  * Inclusion of a major new formalism in the model.
  * Addition of a plurality of pathways.
* A `patch` release is the most common one and is done when only few things have changed in the model, or there are only changes that have to do with format, such as:
  * Adding a single new annotation field.
  * Fixing some chemical formulas/charges.
  * Updating toolboxes.
  * Re-organization of data
  * Refactoring of code.

When releasing, please follow these steps:
  1. Create a pull request from `devel` to `master`, indicating all new features/fixes/etc. and referencing every previous pull request included (examples [here](https://github.com/SysBioChalmers/Sco-GEM/releases)). Tip: if any [issue](https://github.com/SysBioChalmers/Sco-GEM/issues) gets solved in the release, write in the pull request description "Closes #X", where "X" is the issue number. That way the issue will be automatically closed after merge.
  2. Once the pull request is reviewed and accepted, merge to `devel` to `master`.
  3. Switch locally to `master`, pull changes and update `history.md`, by putting at the top the same description of the corresponding pull request from step 1.
  4. Bump version with `code/increaseVersion.py`. Run as `python increaseVersion.py 'bumpType'` where 'bumpType' is either 'major', 'minor', or 'patch'. **NOTE:** The function will error if unexpected changes are occurring. If this happens, probably step 1 was done incorrectly. To fix it, commit in `devel` any necessary changes and make a new pull request.
  5. Commit changes from steps 3 and 4 with the message `chore: version X.Y.Z`, and push to the remote.
  6. Make the new release at GitHub [here](https://github.com/SysBioChalmers/Sco-GEM/releases/new), using the proper tag "vX.Y.Z" and with the same description as the corresponding pull request from step 1.
  7. Review the [Zenodo](https://zenodo.org) release: every new release from Github (step 6) automatically triggers a new release in Zenodo. However, to be sure check that the new release is made available [here](https://zenodo.org/badge/latestdoi/145685631). Note that it might take some minutes for the Zenodo release to appear after you create the release in Github.

## Previous work

``sco-GEM`` leverages the hard labour that has previously been performed in the development of genome-scale models of *Streptomyces coelicolor*, as published in the following papers:

- Borodina I, Krabben P, Nielsen J. _Genome Res_. **2005**;15: 820–9. [doi](http://doi.org/10.1101/gr.3364705)
- Alam MT, Merlo ME, Hodgson DA, Wellington EMH, Takano E, Breitling R. _BMC Genomics_. **2010**;11: 202. [doi](http://doi.org/10.1186/1471-2164-11-202)
- Kim M, Sang Yi J, Kim J, Kim J-N, Kim MW, Kim B-G. Biotechnol J. **2014**;9: 1185–94. [doi](http://doi.org/doi:10.1002/biot.201300539)
- Amara A, Takano E, Breitling R. _BMC Genomics_. **2018**;19: 519. [doi](https://doi.org/10.1186/s12864-018-4905-5)
- Wang H, Marcišauskas S, Sánchez BJ, Domenzain I, Hermansson D, Agren R, Nielsen J, Kerkhoven EJ. _PLOS Comput Biol_. **2018**;14: e1006541. [doi](http://doi.org/10.1371/journal.pcbi.1006541)
- Sulheim S, Kumelj T, Wentzel A, Almaas E. _Biotech J._ **2018**;14: 1800180. [doi](https://doi.org/10.1002/biot.201800180)
- Sulheim S, Kumelj T, van Dissel D, Salehzadeh-Yazdi A, Du C, van Wezel GP, Nieselt K, Almaas E, Wentzel A, Kerkhoven EJ. _iScience_ **2020**;23: 101525. [doi](https://doi.org/10.1016/j.isci.2020.101525)

## Acknowledgments

These contribution guidelines were written based on the contribution guidelines of [SysBioChalmers/yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM/.github/CONTRIBUTING.md).

