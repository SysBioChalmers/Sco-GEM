# -*- coding: utf-8 -*-
"""
This script does:
Increase the version

Example
=======
The script can either be run directly from the command line::

    $ python export.py Sco-GEM.xml --folder model --formats xml yml txt

or using the module function::
    from export import export
    export(Sco_GEM, formats = ["xml", "yml", "txt"])

"""
import cobra
import export
import argparse
import sys
import re
from pathlib import Path

REPO_MAIN_FOLDER = Path(__file__).resolve().parent.parent

def check_git_branch():
    print("Checking that current branch is master...")
    branch_name = subprocess.run(["git","branch","--show-current"], stdout = subprocess.PIPE).stdout.decode("utf-8")
    if branch_name != "master":
        sys.exit("The local git branch is '{0}'. This function is only to increase the version of Sco-GEM in the 'master' branch.".format(branch_name))

def increase_version(model,bumpType):
    f = open(REPO_MAIN_FOLDER + "version.txt", "rt")
    oldVersion = f.read()
    f.close()
    
    oldVersion = model.id
    newVersion = oldVersion.split(".")
    newVersion = list(map(int, newVersion))
    if bumpType == "major":
        newVersion[0] += 1
        newVersion[1] = 0
        newVersion[2] = 0
    elif bumpType == "minor":
        newVersion[1] += 1
        newVersion[2] = 0
    elif bumpType == "patch":
        newVersion[2] += 1
    else:
        sys.exit("bumpType should be either 'major', 'minor' or 'patch'")
    
    newVersion = list(map(str, newVersion))
    s = "_"
    newVersionLabel = s.join(newVersion) # To format the model ID with underscores instead of dots
    model.id = "Sco_GEM_v" + s.join(newVersion)
    s = "."
    newVersion = s.join(newVersion) # To format the version with dots for printing and parsing
    print("Increasing from version {0} to version {1}".format(oldVersion, newVersion))
    return newVersion

def write_version_file(newVersion):
    print("Update version.txt...")
    with open(REPO_MAIN_FOLDER + "version.txt", "w") as filetowrite:
        filetowrite.write(newVersion)

def check_history(newVersion):
    print("Check that history.md details the latest changes...")
    file = open(REPO_MAIN_FOLDER + "history.md","r")
    history = file.read()
    file.close
    if "Sco-GEM v" + newVersion + ":" not in history:
        sys.exit("history.md does not yet contain the latest changes. Copy these from the version-related pull-request (devel to master) into history.md, mentioning the new version.")

def update_modelstats(model):
    # TODO: grab memote score from Travis run from last PR to master, and update in README.md

    print("Update model stats in README.md...")
    # New string containing model stats
    new_string = r"\1 " + str(len(model.reactions)) + " | " + str(len(model.metabolites)) + " | " + str(len(model.genes)) + r" \2"
    
    f = open(REPO_MAIN_FOLDER + "README.md", "rt")
    data = f.read()
    data = re.sub(r"(coelicolor_ A3\(2\) \| iKS1317 \|) \d+ \| \d+ \| \d+ (\| .*\|)", new_string, data) # regex with new string, do not touch memote score
    f.close()
    
    f = open(REPO_MAIN_FOLDER + "README.md", "wt")
    f.write(data)
    f.close()

if __name__ == "__main__":
    # Check that branch is "master"
    check_git_branch()

    # Parse argument, which type of version increase    
    parser = argparse.ArgumentParser(description= "Increase version of Sco-GEM model")
    parser.add_argument("bumpType", help = "string of either 'major', 'minor' or 'patch', indicating the type of version increase")
    args = parser.parse_args()
    bumpType = args.bumpType

    # Load existing model
    try:
        print("Loading model file...")
        model = cobra.io.read_sbml_model(REPO_MAIN_FOLDER / "model/Sco-GEM.xml")
    except OSError as e:
        print("Cannot find model/Sco-GEM in the repository")

    # Modify model to increase version
    newVersion = increase_version(model,bumpType)

    # Check if history file was updated
    check_history(newVersion)

    # Write model file
    export.export(model, formats = ["xml", "yml"], write_requirements = 1, objective = "BIOMASS_SCO_tRNA")

    # Update version.txt file
    # Only do this AFTER writing the model files, and not as part of increase_version routine, as version.txt should only be changed if the model has successfully exported.
    write_version_file(newVersion)

    # Update README.md with model statistics
    update_modelstats(model)


