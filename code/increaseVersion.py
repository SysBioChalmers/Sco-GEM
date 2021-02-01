# -*- coding: utf-8 -*-
"""
This script prepares the model and auxiliary files to make a new release, by
going through the following steps:

1. Check that the current git branch is 'master'
2. Check which type of version bump is desired ('major', 'minor' or 'patch')
3. Load the model from the current branch
4. Define new version number
5. Check that 'history.md' details the changes for the new version
6. Change the version number in the model file and export xml and yml formats
7. Write the new version number in version.txt
8. Update 'README.md' to contain updated model stats (number of reactions,
   metabolites and genes)

Example
=======
The script is run directly from the command line while in the master branch:

    $ python increaseVersion 'bumpType'

where 'bumpType' is either:
    'major'     e.g. increase from version 1.2.3 to 2.0.0
    'minor'     e.g. increase from version 1.2.3 to 1.3.0
    'patch'     e.g. increase from version 1.2.3 to 1.2.4

See contributing guidelines for direction on when which bumpType is appropriate.

"""
import cobra
import export
import argparse
import sys
import re
from dotenv import find_dotenv
import subprocess

# find .env + define paths
REPO_PATH = find_dotenv()
REPO_PATH = REPO_PATH[:-5]
MODEL_PATH = f"{REPO_PATH}/model"

def check_git_branch():
    # Make sure that current branch is master branch, otherwise exit.
    print("Checking that current branch is master...")
    branch_name = subprocess.run(["git","branch","--show-current"], stdout = subprocess.PIPE).stdout.decode("utf-8")
    if branch_name != "master":
        sys.exit("The local git branch is '{0}'. This function is only to increase the version of Sco-GEM in the 'master' branch.".format(branch_name))

def increase_version(model,bumpType):
    # Determine and print new version number.
    f = open(REPO_PATH + "/version.txt", "rt")
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
    # Write /version.txt with the new version number
    print("Update version.txt...")
    with open(REPO_PATH + "/version.txt", "w") as filetowrite:
        filetowrite.write(newVersion)

def check_history(newVersion):
    # Make sure that the history.md contains a line with the new version number,
    # suggesting that history.md has been updated.
    print("Check that history.md details the latest changes...")
    file = open(REPO_PATH + "/history.md","r")
    history = file.read()
    file.close
    if "Sco-GEM v" + newVersion + ":" not in history:
        sys.exit("history.md does not yet contain the latest changes. Copy these from the version-related pull-request (devel to master) into history.md, mentioning the new version.")

def update_modelstats(model):
    # Update some stats from the new model in README.md
    # TODO: grab memote score from Travis run from last PR to master, and also update here

    print("Update model stats in README.md...")
    # New string containing model stats
    new_string = r"\1 " + str(len(model.reactions)) + " | " + str(len(model.metabolites)) + " | " + str(len(model.genes)) + r" \2"
    
    f = open(REPO_PATH + "/README.md", "rt")
    data = f.read()
    data = re.sub(r"(coelicolor_ A3\(2\) \| iKS1317 \|) \d+ \| \d+ \| \d+ (\| .*\|)", new_string, data) # regex with new string, do not touch memote score
    f.close()
    
    f = open(REPO_PATH + "/README.md", "wt")
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
        model = cobra.io.read_sbml_model(MODEL_PATH + "/Sco-GEM.xml")
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
