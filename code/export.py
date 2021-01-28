# -*- coding: utf-8 -*-
"""
This script does:
1) Act as a wrapper for export functions in cobrapy
2) Ensure that the model entities are sorted so that the changes in model files are minimal

Export functions
========================

keys          | filetype
-------------------------
xml, sbml     | SBML 3 v1, fbc v2 (standard)
yml, yaml     | yaml-file
json          | json-file
mat, matlab   | matlab-file
txt, csv, tsv | todo

# It is possible to get different SBML-versions by using the *sbml_level* and *sbml_version* 
# kwargs in export-function


Order of model entities
=======================
The model entities are presented in the following order in the model files,
# sorted alphabetically within each category unless explicitly stated *unsorted*

1. List of units (unsorted)
2. List of objectives (unsorted)
3. List of parameters (unsorted)
4. List of metabolites
    a. Biomass (SBO:0000649)
    b. Others (SBO:0000247)
5. List of reactions
    a. Biomass  (SBO:0000629)
    b. ATP mainteinance (SBO:0000630)
    c. Encapsulating  (SBO:0000395)
    d. Pseudo (SBO:0000631)
    e. Biochemical [normal] (SBO:0000176)
    f. Transport (SBO:0000655)
    g. Demand (SBO:0000628)
    h. Exchange (SBO:0000627)
    i. Sink (SBO:0000632)
    j. Other 
    k. Non-sbo  

Example
=======
The script can either be run directly from the command line::

    $ python export.py Sco-GEM.xml --folder model --formats xml yml txt

or using the module function::
    from export import export
    export(Sco_GEM, formats = ["xml", "yml", "txt"])

"""
import cobra
from collections import OrderedDict
from os import remove
from sys import exit
import argparse
import subprocess
from dotenv import find_dotenv
from re import match

# find .env + define paths
REPO_PATH = find_dotenv()
REPO_PATH = REPO_PATH[:-5]
MODEL_PATH = f"{REPO_PATH}/model"

def export(Sco_GEM, folder = "model", name = "Sco-GEM", formats = ["xml", "yml", "txt"], 
           write_requirements = True, objective = None):
    """
    Sort the model entities and export the model in given formats

    Parameters
    ----------
    Sco_GEM: cobra.core.model.Model
        The model file
    folder: str, optional
        The subfolder in the repository in which all the model files are stored, relative to the root dir of the repository (default: "model")
    name: str, optional
        Model file name (default: "Sco-GEM")
    formats: list of str, optional
        The formats the model should be exported in (default: ["xml", "yml", "txt"])
    use_fbc_package: bool, optional
        If the fbc-package should be used to store the sbml-file (default: True)
    objective: str, optional
        The id of the objective function, default BIOMASS_SCO

    """
    if objective:
        set_objective(Sco_GEM, objective)
    sort_model(Sco_GEM)

    # write xml
    if ("xml" in formats) or ("sbml" in formats):
        model_fn_xml = MODEL_PATH + "/{0}.xml".format(name)
        print("Writing {0}".format(str(model_fn_xml)))
        cobra.io.write_sbml_model(Sco_GEM, str(model_fn_xml))
    
    if ("yml" in formats) or ("yaml" in formats):
        model_fn_yml = MODEL_PATH + "/{0}.yml".format(name)
        print("Writing {0}".format(str(model_fn_yml)))
        cobra.io.save_yaml_model(Sco_GEM, str(model_fn_yml))

    if "json" in formats:
        model_fn_json = MODEL_PATH + "/{0}.json".format(name)
        print("Writing {0}".format(str(model_fn_json)))
        cobra.io.save_json_model(Sco_GEM, str(model_fn_json))
        
    if "txt" in formats:
        print("Can't print txt-files yet. Skipping")
        # raise NotImplementedError
        # model_fn_txt = MODEL_PATH + "/{0}.txt".format(name)
        # print("Writing {0}".format(str(model_fn_txt)))
        # cobra.io.write_sbml_model(Sco_GEM, str(model_fn_txt))

    if "mat" in formats:
        model_fn_mat = MODEL_PATH + "/{0}.mat".format(name)
        print("Writing {0}".format(str(model_fn_mat)))
        cobra.io.save_matlab_model(Sco_GEM, str(model_fn_mat))

    if write_requirements:
        write_requirements_fun(REPO_PATH)

def set_objective(Sco_GEM, r_id):
    Sco_GEM.reactions.BIOMASS_SCO.objective_coefficient = 0
    Sco_GEM.reactions.get_by_id(r_id).objective_coefficient = 1


def sort_model(Sco_GEM):
    # Sort metabolites
    Sco_GEM.metabolites.sort()
    
    # Sort compartments
    # ToDo: This is currently not working
    # Sco_GEM.compartments = sort_dict(Sco_GEM.compartments)

    # Sort genes
    Sco_GEM.genes.sort()


    # Sort reactions
    reaction_order_dict = sort_reactions(Sco_GEM)
    Sco_GEM.reactions.sort(key = lambda r: reaction_order_dict[r.id])

def sort_reactions(Sco_GEM):

    REACTION_SBO_TERMS = {
    "exchange reaction":     "SBO:0000627",
    "demand reaction":       "SBO:0000628",
    "biomass production":    "SBO:0000629",
    "ATP mainteinance":      "SBO:0000630",
    "encapsulating process": "SBO:0000395",
    "transport reaction":    "SBO:0000655",
    "biochemical reaction":  "SBO:0000176",
    "sink reaction":         "SBO:0000632",
    "pseudoreaction":        "SBO:0000631",
}
    reaction_order_dict = {}
    i = 0

    biomass_reactions = []
    pseudo_reactions = []
    normal_reactions = []
    transport_reactions = []
    demand_reactions = []
    exchange_reactions = []
    sink_reactions = []
    non_sbo_reactions = []
    atpm_reactions = []
    encapsulating_reactions = []
    other_reactions = []
    for r in Sco_GEM.reactions:
        try:
            sbo = r.annotation["sbo"]
        except KeyError:
            non_sbo_reactions.append(r.id)
            continue
        if sbo == REACTION_SBO_TERMS["exchange reaction"]:
            exchange_reactions.append(r.id)
        elif sbo == REACTION_SBO_TERMS["demand reaction"]:
            demand_reactions.append(r.id)
        elif sbo == REACTION_SBO_TERMS["biomass production"]:
            biomass_reactions.append(r.id)
        elif sbo == REACTION_SBO_TERMS["ATP mainteinance"]:
            atpm_reactions.append(r.id)
        elif sbo == REACTION_SBO_TERMS["encapsulating process"]:
            encapsulating_reactions.append(r.id)
        elif sbo == REACTION_SBO_TERMS["transport reaction"]:
            transport_reactions.append(r.id)
        elif sbo == REACTION_SBO_TERMS["biochemical reaction"]:
            normal_reactions.append(r.id)
        elif sbo == REACTION_SBO_TERMS["sink reaction"]:
            sink_reactions.append(r.id)
        elif sbo == REACTION_SBO_TERMS["pseudoreaction"]:
            pseudo_reactions.append(r.id)
        else:
            other_reactions.append(r.id)

    all_reaction_lists = [biomass_reactions, atpm_reactions, encapsulating_reactions, pseudo_reactions, 
                          normal_reactions, transport_reactions, demand_reactions, exchange_reactions, sink_reactions,
                          other_reactions, non_sbo_reactions]
    for reaction_list in all_reaction_lists:
        for r_id in sorted(reaction_list):
            reaction_order_dict[r_id] = i
            i += 1
    return reaction_order_dict


def write_requirements_fun(directory, force = True):
    directory = str(directory)
    if force:
        subprocess.run(["pipreqs", "--force", "."], cwd = directory)
    else:
        subprocess.run(["pipreqs", "."], cwd = directory)

def sort_dict(model_dict, by = "key"):
    sorted_dict = OrderedDict()
    if by == "key":
        for key in sorted(model_dict.keys()):
            sorted_dict[key] = model_dict[key]
    else:
        for key in sorted(model_dict, key = model_dict.get):
            sorted_dict[key] = model_dict[key]
    return sorted_dict

def get_earlier_model_unversioned(version):
    # version is either 'master' for the latest model on master branch, or a
    # tag that corresponds to a particular release, e.g. 'v1.2.3'.
    if version == 'master':
        git_result = subprocess.run(["git","show","master:ModelFiles/xml/Sco-GEM.xml"], stdout = open("_earlierModel.xml", "w")) # change "ModelFiles/xml" to "model" before PR to master
    elif match(r"v\d+.\d+.\d",version):
        tagpath = "refs/tags/" + version + ":ModelFiles/xml/Sco-GEM.xml" # change "ModelFiles/xml" to "model" before PR to master
        git_result = subprocess.run(["git","show",tagpath], stdout = open("_earlierModel.xml", "w"))
    else:
        exit("Unclear which earlier model version should be loaded.")
    
    model = cobra.io.read_sbml_model("_earlierModel.xml")
    remove("_earlierModel.xml")
    model.id = "Sco_GEM" # Remove versioning in model ID
    return model

if __name__ == '__main__':
    if 0:
        Sco_GEM_FN = MODEL_PATH + "Sco-GEM.xml"
        Sco_GEM = cobra.io.read_sbml_model(Sco_GEM_FN)
        export(Sco_GEM)
    else:
        parser = argparse.ArgumentParser(description= "Export the Sco-GEM model")
        parser.add_argument("model_path", help = "Path to SBML (xml) model file")
        parser.add_argument("--formats", nargs='+', help = "The different file formats to export the model in", default = ["xml", "yml", "txt"])
        parser.add_argument("--name", type = str, help = "General name of exported model files", default = "Sco-GEM")
        parser.add_argument("--folder", type = str, help = "The general subfolder relative to the repo's root directory to store the model files in", default = "model")
        parser.add_argument("--objective", type = str, help = "The reaction ID of the objective function", default = "BIOMASS_SCO_tRNA")
     
        args = parser.parse_args()
        cwd = Path.cwd()

        model_path = str(cwd / args.model_path)
        try:
            Sco_GEM = cobra.io.read_sbml_model(model_path)
        except OSError as e:
            print("Cant find model file: {0}".format(model_path))
            raise

        else:
            if args.sbml_version != 1:
                print("Currently, only SBML-version 1 is supported. sbml_version parameter ignored.")
            sbml_version = 1
                
            export(Sco_GEM, args.folder, args.name, args.formats)
