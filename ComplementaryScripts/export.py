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
    f. transport (SBO:0000655)
    g. demand (SBO:0000628)
    h. exchange (SBO:0000627)
    i. sink (SBO:0000632)
    j. other 
    k. non-sbo  

Example
=======
The script can either be run directly from the command line::

    $ python export.py scoGEM.xml -- folder ModelFiles --formats xml yml txt

or using the module function::
    from export import export
    export(scoGEM, formats = ["xml", "yml", "txt"])


"""
from consensusModel.fix_SBO_terms import REACTION_SBO_TERMS

import cobra
from collections import OrderedDict
from pathlib import Path
import argparse

REPO_MAIN_FOLDER = Path(__file__).resolve().parent.parent



def export(scoGEM, folder = "ModelFiles", name = "scoGEM", formats = ["xml", "yml", "txt"], sbml_level = 3, sbml_version = 1, use_fbc_package = True):
    """
    Sort the model entities and export the model in given formats

    Parameters
    ----------
    scoGEM: cobra.core.model.Model
        The model file
    folder: str, optional
        The folder in which all the model files are stored (default: "ModelFiles")
    name: str, optional
        Model file name (default: "scoGEM")
    formats: list of str, optional
        The formats the model should be exported in (default: ["xml", "yml", "txt"])
    sbml_level: int, optional
        The sbml level of the xml file (default: 3)
    sbml_version: int, optional
        The version of the sbml level (default: 1)
    use_fbc_package: bool, optional
        If the fbc-package should be used to store the sbml-file (default: True)

    """
    sort_model(scoGEM)

    main_folder_path = REPO_MAIN_FOLDER / folder
    # write xml
    if ("xml" in formats) or ("sbml" in formats):
        if sbml_level == 3:
            model_fn_xml = main_folder_path / "xml" / "{0}.xml".format(name)
            check_folder(model_fn_xml)
            print("Writing {0}".format(str(model_fn_xml)))
            cobra.io.write_sbml_model(scoGEM, str(model_fn_xml), use_fbc_package = use_fbc_package)
        else:
            model_fn_xml = main_folder_path / "xml" / "{0}_lvl2.xml".format(name)
            check_folder(model_fn_xml)
            print("Writing {0}".format(str(model_fn_xml)))
            cobra.io.write_legacy_sbml(scoGEM, str(model_fn_xml), use_fbc_package = use_fbc_package)
    
    if ("yml" in formats) or ("yaml" in formats):
        model_fn_yml = main_folder_path / "yml" / "{0}.yml".format(name)
        check_folder(model_fn_yml)
        print("Writing {0}".format(str(model_fn_yml)))
        cobra.io.save_yaml_model(scoGEM, str(model_fn_yml))

    if "json" in formats:
        model_fn_json = main_folder_path / "json" / "{0}.json".format(name)
        check_folder(model_fn_json)
        print("Writing {0}".format(str(model_fn_json)))
        cobra.io.save_json_model(scoGEM, str(model_fn_json))
        
    if "txt" in formats:
        print("Can't print txt-files yet. Skipping")
        # raise NotImplementedError
        # model_fn_txt = main_folder_path / "txt" / "{0}.txt".format(name)
        # check_folder(model_fn_txt)
        # print("Writing {0}".format(str(model_fn_txt)))
        # cobra.io.write_sbml_model(scoGEM, str(model_fn_txt))

    if "mat" in formats:
        model_fn_mat = main_folder_path / "mat" / "{0}.mat".format(name)
        check_folder(model_fn_mat)
        print("Writing {0}".format(str(model_fn_mat)))
        cobra.io.save_matlab_model(scoGEM, str(model_fn_mat))


def check_folder(file_path):
    if not file_path.parent.is_dir():
        file_path.parent.mkdir(parents = True, exist_ok = True)


def sort_model(scoGEM):
    # Sort metabolites
    scoGEM.metabolites.sort()
    
    # Sort compartments
    # ToDo: This is currently not working
    # scoGEM.compartments = sort_dict(scoGEM.compartments)

    # Sort genes
    scoGEM.genes.sort()


    # Sort reactions
    reaction_order_dict = sort_reactions(scoGEM)
    scoGEM.reactions.sort(key = lambda r: reaction_order_dict[r.id])



def sort_reactions(scoGEM):
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
    for r in scoGEM.reactions:
        try:
            sbo = r.annotation["SBO"]
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





def sort_dict(model_dict, by = "key"):
    sorted_dict = OrderedDict()
    if by == "key":
        for key in sorted(model_dict.keys()):
            sorted_dict[key] = model_dict[key]
    else:
        for key in sorted(model_dict, key = model_dict.get):
            sorted_dict[key] = model_dict[key]
    return sorted_dict

if __name__ == '__main__':
    if 0:
        scoGEM_FN = "../ModelFiles/xml/scoGEM.xml"
        scoGEM = cobra.io.read_sbml_model(scoGEM_FN)
        export(scoGEM)
    else:
        parser = argparse.ArgumentParser(description= "Export the scoGEM model")
        parser.add_argument("model_path", help = "Path to SBML (xml) model file")
        parser.add_argument("--formats", nargs='+', help = "The different file formats to export the model in", default = ["xml", "yml", "txt"])
        parser.add_argument("--name", type = str, help = "General name of exported model files", default = "scoGEM")
        parser.add_argument("--folder", type = str, help = "General name folder to store model files", default = "ModelFiles")
        parser.add_argument("--sbml_level", type = str, help = "SBML level", default = 3)
        parser.add_argument("--sbml_version", type = str, help = "SBML version", default = 1)
        parser.add_argument("--use_fbc_package", type = bool, help = "Use the fbc package", default = True)

        args = parser.parse_args()
        cwd = Path.cwd()

        model_path = str(cwd / args.model_path)
        try:
            scoGEM = cobra.io.read_sbml_model(model_path)
        except OSError as e:
            print("Cant find model file: {0}".format(model_path))
            raise

        else:
            if args.sbml_version != 1:
                print("Currently, only SBML-version 1 is supported. sbml_version parameter ignored.")
            sbml_version = 1
                
            export(scoGEM, args.folder, args.name, args.formats, args.sbml_level, sbml_version, args.use_fbc_package)