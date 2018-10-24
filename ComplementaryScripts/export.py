# -*- coding: utf-8 -*-
"""
Ensure that the model is exported in a given order
The following order is used for metabolites, genes and reactions:
TODO
Reactions
---------
1. biomass reactions
2. atpm reactions
3. encapsulating reactions
4. pseudo reactions
5. normal reactions
6. transport reactions
7. demand reactions
8. exchange reactions
9. sink reactions
10. other reactions
11. non-sbo reactions

"""
from consensusModel.fix_SBO_terms import REACTION_SBO_TERMS

import cobra
from collections import OrderedDict
from pathlib import Path

REPO_MAIN_FOLDER = Path(__file__).resolve().parent.parent



def export(scoGEM, folder = "ModelFiles", name = "scoGEM", filetypes = ["xml", "yml", "txt", "json", "mat"]):
    sort_model(scoGEM)

    main_folder_path = REPO_MAIN_FOLDER / folder
    # write xml
    if ("xml" in filetypes) or ("sbml" in filetypes):
        model_fn_xml = main_folder_path / "xml" / "{0}.xml".format(name)
        check_folder(model_fn_xml)
        cobra.io.write_sbml_model(scoGEM, model_fn_xml)
    
    if ("yml" in filetypes) or ("yaml" in filetypes):
        model_fn_yml = main_folder_path / "yml" / "{0}.yml".format(name)
        check_folder(model_fn_yml)
        cobra.io.save_yaml_model(scoGEM, model_fn_yml)

    if "json" in filetypes:
        model_fn_json = main_folder_path / "json" / "{0}.json".format(name)
        check_folder(model_fn_json)
        cobra.io.save_json_model(scoGEM, model_fn_json)
        
    if "txt" in filetypes:
        raise NotImplementedError
        # model_fn_txt = main_folder_path / "txt" / "{0}.txt".format(name)
        # check_folder(model_fn_txt)
        # cobra.io.write_sbml_model(scoGEM, model_fn_txt)

    if "mat" in filetypes:
        model_fn_mat = main_folder_path / "mat" / "{0}.mat".format(name)
        check_folder(model_fn_mat)
        cobra.io.save_matlab_model(scoGEM, model_fn_mat)
    

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
    scoGEM_FN = "../ModelFiles/xml/scoGEM.xml"
    scoGEM = cobra.io.read_sbml_model(scoGEM_FN)
    export(scoGEM)