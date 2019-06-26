# -*- coding: utf-8 -*-
"""
Author: Snorre Sulheim
Created: 28.02.2019
Updated: 12.03.2019
Email: snorre.sulheim@sintef.no


# Description
The main function of this script is to update the annotation of reactions to different subsystems and pathways. 
The curated annotation of subsystems and pathways is found in "ComplementaryData/curation/pathway_and_subsystem/subsystem_curated.csv".
The curation is performed by using KEGG annotations, BioCyc annotations and subsystem annotations from iMK1208. Every reaction is now annotated to
*only one* subsystem and *only one* pathway. For the reactions with none or multiple pathway annotations in KEGG or BioCyc, we have used adjacent
reactions to determine the annotation. All reactions are given a subsystem, but some reactions are missing a pathway-annotation.

In addition to the main function *update_subsystem_annotation.pt*, there are several functions which has performs the following tasks:
- get_pathways_from_biocyc: Extract pathway annotations from the BioCyc database online
- get_pathways_from_KEGG: Extract pathway annotations from the KEGG database online
- export_reaction_subsystem_and_pathways: Create a csv file with the current reaction annotations (the file used as a basis for curation)


Information about the metacyc pathway hierarchy
https://ask.pathwaytools.com/question/17981/accessing-the-pathway-hierarchy/
"""

import cobra
import pandas as pd 
from pathlib import Path
import json
import logging
from collections import defaultdict



REPO_DIR = Path(__file__).parent.parent.parent
KEGG_PATHWAYS_FN = str(REPO_DIR/"ComplementaryData"/"curation"/"pathway_and_subsystem"/"kegg_pathways.json")

def update_subsystem_annotations(model, csv_fn, remove_old_subsystem_annotations = True):
    """
    This is the main function of this script, and also the one which should be called 
    from the *reconstruct_scoGEM.py* script. Removes old pathway and subsystem annotations by default.Then new annotations based on the curated columns in the 
    spreadsheet subsystem_curated.csv.
    """
    df = pd.read_csv(csv_fn, sep = ";", usecols = ["Reaction ID", "curated pathway", "curated subsystem"], index_col = 0)

    for r in model.reactions:
        # Remove old annotations
        if remove_old_subsystem_annotations:
            try:
                r.annotation.pop("kegg.pathway")
            except KeyError:
                pass
            try:
                r.annotation.pop("kegg.subsystem")
            except KeyError:
                pass

        new_pathway_annotation = df.loc[r.id, "curated pathway"]     
        new_subsystem_annotation = df.loc[r.id, "curated subsystem"]
        
        _add_annotation(r, "pathway", new_pathway_annotation)
        _add_annotation(r, "subsystem", new_subsystem_annotation, ",")
        
def _add_annotation(r, key, value, delimiter  = ";"):
    if isinstance(value, (float, int)):
        return False

    if isinstance(value, str):
        if len(value):
            value = value.split(delimiter)
        else:
            return False

    if isinstance(value, list):
        value = [x.strip() for x in value]
        if len(value) == 1:
            if len(value[0]):
                r.annotation[key] = value[0]
        elif len(value) > 1:
            print("Multiple annotations for {0}: {1}".format(r.id, value))
            r.annotation[key] = value
        else:
            return False
    else:
        return False
    logging.info("Added annotatin to {0}: {1}".format(r.id, r.annotation[key]))
    return True



def get_pathways_from_biocyc(model, csv_fn, biocyc_subsystem_fn, db = "meta", add_annotations_to_model = True):
    """
    ** Warning: Must be run using python 2.7 and PathwayTools running in the background **

    This function use pythoncyc to extract pathway information from the BioCyc database based 
    on the BioCyc annotation of each reaction. The result is a table the rows are reactions (IDs)
    and the columns are the biocyc pathway annotations. Because BioCyC use very small pathways we 
    use the parent pathways as annotations.
    Some key steps are required to run this function:
    - PathwayTools must be running in the background ()
    - You need the pythoncyc package (https://github.com/latendre/PythonCyc), 
      more info at https://bioinformatics.ai.sri.com/ptools/pythoncyc.html
    - Pythoncyc only works with python 2.7

    # Parameters
    - model: SBML-model (imported with cobrapy)
    - csv_fn: Where to store the created csv-file
    - biocyc_subsystem_fn: This is in general the All_pathways_of_MetaCyC.txt file, but can be 
      replaced by similar csv-files.
    - db: Which db in BioCyc to use.
    - add_annotations_to_model: A flag used to turn on/off writing the biocyc annotations to the model reactions

    """
    import sys
    assert sys.version_info[0] < 3, ("Can't use PythonCyc with python 3")
    import pythoncyc # Add this import here, so it is only imported if used

    df_subsystem = pd.read_csv(biocyc_subsystem_fn, sep = "\t", index_col = 0)
    biocyc_db = pythoncyc.select_organism(db)
    pathway_list = []

    for r in model.reactions[::2]:
        print(r.id, end = "\t")
        try:
            biocyc = r.annotation["biocyc"]
        except KeyError:
            print()
            continue

        # Fix erroneous annotations
        if biocyc[:5] == "META:":
            biocyc = biocyc[5:]

        r_db = biocyc_db[biocyc]

        try:
            pathways = r_db["in_pathway"]
        except TypeError:
            print(biocyc, " is not in sco-db")
            continue

        if isinstance(pathways, list):
            sub1_list = []
            sub2_list = []
            for pathway in pathways:
                print(pathway, end = ", ")
                pwy = pathway.replace("|", "")
                try:
                    sub1 = df_subsystem.loc[pwy, "Subsystem 1"].split("//")[0].strip()
                    sub2 = df_subsystem.loc[pwy, "Subsystem 2"].split("//")[0].strip()
                except KeyError:
                    pass
                else:
                    sub1_list.append(sub1)
                    sub2_list.append(sub2)

            pathway_list.append([r.id, ";".join(pathway), ";".join(list(set(sub1_list))), ";".join(list(set(sub2_list)))])
            
            if len(sub1_list) and add_annotations_to_model:
                r.annotation["biocyc.subsystem1"] = list(set(sub1_list))
                r.annotation["biocyc.subsystem2"] = list(set(sub2_list))

            print(sub1_list, sub2_list)
        else:
            print("No pathways given for ", biocyc)


    df = pd.DataFrame(pathway_list)
    df.columns = ["Reaction ID", "Pathway", "Subsystem 1", "Subsystem 2"]
    df.to_csv(csv_fn)
    return model


def get_pathways_from_KEGG(model, update_existing = False):
    """
    This function extracts pathway and subsystem information from KEGG by using the KEGG annotation of each reaction.
    The pathways we use are the ones given here: https://www.genome.jp/kegg/pathway.html, 
    under heading 1.: Metabolism. However we don't use the *1.0 Global and overview maps* or 
    *1.12 Chemical structure and transformation maps*, because they don't 
    represent metabolic subsystems. What we here refer to as *subsustems* are the subheadings under Metabolism, i.e.:
    - Carbohydrate metabolism
    - Energy metabolism
    - Lipid metabolism
    - Nucleotide metabolism
    - Amino acid metabolism
    - Metabolism of other amino acids
    - Glycan biosynthesis and metabolism
    - Metabolism of cofactors and vitamins
    - Metabolism of terpenoids and polyketides
    - Biosynthesis of other secondary metabolites
    - Xenobiotics biodegradation and metabolism

    """
    from bioservices.kegg import KEGG
    kegg = KEGG()
    kegg_dict, kegg_overview_maps = _get_KEGG_pathways()
    inverse_pathway_dict = _get_inverse_pathway_dict(kegg_dict)

    for reaction in model.reactions:
        
        # Skip reactions which already have an kegg.pathway annoatation
        # if update_existing = False
        if not update_existing:
            try:
                reaction.annotation["kegg.pathway"]
            except KeyError:
                pass
            else:
                # Skip this one
                continue

        try:
            kegg_id = reaction.annotation["kegg.reaction"]
        except KeyError:
            continue

        kegg_info = kegg.get(kegg_id, parse = True)

        try:
            full_kegg_pathways = kegg_info["PATHWAY"].values()
        except:
            continue

        kegg_pathways = [x for x in full_kegg_pathways if not x in kegg_overview_maps]
        
        try:
            subsystem = list(set([inverse_pathway_dict[x] for x in kegg_pathways]))
        except:
            print("Error!: ", reaction.id, kegg_pathways)
            continue
        
        print("KEGG Subsystem ", reaction.id, subsystem)

        reaction.annotation["kegg.pathway"] = kegg_pathways
        reaction.annotation["kegg.subsystem"] = subsystem
    return model

def _get_KEGG_pathways():
    with open(KEGG_PATHWAYS_FN, "r") as f:
        kegg_dict = json.load(f)

    kegg_overview_maps = kegg_dict.pop("KEGG overview maps")
    return kegg_dict, kegg_overview_maps

def _get_inverse_pathway_dict(kegg_dict):
    new_dict = {}
    for k, v in kegg_dict.items():
        for v_i in v:
            new_dict[v_i] = k
    return new_dict

def export_reaction_subsystem_and_pathways(model, csv_fn):
    """
    Use pandas to write a csv-file which can be used to curate the subsystem
    annotation. The csv-file is written with the following columns:
    
    Reaction ID, Reaction name, KEGG ID, Biocyc-annotation, 
    KEGG Subsystem, KEGG pathway, Subsystem
    """
    annotation_keys = ["kegg.reaction", "biocyc", "kegg.pathway", "kegg.subsystem", 
                       "biocyc.subsystem1", "biocyc.subsystem2", "subsystem"]
    reactions_list = []
    for r in model.reactions:
        r_list = [r.id, r.name]

        for key in annotation_keys:
            try:
                ann = r.annotation[key]
            except KeyError:
                r_list.append(None)
            else:
                if isinstance(ann, str):
                    r_list.append(ann)
                else:
                    r_list.append(", ".join(ann))

        reactions_list.append(r_list)
    df = pd.DataFrame(reactions_list)
    df.columns = ["Reaction ID", "Reaction name"] +  annotation_keys
    # Add empty column
    df["curated pathway"] = None 
    df["curated subsystem"] = None 
    print(df.head())
    df.to_csv(csv_fn, sep = ";")


def print_subsystem_summary(model, key = "subsystem"):
    subsystem_total = defaultdict(int)
    subsystem_other = defaultdict(int)
    subsystem_sco4 = defaultdict(int)
    subsystem_iAA1259 = defaultdict(int)
    subsystem_iKS1317 = defaultdict(int)

    for r in model.reactions:
        subsystem = r.annotation[key]
        try:
            origin = r.annotation["origin"]
        except KeyError:
            print(r)
            origin = "missing"
        if origin == "Sco4":
            subsystem_sco4[subsystem] += 1
        elif origin == "iAA1259":
            subsystem_iAA1259[subsystem] += 1
        elif origin == "missing":
            subsystem_other[subsystem] += 1
        else:
            subsystem_iKS1317[subsystem] += 1
        subsystem_total[subsystem] += 1

    df = pd.DataFrame([subsystem_total, subsystem_iKS1317, subsystem_sco4, subsystem_iAA1259, subsystem_other]).T
    df.columns = ["Total", "iKS1317", "Sco4", "iAA1259", "Other"]
    print(df)


def export_gene_pathway_list(model):
    gene_pathway_list = []
    maxlen = 0
    for gene in model.genes:
        if gene.id == "s0001":
            # These are spontaneous reactions
            continue

        pathway_list = []
        for r in gene.reactions:
            try:
                pathway = r.annotation["pathway"]
            except:
                continue
            else:
                pathway_list.append(pathway)
        pathway_list = list(set(pathway_list))
        if len(pathway_list):
            for pathway in pathway_list:
                gene_pathway_list.append([gene.id, pathway])
        else:
            gene_pathway_list.append([gene.id, None])

    df = pd.DataFrame(gene_pathway_list, columns = ["Gene", "Pathway"])
    df.to_csv(str(REPO_DIR / "model_gene_pathway_table.tsv"), sep = "\t")




if __name__ == '__main__':
    model_fn = REPO_DIR / "ModelFiles" / "xml" / "scoGEM.xml"
    model = cobra.io.read_sbml_model(str(model_fn))

    if 0:
        # Create file used for subsystem curation
        biocyc_pwy_fn = str(REPO_DIR / "ComplementaryData" / "curation" / "pathway_and_subsystem" / "reaction_biocyc_pathway.csv")
        biocyc_subsystem_fn = str(REPO_DIR / "ComplementaryData" / "curation" / "pathway_and_subsystem" / "All_pathways_of_MetaCyc.txt")
        csv_fn = str(REPO_DIR / "ComplementaryData" / "curation" / "pathway_and_subsystem" / "subsystem.csv")
        model = get_pathways_from_biocyc(model, biocyc_pwy_fn, biocyc_subsystem_fn)
        model = get_pathways_from_KEGG(model)
        export_reaction_subsystem_and_pathways(model, csv_fn)

    if 0:
        import sys
        sys.path.append("C:/Users/snorres/git/scoGEM/ComplementaryScripts")
        import export
        # update the subsystem annotations based on the curated csv-file"
        subsystem_curated_csv = str(REPO_DIR / "ComplementaryData" / "curation" / "pathway_and_subsystem" / "subsystem.txt")
        update_subsystem_annotations(model, subsystem_curated_csv)
        # export.export(model, formats = ["xml", "yml"])

    if 1:
        # Print subsystem numbers
        print_subsystem_summary(model)
    if 0:
        export_gene_pathway_list(model)

