#!/usr/bin/env python
"""
Author: Snorre Sulheim
Date: 11.2018

# Description
Map model metabolites and reactions to metanetx identifiers using KEGG, metacyc and BiGG annotations.
"""
import cobra
import pandas as pd
from pathlib import Path
import libchebipy

def map_model_metabolites(model, metanetx_fn):
    df = pd.read_csv(metanetx_fn, header = None, sep = "\t", comment = "#")
    df.columns = ["db:id", "metanetx", "reason", "name"]

    kegg_df = df[df["db:id"].str.contains("kegg:")]
    bigg_df = df[df["db:id"].str.contains("bigg:")]
    metacyc_df = df[df["db:id"].str.contains("metacyc:")]
    # Deprecated?
    del df
    
    new_df_list = []
    for m in model.metabolites:
        # print(m.annotation)
        mnx_annotations = []
        bigg_match = bigg_df.loc[bigg_df["db:id"] == "bigg:{}".format(m.id[:-2]), :]
        mnx_annotations += list(bigg_match["metanetx"])

        try:
            kegg_id = m.annotation["kegg.compound"]
        except KeyError:
            kegg_match = None
        else:
            kegg_match = kegg_df.loc[kegg_df["db:id"] == "kegg:{}".format(kegg_id), :]
            mnx_annotations += list(kegg_match["metanetx"])

        # metacyc
        if len(mnx_annotations) == 0:

            try:
                metacyc_id = m.annotation["biocyc"]
            except KeyError:
                metacyc_match = None
                print(m.id, None)
            else:
                metacyc_match = metacyc_df.loc[metacyc_df["db:id"] == "metacyc:{}".format(metacyc_id), :]
                mnx_annotations += list(metacyc_match["metanetx"])
                print(m.id, list(metacyc_match["metanetx"]))

        # mnx_annotations = list(set(mnx_annotations))

        if "metanetx.chemical" in list(m.annotation.keys()):
            mnx_annot = m.annotation["metanetx.chemical"]
        else:
            mnx_annot = None

        if not mnx_annot and len(mnx_annotations):
            new_df_list.append([m.id]+ mnx_annotations)
        elif (len(mnx_annotations) == 1) and mnx_annot != mnx_annotations[0]:
            new_df_list.append([m.id]+ mnx_annotations)
        elif (len(mnx_annotations) == 2) and isinstance(mnx_annot, str):
            new_df_list.append([m.id]+ mnx_annotations)
        elif (len(mnx_annotations) == 2) and isinstance(mnx_annot, list):
            print("#######", m.id, mnx_annotations, mnx_annot)
        else:
            pass
        # print(new_df_list)

    new_df = pd.DataFrame(new_df_list, columns = ["Met ID", "MNX 1", "MNX 2"])
    new_df.to_csv("../../ComplementaryData/curation/metanetx_to_change.csv", index_label = "index")
    # print(new_df)

def map_model_reactions(model, metanetx_fn):
    df = pd.read_csv(metanetx_fn, header = None, sep = "\t", comment = "#")
    df.columns = ["db:id", "metanetx", "reaction string"]
    
    kegg_df = df[df["db:id"].str.contains("kegg:")]
    bigg_df = df[df["db:id"].str.contains("bigg:")]
    metacyc_df = df[df["db:id"].str.contains("metacyc:")]    
    deprecated_df = df[df["db:id"].str.contains("deprecated:")]    
    del df

    new_df_list = []
    for r in model.reactions:
        mnx_annotations = []

        # BiGG
        bigg_match = bigg_df.loc[bigg_df["db:id"] == "bigg:{}".format(r.id), :]
        mnx_annotations += list(bigg_match["metanetx"])

        # KEGG
        try:
            kegg_id = r.annotation["kegg.reaction"]
        except KeyError:
            kegg_match = None
        else:
            kegg_match = kegg_df.loc[kegg_df["db:id"] == "kegg:{}".format(kegg_id), :]
            mnx_annotations += list(kegg_match["metanetx"])

        # Metacyc
        try:
            metacyc_id = r.annotation["biocyc"]
        except KeyError:
            metacyc_match = None
        else:
            metacyc_match = metacyc_df.loc[metacyc_df["db:id"] == "metacyc:{}".format(metacyc_id), :]
            mnx_annotations += list(metacyc_match["metanetx"])
            # print(r.id, list(metacyc_match["metanetx"]))

        # Remove duplicates
        mnx_annotations = list(set(mnx_annotations))


        # print("{0:<3} {1:<100} {2:<20} {3:<50} {4}".format(i, str(r), str(mnx_annot), ", ".join(mnx_annotations), origin))

        if len(mnx_annotations):
            new_df_list.append([r.id] + mnx_annotations)

    new_df = pd.DataFrame(new_df_list, columns = ["Reaction ID", "MNX 1", "MNX 2", "MNX 3"])
    new_df.to_csv("../../ComplementaryData/curation/metanetx_reaction_annotations_to_change.csv", index_label = "index")


def map_metabolites_to_chebi(scoGEM, metanetx_fn):
    df = pd.read_csv(metanetx_fn, header = None, sep = "\t", comment = "#")
    df.columns = ["db:id", "metanetx", "reason", "name"]

    chebi_df = df[df["db:id"].str.contains("chebi:")]
    del df

    new_df_list = []
    for m in model.metabolites:
        try:
            mnx_annot = m.annotation["metanetx.chemical"]
        except KeyError:
            print("No metanetx annotation for {0}, {1}".format(m.id, ["{0}:{1}".format(key, value) for key, value in m.annotation.items()]))
            continue

        mnx_annot = as_list(mnx_annot)

        chebi_ids = []
        for mnx_i in mnx_annot:
            mnx_match = chebi_df.loc[chebi_df["metanetx"] == mnx_i]
            chebi_ids += list(mnx_match["db:id"].values)
        chebi_ids = [x.upper() for x in chebi_ids]

        parent_chebis = []
        for chebi_id in list(set(chebi_ids)):
            lib_data = libchebipy.ChebiEntity(chebi_id.upper())
            parent = lib_data.get_parent_id()
            if parent:
                parent_chebis.append(parent)
            else:
                parent_chebis.append(chebi_id.upper())
        parent_chebis = list(set(parent_chebis))
        try:
            current_chebi_list = as_list(m.annotation["chebi"])
        except:
            current_chebi_list = [None]
            in_new_chebis = False
        else:
            in_new_chebis = True   
            for current_chebi in current_chebi_list:
                if not current_chebi in chebi_ids:
                    in_new_chebis = False
                    print("{2}: {0} is not in the new set {1}".format(current_chebi, parent_chebis, m.id))
        
        new_df_list.append([m.id, parent_chebis, current_chebi_list, in_new_chebis])
    new_df = pd.DataFrame(new_df_list, columns = ["Met ID", "New chebi annotation", "Current chebi annotation", "Old chebi in new (including secondary chebis)"])
    new_df.to_csv("../../ComplementaryData/curation/chebi_annotation.csv", index = False)

def as_list(param):
    if isinstance(param, list):
        return param
    else:
        return [param]

def apply_metanetx_mapping(scoGEM, met_to_metanetx_fn):
    """
    Depreceated: moved to fix_issue33_annotation_bugs.py
    """
    df = pd.read_csv(met_to_metanetx_fn, index_col = 0)
    for i, row in df.iterrows():
        m_id = row[0]

        m = scoGEM.metabolites.get_by_id(m_id)
        try:
            old_anno = m.annotation["metanetx.chemical"]
        except KeyError:
            old_anno = None

        if isinstance(row[2], str):
            m.annotation["metanetx.chemical"] = row[1:2]
        elif isinstance(row[1], str):
            m.annotation["metanetx.chemical"] = row[1]
        else:
            continue
        logging.info("Changed metanetx.chemical annotation of metabolite {0} from {1} to {2}".format(
                      m.id, old_anno, m.annotation["metanetx.chemical"]))




if __name__ == '__main__':
    repo_path = Path(__file__).parent.parent.parent
    model_fn = repo_path / "ModelFiles" / "xml" / "scoGEM.xml"
    model = cobra.io.read_sbml_model(str(model_fn))
    # map_model_metabolites(model, metanetx_fn)
    fn = repo_path / "ComplementaryData" / "curation" /"metanetx_to_change.csv"
    if 0:
        apply_metanetx_mapping(model, fn)

    if 0:
        metanetx_fn = repo_path / "ComplementaryData" / "curation" / "metanetx_chem_xref.tsv"
        map_metabolites_to_chebi(model, metanetx_fn)

    if 1:
        metanetx_fn = repo_path / "ComplementaryData" / "curation" / "metanetx_reac_xref.tsv"
        map_model_reactions(model, metanetx_fn)
