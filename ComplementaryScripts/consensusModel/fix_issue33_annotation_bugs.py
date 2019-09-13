#!/usr/bin/env python
"""
Author: Snorre Sulheim

# Description
This file fixes a variety of bugs and smal improvements related to annotations. See issue 33 for more info.
"""
import cobra
import logging
import pandas as pd

def fix(scoGEM):
    fix_annotations(scoGEM)
    fix_demand_biocyc_names(scoGEM)
    fix_misc(scoGEM)
    annotate_germicidin_pathway(scoGEM)
    fix_wrong_chebi_mapping(scoGEM)
    remove_biocyc_annotaions_from_exchanges(scoGEM)

def fix_c_c_in_metabolite_ids(scoGEM):
    for m in scoGEM.metabolites:
        if m.id[-4:] == "_c_c":
            old_id = m.id
            m.id = m.id[:-4] + "_c"
            logging.info("Changed id of metabolite {0} to {1}".format(old_id, m.id))


def fix_annotations(scoGEM):

    for r in scoGEM.reactions:
        try:
            kegg_annotation = r.annotation["kegg.reaction"]
        except KeyError:
            pass
        else:
            if not len(kegg_annotation):
                r.annotation.pop("kegg.reaction")
                logging.info("Removed empty KEGG annotation from reaction".format(r.id))
        
        # Remove empty EC-code annotations
        try:
            ec_annotation = r.annotation["ec-code"]
        except KeyError:
            pass
        else:
            if not len(ec_annotation):
                r.annotation.pop("ec-code")
                logging.info("Removed empty EC-code annotation from metabolite {0}".format(r.id))
        
        # Remove appended _nadh from biocyc annotations
        try:
            biocyc_annotation = r.annotation["biocyc"]
        except KeyError:
            pass
        else:
            if biocyc_annotation.endswith("_NADH"):
                r.annotation["biocyc"] = biocyc_annotation.replace("_NADH", "")
                logging.info("Changed BiGG annotation of {0} from {1} to {2}".format(r.id, biocyc_annotation, r.annotation["biocyc"]))
            elif biocyc_annotation.endswith("_NADPH"):
                r.annotation["biocyc"] = biocyc_annotation.replace("_NADH", "")
                logging.info("Changed BiGG annotation of {0} from {1} to {2}".format(r.id, biocyc_annotation, r.annotation["biocyc"]))

        
    for m in scoGEM.metabolites:
        # Remove appended _c and _e from last part of BiGG annotations
        try:
            bigg_annotation = m.annotation["bigg.metabolite"]
        except KeyError:
            pass
        else:
            if bigg_annotation.endswith("_c") or bigg_annotation.endswith("_e"):
                m.annotation["bigg.metabolite"] = bigg_annotation[:-2]
                logging.info("Changed BiGG annotation of {0} from {1} to {2}".format(m.id, bigg_annotation, bigg_annotation[:-2]))
        
        # Remove empty KEGG annotations
        try:
            kegg_annotation = m.annotation["kegg.compound"]
        except KeyError:
            pass
        else:
            if not len(kegg_annotation):
                m.annotation.pop("kegg.compound")
                logging.info("Removed empty KEGG annotation from metabolite {0}".format(m.id))

       

        # Remove empty metanetx annotations
        try:
            mnx_annotation = m.annotation["metanetx.chemical"]
        except KeyError:
            pass
        else:
            if not len(mnx_annotation):
                m.annotation.pop("metanetx.chemical")
                logging.info("Removed empty MetaNetX annotation from metabolite {0}".format(m.id))

        # Remove appended _c from biocyc annotations
        try:
            biocyc_annotation = m.annotation["biocyc"]
        except KeyError:
            pass
        else:
            if biocyc_annotation.endswith("_c") or biocyc_annotation.endswith("_e"):
                m.annotation["biocyc"] = biocyc_annotation[:-2]
                logging.info("Changed BiGG annotation of {0} from {1} to {2}".format(m.id, biocyc_annotation, biocyc_annotation[:-2]))

def fix_misc(scoGEM):
    # Annotate metabolite 44pcpopd_c to KEGG C21770
    m = scoGEM.metabolites.get_by_id("44dpcopd_c")
    m.annotation["kegg.compound"] = "C21770"
    logging.info("Annotated metabolite 44pcpopd_c to KEGG ID C21770")

    scoGEM.metabolites.get_by_id("4gglutbut_c").annotation["kegg.compound"] = "C15700"
    scoGEM.metabolites.get_by_id("4gglutbut_c").annotation["metanetx.chemical"] = "MNXM1378"
   
    scoGEM.metabolites.get_by_id("xylan_e").annotation["kegg.compound"] = "C00707"

   


   
def annotate_germicidin_pathway(scoGEM):
    # Give biocyc annotation to germicidin
    scoGEM.metabolites.get_by_id("germicidinA_c").annotation["biocyc"] = "CPD1UA-13"
    scoGEM.metabolites.get_by_id("germicidinB_c").annotation["biocyc"] = "CPD1UA-14"
    scoGEM.metabolites.get_by_id("germicidinC_c").annotation["biocyc"] = "CPD1UA-15"
    scoGEM.metabolites.get_by_id("germicidinD_c").annotation["biocyc"] = "CPD1UA-16"
    
    scoGEM.reactions.get_by_id("GERMA").annotation["biocyc"] = "RXN1UA-25"
    scoGEM.reactions.get_by_id("GERMB").annotation["biocyc"] = "RXN1UA-20"
    scoGEM.reactions.get_by_id("GERMC").annotation["biocyc"] = "RXN1UA-21"
    scoGEM.reactions.get_by_id("GERMD").annotation["biocyc"] = "RXN1UA-24"



def fix_metanetx_reaction_annotations(scoGEM, metanetx_fn):
    """
    This csv file used for the mapping is created using the function 
    'map_model_reactions' in the map_to_metanetx.py script.
    """
    # Read in dataframe
    df = pd.read_csv(metanetx_fn, index_col = 0)
    
    change_dict = {}
    for i, row in df.iterrows():
        mnx_annotation = [row[key] for key in ["MNX 1", "MNX 2", "MNX 3"] if pd.notna(row[key])]
        change_dict[row["Reaction ID"]] = mnx_annotation

    # Remove all metanetx annotations
    i = 0
    for r in scoGEM.reactions:
        try:
            old_annotation = r.annotation.pop("metanetx.reaction")
            # If the origin is not "Sco4" I don't trust these metantex-annotations, and neither the rhea
        except KeyError:
            old_annotation = None
        
        try:
            new_annotation = change_dict[r.id]
        except KeyError:
            new_annotation = None
        else:
            r.annotation["metanetx.reaction"] = new_annotation
            i += 1

        if new_annotation:
            logging.info("{0}: Changed metanetx annotation from {1} to {2}".format(r.id, old_annotation, new_annotation))
        elif old_annotation:
            logging.info("Removed metanetx annotation from {0}".format(r.id))
        else:
            logging.info("No metanetx annotation for {0}".format(r.id))

    logging.info("{0} of {1} reactions were given metanetx identifiers".format(i, len(scoGEM.reactions)))
    


def remove_biocyc_annotaions_from_exchanges(scoGEM):
    # Remove biocyc annotations from exchanges
    for r in scoGEM.exchanges:
        try:
            r.annotation.pop("biocyc")
        except KeyError:
            pass
        else:
            logging.info("Removed biocyc annotation from {0}".format(r.id))
        


def fix_metanetx_metabolite_annotations(scoGEM, met_to_metanetx_fn):
    """
    This csv file used for the mapping is created using the function 
    'map_model_metabolites' in the map_to_metanetx.py script.
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
            m.annotation["metanetx.chemical"] = [row[1], row[2]]
        elif isinstance(row[1], str):
            m.annotation["metanetx.chemical"] = row[1]
        else:
            continue
        logging.info("Changed metanetx.chemical annotation of metabolite {0} from {1} to {2}".format(
                      m.id, old_anno, m.annotation["metanetx.chemical"]))

def apply_new_chebi_annotations(scoGEM, chebi_annotation_fn):
    """
    This csv file used for the mapping is created using the function 
    'map_metabolites_to_chebi' in the map_to_metanetx.py script.
    """
    df = pd.read_csv(chebi_annotation_fn, index_col = None, converters = {"New chebi annotation": lambda x: x.strip("[]").replace("'", "").split(", ")})
    # df["New chebi annotation"] = df["New chebi annotation"].apply(lambda x: x[1:-1].replace("'", "").split(", "))
    for i, row in df.iterrows():
        m_id = row["Met ID"]
        new_annotation = row["New chebi annotation"]#.strip("[]").replace("'","").split(", ")

        m = scoGEM.metabolites.get_by_id(m_id)
        if len(new_annotation):
            if not len(new_annotation[0]):
                pass
            else:
                m.annotation["chebi"] = new_annotation
                logging.info("Changed chebi annotation of metabolite {0} from {1} to {2}".format(
                      m.id, row[2], new_annotation))

def fix_wrong_chebi_mapping(scoGEM):
    for m in scoGEM.metabolites:
        try:
            origin = m.annotation["origin"]
        except KeyError:
            continue

        try:
            chebi = m.annotation["chebi"]
        except KeyError:
            continue
        if origin == "Sco4":
            m.annotation.pop("chebi")
        else:
            if isinstance(chebi, str):
                chebi = [chebi]
            new_chebis = []
            for cheb in chebi:
                if not "CHEBI" in cheb:
                    new_chebis.append("CHEBI:{0}".format(cheb.split(":")[-1]))
                else:
                    new_chebis.append(cheb)
            m.annotation["chebi"] = new_chebis



def fix_demand_biocyc_names(scoGEM):
    reaction_ids = ["DM_acetone_c", "DM_ahop_c", "DM_2mborn_c", "DM_33biflav_c", "DM_38biflav_c"]
    for r_id in reaction_ids:
        r = scoGEM.reactions.get_by_id(r_id)
        # Remove "_c" from last past of biocyc annotation
        r.annotation["biocyc"] = r.annotation["biocyc"][:-2]

        # Change name
        met = list(r.metabolites.keys())[0]
        r.name = "Sink needed to allow {0} to leave system".format(met.name)
        logging.info("Changed name and fixed biocyc annotation of reaction {0}".format(r_id))

def fix_mmy_bug(scoGEM):
    """
    Discarded from reconstruction, because it was a bug from the fik_iKS1317_issues.py script, and
    that bug was sorted out there. 
    """
    # Change ID bug with acetoacetyl-mmyA and malonyl-mmyA
    m_ac = scoGEM.metabolites.get_by_id("malACPmmy_c")
    m_ac.id = "actACPmmy_c"
    m_ac.annotation.pop("bigg.metabolite")
    
    m_mal = scoGEM.metabolites.get_by_id("temp")
    m_mal.id = "malACPmmy_c"
    m_mal.annotation.pop("bigg.metabolite")
    logging.info("Fixed a bug with the IDs of acetoacetyl-mmyA and malonyl-mmyA. IDs were swapped")

if __name__ == '__main__':
    model = cobra.io.read_sbml_model("../../ModelFiles/xml/scoGEM.xml")
    # fix_annotations(model)
    # fix_demand_biocyc_names(model)
    # apply_new_chebi_annotations(model, "../../ComplementaryData/curation/chebi_annotation.csv")  
    fix_metanetx_reaction_annotations(model, "../../ComplementaryData/curation/metanetx_reaction_annotations_to_change.csv")