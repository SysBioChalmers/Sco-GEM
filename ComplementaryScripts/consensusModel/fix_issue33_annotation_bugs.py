import cobra
import logging
import pandas as pd

def fix(scoGEM):
    fix_annotations(scoGEM)
    fix_mmy_bug(scoGEM)
    fix_demand_biocyc_names(scoGEM)
    fix_misc(scoGEM)
    # fix_metanetx_annotations(scoGEM, met_to_metanetx_fn)

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
   
def fix_metanetx_annotations(scoGEM, met_to_metanetx_fn):
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
# def fix_metanetx(scoGEM): # fix metanetx annotation of atp and adp
#     scoGEM.metabolites.h2o_c.annotation["metanetx.chemical"] = "MNXM2"
#     scoGEM.metabolites.h2o_e.annotation["metanetx.chemical"] = "MNXM2"
#     scoGEM.metabolites.atp_c.annotation["metanetx.chemical"] = "MNXM3"
#     scoGEM.metabolites.o2_c.annotation["metanetx.chemical"] = "MNXM4"
#     scoGEM.metabolites.o2_e.annotation["metanetx.chemical"] = "MNXM4"
#     scoGEM.metabolites.adp_c.annotation["metanetx.chemical"] = "MNXM7"
#     scoGEM.metabolites.nad_c.annotation["metanetx.chemical"] = "MNXM8"
#     scoGEM.metabolites.pi_c.annotation["metanetx.chemical"] = "MNXM9"
#     {"etha": "MNXM218",
#      "ppi" : ""}

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
    fix_demand_biocyc_names(model)