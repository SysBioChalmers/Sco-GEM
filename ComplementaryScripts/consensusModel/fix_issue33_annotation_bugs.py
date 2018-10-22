import cobra
import logging

def fix(scoGEM):
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
    fix(model)