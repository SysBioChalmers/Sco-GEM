import cobra
import pandas as pd
import numpy as np


sco4_fn = r"C:\\Users\\snorres\\git\\scoGEM\\ComplementaryData\\models\\Sco4.xml"
iKS1317_fn = r"C:\\Users\\snorres\\git\\scoGEM\\ComplementaryData\\models\\iKS1317.xml"
iAA1259_fn = r"C:\\Users\\snorres\\git\\scoGEM\\ComplementaryData\\models\\iAA1259.xml"
iMK1208_fn = r"C:\\Users\\snorres\\OneDrive - SINTEF\\SINTEF projects\\INBioPharm\\SCM\\sbml models\\Kim2014\\kim_with_kegg.xml"
CHANGED_REACTIONS_FN = r"C:\Users\snorres\OneDrive - SINTEF\SINTEF projects\INBioPharm\SCM\sbml models\Snorre2017\changed_reaction_ids.csv"
CHANGED_METABOLITES_FN = r"C:\Users\snorres\OneDrive - SINTEF\SINTEF projects\INBioPharm\SCM\sbml models\Snorre2017\changed_metabolite_ids.csv"



def get_models():
    sco4 = cobra.io.read_sbml_model(sco4_fn)
    iKS1317 = cobra.io.read_sbml_model(iKS1317_fn)
    iAA1259 = cobra.io.read_sbml_model(iAA1259_fn)
    return sco4, iKS1317, iAA1259


def list_genes_reactions_metabolites():
    models = get_models()
    for model in models:
        print("{0}: {1} reactions, {2} metabolites, {3} genes".format(len(model.reaction), len(model.metabolites), len(model.genes)))



def compare_genes():
    sco4, iKS1317, iAA1259 = get_models()
    all_genes = []
    specific_genes = []
    for m in [iKS1317, sco4, iAA1259]:
        gene_ids = []
        for g in m.genes:
            gene_ids.append(g.id)
            all_genes.append(g.id)
        specific_genes.append(gene_ids)
    all_genes = list(set(all_genes))

    has_genes = []

    for gene_id in all_genes:
        has_this_gene = []
        for lst in specific_genes:
            if gene_id in lst:
                has_this_gene.append(1)
            else:
                has_this_gene.append(0)
        has_genes.append(has_this_gene)

    df = pd.DataFrame(has_genes)
    df.index = all_genes
    df.columns = ["iKS1317", "sco4", "iAA1259"]
    df.to_csv("gene_comparison.csv")

    # Unique genes
    additional_genes_in_sco4 = [x for x in specific_genes[1] if not x in specific_genes[0]]
    additional_genes_in_iAA1259 = [x for x in specific_genes[2] if not x in specific_genes[0]]


    print(additional_genes_in_sco4)
    print(additional_genes_in_iAA1259)
    print([x for x in additional_genes_in_sco4 if x in additional_genes_in_iAA1259])

def rename_iKS1317(model_iKS1317):
    change_reaction_ids = pd.read_csv(CHANGED_REACTIONS_FN, header = None, sep = ";", index_col = 0).to_dict()[1]
    change_metabolite_ids = pd.read_csv(CHANGED_METABOLITES_FN, header = None, sep = ",", index_col = 0).to_dict()[1]

    for key, value in change_reaction_ids.items():
        # print(key, value)
        try:
            r = model_iKS1317.reactions.get_by_id(value)
        except KeyError:
            print("Missing reaction: {0}, old id: {1}".format(value, key))
        else:
            r.id = key
    for key, value in change_metabolite_ids.items():
        # print(key, value)
        try:
            m = model_iKS1317.metabolites.get_by_id(value)
        except KeyError:
            print("Missing metabolite: {0}, old id: {1}".format(value, key))
        else:
            m.id = key

    model_iKS1317.reactions.EX_glyc_R_e.id = "EX_glyc__R_e"
    model_iKS1317.reactions.EX_pnto_R_e.id = "EX_pnto__R_e"

def fix_sco4(model_sco4):
    for r in model_sco4.reactions:
        r.id = r.id.replace("__45__", "__").replace("__40__e__41__", "_e").replace("__40__c__41__", "_c")
    for m in model_sco4.metabolites:
        m.id = m.id.replace("__45__", "-").replace("__40__e__41__", "_e").replace("__40__c__41__", "_c")
    # model_sco4.reactions.EX_glc_e.id = "EX_glc__D_e"

def fix_iAA1259(model):
    for r in model.reactions:
        r.id = r.id.replace("_LPAREN_e_RPAREN_", "_e").replace("_LPAREN_c_RPAREN_", "_c").replace("_DASH_", "__")
    for m in model.metabolites:
        m.id = m.id.replace("_LPAREN_e_RPAREN_", "_e").replace("_LPAREN_c_RPAREN_", "_c").replace("_DASH_", "-")
    

def fix_iKS1317(model):
    rename_iKS1317(model)
    for m in model.metabolites:
        m.id = m.id.replace("__", "-")
        for key in "RBDLMST52C4":
            m.id = m.id.replace("_{0}".format(key), "-{0}".format(key))
  
                    

def compare_kim_content_of_models2():
    model_sco4 = cobra.io.read_sbml_model(sco4_fn)
    fix_sco4(model_sco4)
    model_iKS1317 = cobra.io.read_sbml_model(iKS1317_fn)
    fix_iKS1317(model_iKS1317)
    
    model_iAA1259 = cobra.io.read_sbml_model(iAA1259_fn)
    fix_iAA1259(model_iAA1259)
    model_iMK1208 = cobra.io.read_sbml_model(iMK1208_fn)
    model_iMK1208.reactions.EX_glc__D_e.id = "EX_glc_e"
    

    reaction_list = []
    for kim_reaction in model_iMK1208.reactions:
        temp_list = [kim_reaction.id]
        for model in [model_iKS1317, model_sco4, model_iAA1259]:
            try:
                model_reaction = model.reactions.get_by_id(kim_reaction.id)
            except KeyError:
                temp_list.append("Missing")
            else:
                temp_list.append(compare_reactions(kim_reaction, model_reaction))

        reaction_list.append(temp_list)    
    df = pd.DataFrame(reaction_list, columns =  ["iMK1208 ID", "iKS1317", "Sco4", "iAA1259"])
    df.to_csv("../../ComplementaryData/compare_iMK1208_content.csv", index = False, sep = ";")
    

def compare_reactions(reaction, reaction2):
    differences = []
    # Bounds
    if reaction.bounds != reaction2.bounds:
        differences.append("bounds")
    
    # Annotations
    # EC - Code
    # Kegg / Biocyc


    # Genes
    # g1 string = reaction.gene_reaction_rule.replace("(", "").replace(")","").strip(" ")
    # g2 = reaction2.gene_reaction_rule.replace("(", "").replace(")","").strip(" ")
    g1 = set([g.id for g in reaction.genes])
    g2 = set([g.id for g in reaction2.genes])
    if g1 != g2:
        print("##########   {0}   ############".format(reaction.id))
        print(g1)
        print(g2)
        differences.append("genes")

    # Metabolites
    mets_1 = set([(m.id, i) for m, i in reaction.metabolites.items()])
    mets_2 = set([(m.id, i) for m, i in reaction2.metabolites.items()])
    mets_inv = set([(x, -i) for x, i in mets_2])

    if (mets_1 == mets_2) or (reaction.reversibility and mets_1 == mets_inv):
        pass
    else:
        if not len(differences):
            print("##########   {0}   ############".format(reaction.id))
        print(reaction.reaction)
        print(reaction2.reaction)
        differences.append("metabolites")

    if not len(differences):
        return "Identical"
    else:
        return "Unequal {0}".format(" and ".join(differences))
        


if __name__ == '__main__':
    # sco4, iKS1317, iAA1259 = get_models()
    # compare_genes()
    # list_genes_reactions_metabolites()
    compare_kim_content_of_models2()
    # compare_kim_content_of_models()