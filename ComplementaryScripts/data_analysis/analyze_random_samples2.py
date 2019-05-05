# -*- coding: utf-8 -*-
"""
This file analyze the random samples from the enzyme-constrained models 
created by Eduard Kerkhoven by using the Raven Toolbox. 

Author: Snorre Sulheim
Created: 12.03.2019
Email: snorre.sulheim@sintef.no
"""
import pandas as pd
from pathlib import Path
import cobra
from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns

pd.set_option("display.max_rows", 150)
pd.set_option("display.max_columns", 101)

matplotlib.rcParams.update({'font.size': 18})


SELECTED_PATHWAYS_M145 = ["Inositol phosphate metabolism", "Glycolysis/Gluconeogenesis", "Pyruvate metabolism", "Oxidative phosphorylation", "Glycine, serine and threonine metabolism",
                          "Pentose and glucuronate interconversions", "Citric Acid Cycle", "Alanine, aspartate and glutamate metabolism", "Pentose phosphate pathway", "Glyoxylate and dicarboxylate metabolism",
                          "Undecylprodigiosin Biosynthesis", "Calcium-Dependent Antibiotics Biosynthesis", "Actinorhodin Biosynthesis", "Cpk biosynthesis"]#, 
                          # "Fatty acid biosynthesis", "Glycerophospholipid metabolism"]

SELECTED_PATHWAYS_M1152 = ["Inositol phosphate metabolism", "Glycolysis/Gluconeogenesis", "Pyruvate metabolism", "Oxidative phosphorylation", "Glycine, serine and threonine metabolism",
                          "Pentose and glucuronate interconversions", "Citric Acid Cycle", "Alanine, aspartate and glutamate metabolism", "Pentose phosphate pathway", "Glyoxylate and dicarboxylate metabolism"]
# Note, if uptake is used to normalize, glyoxylate is no longer in the top 10 (replaced by cpk, then propanoate and butanoate)
                          

def pathway_analysis_plot(model_fn, random_samples_fn, selected_pwys, row_order = None, labels = None, mask_rows = None):
    model = cobra.io.read_sbml_model(model_fn)
    M145_df, M1152_df = get_random_samples(random_samples_fn, model, key = "pathway")
    
    #M145    
    M145_sub_df  = M145_df.groupby("Subsystem").sum()
    M145_sub_df = M145_sub_df.loc[selected_pwys, :]
    M145_sub_df.index.name = None
    # M145_sub_df.columns = [x.replace("_", " - ")+"h" for x in M145_sub_df.columns.values]
    M145_sub_df.columns = [x.split("_")[1] for x in M145_sub_df.columns.values]

    # M1152
    M1152_sub_df  = M1152_df.groupby("Subsystem").sum()
    M1152_sub_df = M1152_sub_df.loc[selected_pwys, :]
    M1152_sub_df.index.name = None
    # M1152_sub_df.columns = [x.replace("_", " - ")+"h" for x in M1152_sub_df.columns.values]
    M1152_sub_df.columns = [x.split("_")[1] for x in M1152_sub_df.columns.values]
    

    if mask_rows:
        M1152_sub_df.iloc[mask_rows, :] = None
    
    if row_order:
        M1152_sub_df = M1152_sub_df.iloc[row_order, :]
    
    
    # Plot Settings
    # fig, [ax_M145, ax_M1152] = plt.subplots(1, 2, sharey = True)
    figsize = (14, 14)
    cmap = sns.cm.vlag
    cmap.set_bad("gray", alpha = 0)

    g_M145 = sns.clustermap(M145_sub_df, cmap = cmap, z_score = 0, col_cluster = False, vmin = -2, vmax = 2, figsize = figsize)
    # plt.setp(g_M145.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.subplots_adjust(left = 0.03, right = 0.5)
    plt.savefig("C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Eduard random sampling/M145/all_pathways_uptake.svg")
    # plt.savefig("C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Eduard random sampling/M145/all_pathways.png")
    # plt.show()
    plt.close()

    # Set hatches
    g_M1152 = sns.clustermap(M1152_sub_df, cmap = cmap, z_score = 0, col_cluster = False, row_cluster = False, vmin = -2, vmax = 2, figsize = figsize)
    g_M1152.ax_heatmap.patch.set(hatch='//', edgecolor='black')
    # plt.setp(g_M1152.ax_heatmap.xaxis.get_majorticklabels(), rotation=45)
    plt.subplots_adjust(left = 0.03, right = 0.5)
    # plt.savefig("C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Eduard random sampling/M1152/all_pathways.png")
    plt.savefig("C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Eduard random sampling/M1152/all_pathways_uptake.svg")
    # plt.show()
    plt.close()
    

def pathway_analysis(model_fn, random_samples_fn, selected_pwys, strain = "M145", row_order = None, labels = None, mask_rows = None, store = True):
    model = cobra.io.read_sbml_model(model_fn)
    M145_df, M1152_df = get_random_samples(random_samples_fn, model, key = "pathway")
    if strain == "M145":
        sub_df  = M145_df.groupby("Subsystem").sum()
    else:
        sub_df  = M1152_df.groupby("Subsystem").sum()

    # print(sub_df)
    # print(sub_df.std(axis = 1))

    std_df = sub_df.std(axis = 1).sort_values(ascending = False)
    print(std_df)
    if store:
        store_df = sub_df.copy()
        store_df["std"] = sub_df.std(axis = 1)
        norm_type = random_samples_fn.rsplit("_")[1].split(".")[0]
        csv_fn = """C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/
                    scoGEM/random sampling/Eduard random sampling/{0}/sub_df_{1}.svg""".format(strain, norm_type)
        store_df.to_csv(csv_fn)
    

    sub_df = sub_df.loc[selected_pwys, :]
    sub_df.index.names = ["Pathways"]

    if labels:
        print(labels)
        sub_df = sub_df.loc[:, labels]

    if mask_rows:
        sub_df.iloc[mask_rows, :] = None

    
    # Settings
    figsize = (12, 12)
    cmap = sns.cm.vlag
    cmap.set_bad("gray", alpha = 0)
    if row_order:
        sub_df = sub_df.iloc[row_order, :]
        g = sns.clustermap(sub_df, cmap = cmap, z_score = 0, col_cluster = False, row_cluster = False, figsize = figsize, vmin = -2, vmax = 2)
    else:
        g = sns.clustermap(sub_df, cmap = cmap, z_score = 0, col_cluster = False, figsize = figsize, vmin = -2, vmax = 2)
        print("Row order:", g.dendrogram_row.reordered_ind)


    # Set hatches
    g.ax_heatmap.patch.set(hatch='//', edgecolor='black')
    plt.subplots_adjust(left = 0.03, right = 0.85)
    plt.show()

def subsystem_analysis(model_fn, random_samples_fn, strain = "M145", row_order = None, labels = None, store = True):
    model = cobra.io.read_sbml_model(model_fn)
    M145_df, M1152_df = get_random_samples(random_samples_fn, model)
    if strain == "M145":
        sub_df  = M145_df.groupby("Subsystem").sum()
    else:
        sub_df  = M1152_df.groupby("Subsystem").sum()


    # print(sub_df.columns)
    # print(sub_df.std(axis = 1))
    # Store the dataframe and the std
    std_df = sub_df.std(axis = 1).sort_values(ascending = False)
    print(std_df)
    if store:
        store_df = sub_df.copy()
        store_df["std"] = sub_df.std(axis = 1)
        norm_type = random_samples_fn.rsplit("_")[1].split(".")[0]
        csv_fn = """C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/
                    scoGEM/random sampling/Eduard random sampling/{0}/sub_df_{1}.svg""".format(strain, norm_type)
        store_df.to_csv(csv_fn)
    
    if labels:
        print(labels)
        sub_df = sub_df.loc[:, labels]
    
    # Settings
    figsize = (24, 12)
    if row_order:
        sub_df = sub_df.iloc[row_order, :]
        g = sns.clustermap(sub_df, cmap = "vlag", z_score = 0, col_cluster = False, row_cluster = False, figsize = figsize, vmin = -2, vmax = 2)
    else:
        g = sns.clustermap(sub_df, cmap = "vlag", z_score = 0, col_cluster = False, figsize = figsize, vmin = -2, vmax = 2)
        print("Row order:", g.dendrogram_row.reordered_ind)
    plt.subplots_adjust(left = 0.03, right = 0.85)
    plt.show()

def get_random_samples(fn, model, key = "subsystem"):
    df = pd.read_csv(fn, index_col = 0, header = 0)
    df_data = df.iloc[:, 2:].abs()
    subsystem_df = add_subsystem_to_df(df_data, model, key)

    M1152_columns = [x for x in df.columns.values if "M1152" in x] + ["Subsystem"]
    M1152_df = subsystem_df[M1152_columns]
    M145_columns = [x for x in df.columns.values if "M145" in x] + ["Subsystem"]
    M145_df = subsystem_df[M145_columns]
    return M145_df, M1152_df


def add_subsystem_to_df(df, model, key = "subsystem"):
    """
    df is a dataframe where the indexes are reaction ids and the columns are different timepoints / strains
    the values are the different fluxes
    """
    df["Subsystem"] = None
    keep_columns = list(df.columns.values)
    add_rows = []
    df["Keep"] = True
    for r_id, row in df.iterrows():
        try:
            r = model.reactions.get_by_id(r_id)
        except KeyError:
            df.loc[r_id, "Keep"] = False
            # print(r_id)
            continue
        try:
            subsystem = r.annotation[key]
        except KeyError:
            df.loc[r_id, "Keep"] = False
            print("Missing subsystem: {0}".format(r_id))
            continue

        if isinstance(subsystem, str) and len(subsystem):
            df.loc[r_id, "Subsystem"] = subsystem
            # print(r_id, subsystem)
        else:
            print("Multiple subsystem annotations for: ", r.id)
            df.loc[r_id, "Keep"] = False

    df = df.loc[df["Keep"], :]
    df_to_keep = df.loc[:, keep_columns]
    # df_to_keep.loc[:, keep_columns[:-1]] = df_to_keep.loc[:, keep_columns[:-1]].abs()/df_to_keep.loc[:, keep_columns[:-1]].abs().sum()
    return df_to_keep

def compare_pathway_ratios():
    pass

if __name__ == '__main__':
    co2_normalized_random_samples_fn = "C:/Users/snorres/Google Drive/scoGEM community model/Supporting information/Model/ec-RandSampComb_AllSamples_co2.csv"
    uptake_normalized_random_samples_fn = "C:/Users/snorres/Google Drive/scoGEM community model/Supporting information/Model/ec-RandSampComb_AllSamples_uptake.csv"
    non_normalized_random_samples_fn = "C:/Users/snorres/Google Drive/scoGEM community model/Supporting information/Model/ec-RandSampComb_AllSamples_none.csv"

    model_fn = "../../ModelFiles/xml/scoGEM.xml"
    if 0:
        labels = ["M145_21", "M145_29", "M145_33", "M145_37", "M145_41", "M145_45", "M145_49", "M145_53", "M145_57"]  #
        subsystem_analysis(model_fn, co2_normalized_random_samples_fn, strain = "M145", labels = labels)
    if 0:
        labels = ["M1152_33", "M1152_41", "M1152_45", "M1152_49", "M1152_53", "M1152_57", "M1152_61", "M1152_65"] # 
        # row_order = [2, 8, 10, 6, 13, 4, 14, 1, 9, 5, 7, 0, 12, 3, 11]
        # row_order = [6, 13, 2, 10, 0, 5, 7, 9, 14, 1, 4, 8, 11, 3, 12]
        row_order = None
        subsystem_analysis(model_fn, co2_normalized_random_samples_fn, strain = "M1152", labels = labels, row_order = row_order)

    if 0:
        labels = ["M145_21", "M145_29", "M145_33", "M145_37", "M145_41", "M145_45", "M145_49", "M145_53", "M145_57"]  #,
        subsystem_analysis(model_fn, uptake_normalized_random_samples_fn, strain = "M145", labels = labels)
    if 0:
        labels = ["M1152_33", "M1152_41", "M1152_45", "M1152_49", "M1152_53", "M1152_57", "M1152_61", "M1152_65"] #  
        # row_order = [3, 0, 14, 4, 1, 9, 5, 7, 8, 12, 6, 13, 2, 10, 11] # all x33
        row_order = [3, 12, 8, 0, 14, 4, 1, 7, 5, 9, 6, 13, 2, 10, 11]
        subsystem_analysis(model_fn, uptake_normalized_random_samples_fn, strain = "M1152", labels = labels, row_order = row_order)


    if 1:
        pathway_analysis(model_fn, co2_normalized_random_samples_fn, SELECTED_PATHWAYS_M145, strain = "M145")
        # row_order = [7, 6, 2, 9, 12, 3, 10, 13, 11, 0, 1, 5, 4, 8]
        # mask_rows = [10,11,12,13]
        # pathway_analysis(model_fn, co2_normalized_random_samples_fn, SELECTED_PATHWAYS_M145, strain = "M1152", row_order = row_order, mask_rows = mask_rows)

        # pathway_analysis(model_fn, uptake_normalized_random_samples_fn, SELECTED_PATHWAYS_M145, strain = "M145")
        # pathway_analysis(model_fn, uptake_normalized_random_samples_fn, SELECTED_PATHWAYS_M1152, strain = "M1152")
        # pathway_analysis(model_fn, uptake_normalized_random_samples_fn, SELECTED_PATHWAYS, strain = "M145")
        # pathway_analysis(model_fn, non_normalized_random_samples_fn, SELECTED_PATHWAYS, strain = "M145")
    if 0:
        row_order = [7, 6, 2, 9, 12, 3, 10, 13, 11, 0, 1, 5, 4, 8]
        mask_rows = [10,11,12,13]
        pathway_analysis_plot(model_fn, co2_normalized_random_samples_fn, SELECTED_PATHWAYS_M145, row_order = row_order, mask_rows = mask_rows)
        # pathway_analysis_plot(model_fn, uptake_normalized_random_samples_fn, SELECTED_PATHWAYS_M145, row_order = row_order, mask_rows = mask_rows)
