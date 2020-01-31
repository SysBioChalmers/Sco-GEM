# -*- coding: utf-8 -*-
"""
This file analyze the random samples from the enzyme-constrained models 
created by Eduard Kerkhoven by using the Raven Toolbox.

This file is more a collection of snippets used to create figures for the manuscript than anything else. It is poorly commented / documented.


Author: Snorre Sulheim
Created: 12.03.2019
Modified: 30.09.2019
Email: snorre.sulheim@sintef.no
"""
import pandas as pd
from pathlib import Path
import cobra
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns

pd.set_option("display.max_rows", 150)
pd.set_option("display.max_columns", 101)


# MAT
matplotlib.rcParams.update({'font.size': 10, 'legend.loc':'upper right'})


# Definition of different groups of pathways used om plotting

SELECTED_PATHWAYS_M145 = ["Glycolysis/Gluconeogenesis", "Pyruvate metabolism", "Oxidative phosphorylation", "Glycine, serine and threonine metabolism",
                          "Citric Acid Cycle", "Alanine, aspartate and glutamate metabolism", "Pentose phosphate pathway",
                          "Fatty acid biosynthesis", "Valine, leucine and isoleucine degradation", "Nucleotide biosynthesis",
                          "Undecylprodigiosin Biosynthesis", "Calcium-Dependent Antibiotics Biosynthesis", "Actinorhodin Biosynthesis", "Coelimycin biosynthesis"]#, 
                          # "Fatty acid biosynthesis", "Glycerophospholipid metabolism"]

SELECTED_PATHWAYS_M1152 = ["Glycolysis/Gluconeogenesis", "Pyruvate metabolism", "Oxidative phosphorylation", "Glycine, serine and threonine metabolism",
                          "Citric Acid Cycle", "Alanine, aspartate and glutamate metabolism", "Pentose phosphate pathway", "Glyoxylate and dicarboxylate metabolism",
                          "Fatty acid biosynthesis", "Valine, leucine and isoleucine degradation", "Nucleotide biosynthesis"]

NUCLEOTIDE_METABOLISM = ["Purine metabolism", "Pyrimidine metabolism"]                          

CENTRAL_CARBON =  ["Glycolysis/Gluconeogenesis", 
                  "Citric Acid Cycle",
                  "Pentose phosphate pathway",
                  "Fructose and mannose metabolism",
                  "Alternate Carbon metabolism",
                  "Pyruvate metabolism",
                  "Pentose and glucuronate interconversions",
                  "Glycogen metabolism",
                  "Peptidoglycan biosynthesis",
                  "Murein Recycling",
                  "Murein Biosynthesis",
                  "Glyoxylate and dicarboxylate metabolism",
                  "Starch and sucrose metabolism",
                  "Butanoate metabolism",
                  "Propanoate metabolism",
                  "Carbon fixation pathways in prokaryotes",
                  "Amino sugar and nucleotide sugar metabolism",
                  "C5-Branched dibasic acid metabolism",
                  "Inositol phosphate metabolism",
                ]

AMINO_ACIDS = ["Alanine, aspartate and glutamate metabolism", 
               "Glycine, serine and threonine metabolism",
               "Valine, leucine and isoleucine degradation",
               "Phenylalanine, tyrosine and tryptophan biosynthesis",
               "Valine, leucine and isoleucine biosynthesis", 
               "Arginine and proline metabolism", 
               "Cysteine and methionine metabolism",
               "Lysine degradation",
               "Arginine biosynthesis",
               "Cysteine metabolism",
               "Lysine biosynthesis",
               "Histidine metabolism",
               "Phenylalanine metabolism",
               "Methionine metabolism",
               "Thiamine metabolism",
               "Tyrosine metabolism",
               "Tryptophan metabolism",
               "D-Alanine metabolism",
               "Ascorbate and aldarate metabolism",
               "D-Glutamine and D-glutamate metabolism",
               "Nitrogen metabolism"]

LIPIDS =  ["Fatty acid biosynthesis",
          "Glycerolipid metabolism",
          "Fatty acid degradation",
          "Glycerophospholipid metabolism"]

VITAMINS_COFACTORS =  ["Vitamin B6 metabolism",
                      "Pantothenate and CoA biosynthesis",
                      "beta-Alanine metabolism",
                      "Nicotinate and nicotinamide metabolism",
                      "Folate biosynthesis",
                      "Porphyrin and chlorophyll metabolism",
                      "Riboflavin metabolism",
                      "Folate metabolism",
                      "One carbon pool by folate",
                      "Lipoic acid metabolism",
                      "Biotin metabolism"]

TOXIC_COMPOUND_DEGRADATION = ["Atrazine degradation",
                              "Benzoate degradation",
                              "Aminobenzoate degradation",
                              ]

ALL_BGC =  ["Terpenoid backbone biosynthesis",
           "Desferrioxamine biosynthesis",
           "Germicidin Biosynthesis",
           "Sesquiterpenoid and triterpenoid biosynthesis",
           "Ubiquinone and other terpenoid-quinone biosynthesis",
           "Albaflavenol biosynthesis",
           "Geosmin synthase",
           "Biflaviolin synthase",
           "2-methylisoborneol synthase"]

REMOVED_BGC = ["Undecylprodigiosin Biosynthesis",
               "Actinorhodin Biosynthesis",
               "Calcium-Dependent Antibiotics Biosynthesis",
               "Coelimycin biosynthesis"]

OXIDATIVE_STRESS = ["Mycothiol metabolism",
                    "Catalase",
                    "Superoxide dismutase"]


IGNORE = ["D-Arginine and D-ornithine metabolism",
          "Toluene degradation",
          "Polyketide sugar unit biosynthesis",
          "Galactose metabolism",
          "Biosynthesis of type II polyketide backbone",
          "Glutathione metabolism"
          ]

OTHER = ["Nucleotide biosynthesis",
         "Methane metabolism",
         "Sulfur metabolism",
         "Oxidative phosphorylation",
         "NADH repair",
         "tRNA Charging",
         ]

def pathway_analysis_plot(model_fn, random_samples_fn, selected_pwys, row_order = None, labels = None, mask_rows = None, sep = "\t", key = "pathway", absolute_values = True):
    """
    Plot random sampling heatmap for both M145 and M1152 in idividual plots
    """
    model = cobra.io.read_sbml_model(model_fn)
    df = get_random_samples(random_samples_fn, model, key = key, sep = sep, column_key = "MEAN", absolute_values = absolute_values)
    sub_df = df.groupby(key).sum()

    # Remove all zero rows
    sub_df = sub_df.loc[(sub_df != 0).any(axis = 1), :]

    # Standardize
    standardized_df = sub_df.subtract(sub_df.mean(axis = 1), axis = 0).divide(sub_df.std(axis = 1), axis = 0)
    
    #M145
    M145_columns = [x for x in standardized_df.columns if "M145" in x]
    M145_sub_df = standardized_df.loc[selected_pwys, M145_columns]
    M145_sub_df.index.name = None
    # M145_sub_df.columns = [x.replace("_", " - ")+"h" for x in M145_sub_df.columns.values]
    M145_sub_df.columns = [x.split("_")[-1] for x in M145_sub_df.columns.values]

    # M1152
    M1152_columns = [x for x in standardized_df.columns if "M1152" in x]
    M1152_sub_df = standardized_df.loc[selected_pwys, M1152_columns]
    M1152_sub_df.index.name = None
    # M1152_sub_df.columns = [x.replace("_", " - ")+"h" for x in M1152_sub_df.columns.values]
    M1152_sub_df.columns = [x.split("_")[-1] for x in M1152_sub_df.columns.values]
    

    if mask_rows:
        M1152_sub_df.iloc[mask_rows, :] = None
    

    if absolute_values:
        abs_string = "abs"
    else:
        abs_string = "real"

    # Plot Settings
    # fig, [ax_M145, ax_M1152] = plt.subplots(1, 2, sharey = True)
    figsize = (14, 14)
    cmap = sns.cm.vlag
    cmap.set_bad("gray", alpha = 0)
    print(M145_sub_df)
    g_M145 = sns.clustermap(M145_sub_df, cmap = cmap, col_cluster = False, vmin = -1.3, vmax = 1.3, figsize = figsize)
    print("Row order:", g_M145.dendrogram_row.reordered_ind)
    # plt.setp(g_M145.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.subplots_adjust(left = 0.03, right = 0.5)

    plt.savefig("C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Random sampling july/M145/all_pathways_uptake_{0}.svg".format(abs_string))
    # plt.savefig("C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Eduard random sampling/M145/all_pathways.png")
    # plt.show()
    plt.close()

    if row_order:
        M1152_sub_df = M1152_sub_df.iloc[g_M145.dendrogram_row.reordered_ind, :]
    # Set hatches
    g_M1152 = sns.clustermap(M1152_sub_df, cmap = cmap, col_cluster = False, row_cluster = False, vmin = -1.3, vmax = 1.3, figsize = figsize)
    g_M1152.ax_heatmap.patch.set(hatch='//', edgecolor='black')
    # plt.setp(g_M1152.ax_heatmap.xaxis.get_majorticklabels(), rotation=45)
    plt.subplots_adjust(left = 0.03, right = 0.5)
    # plt.savefig("C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Eduard random sampling/M1152/all_pathways.png")
    plt.savefig("C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Random sampling july/M1152/all_pathways_uptake_{0}.svg".format(abs_string))
    # plt.show()
    plt.close()
    

def pathway_analysis(model_fn, random_samples_fn, selected_pwys, strain = "M145", row_order = None, labels = None, 
                     mask_rows = None, store = True, key = "pathway", sep = ",", absolute_values = True,
                     show_cbar = False, mask_cells = None):
    """
    Plot pathway heatmap for either M145, M1152 or both in one and same plot
    """
    model = cobra.io.read_sbml_model(model_fn)
    df = get_random_samples(random_samples_fn, model, key = key, sep = sep, column_key = "MEAN", absolute_values = absolute_values)
    sub_df = df.groupby(key).sum()

    # Remove all zero rows
    sub_df = sub_df.loc[(sub_df.abs() > 1e-8).any(axis = 1), :]
    print(sub_df.loc[(sub_df.abs() < 1e-8).any(axis = 1), :].abs().max(axis = 1).sort_values(ascending = False))
    # Standardize
    standardized_df = sub_df.subtract(sub_df.mean(axis = 1), axis = 0).divide(sub_df.std(axis = 1), axis = 0)
   
    # print(sub_df.mean(axis = 1).sort_values(ascending = False))
    std_df = sub_df.iloc[:, 1:].std(axis = 1).sort_values(ascending = False)
    # rstd_df = sub_df.std(axis = 1).divide(sub_df.mean(axis = 1), axis = 0).sort_values(ascending = False)
    print(std_df)
    # print(rstd_df)

    if store:
        store_df = sub_df.copy()
        store_df["std"] = sub_df.std(axis = 1)
        norm_type = random_samples_fn.rsplit("_")[-1].split(".")[0]
        csv_fn = "C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/{0}/sub_df_{1}.svg".format(strain, norm_type)
        store_df.to_csv(csv_fn)
    
    strain_columns = [x for x in standardized_df.columns if strain in x]
    strain_df = standardized_df.loc[:, strain_columns]

    # Select top X
    # selected_pwys = list(std_df.index)[:70]
    if not selected_pwys:
        selected_pwys = strain_df.index
    selected_df = strain_df.loc[selected_pwys, :]
    selected_df.index.names = ["Pathways"]


    if labels:
        print(labels)
        selected_df = selected_df.loc[:, labels]

    # if isinstance(mask_cells, np.ndarray):
    #     print(selected_df.shape, mask_cells.shape)
    #     selected_df[mask_cells] = None

    selected_df.columns = ["-".join(x.split("_")[1:]) for x in selected_df.columns.values]
    
    # Settings
    figsize = (10, len(selected_pwys)/3)
    cmap = sns.cm.vlag
    cmap.set_bad("gray", alpha = 0)
    print(selected_df)
    if row_order:
        selected_df = selected_df.iloc[row_order, :]
        g = sns.clustermap(selected_df, cmap = cmap, col_cluster = False, row_cluster = False, figsize = figsize, vmin = -1.2, vmax = 1.2, yticklabels = True)
    else:
        g = sns.clustermap(selected_df, cmap = cmap, col_cluster = False, figsize = figsize, vmin = -1.2, vmax = 1.2, yticklabels = True, mask=mask_cells)
        print("Row order:", g.dendrogram_row.reordered_ind)
        
    # Remove labels
    ax = g.ax_heatmap
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Show colorbar
    g.cax.set_visible(show_cbar)

    # Set hatches
    g.ax_heatmap.patch.set(hatch='//', edgecolor='black')
    # .xticks(rotation=45)
    #plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha = "center")
    plt.subplots_adjust(left = 0.03, right = 0.60, bottom = 0.27, top = 0.97)
    plt.show()


def get_random_samples(fn, model, key = "subsystem", sep = "\t", column_key = "MEAN", skipcolumns = 0, absolute_values = True):
    """
    Reads the random samples and performs some data preprocessing. 
    """
    df = pd.read_csv(fn, index_col = 0, header = 0, sep = sep)
    df_data = df.iloc[:, skipcolumns:]
    
    # A few reactions are changed sign of to make the total flux even out with a correlated "looped" reaction
    # PGM and PGK is just defined in the opposite direction in the model

    df_data.loc[["SUCOAS", "PGM", "PGK"], :] = df_data.loc[["SUCOAS", "PGM", "PGK"], :].abs()
    df_data.loc["HACD1", :] *=-1

    
    if absolute_values:
        df_data = df_data.abs()
    if column_key:
        columns = [x for x in df.columns.values if column_key in x]
        df_data = df_data.loc[:, columns]
        

    subsystem_df = add_subsystem_to_df(df_data, model, key)
    return subsystem_df

def get_random_samples_raw(fn, sep = "\t", column_key = "MEAN", skipcolumns = 0):
    """
    Reads the random sampling data without the preprocessing
    """
    df = pd.read_csv(fn, index_col = 0, header = 0, sep = sep)
    df_data = df.iloc[:, skipcolumns:]

    if column_key:
        columns = [x for x in df.columns.values if column_key in x]
        df_data = df_data.loc[:, columns]
    return df_data



def get_random_samples_and_split(fn, model, key = "subsystem"):
    df = pd.read_csv(fn, index_col = 0, header = 0, sep = "\t")
    df_data = df.iloc[:, 2:].abs()
    subsystem_df = add_subsystem_to_df(df_data, model, key)

    M1152_columns = [x for x in df.columns.values if "MEAN_M1152" in x] + [key]
    M1152_df = subsystem_df[M1152_columns]
    M145_columns = [x for x in df.columns.values if "MEAN_M145" in x] + [key]
    M145_df = subsystem_df[M145_columns]
    return M145_df, M1152_df


def add_subsystem_to_df(df, model, key = "subsystem"):
    """
    df is a dataframe where the indexes are reaction ids and the columns are different timepoints / strains
    the values are the different fluxes
    """
    df[key] = None
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
            if "Valine, leucine and isoleucine" in subsystem:
                # NOTE FOR REVIEW: Undid this merging 
                # subsystem = "Valine, leucine and isoleucine metabolism"
                pass
            # elif subsystem in OTHER_AMINO_ACIDS:
            #     subsystem = "Metabolism of other amino acids"

            elif subsystem == "FERI metabolism":
                # FERI metabolism is only one reaction, FNOR, which recycles NADP with ferredoxin used by AKGDH2 reaction
                subsystem = "Citric Acid Cycle"
            elif subsystem in NUCLEOTIDE_METABOLISM:
                subsystem = "Nucleotide biosynthesis"
            elif r_id == "CDAS11":
                subsystem = "Fatty acid degradation"
            elif r.id == "MALTHIK":
                subsystem = "Glyoxylate and dicarboxylate metabolism"
            elif r.id == "HEX1":
                subsystem = "Glycolysis/Gluconeogenesis"
            elif r.id == "HACD1":
                subsystem = "Butanoate metabolism"



            df.loc[r_id, key] = subsystem
            # print(r_id, subsystem)
        else:
            print("Multiple subsystem annotations for: ", r.id)
            df.loc[r_id, "Keep"] = False

    df = df.loc[df["Keep"], :]
    df_to_keep = df.loc[:, keep_columns]
    # df_to_keep.loc[:, keep_columns[:-1]] = df_to_keep.loc[:, keep_columns[:-1]].abs()/df_to_keep.loc[:, keep_columns[:-1]].abs().sum()
    return df_to_keep

def make_contribution_subsystems_plots(model_fn, random_samples_fn, sep = "\t", key = "pathway", absolute_values = False):
    """
    Used for data QC. Plots histogram to see which reactions that have the largest contribution to each pathway
    """
    model = cobra.io.read_sbml_model(model_fn)
    df = get_random_samples(random_samples_fn, model, key = key, sep = sep, column_key = "MEAN", absolute_values = absolute_values)
    
    if absolute_values:
        abs_string = "abs"
    else:
        abs_string = "real"


    for pathway, pwy_df in df.groupby(key):
        fig, ax = plt.subplots(1)
        pwy_df.drop(key, axis = 1, inplace = True)
        # pwy_df = pwy_df.abs()
        pwy_df.plot(kind = "bar", logy = absolute_values, title = pathway, figsize = (20, 14), ax = ax).legend(loc = 1)
        ax.set_ylim(1.1*pwy_df.min().min(), 1.1*pwy_df.max().max())
        pathway = pathway.replace("/", "-")
        fn = "C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Random sampling july/pathway contribution plot/{0}_{1}.svg".format(pathway, abs_string)
        fig.savefig(fn)
        plt.close()


def analyze_accoa_consumption(model_fn, random_samples_fn, sep = "\t", subplot = True):
    metabolite_id = "accoa_c"
    model = cobra.io.read_sbml_model(model_fn)
    df = get_random_samples_raw(random_samples_fn, sep = sep, column_key = "MEAN")
    
    metabolite = model.metabolites.get_by_id(metabolite_id)

    reaction_ids = []
    reaction_multiplier = []
    for r in metabolite.reactions:
        reaction_ids.append(r.id)
        reaction_multiplier.append(r.get_coefficient(metabolite_id))

    df_selected = (df.loc[reaction_ids, :].T*reaction_multiplier).T
    strain_annotation = [x.split("_")[1] for x in df_selected.columns.values]
    timepoint_annotation = [1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8]
    print(len(strain_annotation), len(timepoint_annotation))

    # Sum ACTCOALIG anf PTAr
    print(df_selected)

    df_all_0 = (df_selected.abs() < 1e-4).all(axis = 1)
    df_selected = df_selected.loc[~df_all_0, :].T

    # ACTOALIG and PTAr is creating a loop
    # df_selected["ACTCOALIG - PTAr"] = df_selected["ACTCOALIG"] + df_selected["PTAr"]
    # df_selected.drop(["PTAr", "ACTCOALIG"], axis = 1, inplace = True)


    df_selected["Strain"] = strain_annotation
    df_selected["Time point"] = timepoint_annotation
    df_long_form = pd.melt(df_selected, id_vars = ["Strain", "Time point"], value_vars = df_selected.columns.values[:-2], value_name = "Flux")

    # Production
    df_production = df_long_form.loc[df_long_form["Flux"] > 1e-8, :]
    df_consumption = df_long_form.loc[df_long_form["Flux"] < -1e-8, :]
    df_consumption["Flux"] *= -1

    if subplot:
        fig, [ax1, ax2] = plt.subplots(2)
    else:
        fig1, ax1 = plt.subplots(1)
        fig2, ax2 = plt.subplots(1)

    sns.lineplot(data = df_production, x = "Time point", y = "Flux", hue = "rxns", style = "Strain", ax = ax1)
    ax1.set_title("Acetyl-CoA production")
    ax1.set_yscale("log")

    sns.lineplot(data = df_consumption, x = "Time point", y = "Flux", hue = "rxns", style = "Strain", ax = ax2)
    ax2.set_title("Acetyl-CoA consumption")
    ax2.set_yscale("log")
    plt.show()


def analyze_malcoa_consumption(model_fn, random_samples_fn, sep = "\t", subplot = True):
    """
    Make two plots: One of the reactions producing malonyl-CoA and on for the reactions consuming Malonyl-CoA
    """
    metabolite_id = "malcoa_c"
    model = cobra.io.read_sbml_model(model_fn)
    df = get_random_samples_raw(random_samples_fn, sep = sep, column_key = "MEAN")
    
    metabolite = model.metabolites.get_by_id(metabolite_id)

    reaction_ids = []
    reaction_multiplier = []
    for r in metabolite.reactions:
        reaction_ids.append(r.id)
        reaction_multiplier.append(r.get_coefficient(metabolite_id))

    df_selected = (df.loc[reaction_ids, :].T*reaction_multiplier).T
    strain_annotation = [x.split("_")[1] for x in df_selected.columns.values]
    timepoint_annotation = [1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8]
    print(len(strain_annotation), len(timepoint_annotation))



    df_all_0 = (df_selected.abs() < 1e-10).all(axis = 1)
    df_selected = df_selected.loc[~df_all_0, :].T

    # Sum RED, ACT
    red_columns = [x for x in df_selected.columns if x[:4] == "REDS"]
    df_selected["RED"] = df_selected.loc[:, red_columns].sum(axis = 1)
    df_selected.drop(red_columns, axis = 1, inplace = True)

    cpk_columns = [x for x in df_selected.columns if x[:4] == "CPKS"]
    df_selected["CPK"] = df_selected.loc[:, cpk_columns].sum(axis = 1)
    df_selected.drop(cpk_columns, axis = 1, inplace = True)

    cda_columns = [x for x in df_selected.columns if x[:4] == "CDAS"]
    df_selected["CDA"] = df_selected.loc[:, cda_columns].sum(axis = 1)
    df_selected.drop(cda_columns, axis = 1, inplace = True)

    act_columns = [x for x in df_selected.columns if x[:4] == "ACTS"]
    df_selected["ACT"] = df_selected.loc[:, act_columns].sum(axis = 1)
    df_selected.drop(act_columns, axis = 1, inplace = True)

    # SUM ACCOAC and MCOATA
    # df_selected["ACCOAC + MCOATA + ACCOAC_1"] = df_selected["ACCOAC"] + df_selected["MCOATA"] + df_selected["ACCOAC"]
    # df_selected.drop(["ACCOAC", "MCOATA"], axis = 1, inplace = True)



    df_selected["Strain"] = strain_annotation
    df_selected["Time point"] = timepoint_annotation
    print(df_selected)
    df_selected = df_selected.reindex(["ACCOAC", "ACCOAC_1", "ACT", "MCOATA", "CPK", "RED", "CDA", "THYDNAPS", "Strain", "Time point"], axis = 1)
    df_long_form = pd.melt(df_selected, id_vars = ["Strain", "Time point"], value_vars = df_selected.columns.values[:-2], value_name = "Normalized flux")
    
    # Production
    df_production = df_long_form.loc[df_long_form["Normalized flux"] > 1e-8, :]
    df_consumption = df_long_form.loc[df_long_form["Normalized flux"] < -1e-8, :]
    df_consumption["Normalized flux"] *= -1
    print(df_consumption)

    if subplot:
        fig, [ax1, ax2] = plt.subplots(1, 2, sharey = True)
    else:
        fig1, ax1 = plt.subplots(1)
        fig2, ax2 = plt.subplots(1)

    sns.lineplot(data = df_production, x = "Time point", y = "Normalized flux", hue = "rxns", style = "Strain", ax = ax1, lw = 2)
    ax1.set_title("Malonyl-CoA production")
    ax1.set_yscale("log")
    ax1.legend(loc='lower left')#, bbox_to_anchor=(1.01, 0.5), ncol=1)
    # ax1.spines['right'].set_visible(False)
    # ax1.spines['top'].set_visible(False)

    sns.lineplot(data = df_consumption, x = "Time point", y = "Normalized flux", hue = "rxns", style = "Strain", ax = ax2, lw = 2)
    ax2.set_title("Malonyl-CoA consumption")
    ax2.set_yscale("log")
    ax2.legend(loc='lower left')
    # ax2.spines['right'].set_visible(False)
    # ax2.spines['top'].set_visible(False)
    sns.despine()

    plt.show()

    
    print(df_consumption[df_consumption["rxns"]!="MCOATA"].groupby(["Strain", "Time point"]).sum())
    print(df_consumption[df_consumption["rxns"]=="MCOATA"].groupby(["Strain", "Time point"]).sum())


def analyze_met_consumption(model_fn, random_samples_fn, metabolite_id = "nadh_c", metabolite_name = "NADH", sep = "\t", lim = 1e-3):
    """
    Make two plots: One of the reactions producing the chosen metabolite and one for the reactions consuming the metabolite
    """
    model = cobra.io.read_sbml_model(model_fn)
    df = get_random_samples_raw(random_samples_fn, sep = sep, column_key = "MEAN")
    
    metabolite = model.metabolites.get_by_id(metabolite_id)

    reaction_ids = []
    reaction_multiplier = []
    for r in metabolite.reactions:
        reaction_ids.append(r.id)
        reaction_multiplier.append(r.get_coefficient(metabolite_id))

    df_selected = (df.loc[reaction_ids, :].T*reaction_multiplier).T
    strain_annotation = [x.split("_")[1] for x in df_selected.columns.values]
    timepoint_annotation = [1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8]
    print(strain_annotation, timepoint_annotation)



    df_all_0 = (df_selected.abs() < lim).all(axis = 1)
    df_selected = df_selected.loc[~df_all_0, :].T


    df_selected["Strain"] = strain_annotation
    df_selected["Time point"] = timepoint_annotation
    df_long_form = pd.melt(df_selected, id_vars = ["Strain", "Time point"], value_vars = df_selected.columns.values[:-2], value_name = "Flux")

    print(df_long_form)
    
    # Production
    df_production = df_long_form[df_long_form["Flux"] > 1e-8]
    df_consumption = df_long_form[df_long_form["Flux"] < -1e-8]
    df_consumption["Flux"] *= -1

    print(df_production.groupby(["Strain", "Time point"]).sum())


    fig, [ax1, ax2] = plt.subplots(2)

    sns.lineplot(data = df_production, x = "Time point", y = "Flux", hue = "rxns", style = "Strain", ax = ax1, style_order = ["M145", "M1152"])
    ax1.set_title("{0} production".format(metabolite_name))
    ax1.set_yscale("log")

    sns.lineplot(data = df_consumption, x = "Time point", y = "Flux", hue = "rxns", style = "Strain", ax = ax2, style_order = ["M145", "M1152"])
    ax2.set_title("{0} consumption".format(metabolite_name))
    ax2.set_yscale("log")
    plt.show()

def plot_key_metabolic_reactions(random_samples_fn, title = "", sep = ","):
    """
    Plot the mean predicted flux through the following reactions selected by Tjasa Kumelj
    - PFK
    - G6PDH2r
    - TAGS140
    - PDH
    - MCOATA
    - ACTS1
    - ICDHyr
    - AKGDH2
    - CYO2a
    - CYO2b
    - NADH17b
    - GLUN
    - GLUDxi
    """
    key_reactions = ["TAGS140","MCOATA","PFK", "PDH","ACTS1","AKGDH2","ICDHyr","CYO2a","CYO2b","G6PDH2r","NADH17b","GLUN","GLUDxi"]
    df = get_random_samples_raw(random_samples_fn, sep = sep, column_key = "MEAN")
    df_selected = df.loc[key_reactions, :]
    df_selected.columns = ["-".join(x.split("_")[1:]) for x in df_selected.columns.values]
    # print(df_selected)
    # Standardize
    print(df_selected.subtract(df_selected.mean(axis = 1), 0))
    df_centered = df_selected.subtract(df_selected.mean(axis = 1), 0)
    df_z_score = df_centered.div(df_centered.std(axis = 1), 0)
    # strain_annotation = [x.split("_")[1] for x in df_selected.columns.values]   
    # timepoint_annotation = [1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8]
    cmap = sns.cm.vlag
    # g = sns.clustermap(df_selected, row_cluster = False, col_cluster = False, z_score = 0, cmap = cmap, vmin = -1.5, vmax = 1.5)
    sns.heatmap(df_z_score, cmap = cmap, vmin = -1.5, vmax = 1.5, cbar_kws={'label': 'Flux Z-score'})
    plt.title(title)
    plt.subplots_adjust(bottom = 0.2)
    save_key = random_samples_fn.split("/")[-1].split(".")[0]
    plt.savefig("C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Random sampling july/key_metabolic_reactions_{0}.svg".format(save_key))    
    # plt.show()    
    plt.close()

def plot_selected_reactions(random_samples_fn, reactions, title = "", sep = ",", strain = "M145"):
    """
    Plot the flux for a list of chosen reactions
    """
    df = get_random_samples_raw(random_samples_fn, sep = sep, column_key = "MEAN")
    df_selected = df.loc[reactions, :]
    df_selected.columns = ["-".join(x.split("_")[1:]) for x in df_selected.columns.values]
    if strain:
        selected_columns = [x for x in df_selected.columns if strain in x]
        df_selected = df_selected.loc[:, selected_columns]
        df_selected.columns = [x.split("-")[1] for x in df_selected.columns]

    print(df_selected)
    
    ax = df_selected.T.plot(alpha = 0.9, logy = False)
    ax.set_ylabel("Normalized mean flux")
    ax.set_xlabel("Hours")
    plt.show()

def print_genes_in_subsystem(model_fn, pathway):
    model = cobra.io.read_sbml_model(model_fn)

    pathway_genes = []
    for r in model.reactions:
        try:
            p_r = r.annotation["pathway"]
        except:
            continue
        if p_r == pathway:
            pathway_genes += [g.id for g in r.genes]

    print("Genes related to {0};".format(pathway))
    print(list(set(pathway_genes)))

def plot_all_reactions_for_metabolite(model_fn, metabolite_id, random_samples_fn, sep = "\t"):
    model = cobra.io.read_sbml_model(model_fn)
    df = get_random_samples_raw(random_samples_fn, sep = sep, column_key = "MEAN")

    metabolite = model.metabolites.get_by_id(metabolite_id)
    reactions = [r.id for r in metabolite.reactions]

    coeff = [r.metabolites[metabolite] for r in metabolite.reactions]
    df_scaled = df.loc[reactions, :].multiply(coeff, axis = 0)
    # print(df_scaled)    

    std_df = df_scaled.std(axis = 1).sort_values(ascending = False)
    print(std_df)
    # print(df_scaled[reactions, :].head())
    df_selected = df_scaled.loc[std_df > 1e-3, :]

    # Row colors
    lut = dict(zip([1,-1], "gr"))
    row_sign = [np.sign(r.metabolites[metabolite]) for r in metabolite.reactions if r.id in df_selected.index]
    row_colors = [lut[x] for x in row_sign]

    g = sns.clustermap(df_selected, z_score = None, cmap = sns.cm.vlag, col_cluster = False, vmin = -0.1, vmax = 0.1)#row_colors = row_colors)
    plt.show()


if __name__ == '__main__':
    co2_normalized_random_samples_fn = "C:/Users/snorres/Google Drive/scoGEM community model/Supporting information/Model/randomsampling_july/ec-RandSampComb_proteomics_CO2norm.tsv"
    gluglc_normalized_random_samples_fn = "C:/Users/snorres/Google Drive/scoGEM community model/Supporting information/Model/randomsampling_july/ec-RandSampComb_proteomics_GlcGluNorm.tsv"
    growth_normalized_random_samples_fn = "C:/Users/snorres/Google Drive/scoGEM community model/Supporting information/Model/randomsampling_july/ec-RandSampComb_proteomics_growthNorm.tsv"
    all_flux_normalized_random_samples_fn = "C:/Users/snorres/Google Drive/scoGEM community model/Supporting information/Model/randomsampling_july/ec-RandSampComb_proteomics_AllFluxNorm.tsv"
    
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


    if 0:
        matplotlib.rcParams.update({'font.size': 6, 'legend.loc':'upper right'})
        row_order = [4, 58, 19, 40, 26, 66, 23, 24, 69, 72, 46, 79, 18, 54, 14, 27, 7, 5, 67, 53, 29, 21, 74, 51, 48, 8, 52, 42, 49, 63, 70, 59, 81, 76, 44, 47, 61, 22, 30, 35, 56, 64, 80, 32, 45, 38, 65, 62, 43, 13, 28, 9, 71, 73, 31, 57, 36, 37, 15, 1, 78, 0, 12, 34, 17, 77, 20, 6, 11, 2, 55, 16, 25, 60, 75, 3, 10, 41, 50, 39, 33, 68]
        # subplots adjust left = 0.03, right = 0.35, bottom = 0.27, top = 0.97
        pathway_analysis(model_fn, co2_normalized_random_samples_fn, SELECTED_PATHWAYS_M145, strain = "MEAN", sep = "\t", store = False, absolute_values = False)
        #pathway_analysis(model_fn, gluglc_normalized_random_samples_fn, SELECTED_PATHWAYS_M145, strain = "MEAN", sep = "\t", store = False, absolute_values = False, row_order = row_order)
        
        # row_order = [7, 6, 2, 9, 12, 3, 10, 13, 11, 0, 1, 5, 4, 8]
        # mask_rows = [10,11,12,13]
        # pathway_analysis(model_fn, co2_normalized_random_samples_fn, SELECTED_PATHWAYS_M145, strain = "M1152", row_order = row_order, mask_rows = mask_rows)

    if 0:
        row_order = [12, 13, 10, 11, 3, 2, 8, 6, 7, 9, 5, 4, 0, 1]
        mask_rows = [10,11,12,13]
        pathway_analysis_plot(model_fn, co2_normalized_random_samples_fn, SELECTED_PATHWAYS_M145, row_order = True, mask_rows = mask_rows, absolute_values = False)
        # pathway_analysis_plot(model_fn, uptake_normalized_random_samples_fn, SELECTED_PATHWAYS_M145, row_order = row_order, mask_rows = mask_rows)
    if 0:
        analyze_accoa_consumption(model_fn, co2_normalized_random_samples_fn, subplot = False)
    if 1:
        analyze_malcoa_consumption(model_fn, co2_normalized_random_samples_fn, subplot = True)

    if 0:
        analyze_pyruvate_consumption(model_fn, co2_normalized_random_samples_fn)
    if 0:
        analyze_glutamate_consumption(model_fn, co2_normalized_random_samples_fn) 
    if 0:
        analyze_met_consumption(model_fn, co2_normalized_random_samples_fn, "akg_c", "Glutamate", lim = 1e-2)
    if 0:
        make_contribution_subsystems_plots(model_fn, co2_normalized_random_samples_fn, absolute_values = False)
    if 0:
        plot_key_metabolic_reactions(co2_normalized_random_samples_fn, "CO2-normalized ", "\t")
        plot_key_metabolic_reactions(gluglc_normalized_random_samples_fn, "Glucose and glutamate carbon uptake normalized","\t")
        plot_key_metabolic_reactions(growth_normalized_random_samples_fn, "Growth rate normalized", "\t")
        plot_key_metabolic_reactions(all_flux_normalized_random_samples_fn, "Normalized by sum of fluxes", "\t")

    if 0:
        # plot_selected_reactions(co2_normalized_random_samples_fn, ["CS", "ICDHyr"], "Aceyl-CoA consumption into TCA cycle", "\t")
        # plot_selected_reactions(co2_normalized_random_samples_fn, ["FBA", "PFK", "ENO", "G6PDH2r"], "Aceyl-CoA consumption into TCA cycle", "\t", strain = "M145")
        plot_selected_reactions(co2_normalized_random_samples_fn, ["ILETA", "LEUTA", "VALTA", "ILEDHr", "VALDHr", "LLEUDr"], "Branched-chain amino acids", "\t", strain = "M145")
    if 0:
        # Plot Supplementary figure 3 as several subpanels
        matplotlib.rcParams.update({'font.size': 10, 'legend.loc':'upper right'})

        # Glycolysis
        #pathway_analysis(model_fn, co2_normalized_random_samples_fn, CENTRAL_CARBON, strain = "MEAN", sep = "\t", store = False, absolute_values = False)

        # Amino acids
        # pathway_analysis(model_fn, co2_normalized_random_samples_fn, AMINO_ACIDS, strain = "MEAN", sep = "\t", store = False, absolute_values = False)

        
        # for key in [OXIDATIVE_STRESS]:
        #     # Lipids
        #     pathway_analysis(model_fn, co2_normalized_random_samples_fn, key, strain = "MEAN", sep = "\t", store = False, absolute_values = False, show_cbar = False)

        # BGC
        mask_cells = np.zeros((len(REMOVED_BGC+ALL_BGC), 17), dtype = bool)
        mask_cells[0:4, 9:] = 1
        # mask_cells = mask_cells.T
        pathway_analysis(model_fn, co2_normalized_random_samples_fn, REMOVED_BGC + ALL_BGC, strain = "MEAN", sep = "\t", store = False, absolute_values = False, mask_cells = mask_cells)
    if 0:
        # print_genes_in_subsystem(model_fn, "Alanine, aspartate and glutamate metabolism")
        print_genes_in_subsystem(model_fn, "D-Glutamine and D-glutamate metabolism")
    if 1:
        plot_all_reactions_for_metabolite(model_fn, "glu__L_c", co2_normalized_random_samples_fn)