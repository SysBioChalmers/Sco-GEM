# -*- coding: utf-8 -*-
"""
Author: Snorre Sulheim
Date: 09.04.2019

# Description
Create plots of how the flux changes during the timecourse and the corresponding genes

## Plot 1
Plot the sampled flux for a given reaction for all timesteps for both strains
function: plot_all_flux()

##



"""

import cobra
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
import re


pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)

M145_FERMENTORS = ["F516", "F517", "F518"]
M1152_FERMENTORS = ["F519", "F521", "F522"]

# Delayed reactions in M1152
DELAYED_M1152 = ["MCOATA", "TAGS140"]

# Different reactions in M1152
DIFFERENT_M1152 = ["PDH", "G6PDH2r", "PGL", "CYO2a", "ACONTa", "GLUDxi"]

# Genes in pho-regulon
DELAYED_GENES = ["SCO1196", "SCO1968", "SCO2286", "SCO1906"]
NOT_RESPONDING_GENES = ["SCO1565", "SCO4229"]
EXPORTED_GDPD_GENES = ["SCO1565", "SCO1968", "SCO7550"]
CYTOSOLIC_GDPD_GENES = ["SCO1090", "SCO1419", "SCO3976", "SCO5661"]
TAT_SECRETED_GENES = ["SCO0736", "SCO1172", "SCO1196", "SCO1356", "SCO1432", "SCO1565", "SCO1590", "SCO1639", "SCO1906", "SCO2068", "SCO2286", "SCO2758", "SCO2780", "SCO2786", "SCO3484", "SCO3790", "SCO4672", "SCO6052", "SCO6198", "SCO6272", "SCO6281", "SCO6457", "SCO6580", "SCO6594", "SCO6691", "SCO7631", "SCO7677"]
CYO2a_GENES = ["SCO2150", "SCO2149", "SCO7236", "SCO2148", "SCO7120"]
ATP_SYNTHASE_GENES = ["SCO5368", "SCO5370", "SCO5373", "SCO5374"]
ALL_ATP_SYNTHASE_GENES = ["SCO5366", "SCO5367", "SCO5368", "SCO5369", "SCO5370", "SCO5371", "SCO5372", "SCO5373","SCO5374"]
INOSITOL_DEHYDROGENASE_GENES = ["SCO6255", "SCO6984", "SCO7254", "SCO1527"]
SPODM = ["SCO0999", "SCO2633", "SCO5254"]
CAT = ["SCO0379", "SCO0560", "SCO0666", "SCO2529", "SCO6204", "SCO7590"]

# Stress genes as asked for by reviewers
# From paper
STRESS_GENES_1 = ["SCO5126","SCO6635","SCO0596","SCO0167","SCO0200","SCO7299","SCO5031","SCO0560","SCO0999","SCO2633","SCO3890","SCO0641","SCO2367","SCO2368","SCO3767","SCO4277","SCO5806"] # From https://pubs.acs.org/doi/10.1021/acs.jproteome.7b00163

# From the identified differentially expressed genes
STRESS_GENES_2 = ["SCO2633", "SCO4834", "SCO4835", "SCO0185", "SCO0186", "SCO0187","SCO0188","SCO7416", "SCO7417",
                  "SCO7418", "SCO7419", "SCO7420", "SCO7421", "SCO7422", "SCO2885"]

# Glutamate related genes
ALL_GLUTAMATE_GENES = ['SCO0386', 'SCO0401', 'SCO0603', 'SCO0910', 'SCO0985', 'SCO1018', 'SCO1042', 'SCO1204', 'SCO1254', 'SCO1483', 'SCO1484', 'SCO1487', 'SCO1522', 'SCO1570', 'SCO1578', 'SCO1579', 'SCO1613', 'SCO1773', 'SCO1776', 'SCO1977', 'SCO2025', 'SCO2026', 'SCO2051', 'SCO2086', 'SCO2089', 'SCO2117', 'SCO2198', 'SCO2210', 'SCO2234', 'SCO2238', 'SCO2241', 'SCO2587', 'SCO2614', 'SCO2664', 'SCO2789', 'SCO2829', 'SCO2830', 'SCO2831', 'SCO2951', 'SCO2999', 'SCO3213', 'SCO3382', 'SCO3411', 'SCO3416', 'SCO3435', 'SCO3615', 'SCO3629', 'SCO3658', 'SCO3807', 'SCO3851', 'SCO4078', 'SCO4086', 'SCO4089', 'SCO4162', 'SCO4469', 'SCO4498', 'SCO4645', 'SCO4683', 'SCO4740', 'SCO4780', 'SCO4785', 'SCO4911', 'SCO4984', 'SCO5042', 'SCO5150', 'SCO5520', 'SCO5655', 'SCO5676', 'SCO5774', 'SCO5775', 'SCO5776', 'SCO5777', 'SCO5778', 'SCO5779', 'SCO6060', 'SCO6222', 'SCO6412', 'SCO6702', 'SCO6732', 'SCO6789', 'SCO6962', 'SCO7034', 'SCO7035', 'SCO7036', 'SCO7049', 'SCO7152']

TARGET_GLUTAMATE_GENES = ["SCO4159","SCO2213","SCO2198","SCO2210","SCO5584","SCO5585","SCO2234","SCO5774","SCO5775",
                          "SCO5776","SCO5777","SCO5778","SCO5779","SCO4683","SCO2999","SCO1977","SCO2025","SCO2026",
                          "SCO3416","SCO1613","SCO2241","SCO6962","SCO4078","SCO7049","SCO5520","SCO6412","SCO5676",
                          "SCO7034","SCO1018"]

def convert_sample_number_to_timepoint(row):
    if row["Strain"] == "M145":
        timepoint_dict = {"P6":1, "P14":2, "P18":3, "P22":4, "P26":5, "P30":6, "P34":7, "P38":8, "P42":9}
    else:
        timepoint_dict = {"P18":1, "P26":2, "P30":3, "P34":4, "P38":5, "P42":6, "P46":7, "P50":8}
    return timepoint_dict[row["Sample number"]]

def convert_to_timepoint2(row):
    if row["Strain"] == "M145":
        timepoint_dict = {21:1, 29:2, 33:3, 37:4, 41:5, 45:6, 49:7, 53:8, 57:9}
    else:
        timepoint_dict = {33:1, 41:2, 45:3, 49:4, 53:5, 57:6, 61:7, 65:8, 69:9}
    return timepoint_dict[int(row["Hours"])]

def convert_to_timepoint(strain, hour):
    if strain == "M145":
        timepoint_dict = {21:1, 29:2, 33:3, 37:4, 41:5, 45:6, 49:7, 53:8, 57:9}
    else:
        timepoint_dict = {33:1, 41:2, 45:3, 49:4, 53:5, 57:6, 61:7, 65:8, 69:9}
    return timepoint_dict[int(hour)]

def get_all_flux(folder, reactions):
    """
    DEPRECATED
    """
    folder = Path(folder)
    df_dict = {}
    for i, file in enumerate(folder.glob("*.csv")):
        _, strain, hours = file.stem.split("_")
        print(strain, hours)
        df = pd.read_csv(folder / file, sep = "\t", header = None, index_col = 0)
        df = df.loc[df.index.str.contains("|".join(reactions)), :]
        df = revert_GECKO_reactions(df)
        # print(df.head())
        df_dict[(strain, hours)] = df.loc[reactions, :]
        
    for r in reactions:
        df_list = []
        for (strain, hours), df in df_dict.items():
            new_df = pd.DataFrame()
            new_df[r] = df.T[r]
            new_df["Strain"] = strain
            new_df["Hours"]  = hours
            new_df["Time point"] = convert_to_timepoint(strain, hours)
            df_list.append(new_df)
        r_df = pd.concat(df_list, axis = 0)
        fn = folder / "reactions"/ "{0}.csv".format(r)
        r_df.to_csv(fn, index = False)

def plot_all_flux(folder, reaction):
    df = pd.read_csv(Path(folder) / "reactions" / "{0}.csv".format(reaction))
    print(df)
    sns.violinplot(x = "Time point", y = reaction, hue = "Strain", data = df, inner = "stick", split = False, cut = 0, scale = "count")
    # sns.stripplot(x = "Time point", y = reaction, hue = "Strain", data = df, inner = None, split = False)
    # sns.swarmplot(x = "Time point", y = reaction, hue = "Strain", data = df)
    # sns.boxplot(x = "Time point", y = reaction, hue = "Strain", data = df)

    plt.show()

def get_gene_data(gene_expression_csv, gene_ids):
    # The csv file contains normalized rna-seq data for each of the three biological replicates for each strain
    df_gene_fermentors = pd.read_csv(gene_expression_csv, sep = "\t")
    df_gene_fermentors.set_index("Identifier", inplace = True)
    #print(df_gene_fermentors.T)

    df_gene = df_gene_fermentors.T

    # Set strain based on index
    df_gene.loc[df_gene.index.str.contains("|".join(M145_FERMENTORS)), "Strain"] = "M145"    
    df_gene.loc[df_gene.index.str.contains("|".join(M1152_FERMENTORS)), "Strain"] = "M1152"    
    
    # print(df_gene.head())
    # Set timepoint
    df_gene.loc[:, "Hours"] = df_gene.index.str[5:7]
    
    df_selected_genes = df_gene.loc[:, gene_ids+["Strain", "Hours"]]
    df_selected_genes["Time point"] = df_selected_genes.apply(convert_to_timepoint2, axis = 1)

    return df_selected_genes

def get_single_flux(mean_flux_fn, sep = "\t"):
    df_flux = pd.read_csv(mean_flux_fn, index_col = 0, sep = sep)
    print(df_flux)
    #Remove the gene annotation and name
    df_flux = df_flux.iloc[:, 2:].T
    # df_flux.loc[df_flux.index.str.contains("M145"), "Strain"] = "M145"
    # df_flux.loc[df_flux.index.str.contains("M1152"), "Strain"] = "M1152"
    # df_flux[["Strain", "Hours"]] = 
    df_flux.index = df_flux.index.str.split("_", expand = True)
    df_flux.index.set_names(["Strain", "Hours"], inplace = True)
    df_flux.reset_index(inplace = True)
    return df_flux 


def plot_gene_and_single_flux(mean_flux_fn, gene_expression_csv):
    df_flux = get_single_flux(mean_flux_fn)
    df_gene = get_gene_data(gene_expression_csv)
    
def plot_reactions(mean_flux_fn, reaction_ids):
    df_flux = get_single_flux(mean_flux_fn)
    # df_flux = df_flux.loc[:, reaction_ids]
    #
    print(df_flux.head())
    for r in reaction_ids:
        sns.lineplot(x = "Hours", y = r, hue = "Strain", data = df_flux)
        plt.show()



def revert_GECKO_reactions(full_flux_frame):
    """
    Sum all the different parallell reactions created by the GECKO-method
    - irreversible to reversible (REV)
    - arm_reactions?
    - _Nox reactions
    """
    all_reaction_ids = list(full_flux_frame.index)
    all_columns = list(full_flux_frame.columns)
    full_flux_frame.loc[:, "keep"] = False
    full_flux_frame.loc[:, "base id"] = list(full_flux_frame.index)
    # get arm reactions
    arm_ids = []
    non_arm_ids  = []

    for r_id in all_reaction_ids:
        if r_id[:4] == "arm_":
            full_flux_frame.loc[r_id, "keep"] = True
            if r_id[-4:] == "_REV":
                full_flux_frame.loc[r_id, "base id"] = r_id[4:-4]
                full_flux_frame.loc[r_id, all_columns] *= -1
            else:
                full_flux_frame.loc[r_id, "base id"] = r_id[4:]
            arm_ids.append(r_id[4:])
        else:
            non_arm_ids.append(r_id)


    for r_id in non_arm_ids:
        no_stripped_id = re.sub(r"No\d+\Z", "", r_id)
        if not no_stripped_id in arm_ids:
            full_flux_frame.loc[r_id, "keep"] = True
        
            if no_stripped_id[-4:] == "_REV":
                full_flux_frame.loc[r_id, "base id"] = no_stripped_id[:-4]
                full_flux_frame.loc[r_id, all_columns] *= -1
            else:
                full_flux_frame.loc[r_id, "base id"] = no_stripped_id
        
    full_flux_frame.set_index("base id", inplace = True)
    flux_frame = full_flux_frame.loc[full_flux_frame["keep"], all_columns]
    flux_frame = flux_frame.groupby(flux_frame.index, axis = 0).sum()
    # print(flux_frame)
    # print(flux_frame_sum)
    return flux_frame

def plot_genes_transcriptome_proteome(proteome_norm_tsv, gene_expression_csv, gene_ids):
    """
    Plot the timecourse for a list of genes, one gene ni each plot for both
    transcriptome and proteome data
    """
    df_transcr = get_gene_data(gene_expression_csv, gene_ids)
    df_prot = get_proteome(proteome_norm_tsv, gene_ids)


    for i, gene in enumerate(gene_ids):
        # ax = axes[i]
        fig, ax = plt.subplots(1, figsize  = (8, 6))
        ax = sns.lineplot(x = "Time point", y = gene, hue = "Strain", data = df_transcr, legend = False)#, ax = ax)
        ax2 = ax.twinx()
        ax2 = sns.lineplot(x = "Time point", y = gene, hue = "Strain", data = df_prot, err_style = "bars", ax = ax2,
                           legend = False, alpha = 0.8)

        for i in range(len(ax2.lines)):
            # print(ax2.lines[i].__dict__)
            # if ax2.lines[i]._linestyle is not None:
            ax2.lines[i].set_linestyle("--")

        lines = ax.lines + ax2.lines[::2]
        labels = ["M145 transcriptome", "M1152 transcriptome", "M145 proteome", "M1152 proteome"]
        ax.legend(lines, labels, bbox_to_anchor=(1,0), loc="lower right", 
                bbox_transform=fig.transFigure, ncol=4)

        ax.set_ylabel("Normalized transcriptome data")
        ax2.set_ylabel("Normalized proteome data")
        ax.set_title(gene)
        plt.subplots_adjust(right = 0.85, bottom = 0.2)
        plt.show()

def plot_genes(gene_expression_csv, gene_ids):
    df_selected_genes = get_gene_data(gene_expression_csv, gene_ids)
    
    df_melt = df_selected_genes.melt(["Time point", "Strain"], value_vars = gene_ids, var_name = "Gene", value_name = "Log2 normalized count")
    #print(df_melt.head())
    sns.lineplot(x = "Time point", y = "Log2 normalized count", hue = "Gene", style = "Strain", data = df_melt, err_style = "bars")
    plt.show()

def plot_genes_clustermap(gene_expression_csv, gene_ids, z_score = 0, change_min = 1, max_min = 9):
    df_selected_genes = get_gene_data(gene_expression_csv, gene_ids)
    #print(df_selected_genes.head())
    df_mean = df_selected_genes.groupby(["Strain", "Hours"]).mean()
    df_mean.reset_index(inplace  = True)
    df_mean.index = df_mean["Strain"] + "-" + df_mean["Hours"].map(str)+"h"
    #print(df_mean.index)
    print(df_mean.index)
    df_mean.sort_values(by = ["Strain", "Hours"], ascending = [False, True], inplace = True)
    print(df_mean.index)

    df_change = df_mean.loc[:, gene_ids].max(axis = 0) - df_mean.loc[:, gene_ids].min(axis = 0)
    
    df_idx = (df_change > change_min) & (df_mean.loc[:, gene_ids].max(axis = 0) >max_min)

    df_gene = df_mean.loc[:, gene_ids]
    df_selected = df_gene.loc[:, df_idx]
    # df_selected.sort_index(inplace = True, ascending = False)

    g = sns.clustermap(df_selected.T, col_cluster = False, cmap  = "Reds", 
                        z_score = z_score, figsize = (8, 10), cbar_kws={'label': 'log2 RNA-seq'})
    ax = g.ax_heatmap
    ax.set_xlabel("")
    ax.set_ylabel("")
    #plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
    # top=0.92,bottom=0.255,left=0.14,right=0.755,hspace=0.19,wspace=0.2
    plt.show()

def get_proteome(proteome_norm_tsv, gene_ids):
    df = pd.read_csv(proteome_norm_tsv, sep = "\t", index_col = 0).T
    df_selected_genes = df.loc[:, gene_ids] # Select only a few genes
    df_selected_genes.index = df_selected_genes.index.str.split("_", expand = True) # Split index into strain, f-number and sample time
    df_selected_genes.index.set_names(["Strain", "F number", "Sample number"], inplace = True) # Set index names
    df_selected_genes.reset_index(inplace = True) # Convert multiindex frame to df with three extra colums
    print(df_selected_genes)
    df_selected_genes["Time point"] = df_selected_genes.apply(convert_sample_number_to_timepoint, axis = 1)
    return df_selected_genes

def plot_proteome(proteome_norm_tsv, gene_ids):
    df_prot = get_proteome(proteome_norm_tsv, gene_ids)
    # sns.lineplot(data = df_prot)
    # plt.show()
    df_melt = df_prot.melt(["Time point", "Strain"], value_vars = gene_ids, var_name = "Protein", value_name = "Normalized counts")
    sns.lineplot(x = "Time point", y = "Normalized counts", hue = "Protein", style = "Strain", data = df_melt, err_style = "bars")
    plt.show()

def plot_germicidin():
    germicidin_fn = "../../ComplementaryData/data/germicidin.csv"
    df = pd.read_csv(germicidin_fn, header = 0, sep = ",")
    print(df.head())
    ax = sns.lineplot(x = "Time [h]", y = "Concentration [ng/ml]", hue = "Strain", style = "Compound", data = df, palette = ["b", "r", "g"], err_style = "band", markers=True, ci = 'sd')


    ax.axvline(x = 47, ymin = 0, ymax = 1, c = "r", ls = "--")
    ax.axvline(x = 35, ymin = 0, ymax = 1, c = "b", ls = "--")
    ax.axvline(x = 38, ymin = 0, ymax = 1, c = "g", ls = "--")
    plt.show()

def plot_histogram_carbon_nitrogen_uptake_and_secretion(mean_flux_fn, cols = ["M145_29", "M1152_41"], 
                        in_reactions = ["EX_glc__D_e", "EX_glu__L_e"], out_reactions = ["EX_nh4_e", "EX_ac_e"]):
    df_flux = pd.read_csv(mean_flux_fn, index_col = 0, sep = "\t")
    
    reactions = in_reactions + out_reactions
    df_mean = df_flux.loc[reactions, ["MEAN_{0}".format(x) for x in cols]]
    df_mean.loc[in_reactions, :] *= -1

    df_mean.reset_index(inplace = True)
    df_mean.columns = ["Reaction ID"] + cols
    # df_mean.plot(kind = "bar", use_index = True)
    # plt.show()
    
    df_melt = df_mean.melt(["Reaction ID"], value_vars = cols, var_name = "Strain", value_name = "CO2 normalized flux")
    sns.barplot(data = df_melt, x = "Reaction ID", y = "CO2 normalized flux", hue = "Strain")
    plt.show()


def read_differentially_expressed_genes(fn):
    df = pd.read_csv(fn, sep = ",", header = 0)
    return list(df["Gene (SCO-number)"])

if __name__ == '__main__':
    mean_flux_co2_scaled = "C:/Users/snorres/Google Drive/scoGEM community model/Supporting information/Model/randomsampling_july/ec-RandSampComb_proteomics_CO2norm.tsv"
    all_random_samples_folder = "C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Eduard random sampling"
    gene_expression_csv = "C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/manuscript/m145-m1152-normalized_allfermenters.tsv"
    proteome_norm_tsv = "C:/Users/snorres/Google Drive/scoGEM community model/Supporting information/Correlation_ProtRNA_Ed/proteomeNorm.tsv"
    differentially_expressed_genes_csv = "C:/Users/snorres/git/Sco-GEM/ComplementaryData/data/differentially_expressed_genes.csv"


    # plot_gene_and_single_flux(mean_flux_co2_scaled, gene_expression_csv)
    # plot_reactions(mean_flux_co2_scaled, reaction_ids = ["EX_glc__D_e", "EX_glu__L_e"])

    if 0:
        # Plot no 1
        get_all_flux(all_random_samples_folder, DIFFERENT_M1152+DELAYED_M1152)
        for r in DIFFERENT_M1152+DELAYED_M1152:
            plot_all_flux(all_random_samples_folder, r)

    if 1:
        plot_germicidin()
    if 0:
        # plot no 2
        gene_ids = ["SCO6984", "SCO1527"]
        # plot_genes(gene_expression_csv, gene_ids)
        # plot_genes(gene_expression_csv, ALL_ATP_SYNTHASE_GENES)
        plot_genes(gene_expression_csv, SPODM)
        plot_genes(gene_expression_csv, CAT)


    if 0:
        #  Plot proteome (no 3)
        gene_ids = ["SCO1196", "SCO1968", "SCO2286"]
        plot_proteome(proteome_norm_tsv, gene_ids)


    if 0:
        #  Plot proteome (no 4)
        gene_ids = ["SCO1565", "SCO4229"]
        # plot_proteome(proteome_norm_tsv, gene_ids)
        # plot_proteome(proteome_norm_tsv, ALL_ATP_SYNTHASE_GENES)
        plot_proteome(proteome_norm_tsv, SPODM)
        plot_proteome(proteome_norm_tsv, CAT)

    if 0:
        # Plot proteome and transcriptome in same plot
        
        # plot_delayed_genes_transcriptome_proteome(proteome_norm_tsv, gene_expression_csv, CYO2a_GENES)
        plot_genes_transcriptome_proteome(proteome_norm_tsv, gene_expression_csv, ATP_SYNTHASE_GENES)
        # plot_delayed_genes_transcriptome_proteome(proteome_norm_tsv, gene_expression_csv, DELAYED_GENES)
        # plot_genes_transcriptome_proteome(proteome_norm_tsv, gene_expression_csv, INOSITOL_DEHYDROGENASE_GENES)
        # plot_genes_transcriptome_proteome(proteome_norm_tsv, gene_expression_csv, ["SCO1565", "SCO4229"])
    if 0:
        get_all_flux(all_random_samples_folder, ["CS", "ACONTa"])
        plot_all_flux(all_random_samples_folder, "CS")

    if 0:
        plot_histogram_carbon_nitrogen_uptake_and_secretion(mean_flux_co2_scaled)

    if 0:
        # Plot gene expression and proteome abundance of genes related to oxidative stress
        genes = list(set(STRESS_GENES_2 + SPODM + CAT))
        genes.sort()
        
        differentially_expressed_genes = read_differentially_expressed_genes(differentially_expressed_genes_csv)
        diff_genes = []
        for g in genes:
            if not g in differentially_expressed_genes:
                print("{0} is not differentially expressed".format(g))
            else:
                diff_genes.append(g)



        plot_genes_clustermap(gene_expression_csv, diff_genes)

    if 0:
        # Plot gene expression  of genes related to glutamate
        genes = list(set(TARGET_GLUTAMATE_GENES))
        genes.sort()
        print(genes)
        differentially_expressed_genes = read_differentially_expressed_genes(differentially_expressed_genes_csv)
        diff_genes = []
        i = 0
        for g in genes:
            if not g in differentially_expressed_genes:
                print("{0} is not differentially expressed".format(g))
                i+=1
            else:
                diff_genes.append(g)
        print(len(genes), i)


        plot_genes_clustermap(gene_expression_csv, genes, z_score = 0, change_min = 0, max_min = 0)
