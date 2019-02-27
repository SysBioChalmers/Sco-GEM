import pandas as pd
from pathlib import Path
import scipy.stats as st
import numpy as np
from matplotlib import pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import cobra
import seaborn as sns
import re


REPO_DIR = Path(__file__).parent.parent.parent
RESULT_FOLDER = REPO_DIR / "ComplementaryData" / "ecmodel" / "random_sampling" 

def get_df(strain = None, timepoint = None, iterations = None, processes = None, filename = None):
    if not filename:
        filename = "randomsampling_{0}_{1}h_{2}_iterations_{3}_processes.csv".format(strain, timepoint, iterations, processes)
    df = pd.read_csv(str(filename), index_col = 0)
    df_agg = df.agg(["mean", "std"]).T
    idx = (df_agg["mean"] == 0) & (df_agg["std"] == 0)
    df_agg_nz = df_agg.loc[~idx]
    df_nz = df.T[~idx]
    n = len(df.index)

    return df_agg, df_agg_nz, df_nz, n

def revert_GECKO_reactions(full_flux_frame, key = "fluxes"):
    """
    Sum all the different parallell reactions created by the GECKO-method
    - irreversible to reversible (REV)
    - arm_reactions?
    - _Nox reactions
    """
    all_reaction_ids = list(full_flux_frame.index)
    full_flux_frame["keep"] = False
    full_flux_frame["base id"] = list(full_flux_frame.index)

    # get arm reactions
    arm_ids = []
    non_arm_ids  = []

    for r_id in all_reaction_ids:
        if r_id[:4] == "arm_":
            full_flux_frame.loc[r_id, "keep"] = True
            if r_id[-4:] == "_REV":
                full_flux_frame.loc[r_id, "base id"] = r_id[4:-4]
                full_flux_frame.loc[r_id, key] *= -1
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
                full_flux_frame.loc[r_id, key] *= -1
            else:
                full_flux_frame.loc[r_id, "base id"] = no_stripped_id
        

    flux_frame = full_flux_frame.loc[full_flux_frame["keep"], [key, "base id"]]
    # print(flux_frame)

    flux_frame_sum = flux_frame.groupby("base id").sum()
    # print(flux_frame_sum)
    return flux_frame_sum

def compare_two(fn1, fn2):
    df1_agg, df1_agg_nz, df1_nz, n1 = get_df(filename = fn1)
    df2_agg, df2_agg_nz, df2_nz, n2 = get_df(filename = fn2)
    
    for r_id in df2_agg.index:
        try:
            m1, s1 = df1_agg.loc[r_id, :]
            m2, s2 = df2_agg.loc[r_id, :]
        except KeyError:
            print(r_id, "is not in both datasets")
            continue
        try:
            t, dof, p = welch_ttest(m1, s1, m2, s2, n1, n2)
        except ZeroDivisionError:
            continue
        print(r_id, p, m1-m2)

def merge_means(fn_list):
    df_list = []
    for fn in fn_list:
        df_agg, df_agg_nz, _, _ = get_df(filename = fn)
        df_list.append(df_agg["mean"])

    df = df_list.pop(0).to_frame()
    # df.sort_values(by = "mean", inplace = True)
    print(df.head())


    for i, df_i in enumerate(df_list):
        print(i)
        print(df_i["NNDPRNo1"])
        df = df.join(df_i.to_frame(), rsuffix = "_{0}".format(i+1), how = "inner")
    print(df.head())
    df.to_csv(str(RESULT_FOLDER / "mean_df.csv"))
    return df


def welch_ttest(m1, s1, m2, s2, n1, n2 = None):
    if not n2:
        n2 = n1
    v1 = s1**2
    v2 = s2**2
    t = (m1 - m2) / np.sqrt(v1 / n1 + v2 / n2)
    dof = (v1 / n1 + v2 / n2)**2 / (v1**2 / (n1**2 * (n1 - 1)) + v2**2 / (n2**2 * (n2 - 1)))
    p = 2 * st.t.cdf(-abs(t), dof)
    return t, dof, p


def plot_means(fn_list):
    # df = merge_means(fn_list)
    # df = df.reset_index()
    df = pd.read_csv(str(RESULT_FOLDER / "mean_df.csv"), index_col = 0)
    # df = df.iloc[:100, :]
    df.sort_values(by = "mean_1", inplace = True)
    df[df<0] = 0
    draw_prot_idx = df.index.str.contains("draw_prot_")
    df_prot = df[draw_prot_idx]
    df = df[~draw_prot_idx]
    print(df_prot)
    df_prot = df_prot.reset_index()
    df = df.reset_index()

    df_prot.plot(lw = 0, marker = "o", alpha = 0.7, markeredgecolor = "k")

    print(df)
    # df.sort_values(by = "mean_1", inplace = True, ascending = False, kind = "mergesort")
    # df.reset_index(inplace = True)
    df.plot(lw = 0, marker = "o", alpha = 0.7, markeredgecolor = "k")
    # df.plot(y = "mean_1", logy = True, lw = 0, marker = "o", alpha = 0.7, markeredgecolor = "k")
    # plt.show()

def aggregate_pca(fn_list, labels, key = "mean", name = "temp", load = False):
    """
    Calculate and displaye a PCA-plot for the given strains / timepoints for the mean/median/min/max values of the random samples
    """
    df = aggregate_random_samples(fn_list, labels, key, name, load)
    
    features = list(df.columns.values)
    # features.remove("label")
    # features.remove("strain no")
    print(len(features))
    # protein reactions
    draw_prot_features = [x for x in features if x[:10] == "draw_prot_"]
    prot_features = [x for x in features if x[:5] == "prot_"]
    all_prot_features = draw_prot_features + prot_features
    non_prot_features = [x for x in features if not x in all_prot_features]
    

    draw_prot_x = df.loc[:, draw_prot_features].values
    prot_x = df.loc[:, prot_features].values
    all_prot_x = df.loc[:, all_prot_features].values
    non_prot_x = df.loc[:, non_prot_features].values
    y = df.loc[:, "label"].values


    pca(draw_prot_x, y, labels, "Draw protein reactions", draw_prot_features)
    pca(prot_x, y, labels, "Protein exchange reactions", prot_features)
    pca(all_prot_x, y, labels, "All protein reactions", all_prot_features)
    pca(non_prot_x, y, labels, "All non-protein reactions", non_prot_features)

def full_pca(fn_list, labels):
    # df = merge_data(fn_list, labels)
    df = pd.read_csv(str(RESULT_FOLDER / "merge_data.csv"), index_col = 0)
    print(df.T.index)
    # print(df)
    # df = revert_GECKO_reactions(df.T, list(df.index))
    print(df.shape)

    # Select strains
    idx = df.loc[:,"label"].isin(labels)
    df = df.loc[idx, :]
    print(df.shape)
    # df[df.isna()] = 0

    features = list(df.columns.values)
    features.remove("label")
    features.remove("strain no")
    print(len(features))
    # protein reactions
    draw_prot_features = [x for x in features if x[:10] == "draw_prot_"]
    prot_features = [x for x in features if x[:5] == "prot_"]
    all_prot_features = draw_prot_features + prot_features
    non_prot_features = [x for x in features if not x in all_prot_features]
    

    draw_prot_x = df.loc[:, draw_prot_features].values
    prot_x = df.loc[:, prot_features].values
    all_prot_x = df.loc[:, all_prot_features].values
    non_prot_x = df.loc[:, non_prot_features].values
    y = df.loc[:, "label"].values


    pca(draw_prot_x, y, labels, "Draw protein reactions", draw_prot_features)
    pca(prot_x, y, labels, "Protein exchange reactions", prot_features)
    pca(all_prot_x, y, labels, "All protein reactions", all_prot_features)
    pca(non_prot_x, y, labels, "All non-protein reactions", non_prot_features)

def pca(x, y, labels, title, features):
    """
    Maxe a PCA plot for the the four first principal components of th data
    """
    print(title, x.shape)
    # Replace nan with 0
    x[np.isnan(x)] = 0

    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components = 4)
    princcomp = pca.fit_transform(x)
    pca_df = pd.DataFrame(data = princcomp, columns = ["PCA 1", "PCA 2", "PCA 3", "PCA 4"])
    pca_df["label"] = y


    # Weights
    weight_df = pd.DataFrame(data = pca.components_, columns = features)
    print(weight_df)
    print(weight_df.shape)


    # pca_df.to_csv(str(RESULT_FOLDER / "pca_df.csv"))
    
    fig, [ax1, ax2] = plt.subplots(1, 2, figsize = (20, 10))
    fig.suptitle(title)
    n = len(labels)
    colors = plt.cm.tab20(np.arange(n))
    for i, label in enumerate(labels):
        keep_idx = pca_df["label"] == label
        ax1.scatter(pca_df.loc[keep_idx, "PCA 1"], pca_df.loc[keep_idx, "PCA 2"], color = colors[i], label = label)
        ax2.scatter(pca_df.loc[keep_idx, "PCA 3"], pca_df.loc[keep_idx, "PCA 4"], color = colors[i], label = label)
    ax1.set_xlabel("PCA 1 ({0:.1f}% explained variance)".format(pca.explained_variance_ratio_[0]*100))
    ax1.set_ylabel("PCA 2 ({0:.1f}% explained variance)".format(pca.explained_variance_ratio_[1]*100))
    ax2.set_xlabel("PCA 3 ({0:.1f}% explained variance)".format(pca.explained_variance_ratio_[2]*100))
    ax2.set_ylabel("PCA 4 ({0:.1f}% explained variance)".format(pca.explained_variance_ratio_[3]*100))
    ax1.legend()
    ax2.legend()
    plt.show()

    fig, axes = plt.subplots(2,2, figsize = (24, 16))
    axes = [ax for sublist in axes for ax in sublist]
    for ax, (index, row) in zip(axes, weight_df.iterrows()):
        largest_n = row.abs().nlargest(100)
        largest_n.plot(kind = "barh", ax = ax)
        ax.set_xlabel("Weight principal component {0}".format(index+1))
    fig.suptitle("{0}: Top 100 weights PCA".format(title))
    plt.show()


def merge_data(fn_list, labels):
    """
    Concatenate data from random samples into one data frame
    """
    df_list = []
    for i, fn in enumerate(fn_list):
        df = pd.read_csv(str(fn), index_col = 0)
        df["strain no"] = i+1
        df["label"] = labels[i]
        df_list.append(df)
    df = pd.concat(df_list, sort = False)
    df.to_csv(str(RESULT_FOLDER / "merge_data.csv"))
    # df[df.isna()] = 0
    return df

def aggregate_random_samples(fn_list, labels, key = "mean", name = "temp", load = False):
    df_list = []
    if load:
        df = pd.read_csv(str(RESULT_FOLDER / "{0}_{1}.csv".format(key, name)), index_col = 0)
    else:
        for i, fn in enumerate(fn_list):
            df = pd.read_csv(str(fn), index_col = 0)
            print(fn, "\n", df.loc[:, ["EX_glc__D_e_REV", "EX_glu__L_e_REV"]])
            
            df[df.isna()] = 0
            agg_df = df.agg([key]).T
            agg_df["strain no"] = i
            
            reverted_df = revert_GECKO_reactions(agg_df, key)[key]
            print(reverted_df.index)
            df_list.append(reverted_df)

        df = pd.concat(df_list, axis  = 1, keys = labels, sort = False)
        df.to_csv(str(RESULT_FOLDER / "{0}_{1}.csv".format(key, name)))
    return df

def QC_of_aggregated_samples(fn_list, labels, key = "mean", name = "temp", load = False):
    df = aggregate_random_samples(fn_list, labels, key, name, load)

    # Find number of large values Large values
    print("Values larger than 10")
    print(df[(df.abs()>10).any(axis=1)])

    # # plot mean values
    # df.plot(kind = "bar")
    # plt.show()

    # Histogram
    df.plot(kind = "hist", bins = 1000, logy = True, stacked = True)
    plt.show()



def analyze_subsystem_max_min_median(fn_list, labels, key = "mean", name = "temp", load = False, row_order = None):
    df = aggregate_random_samples(fn_list, labels, key, name, load)
    df[df.isna()] = 0
    
    # Remove outlier rows
    outlier_row_idx = (df.abs()>10).any(axis=1)
    print("Removed the following rows:", "\n", df[outlier_row_idx])
    df = df[~outlier_row_idx]
    
    print(df.abs().sum())
    carbon_sum = df.loc[["EX_glc__D_e", "EX_glu__L_e"], labels].sum() * -1
    df_abs = add_subsystem_to_df(df.abs()/carbon_sum)
    sub_df = df_abs.groupby("Subsystem").sum()
    print(sub_df)

    if row_order:
        sub_df = sub_df.iloc[row_order, :]
        g = sns.clustermap(sub_df, cmap = "vlag", z_score = 0, col_cluster = False, row_cluster = False)
    else:
        g = sns.clustermap(sub_df, cmap = "vlag", z_score = 0, col_cluster = False)
        print("Row order:", g.dendrogram_row.reordered_ind)

    # sub_df.plot(kind  ="barh")
    plt.show()

def analyze_individual_max_min_median(fn_list, labels, key = "mean", name = "temp", load = False):
    df = aggregate_random_samples(fn_list, labels, key, name, load)
    df[df.isna()] = 0
    # print(df.isna().sum().sum())
    ## Plot clustermap
    # Normalize by dividing by the sum of carbon uptake
    carbon_sum = df.loc[["EX_glc__D_e", "EX_glu__L_e"], labels].sum() * -1
    df.loc[:, labels] = df.loc[:, labels] / carbon_sum
    g = sns.clustermap(df.loc[:, labels], z_score = True, cmap = "vlag", col_cluster = False)
    plt.show()


def add_subsystem_to_df(df):
    """
    df is a dataframe where the indexes are reaction ids and the columns are different timepoints / strains
    the values are the different fluxes
    """
    model_fn = str(REPO_DIR / "ModelFiles"/"xml"/"scoGEM.xml")
    model = cobra.io.read_sbml_model(model_fn)
    df["Subsystem"] = "Missing annotation"
    keep_columns = list(df.columns.values)
    add_rows = []
    df["Keep"] = True
    for r_id, row in df.iterrows():
        try:
            r = model.reactions.get_by_id(r_id)
        except KeyError:
            df.loc[r_id, "Keep"] = False
            print(r_id)
            continue
        try:
            subsystem = r.annotation["kegg.subsystem"]
        except KeyError:
            df.loc[r_id, "Keep"] = False
            print("Missing subsystem: {0}".format(r_id))
            continue

        if isinstance(subsystem, str) and len(subsystem):
            df.loc[r_id, "Subsystem"] = subsystem
        elif isinstance(subsystem, list):
            print(r_id, subsystem)
            df.loc[r_id, "Subsystem"] = subsystem[0]
            for ss in subsystem[1:]:
                print(r_id, subsystem)
                new_row = row.copy()
                new_row["Subsystem"] = ss
                df = df.append(new_row)
    df = df.loc[df["Keep"], :]
    df_to_keep = df.loc[:, keep_columns]
    df_to_keep.loc[:, keep_columns[:-1]] = df_to_keep.loc[:, keep_columns[:-1]].abs()/df_to_keep.loc[:, keep_columns[:-1]].abs().sum()
    return df_to_keep





if __name__ == '__main__':
    fn_list = [
    # RESULT_FOLDER / "randomsample_M145_21h_1000_iterations.csv",
    RESULT_FOLDER / "randomsample_M145_29h_1000_iterations.csv",
    RESULT_FOLDER / "randomsample_M145_33h_1000_iterations.csv",
    RESULT_FOLDER / "randomsample_M145_37h_1000_iterations.csv",
    RESULT_FOLDER / "randomsample_M145_41h_1000_iterations.csv",
    RESULT_FOLDER / "randomsample_M145_45h_1000_iterations.csv",
    RESULT_FOLDER / "randomsample_M145_49h_1000_iterations.csv",
    RESULT_FOLDER / "randomsample_M145_53h_1000_iterations.csv",
    RESULT_FOLDER / "randomsample_M145_57h_1000_iterations.csv",
    # RESULT_FOLDER / "randomsample_M1152_33h_1000_iterations.csv",
    # RESULT_FOLDER / "randomsample_M1152_41h_1000_iterations.csv",
    # RESULT_FOLDER / "randomsample_M1152_45h_1000_iterations.csv",
    # RESULT_FOLDER / "randomsample_M1152_49h_1000_iterations.csv",
    # RESULT_FOLDER / "randomsample_M1152_53h_1000_iterations.csv",
    # RESULT_FOLDER / "randomsample_M1152_57h_1000_iterations.csv",
    # RESULT_FOLDER / "randomsample_M1152_61h_1000_iterations.csv",
    # RESULT_FOLDER / "randomsample_M1152_65h_1000_iterations.csv",
    ]

    labels = ["M145-29", "M145-33", "M145-37", "M145-41", "M145-45", "M145-49", "M145-53", "M145-57"]  #"M145-21", 
    # labels = ["M1152-33", "M1152-41", "M1152-45", "M1152-49", "M1152-53", "M1152-57", "M1152-61", "M1152-65"]
    # labels = ["M1152-41", "M1152-45", "M1152-49", "M1152-53", "M1152-57", "M1152-61", "M1152-65"]
    
    # compare_two(fn, fn2)
    # plot_means([fn1, fn2, fn3, fn4, fn5, fn6])


    # aggregate_pca(fn_list, labels = labels, key = "mean", name = "M1152", load = False)
    # QC_of_aggregated_samples(fn_list, labels = labels, key = "mean", name = "M1152", load = True)
    # M145_row_order = [1, 12, 9, 3, 5, 11, 7, 4, 10, 0, 8, 2, 6, 13]
    # analyze_subsystem_max_min_median(fn_list, labels, key = "mean", name = "M1152_x33", load = True, row_order = M145_row_order)
    analyze_subsystem_max_min_median(fn_list, labels, key = "mean", name = "M145", load = True)

