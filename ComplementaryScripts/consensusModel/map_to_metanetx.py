import cobra
import pandas as pd
from pathlib import Path

def map_model_metabolites(model, metanetx_fn):
    df = pd.read_csv(metanetx_fn, header = None, sep = "\t", comment = "#")
    df.columns = ["db:id", "metanetx", "reason", "name"]

    kegg_df = df[df["db:id"].str.contains("kegg:")]
    bigg_df = df[df["db:id"].str.contains("bigg:")]
    metacyc_df = df[df["db:id"].str.contains("metacyc:")]
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
    new_df.to_csv("metanetx_to_change.csv", index_label = "index")
    # print(new_df)

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
    metanetx_fn = repo_path / "ComplementaryData" / "curation" / "metanetx_chem_xref.tsv"
    model_fn = repo_path / "ModelFiles" / "xml" / "scoGEM.xml"
    model = cobra.io.read_sbml_model(str(model_fn))
    # map_model_metabolites(model, metanetx_fn)
    fn = repo_path / "ComplementaryData" / "curation" /"metanetx_to_change.csv"
    apply_metanetx_mapping(model, fn)