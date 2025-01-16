import os
import sys
import pandas as pd

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from default import DATAPATH, ST_TYPES

pathogens = [name for name in os.listdir(DATAPATH) if os.path.isdir(os.path.join(DATAPATH, name))]


params = ["Total all molecules",
          "Positive percentage All lc",
          "Positive percentage All hc",
          "Total Org molecules",
          "Positive percentage Org lc",
          "Positive percentage Org hc",
          "Total Prot molecules",
          "Positive percentage Prot hc", 
          "Positive percentage Prot lc", 
          "Total assays collected",
          "Positive percentage assays"]      

summary_all = []

for p in pathogens:
    print(p)
    s = pd.read_csv(os.path.join(DATAPATH, p, f'{p}_summary.csv'))
    sp = {}
    total_mols_all = s.loc[s['Concept'] == 'Dataset Original Unique molecules:', 'Cases'].values[0]
    pos_lc_all = s.loc[s['Concept'] == 'Dataset Original, low cutoff positive molecules:', 'Cases'].values[0]
    pos_hc_all = s.loc[s['Concept'] == 'Dataset Original, high cutoff positive molecules:', 'Cases'].values[0]
    sp["Total All molecules"] = total_mols_all
    sp["Positive percentage All lc"] = pos_lc_all/total_mols_all*100
    sp["Positive percentage All hc"] = pos_hc_all/total_mols_all*100
    total_mols_org = s.loc[s['Concept'] == 'Dataset Org All, low cutoff total data size', 'Cases'].values[0]
    sp["Total Org Molecules"] = total_mols_org
    try:
        pos_lc_org = s.loc[s['Concept'] == 'Dataset Org All, low cutoff positive molecules:', 'Cases'].values[0]
        sp["Total Org Molecules"] = total_mols_org
        sp["Positive percentage Org lc"] = pos_lc_org/total_mols_org*100
    except:
        sp["Positive percentage Org lc"] = 0
    try:
        pos_hc_org = s.loc[s['Concept'] == 'Dataset Org All, high cutoff positive molecules:', 'Cases'].values[0]
        sp["Positive percentage Org hc"] = pos_hc_org/total_mols_org*100
    except:
        sp["Positive percentage Org hc"] = 0
    total_mols_prot = s.loc[s['Concept'] == 'Dataset Prot All, low cutoff total data size', 'Cases'].values[0]
    sp["Total Prot Molecules"] = total_mols_prot
    try:
        pos_lc_prot = s.loc[s['Concept'] == 'Dataset Prot All, low cutoff positive molecules:', 'Cases'].values[0]
        sp["Positive percentage Prot lc"] = pos_lc_prot/total_mols_prot*100
    except:
        sp["Positive percentage Prot lc"] = 0
    try:
        pos_hc_prot = s.loc[s['Concept'] == 'Dataset Prot All, high cutoff positive molecules:', 'Cases'].values[0]
        sp["Positive percentage Prot hc"] = pos_hc_prot/total_mols_prot*100
    except:
        sp["Positive percentage Prot hc"] = 0
    assay_count = 0
    if sp["Positive percentage Org lc"] or sp["Positive percentage Org hc"]!= 0:
        assay_count = assay_count+1
    if sp["Positive percentage Prot lc"] or sp["Positive percentage Prot hc"]!= 0:
        assay_count = assay_count+1
    for x in [0,1,2]:
        i=0
        for index, row in s.iterrows():
            if (f"Dataset Org Top_{x}," in row['Concept'] and "positive molecules:" in row['Concept']):
                i = 1
                break
        assay_count = assay_count+i
    for x in [0,1,2]:
        i=0
        for index, row in s.iterrows():
            if (f"Dataset Prot Top_{x}," in row['Concept'] and "positive molecules:" in row['Concept']):
                i = 1
                break
        assay_count = assay_count+i
    for x in ["MIC", 'IZ', "IC50", "Inhibition", "Activity"]:
        i=0
        for index, row in s.iterrows():
            if (f"Dataset Bioassay {x} merged," in row['Concept'] and "positive molecules:" in row['Concept']):
                i = 1
                break
        assay_count = assay_count+i
    
    sp["Total assays collected"] = assay_count

    for x in ["MIC", 'IZ', "IC50", "Inhibition", "Activity"]:
        for index, row in s.iterrows():
            if (f"Dataset Bioassay {x} merged, low cutoff" in row['Concept'] and "positive molecules:" in row['Concept']):
                total_bioa_lc = s.loc[s['Concept'] == f'Dataset Bioassay {x} merged, low cutoff total data size', 'Cases'].values[0]
                pos_bioa_lc = s.loc[s['Concept'] == f'Dataset Bioassay {x} merged, low cutoff positive molecules:', 'Cases'].values[0]
                sp[f"Positive percentage Org {x} lc"] = pos_bioa_lc / total_bioa_lc*100
            if (f"Dataset Bioassay {x} merged, high cutoff" in row['Concept'] and "positive molecules:" in row['Concept']):
                total_bioa_lc = s.loc[s['Concept'] == f'Dataset Bioassay {x} merged, high cutoff total data size', 'Cases'].values[0]
                pos_bioa_lc = s.loc[s['Concept'] == f'Dataset Bioassay {x} merged, high cutoff positive molecules:', 'Cases'].values[0]
                sp[f"Positive percentage Org {x} hc"] = pos_bioa_lc / total_bioa_lc*100

    summary_all += [sp]

# Create DataFrame from the list of dictionaries
df = pd.DataFrame(summary_all, index=pathogens)

# Transpose the DataFrame to have pathogens as columns
df = df.transpose()
df = df.reset_index()
df = df.rename(columns={'index': 'Concept'})

df.to_csv(os.path.join(DATAPATH, "summary_all.csv"), index=False)