
import os
import sys
import pandas as pd
import collections
import random
import subprocess
import time
from tqdm import tqdm
from llama_index.core.llms import ChatMessage, MessageRole
from llama_index.llms.llamafile import Llamafile

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)

from default import CONFIGPATH, TMPDIR

LLM_BIN_FILENAME = "gemma-2-9b-it.Q6_K.llamafile"

llm_bin_path = os.path.abspath(os.path.join(root, "..", "bin", LLM_BIN_FILENAME))

print(f"LLM bin path: {llm_bin_path}")

llm_model = Llamafile()


def run_llamafile():
    print("Running LLM...")
    subprocess.Popen(f'bash {llm_bin_path} --server', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("LLM running...")


def classify_activity_comments(file_path, rewrite=False):
    """
    Classify activity comments as active (1), inactive (0) or inconclusive (None).
    """
    df = pd.read_csv(file_path, dtype=str)
    columns = list(df.columns)
    if "activity_classified" not in columns:
        df["activity_classified"] = ""
        df.to_csv(file_path, index=False)
    if rewrite:
        df["activity_classified"] = ""
        df.to_csv(file_path, index=False)
    texts = df["activity_comment"].tolist()
    classifications_ = df["activity_classified"].tolist()
    classifications = []
    for classification in classifications_:
        if str(classification) == "nan":
            classifications.append("")
        else:
            classifications.append(str(classification))
    i = 0
    current_classifications = [x for x in classifications]
    print(len(texts), len(classifications))
    for comment, classification in tqdm(zip(texts, classifications)):
        comment = str(comment)
        if comment == "nan":
            current_classifications[i] = "0"
            i += 1
            continue
        if comment.isnumeric():
            current_classifications[i] = "0"
            i += 1
            continue
        if comment.startswith("Median N="):
            current_classifications[i] = "0"
            i += 1
            continue
        if "no significant biotransformation" in comment.lower():
            current_classifications[i] = "-1"
            i += 1
            continue
        if ". significant biotransformation" in comment.lower():
            current_classifications[i] = "1"
            i += 1
            continue
        if comment.startswith("inhibitor ["):
            current_classifications[i] = "1"
            i += 1
            continue
        if comment.startswith("substrate ["):
            current_classifications[i] = "0"
            i += 1
            continue
        if "putative" in comment.lower() and "is biotransformed by" in comment.lower():
            current_classifications[i] = "1"
            i += 1
            continue
        if "biotransformation" in comment.lower() and "could not be" in comment.lower():
            current_classifications[i] = "0"
            i += 1
            continue
        if "biotransformation" in comment.lower() and "was proven to be" in comment.lower():
            current_classifications[i] = "1"
            i += 1
            continue
        prompt_text = (
            f"You are a biomedical data curator from ChEMBL. Your task is to classify comments about the activity of compounds.\n"
            f"Activity comments should be classified as 1. This includes terms like 'active', 'significant activity', 'inhibition', 'binding', 'toxic', 'agonist', 'antagonist', 'weakly active', 'positive', 'substantial inhibition', 'significant', etc.\n"
            f"Inactivity or no activity comments should be classified as -1. This includes terms like 'inactive', 'no activity', 'not active', 'no inhibition', 'non-toxic', 'negative', 'no significant', etc.\n"
            f"Inconclusive or irrelevant comments should be classified as 0. This includes empty text, terms like 'inconclusive', 'not measured', etc., terms that simply explain the type of measurement, and terms that seem unrelated to the measured activity.\n"
            f"Simply answer 1, 0 or -1. Do not explain anything. Just give the numerical value.\n"
            f"This is the comment you have to curate:"
            f"- {comment}\n"
        )
        if classification != "":
            i += 1
            continue
        try:
            response = llm_model.chat(
                messages=[ChatMessage(role=MessageRole.USER, content=prompt_text)]
            )
        except:
            print("LLM error. Starting LLM...")
            run_llamafile()
            print("LLM started. Waiting for 10 seconds...")
            time.sleep(10)
            response = llm_model.chat(
                messages=[ChatMessage(role=MessageRole.USER, content=prompt_text)]
            )
        prediction = response.message.content.strip()
        prediction = prediction.replace("<|eot_id|>", "")
        prediction = prediction.replace("<end_of_turn>", "")
        prediction = prediction.strip()
        print(f"Prediction: {prediction} / Comment: {comment}")
        if prediction == "1":
            current_classifications[i] = "1"
        elif prediction == "-1":
            current_classifications[i] = "-1"
        elif prediction == "0":
            current_classifications[i] = "0"
        else:
            current_classifications[i] = "To resolve"
        i += 1
        df["activity_classified"] = current_classifications
        df.to_csv(file_path, index=False)
    df["activity_classified"] = current_classifications
    df.to_csv(file_path, index=False)


def classify_activity_standards_with_direction(file_path, rewrite=False):
    """
    Classify activity standards with the right direction (-1: the lower the better, 1: the higher the better, 0: inconclusive).
    """
    df = pd.read_csv(file_path, dtype=str)
    columns = list(df.columns)
    if "activity_direction" not in columns:
        df["activity_direction"] = ""
        df.to_csv(file_path, index=False)
    if rewrite:
        df["activity_direction"] = ""
        df.to_csv(file_path, index=False)
    current_directions_ = df["activity_direction"].tolist()
    current_directions = []
    for direction in current_directions_:
        if str(direction) == "nan":
            current_directions += [""]
        else:
            current_directions += [str(direction)]
    for i, v in enumerate(df.values):
        std_type = v[1]
        definition = v[2]
        std_unit = v[3]
        if current_directions[i] != "":
            continue
        prompt_text = (
            f"You are a biomedical data curator from ChEMBL. Your task is to classify comments types of activity measures of bioactivity assays.\n"
            f"When the lower the value, the higher the activity (including toxicity), classify as -1. For example, in an IC50 measure, smaller concentrations indicate more potency. Likewise, smaller doses, etc.\n"
            f"When the higher the value, the higher the activity (including toxicity), classify as 1. For example, in a percentage inhibition, larger values indicate more potency. Likewise, high zone of inhibition, etc.\n"
            f"When the direction is inconclusive, classify as 0. This includes assays related to pharmacokinetics in the human body (e.g. plasma protein binding, bioavailability, fraction unbound in plasma, etc.), or physicochemical measurements not related to antipathogen activity data.\n"
            f"Simply answer 1, 0 or -1. Do not explain anything. Just give the numerical value.\n"
            f"This is the activity type and units you have to curate:"
            f"- {std_type} : {definition}. Units: {std_unit}\n"
        )
        response = llm_model.chat(
            messages=[ChatMessage(role=MessageRole.USER, content=prompt_text)]
        )
        prediction = response.message.content.strip()
        prediction = prediction.replace("<|eot_id|>", "")
        prediction = prediction.replace("<end_of_turn>", "")
        prediction = prediction.strip()
        print(f"Prediction: {prediction} / {std_type} : {definition}. Units: {std_unit}")
        if prediction == "1":
            current_directions[i] = "1"
        elif prediction == "-1":
            current_directions[i] = "-1"
        elif prediction == "0":
            current_directions[i] = "0"
        else:
            current_directions[i] = "0"
        df["activity_direction"] = current_directions
        df.to_csv(file_path, index=False)
    df["activity_direction"] = current_directions
    df.to_csv(file_path, index=False)


def get_sample_assay_descriptions_for_standard_units():
    print("Getting sample assay descriptions for standard units...")
    print("Getting assay descriptions...")
    df_ass = pd.read_csv(os.path.join(TMPDIR, "assay_descriptions.csv"), low_memory=False)
    assay_descs = {}
    for v in df_ass[["assay_id", "description"]].values:
        assay_descs[v[0]] = v[1]
    print("Getting activities")
    df_act = pd.read_csv(os.path.join(TMPDIR, "activities.csv"), low_memory=False)
    print("Done reading activities")
    std_unit_assays = collections.defaultdict(list)
    for v in tqdm(df_act[["assay_id", "standard_type", "standard_units"]].values):
        std_unit_assays[(v[1], v[2])] += [v[0]]
    R = []
    for k,v in tqdm(std_unit_assays.items()):
        random.shuffle(v)
        v = [assay_descs[x] for x in v[:10]]
        v = [x for x in v if str(v) != "" and str(v) != "nan"]
        v = v[:3]
        if len(v) == 2:
            v = v + [""]
        elif len(v) == 1:
            v = v + ["", ""]
        else:
            pass
        R += [[k[0], k[1], v[0], v[1], v[2]]]
    df = pd.DataFrame(R, columns = ["standard_type", "standard_unit", "assay_description_1", "assay_description_2", "assay_description_3"])
    df.to_csv(os.path.join(CONFIGPATH, "activity_std_units_with_3_assay_descriptions.csv"), index=False)


def classify_all_activity_standards_with_direction(file_path, rewrite=False):
    """
    Classify all activity standard units with the right direction (-1: the lower the more active, 1: the higher the more active, 0: inconclusive).
    """
    df = pd.read_csv(file_path, dtype=str)
    da = pd.read_csv(os.path.join(CONFIGPATH, "activity_std_units_with_3_assay_descriptions.csv"), dtype=str)
    std_unit_descriptions = {}
    for v in da.values:
        std_unit_descriptions[(v[0], v[1])] = [v[2], v[3], v[4]]
    columns = list(df.columns)
    if "activity_direction" not in columns:
        df["activity_direction"] = ""
        df.to_csv(file_path, index=False)
    if rewrite:
        df["activity_direction"] = ""
        df.to_csv(file_path, index=False)
    current_directions_ = df["activity_direction"].tolist()
    current_directions = []
    for direction in current_directions_:
        if str(direction) == "nan":
            current_directions += [""]
        else:
            current_directions += [str(direction)]
    for i, v in tqdm(enumerate(df.values)):
        std_type = v[0]
        std_unit = v[1]
        descriptions = std_unit_descriptions[(std_type, std_unit)]
        if current_directions[i] != "":
            continue
        prompt_text = (
            f"You are a biomedical data curator from ChEMBL. Your task is types of activity measures of bioactivity assays.\n"
            f"You will be provided with an activity type and units, and three assay descriptions where those units are used. You can use these assay descriptions to get hints about the activity type and units. The most important determinant for your decision should be the activity type and unit.\n"
            f"Assign -1 when:\n"
            f"- The lower the value, the higher the activity. That is, smaller values indicate more potency.\n"
            f"- Examples: IC50 (uM, nM, etc.), MIC, MIC90, etc.\n"
            f"Assign 1 when:\n"
            f"- The higher the value, the higher the activity. That is, larger values indicate more potency.\n"
            f"- Examples: % inhibition, activity (percentage or proportion (no unit)), IZ (zone of inhibition, in mm), etc.\n"
            f"Assign 0 when:\n"
            f"- The direction is inconclusive. That is, it is not clear if higher or lower values are preferrable\n"
            f"- Examples: logP, other physicochemical measurements, some pharmacokinetic properties, etc.\n"
            f"Attention: Check the units. Sometimes the same activity type can have different directions depending on the units. For example, -log (minus log, as in pIC50) will reverse the IC50 activity trend from -1 to 1.\n"
            f"Here are some assay descriptions you can use to know more about the activity type unit:\n"
            f"1. {descriptions[0]}\n"
            f"2. {descriptions[1]}\n"
            f"3. {descriptions[2]}\n"
            f"This is the activity type and units you have to classify:"
            f"- Activity type: {std_type}. Unit: {std_unit}\n"
            f"Simply answer 1, 0 or -1. Do not explain anything. Do not give more than one value. Just give the numerical value.\n"
        )
        try:
            response = llm_model.chat(
                messages=[ChatMessage(role=MessageRole.USER, content=prompt_text)]
            )
        except:
            print("LLM error. Starting LLM...")
            run_llamafile()
            print("LLM started. Waiting for 10 seconds...")
            time.sleep(10)
            response = llm_model.chat(
                messages=[ChatMessage(role=MessageRole.USER, content=prompt_text)]
            )
        prediction = response.message.content.strip()
        prediction = prediction.replace("<|eot_id|>", "")
        prediction = prediction.replace("<end_of_turn>", "")
        prediction = prediction.strip()
        print(f"Prediction: {prediction} / {std_type} Units: {std_unit} {descriptions[0]}.")
        if prediction == "1":
            current_directions[i] = "1"
        elif prediction == "-1":
            current_directions[i] = "-1"
        elif prediction == "0":
            current_directions[i] = "0"
        else:
            print("To revise")
            print(prediction)
            current_directions[i] = "To revise"
        df["activity_direction"] = current_directions
        df.to_csv(file_path, index=False)
    df["activity_direction"] = current_directions
    df.to_csv(file_path, index=False)


if __name__ == "__main__":
    print("Starting with activity comments classification...")
    file_path = os.path.join(CONFIGPATH, "activity_comments.csv")
    classify_activity_comments(file_path, rewrite=False)
    print("Activity comments classification done.")
    print("Starting with activity standards classification (direction)...")
    file_path = os.path.join(CONFIGPATH, "activity_stds_lookup.csv")
    classify_activity_standards_with_direction(file_path, rewrite=False)
    print("Activity standards classification done.")
    print("Getting sample assay descriptions for standard units...")
    #get_sample_assay_descriptions_for_standard_units()
    print("Sample assay descriptions for standard units done.")
    print("Starting with all activity standards classification (direction)...")
    file_path = os.path.join(CONFIGPATH, "activity_std_units.csv")
    classify_all_activity_standards_with_direction(file_path, rewrite=False)