import os
import shutil
import pandas as pd


def move_samples_to_one_folder(old_dir: str, new_dir):
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    entries = os.listdir(old_dir)
    # Filter out and return only the directories
    folders = [entry for entry in entries if os.path.isdir(os.path.join(old_dir, entry))]
    phenotype_folders = [entry for entry in folders if (entry.startswith("COSO") and not entry.endswith("results"))]
    samples = []
    for phenotype in phenotype_folders:
        phenotype_folder = os.path.join(old_dir, phenotype)
        files = os.listdir(phenotype_folder)
        samples += [(phenotype_folder, file) for file in files if file.startswith("COSS")]

    for folder, file in samples:
        src_file = os.path.join(folder, file)
        dst_file = os.path.join(new_dir, file)
        shutil.copy2(src_file, dst_file)


def add_sample_name_as_column(output_dir):
    entries = os.listdir(output_dir)
    for entry in entries:
        file = os.path.join(output_dir, entry)
        df = pd.read_csv(file, sep="\t")
        df["FILTER"] = entry[:-4:]
        df.to_csv(file, index=False, sep="\t")
