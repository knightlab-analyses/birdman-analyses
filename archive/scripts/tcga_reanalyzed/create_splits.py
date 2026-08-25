#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga_reanalyzed/%x.%a.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --array=1-20

import biom
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split, StratifiedKFold


def df_to_biom(df):
    tbl = biom.Table(
        df.T.values,
        sample_ids=df.index,
        observation_ids=df.columns
    )
    return tbl


def main():
    PROC_PATH = Path("data/tcga_reanalyzed/processed")
    TBL_FILE = PROC_PATH / "merged_tbl.biom"
    MD_FILE = PROC_PATH / "processed_md.tsv"

    TBL = biom.load_table(TBL_FILE)
    MD = pd.read_table(MD_FILE, sep="\t", index_col=0)

    y = MD["investigation"]
    X = TBL.to_dataframe(dense=True).T.loc[y.index]

    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.8)

    test_tbl = df_to_biom(X_test)
    test_tbl_path = str(PROC_PATH / "merged_tbl.test.biom")
    with biom.util.biom_open(test_tbl_path, "w") as f:
        test_tbl.to_hdf5(f, "split")


    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    for i, (train_index, valid_index) in enumerate(skf.split(X_train, y_train)):
        train_tbl = X.iloc[train_index]
        valid_tbl = X.iloc[valid_index]

        train_tbl = biom.Table(train_tbl.T.values, sample_ids=train_tbl.index,
                               observation_ids=train_tbl.columns)
        valid_tbl = biom.Table(valid_tbl.T.values, sample_ids=valid_tbl.index,
                               observation_ids=valid_tbl.columns)

        split_dir = PROC_PATH / f"split_{i}"
        split_dir.mkdir(exist_ok=True)

        train_tbl_path = str(split_dir / f"merged_tbl.train.{i}.biom")
        valid_tbl_path = str(split_dir / f"merged_tbl.valid.{i}.biom")

        with biom.util.biom_open(train_tbl_path, "w") as f:
            train_tbl.to_hdf5(f, "split")

        with biom.util.biom_open(valid_tbl_path, "w") as f:
            valid_tbl.to_hdf5(f, "split")


if __name__ == "__main__":
    main()
