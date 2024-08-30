import h5py
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

# out_dir = Path("/home2/ayh8/predictions/lcl_subsample/") # wherever your predictions are stored
# procap_fp = "/home2/ayh8/data/lcl/fixed_windows/concat_procap_0.csv.gz" # wherever your experimental data is stored
import gzip
import csv

# Replace 'your_file.csv.gz' with the path to your .csv.gz file
with gzip.open('/fs/cbsubscb17/storage/projects/CLIPNET/data/gse110638/fixed_windows/data_folds/procap/concat_procap_0.csv.gz', mode='rt') as file:
    reader = csv.reader(file)
    for i, line in enumerate(reader):
        if i < 5:
            print(line)
        else:
            break

# procap = pd.read_csv(procap_fp, index_col=0, header=None)

# runs = range(5)
# for r in runs:
#     names = tuple(filter(None, Path(f"n1_run{r}/lines.txt").read_text().split("\n")))
#     out_fp = out_dir.joinpath(f"n1_run{r}_fold_0_predictions.h5")
#     with h5py.File(out_fp, "r") as hf:
#         profile = hf["track"][:]
#         quantity = hf["quantity"][:]
#         prediction = (profile / np.sum(profile, axis=1)[:, None]) * quantity
#     # filter rows by name
#     mask = procap.index.str.startswith(names)
#     procap_sub = np.array(procap)[mask][:, np.r_[250:750, 1250:1750]]
#     prediction_sub = prediction[mask]
#     # calculate performance metrics
#     profile_pearson = pd.DataFrame(procap_sub).corrwith(pd.DataFrame(prediction_sub), axis=1)
#     quantity_pcc = pearsonr(np.log(procap_sub.sum(axis=1) + 1e-3), np.log(prediction_sub.sum(axis=1) + 1e-3))[0]
#     print(f"Run {r}: {profile_pearson.median()}, {quantity_pcc}")