import h5py
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import gzip
import csv
import numpy as np

# Load the NPZ file
data = np.load('/home2/npp8/data/1_subsample_run0/final/concat_procap_0.npz')

# Assuming there's a key in the NPZ file that you want the length of
# You need to know the key (variable name) within the .npz file
key = list(data.keys())[0]  # This gets the first key, adjust based on your needs
data_length = len(data[key])
print(f"Length of '{key}' dataset in 'concat_procap_0.npz': {data_length}")

import h5py

# Open the HDF5 file
with h5py.File('n1_run0_prediction.hdf5', 'r') as hf:
    # Assuming 'track' or 'quantity' is the dataset you want to get the length of
    data_length = len(hf['track'])  # or use 'quantity' depending on which dataset you're interested in
    print(f"Length of 'track' dataset in 'n1_run0_prediction.hdf5': {data_length}")


# with gzip.open('/fs/cbsubscb17/storage/projects/CLIPNET/data/gse110638/fixed_windows/data_folds/procap/concat_procap_0.csv.gz', mode='rt') as file:
#     reader = csv.reader(file)
#     for i, line in enumerate(reader):
#         if i < 5:
#             print(line)
#         else:
#             break
# out_dir = Path("/home2/ayh8/predictions/lcl_subsample/") # wherever your predictions are stored
# procap_fp = "/home2/ayh8/data/lcl/fixed_windows/concat_procap_0.csv.gz" # wherever your experimental data is stored

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