import glob
import pandas as pd

all_folders = glob.glob('Data\*')
all_datasets = {}

# Iterates through all folders and combines the data in the files
for folder in all_folders:
    folder_filenames = glob.glob(folder + '\*.dat')
    times = []
    dataframes = []
    for filename in folder_filenames:
        # Skipped 2 rows because the file has the file name and date in those rows and not data
        time = float(
            pd.read_csv(filename, nrows=2, names=["metadata"])["metadata"][1][4:])
        times.append(time)
        df = pd.read_csv(filename, skiprows=2, sep="\t", names=["energy", "count"])
        dataframes.append(df)
        all_datasets[folder[5:]] = {"elapsed time": sum(times),
                                    "data": pd.concat(dataframes).groupby("energy")["count"].sum().reset_index()}

# Prints out the data in the datasets
print("All Datasets:")
for folder, dataset in all_datasets.items():
    print("Folder:", folder)
    print(dataset)



