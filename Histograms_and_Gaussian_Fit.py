import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from Read_Data import all_datasets

# Gaussian (unbinned) fit based on the Gaussian formula (not applicable to background measurements)
def gaussian_fit(x, a, b, c):
   return a * np.exp(-(x - b) ** 2 / (2 * c ** 2))

# Plots a histogram plot of the datasets with a fitted gaussian curve
def hist_plotter(dataset, left_crop, right_crop):
   energy = list(dataset[1]["data"]["energy"])
   counts = list(dataset[1]["data"]["count"])
   total_count = sum(counts)
   elapsed_time = dataset[1]["elapsed time"]


   # Organizes and crops data
   counted_data = {energy[i]: counts[i] for i in range(len(energy))}
   cropped_data = {energy[i]: counts[i] for i in range(len(energy)) if
                   energy[i] > left_crop and energy[i] < right_crop}
   cropped_energy = list(cropped_data.keys())
   cropped_counts = list(cropped_data.values())

   # Plots histogram
   plt.figure(figsize=(15, 10))
   hits, bin_edges, _ = plt.hist(cropped_data.keys(), alpha=0.5,
                                 bins=cropped_energy,
                                 weights=cropped_data.values(), rwidth=.9)
   bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

   # Filters out background datasets
   split_foldername = [i for i in dataset[0].split("deg")]
   fit_uncertainties = np.zeros(3)
   if "background" not in split_foldername[1][-10:]:
      # Fits gaussian
      run_fit = True
      while run_fit:
         try:
            gaussian_popt, gaussian_pcov = curve_fit(gaussian_fit, bin_centers,
                                                     hits)
            binned_fit = gaussian_fit(bin_centers, *gaussian_popt)
            cropped_counts = binned_fit
            fit_uncertainties = np.sqrt(np.diag(gaussian_pcov))
            plt.plot(bin_centers, cropped_counts, 'g',
                     label="Gaussian unbinned fit")
            run_fit = False
         except RuntimeError:
            run_fit = False

   # Finds peaks
   peak_index, _ = find_peaks(cropped_counts, distance=100)
   peak_energy = [bin_centers[i] for i in peak_index]
   peak_counts = [hits[i] for i in peak_index]
   peaks = list(zip(peak_energy, peak_counts))

   # Plots figure
   plt.errorbar(bin_centers, hits, yerr=np.sqrt(hits), fmt='b.', capsize=3,
                label="Error")
   plt.vlines(peaks, 0, max(cropped_counts), color='r', label="Peak")
   plt.xlim(left_crop, right_crop)
   plt.title(
      dataset[0] + " (elapsed time: ~{:.0f} min)".format(elapsed_time / 60))
   plt.xlabel("Energy [e-19 eV]")
   plt.ylabel("Particle Counts [counts]")
   plt.grid(alpha=0.3)
   plt.legend()
   plt.show()

   return peaks, elapsed_time, total_count, fit_uncertainties

left_crop = 0.75
right_crop = 2.1
hist_data = {}
for dataset in all_datasets.items():
    histogram = hist_plotter(dataset, left_crop, right_crop)
    hist_data[dataset[0]] = {"peaks": histogram[0],
                            "elapsed time": histogram[1],
                            "total count": histogram[2],
                            "fit uncertainties": histogram[3]}