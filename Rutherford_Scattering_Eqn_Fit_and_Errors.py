import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy.optimize import curve_fit
from Histograms_and_Gaussian_Fit import hist_data

def rutherford_fit(x, a):
   x = np.array(x)
   return a * x

def fit_plotter(angles, counts, elapsed_times, title):
   detector_radius = 0.009525
   detector_to_target = 0.0315
   detector_area = np.pi * detector_radius ** 2
   source_aperature = 0.0015875
   # Detector_solid_angle = 4*np.pi*(np.sin(0.5*np.arcsin(2*detector_radius/detector_to_target)))**2
   total_emitted_particles_per_sec = 37000 * (source_aperature / 2) ** 2 / (
              4 * 0.01 ** 2)
   rates = [counts[i] / (elapsed_times[i] * total_emitted_particles_per_sec) for
            i in range(len(counts))]
   rates_err = [
      np.sqrt(counts[i]) / (elapsed_times[i] * total_emitted_particles_per_sec)
      for i in range(len(counts))]

   # Formats points
   points = list(
      zip(angles, [(rates[i], rates_err[i]) for i in range(len(rates))]))
   points.sort()
   angles = [i[0] for i in points]
   rates = [i[1][0] for i in points]
   rates_err = [i[1][1] for i in points]

   # Fits to rutherford scattering eqn
   fit_angles = [(np.sin((angle * np.pi / 180) / 2)) ** (-4) for angle in angles
                 if angle > 5]
   fit_rates = rates[-len(fit_angles):]
   fit_rates_err = rates_err[-len(fit_angles):]
   x = np.linspace(min(fit_angles), max(fit_angles), 100)
   rut_popt, rut_pcov = curve_fit(rutherford_fit, fit_angles, fit_rates)
   fit_uncertainties = np.sqrt(np.diag(rut_pcov))
   binned_fit = rutherford_fit(x, *rut_popt)

   # Chi-squared gof test
   chi_squared = sum([((fit_rates[i] - rutherford_fit(fit_angles, *rut_popt)[
      i]) / (fit_rates[i])) ** 2 for i in range(len(fit_rates))])
   
   # Chi_squared = chisquare(fit_rates, rutherford_fit(fit_angles, *rut_popt), fit_rates_err)
   dof = len(fit_angles) - 1
   chi_squared /= dof
   p_val = chi2.cdf(chi_squared, dof)

   # Plots figure for particle count rate vs scattering angle
   plt.figure(figsize=(15, 10))
   plt.plot(angles, rates, 'o-')
   plt.errorbar(angles, rates, yerr=rates_err, fmt='b.', capsize=3)
   plt.title(title)
   plt.xlabel("Scattering Angle [$^\circ$]")
   plt.ylabel("Particle Count Ratio [$N_o / N_i$]")
   plt.grid(alpha=0.3)
   plt.show()

   # Plots figure with Rutherford fit
   plt.figure(figsize=(15, 10))
   plt.plot(fit_angles, fit_rates, 'o')
   plt.plot(x, binned_fit, 'g')
   plt.errorbar(fit_angles, fit_rates, yerr=fit_rates_err, fmt='b.', capsize=3)
   plt.title(title + "\n(Fitted to the Rutherford Scattering Equation)")
   plt.xlabel("Rescaled Scattering Angle, $\\sin^{-4}(\\theta/2)$")
   plt.ylabel("Particle Count Ratio [$N_o / N_i$]")
   plt.grid(alpha=0.3)
   plt.show()

   return rut_popt, fit_uncertainties, chi_squared, p_val

angles = []
peaks = []
total_particles = []
elapsed_times = []

background_angles = []
background_peaks = []
background_total_particles = []
background_elapsed_times = []

for hist in hist_data.items():
    split_foldername = [i for i in hist[0].split("deg")]
    angle = split_foldername[0]
    if "background" not in split_foldername[1][-10:]:
        angles.append(float(angle))
        peaks.append(hist[1]["peaks"][0][1])
        total_particles.append(hist[1]["total count"])
        elapsed_times.append(hist[1]["elapsed time"])
    else:
        continue

rut_fit, rut_fit_err, rut_chi_squared, rut_p_val = fit_plotter(angles, peaks, elapsed_times, "Peak Particle Count vs Scattering Angle")

C_observed = rut_fit[0]
C_observed_err = rut_fit_err[0]
#Fit_plotter(background_angles, background_peaks, background_elapsed_times, "Background (w/o vacuum)")
print("C_observed = {} \u00B1 {}".format(rut_fit[0],rut_fit_err[0]))
print("\u03c7^2/dof =",rut_chi_squared)
print("p-value =",rut_p_val)

# Calculate expected Rutherford Scattering "constant"
n = 5.9040605784e22
L = 1.3652877e-4
Z = 95
k = 8.9875517923e15
e = 1.602176634e-19
E = 6.0402059102e-9
r = 3.15

C_expected = (n * L * Z**2 * k**2 * e**4)/(4 * r**2 * E**2)
print("C_expected =", C_expected)
def rutherford_equation(C, theta):
   return C * (np.sin(theta / 2 * np.pi / 180)) ** (-4)

# Plots the Rutherford Scattering Equation with the C_observed and C_expected
theta = np.linspace(1, 25, 100)
plt.figure(figsize=(15, 10))
plt.plot(theta, rutherford_equation(C_observed, theta), label="Observed")
plt.plot(theta, rutherford_equation(C_expected, theta), label="Expected")
plt.title("Rutherford Scattering Equation")
plt.xlabel("Scattering Angle [$^\circ$]")
plt.ylabel("Particle Ratio [$No/N_i$] per unit area")
plt.grid(alpha=0.3)
plt.legend()
plt.show()




