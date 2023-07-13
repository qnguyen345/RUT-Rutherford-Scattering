# Rutherford Scattering Data Analysis
The theory and full data analysis will be done in the experimental report.

**Read_Data.py**: Reads the ".dat" files and creates a table with columns "energy" and "count". 
* The code assumes that the Data file is in the same directory as the Read_Data.py file. 
* The Data file has subfiles that contains the ".dat" files like "1deg_01_5min" or "2.5deg_01_5min". 
* Note: The data directly comes from the lab computer.

**Histograms_and_Gaussian_Fit.py**: Makes a histogram plot of the data (Energy vs. Particle Counts) through hist_plotter method. 
* The plot also displays the peak energy of the histogram and a gaussian fit.


**Rutherford_Scattering_Eqn_Fit_and_Errors.py**: Makes a linear fit based on the Rutherford Scattering Formula.
* From the Rutherford Scattering Equation, the particle count ratio is given by some constant (C) * sin^-4(θ/2).
* So, if we scaled the x-axis to sin^-4(θ/2), we can make a linear fit and find the observed constant (C_observed).
* We compare our observed constant with the expected constant calculated from the Rutherford Scattering formula.
* We choose to do the particle count ratio because it's easier to fit and analyze the data.
