<h2 align="center">Random Forest prediction of elemental concentrations from soil NIR spectra</h1>

The measurement of elemental concentrations in soil is an expensive and hazardous process. Near-infrared (NIR) spectroscopy is an alternative analytical method which can be implemented based on calibration models for prediction. For such, Random Forest (RF) machine learning model can be used.

I bring here my own R script for prediction of soil elemental concentrations from NIR spectra using RF, followed by an example dataset containing 178 samples with Al, Fe, Mn, Ti, Th and V concentrations, and the outputs of the script.

The dataset must contain columns of samples id, known elemental concentrations measured by standard methods (for calibration), and the reflectance values for each individual wavelength (spectra). Any column containing sample information must be located between sample id and the first column of elemental concentration.

The script implements the modeling procedure in the following steps:<br>
>1) Importing and loading required packages<br>
>2) Declaring functions and fonts to be used<br>
>3) Importing the dataset and defining the variables which will be used<br>
>4) Cleaning and filtering the dataset<br>
>5) Splitting the dataset into calibration and validation sets<br>
>6) Preprocessing spectra using five different methods (Savitzky-Golay smoothening, Standard Normal Variate, Multiplicative Scatter Correction, Detrend Normalization, and Continuum Removal) <br>
>8) Modeling and prediction procedures condensed in a FOR loop for each individual element in the dataset<br>
>9) Combining the evaluation results in csv tables<br>

The outputs of the script are figures of the predicted versus observed values, tables of calibration results and validation results for each individual prediction and two tables containing all of the calibration and validation results for all predictions. All models are stored as .rda files (R objects) and can be loaded afterwards.
