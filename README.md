<h2 align="center">Random Forest prediction of elemental concentrations from soil NIR spectra</h1>

The measurement of elemental concentrations in soil is an expensive and hazardous process. Near-infrared (NIR) spectroscopy is an alternative analytical method which can be implemented based on calibration models for prediction. For such, Random Forest (RF) machine learning model can be used.

In this repository I bring my own R script for prediction of soil elemental concentrations from NIR spectra using RF, comprising an example dataset and the outputs of the script.

The script implements the modeling procedure in the following steps:<br>
|     Importing and loading required packages<br>
|     Declaring functions and fonts to be used<br>
|     Importing the dataset and defining the variables<br>
|     Cleaning and filtering the dataset<br>
|     Splitting the dataset into calibration and validation sets<br>
|     Preprocessing spectra using five different methods (Savitzky-Golay derivative, Standard Normal Variate, Multiplicative Scatter Correction, Detrend Normalization, and Continuum Removal) <br>
|     Modeling and prediction procedures condensed in a FOR loop for each individual element in the dataset<br>
|     Combining the results in csv tables<br>

The outputs of the script are figures of the predicted versus observed values, figures of the variables importance in the prediction, tables of calibration results and validation results for each individual prediction and two tables containing all of the calibration and validation results for all predictions.
