<h2 align="center">Random Forest prediction of elemental concentrations from soil NIR spectra</h1>

Near-infrared (NIR) spectroscopy is an alternative analytical method which can be implemented for the estimation of elemental concentrations using regression models for prediction. 

In this repository I provide my own R script for prediction of soil elemental concentrations from NIR spectra using Random Forest regression, comprising an example dataset and the outputs of the script.

```
The script implements the modeling procedure in the following steps:

1. Loading required packages
2. Declaring functions
3. Importing the dataset
4. Cleaning and filtering data
5. Splitting the dataset into calibration and validation sample sets
6. Preprocessing spectra using five different methods: Savitzky-Golay derivative, Standard Normal Variate, Multiplicative Scatter Correction, Detrend Normalization, and Continuum Removal
7. Modeling Random Forest Regression predictions in a 'for' loop for each element in the dataset
8. Combining the results and exporting in csv tables
```

The output of the script: tables of calibration results and validation results for each individual prediction and two tables containing all of the calibration and validation results for all predictions, figures of the predicted versus observed values and the variables importance in the prediction.

You can check out a script to generate prediction scores visualization on [this repository](https://github.com/angelomaia/Pred_Scores_Visualization_R).
