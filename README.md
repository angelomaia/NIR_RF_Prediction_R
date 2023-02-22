<h2 align="center">Random Forest prediction of elemental concentrations from soil NIR spectra</h1>

Near-infrared (NIR) spectroscopy is an alternative analytical method which can be implemented using regression models for prediction. 

In this present work, Random Forest (RF) machine learning model was used for prediction.

In this repository I bring my own R script for prediction of soil elemental concentrations from NIR spectra using RF, comprising an example dataset and the outputs of the script.

```
The script implements the modeling procedure in the following steps:

1. Importing and loading required packages
2. Declaring functions and fonts to be used
3. Importing the dataset and defining the variables
4. Cleaning and filtering the dataset
5. Splitting the dataset into calibration and validation sets
6. Preprocessing spectra using five different methods (Savitzky-Golay derivative, Standard Normal Variate, Multiplicative Scatter Correction, Detrend Normalization, and Continuum Removal) 
7. Modeling and prediction procedures condensed in a FOR loop for each individual element in the dataset
8. Combining the results in csv tables
```

The outputs of the script are tables of calibration results and validation results for each individual prediction and two tables containing all of the calibration and validation results for all predictions, figures of the predicted versus observed values, and figures of the variables' importance in the prediction.

You can check out a script to generate prediction scores visualization on [this repository](https://github.com/angelomaia/Pred_Scores_Visualization_R).
