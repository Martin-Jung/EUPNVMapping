
<!-- README.md is generated from README.Rmd. Please edit that file -->

# European potential Natural vegetation mapping

The goal of this exercise is to utilize “machine learning” approaches to
determine the potential natural vegetation (PNV) in Europe. Predictions
are made for six different vegetation classes following the MAES Habitat
Classification system.

This repository contains the analytical scripts from data preparation to
prediction.

Description of scripts:

- 00_Functions.R = Helper functions for the analysis

- 01_ParameterSettings.R = General parameter script to support
  predictions on different machines and settings

- 02_PrepareOccurrenceData.R = The preparation and harmonization of
  input training data for the prediction.

- 03_variablePrep.R = Preparation of input covariates for current and
  future states.

- 04_runIBISpnvHabitats.R = The primary inference and prediction script
  using the [ibis.iSDM](https://iiasa.github.io/ibis.iSDM/) package.

- 05_postprocessPNV.R = Posthoc correction of predictions.

- 07_XXX = Scripts to create output figures.

- 08_preparePublicRelease.R = Script to prepare all predictions for
  upload in Zenodo.

**The data has been uploaded here:**  
Martin, J. (2024). *Current and future European potential vegetation
types* (1.0) \[Data set\]. Zenodo.
<https://doi.org/10.5281/zenodo.13686776>

**A methodological description of the work can be found here:**  
Martin, J. (2024) *Mapping current and future European potential
vegetation to support restoration planning*
<https://doi.org/10.31223/X59H71>
