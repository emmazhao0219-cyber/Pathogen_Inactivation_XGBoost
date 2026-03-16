# Pathogen Inactivation Prediction via Species-Agnostic XGBoost

## 📌 Overview
This repository contains the source code, pre-trained machine learning model, and data processing pipelines for our study on predicting pathogen inactivation in wastewater treatment. 

By leveraging quantitative metagenomics and extracting microdiversity-driven physical/genetic damage features (e.g., $\Delta \text{Breadth}$, $\Delta \text{SNV Density}$), we developed a species-agnostic eXtreme Gradient Boosting (XGBoost) model capable of predicting the Log Removal Value (LRV) of 31 diverse pathogens across multiple disinfection technologies.

## 📂 Repository Structure
* `/data/`: Contains raw inputs (metagenomic features) and processed data.
* `/models/`: Contains the pre-trained `XGBoost_Species_Agnostic_v1.model` and the feature encoder (`.rds`).
* `/scripts/`: Modular R scripts for data processing, training, inference, and visualization.
* `/results/`: Output prediction tables and high-resolution academic figures.

## 🚀 Usage Instructions
**Prerequisites (R packages):**
```R
install.packages(c("tidyverse", "xgboost", "caret", "readxl", "writexl", "ggpubr", "extrafont"))