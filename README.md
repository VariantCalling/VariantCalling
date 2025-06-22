# 🧬 PicklSeq and Haplo Caller Demo Guide

This guide explains how to use the `PickleSeq_and_HaploCaller_WetLab_Demo.ipynb` notebook along with the supporting script `haplo_caller_notebook_script.py`. The tools are used to perform haplotype-based variant calling and evaluate clone similarity using a trained comparator model.

---

## 📁 Files Overview

- **Notebook**: `PickleSeq_and_HaploCaller_WetLab_Demo.ipynb`  
  A step-by-step Jupyter notebook that walks through sequencing, preprocessing, and variant comparison.

- **Script**: `haplo_caller_notebook_script.py`  
  A supporting Python script that loads a trained model and performs inference for clone similarity prediction. This script was created to continue
  supporting Python 3.10 in the Google Colab environment.

---

## 🛠️ User Inputs Required

### 1. In the Notebook

Inputs include: 

sequence_selected:  "CRT" | "DHPS" | "DHFR" | "RBM200"

For CRT/DHPS/DHFR, example dataset, the user selects 1 out of the 10 data available:

barcode = "02" # Options: 00-09

---

### 2. In the Script (`haplo_caller_notebook_script.py`)

The inputs shall be updated in the script too.

## 🧪 Output

The script will output:

- Clone name and its similarity proportion
- The predicted distribution for the unknown sample
- Total elapsed time for the inference run

Sample output:

```
Clone         Proportion
3D7:          35.22%
HB3:          42.13%
DD2:          20.41%
Unknown:      2.24%
Time elapsed (s):  1.23
```


# Original VariantCalling

## Library
### VariantCalling.py
This file contains the infrastructure to simulate the clone sequence as seen in each alignment reading. See file for information.

## Test Scripts
### generating_data_sample.ipynb
This notebook generates random sequence of genes.

### gen_clone_ratio.ipynb
This notebook is created to test the functions developed for the 

It generates one single alignment read with coverage of 100 based on the sample clones provided in clones.txt

## PLOS_CB_latex
This folder contains the latex package for the write-up. As the free version of overleaf doesn't support auto sync with GitHub, this folder will be manually synced and serves as a back-up for the write-up.