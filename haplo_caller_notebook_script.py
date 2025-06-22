import tensorflow as tf
import VariantCalling as vc
import Comparator
import pickle
import numpy as np
import pandas as pd
import importlib
importlib.reload(Comparator)

# User-defined Parameters
sequence_selected = "CRT" #  "CRT" | "DHPS" | "DHFR" | "RBM200"
barcode = "02" # Options: 00-09

# Declaration of parameters based on user input
if sequence_selected == "CRT":
    img_row, seq_length = 101, 178
    clone_names = ["3D7","DD2","7G8"]
    truncation_range = [[0,178]]
    nb_mutations = len(clone_names)
    clone_file = "crt_clones.txt"
    model_name = "Comparator_CRT.keras"
elif sequence_selected == "DHPS":
    img_row, seq_length = 101, 642
    clone_names = ["3D7","DD2","HB3"]
    truncation_range = [[70,130],[600,642]]
    nb_mutations = len(clone_names)
    clone_file = "dhps_clones.txt"
    model_name = "Comparator_DHPS_truncate.keras"
elif sequence_selected == "DHFR":
    img_row, seq_length = 101, 491
    clone_names = ["3D7","7G8","DD2","HB3"]
    truncation_range = [[100,175],[280,350]]
    nb_mutations = len(clone_names)
    clone_file = "dhfr_clones.txt"
    model_name = "Comparator_DHFR_truncate.keras"
elif sequence_selected == "RBM200":
    img_row, seq_length = 101, 200
    clone_names = ["WUHAN","ALPHA","DELTA","OMICRON"]
    truncation_range = [[0,200]]
    nb_mutations = len(clone_names)
    clone_file = "rbm_clones_200.txt"
    model_name = "Comparator_SARS-Cov2.keras"

# Define the output file
in_file = r"/content/FAY33695_pass_barcode" + barcode + ".fastq"
out_file = in_file.replace(".fastq",".pkl")

import time
time_start = time.time()

model_path = r"/content/VariantCalling/comparator_models/" + model_name
model = Comparator.load_comparator(model_path)
tracker, pred_list = Comparator.run_comparator_inference(model,out_file,clone_file)
print("\nClone\t\tProportion")
for i in range(len(tracker)-1):
    print(clone_names[i] + ":\t\t" + str(round(tracker[i]*100/sum(tracker),2)) + "%")
print("Unknown:\t" + str(round(tracker[len(tracker)-1]*100/sum(tracker),2)) + "%")
print(in_file)
time_end = time.time()
print("Time elapsed (s):\t",time_end - time_start)
