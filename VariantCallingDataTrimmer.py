import pandas as pd
import subprocess
import re
import pickle
from datetime import datetime

cigar_num_pattern = r"\d+"
cigar_alpha_pattern = r"[A-Z]"
seq_len = 491
now = datetime.now()

print("Running subprocess")
process = subprocess.run(['bash', './command.sh'], capture_output=True, text=True)
with open('./data/log/_log', 'a') as log_file:
    log_file.write("################\n")
    log_file.write(f"Log for {now}\n")
    log_file.write("################\n")
    log_file.write(process.stdout)
    log_file.write(process.stderr)
    log_file.write("######END#######\n")
print("Process Completed")


data = []
with open('./data/output/sorted_example_alignment.txt', 'r') as f:
    for line in f:
        line_array = line.split('\t')
        dna_seq = line_array[9]
        cigar_str = line_array[5]

        cigar_nums = re.findall(cigar_num_pattern, line_array[5])
        cigar_alphas = re.findall(cigar_alpha_pattern, line_array[5])
        if len(cigar_nums) != 0:
            if cigar_alphas[0] == 'S':
                dna_seq = dna_seq[int(cigar_nums[0]):]
            if len(dna_seq) < seq_len:
                continue

            interested = [line_array[2], cigar_str, dna_seq[:seq_len]]
            data.append(interested)

print("number of alignments", len(data))

# Open a file in binary mode
with open(f"./data/pickle/trimmed_data_{now.strftime(r'%Y%m%d%H%M%S')}.pkl", "wb") as file:
    # Dump the list to the file using pickle
    pickle.dump(data, file)

print("files saved at", f"trimmed_data_{now.strftime(r'%Y%m%d%H%M%S')}.pkl")