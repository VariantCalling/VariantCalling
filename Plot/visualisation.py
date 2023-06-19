# -*- coding = utf-8 -*-
# @time:19/06/2023 09:20
# Author:Yunbo Long
# @File:visualisation.py
# @Software:PyCharm

import matplotlib.pyplot as plt
import numpy as np

import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
# try

NUCLEOTIDES = ['A', 'C', 'G', 'T']

def plot_alignment(alignment):
    if isinstance(alignment, list):
        alignment = np.array(alignment)

    transdict = {nucleotide: i for i, nucleotide in enumerate(NUCLEOTIDES)}
    alignment_ints = np.vectorize(transdict.get)(alignment)

    nucleotide_colors = ['darkgreen', 'darkblue', 'darkred', 'black']
    alignment_colors = [nucleotide_colors[n] for n in alignment_ints.flatten()]
    plt.imshow(alignment_ints, cmap='jet')
    plt.xticks([])
    plt.yticks([])
    plt.grid(visible=False)
    plt.colorbar(ticks=[0, 1, 2, 3], orientation='vertical')
    plt.imshow(np.zeros_like(alignment_ints), cmap='jet', alpha=0.4)
    for i in range(alignment_ints.shape[0]):
        for j in range(alignment_ints.shape[1]):
            plt.text(j, i, alignment[i][j], ha='center', va='center', color=alignment_colors[i * alignment_ints.shape[1] + j])
    plt.show()


# Example usage
alignment_example = [
    ['A', 'C', 'C', 'T'],
    ['A', 'C', 'G', 'T'],
    ['A', 'C', 'C', 'T'],
    ['A', 'C', 'G', 'T']
]
plot_alignment(alignment_example)