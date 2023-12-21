import numpy as np
import matplotlib.pyplot as plt
import pickle

clone_3D7 = "TGTCTTGGTAAATGTGCTCATGTGTTTAAACTTATTTTTAAAGAGATTAAGGATAATATTTTTATTTATATTTTAAGTATTATTTATTTAAGTGTATGTGTAATGAATAAAATTTTTGCTAAAAGAACTTTAAACAAAATTGGTAACTATAGTTTTGTAACATCCGAAACTCACAACT"
clone_DD2 = "TGTCTTGGTAAATGTGCTCATGTGTTTAAACTTATTTTTAAAGAGATTAAGGATAATATTTTTATTTATATTTTAAGTATTATTTATTTAAGTGTATGTGTAATTGAAACAATTTTTGCTAAAAGAACTTTAAACAAAATTGGTAACTATAGTTTTGTAACATCCGAAACTCACAACT"
clone_7G8 = "TGTCTTGGTAAATGTGCTCATGTGTTTAAACTTATTTTTAAAGAGATTAAGGATAATATTTTTATTTATATTTTAAGTATTATTTATTTAAGTGTAAGTGTAATGAATACAATTTTTGCTAAAAGAACTTTAAACAAAATTGGTAACTATAGTTTTGTAACATCCGAAACTCACAACT"

length = len(clone_3D7)

x_values = np.linspace(0, length - 1, length)
std_dev_value = 1

def all_same(*args):
    return all(var == args[0] for var in args[1:])


def gaussian(x, mean, std_dev=1):
    """
    Gaussian (normal) distribution function.

    Parameters:
    - x: Input values.
    - mean: Mean (center) of the distribution.
    - std_dev: Standard deviation (spread or width) of the distribution.

    Returns:
    - Gaussian distribution values corresponding to the input x.
    """
    return 1 / (std_dev * np.sqrt(2 * np.pi)) * np.exp(-(x - mean)**2 / (2 * std_dev**2))


def get_weight_arr(*clones):
    """
    Generate the weight array to apply the weights to the location with likely mutations
    """
    weight_arr = np.zeros(length)
    index = 0
    for chars in zip(*clones):
        if not all_same(*chars):
            weight_window = gaussian(x_values, index, std_dev_value)
            normalised = weight_window / np.max(weight_window)
            plt.plot(normalised, linestyle=':')
            weight_arr += normalised
        index += 1

    weight_arr += 2  # This is needed to increase the importance of other terms
    weight_arr = list(weight_arr)
    return weight_arr

weight_arr = get_weight_arr(clone_3D7, clone_DD2, clone_7G8)

# plt.plot(weight_arr)
# plt.show()

def get_vector(sequence: str, *clones):
    vector = [0 for _ in range(len(clones))]
    
    for i, char in enumerate(sequence):
        for j, clone in enumerate(clones):
            if char != clone[i]:
                vector[j] += weight_arr[i]

    return np.array(vector)

clone_3D7_vector = get_vector(clone_3D7, clone_3D7, clone_DD2, clone_7G8)
clone_DD2_vector = get_vector(clone_DD2, clone_3D7, clone_DD2, clone_7G8)
clone_7G8_vector = get_vector(clone_7G8, clone_3D7, clone_DD2, clone_7G8)


def compare_vector(sequence, *clone_vectors):
    test_vector = get_vector(sequence, clone_3D7, clone_DD2, clone_7G8)
    
    compare_index = np.array([np.linalg.norm(test_vector - clone_vector) for clone_vector in clone_vectors])
    
    clone_type = np.argmin(compare_index)
    confidence_index = 1 / np.sum(compare_index)
    return clone_type, confidence_index, compare_index


def run_calling(file):
    with open(file, 'rb') as f:
        data = pickle.load(f)

    count = [0, 0, 0]
    for datum in data:
        seq = datum[2]
        clone_type, _, _ = compare_vector(seq, clone_3D7_vector, clone_DD2_vector, clone_7G8_vector)
        count[clone_type] += 1

    print(np.array(count) / len(data))

run_calling("pickle/CRT_denoised/trimmed_data_CRT3D71.pkl")
run_calling("pickle/CRT_denoised/trimmed_data_CRT3D72.pkl")

run_calling("pickle/CRT_denoised/trimmed_data_CRTDD21.pkl")
run_calling("pickle/CRT_denoised/trimmed_data_CRTDD22.pkl")

run_calling("pickle/CRT_denoised/trimmed_data_CRT7G81.pkl")
run_calling("pickle/CRT_denoised/trimmed_data_CRT7G82.pkl")

def wc_evaluator_file(pickle_file, *clones):  # should be methods of a class so that the clones can be initialised in the INIT function
    """Takes in the pickle file then output the map ratio of each clone of the pickle file"""
    pass

def wc_evaluator_array(pileup, *clones):
    """Takes in the pileup then output the map ratio of each clone of the pileup"""
    pass