import random
import math
import numpy as np
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm

HEIGHT = 100
WIDTH = 178

class VariantCalling:
    def __init__(self, mutation_labels, mutation_types_names, file_name) -> None:
        self.mutation_labels = mutation_labels
        self.mutation_type_names = mutation_types_names
        self.NUCLEOTIDES = "ACGT"
        self.transdict = {"A":0, "C": 1, "G":2, "T":3}
        self.reverse_transdict = {0: "A", 1: "C", 2: "G", 3: "T"}
        self.nb_clones = len(self.clones)
        self.clones_int = self.char_to_int(self.clones)

class VariantCallingDataMini(VariantCalling):
    """Class for simulated data generation"""
    def __init__(self,
                 mutation_labels={"no_SNP": 0, "heterozygous_SNP": 1, "homozygous_SNP": 2},
                 mutation_types_names={0: "No mutation", 1: "Heterozygous SNP", 2: "Homozygous SNP"},
                 file_name="clones.txt"
                 ) -> None:
        self.generate_random_clones()
        super().__init__(mutation_labels=mutation_labels, mutation_types_names=mutation_types_names,
                         file_name=file_name)
        self.alignments = None
        self.mutation_types = None


    def char_to_int(self, alignments=None) -> np.ndarray:
        """Maps the char ACGT to the corresponding integers"""
        if alignments is None:
            alignments = self.alignments
        return np.vectorize(self.transdict.get)(alignments)



    def _simulate_read(self, clone_type_index: int, alignment_error_prob: float, sequencing_error_prob: float):
        sim_read = []
        copied_clone = list(self.clones[clone_type_index])  # This is so that run-time is o(n) not o(2n)
        pointer = 0
        second_alignment_error_prob = alignment_error_prob / 5  # guess work right here
        third_alignment_error_prob = alignment_error_prob /  8  # guess work right here
        forth_alignment_error_prob = alignment_error_prob /  15  # guess work right here
        while len(sim_read) < len(self.clones[clone_type_index]):
            alignment = random.uniform(0, 1)
            direction = random.choice([1, -1])
            if alignment <=  forth_alignment_error_prob:
                pointer += direction * 4
            elif alignment <= third_alignment_error_prob:
                pointer += direction * 3
            elif alignment <= second_alignment_error_prob:
                pointer += direction * 2
            elif alignment <=  alignment_error_prob:
                pointer += direction

            if pointer < 0 or pointer >= len(self.clones[clone_type_index]):
                sim_read.append(random.choice(['A', 'C', 'G', 'T']))
            else:
                sequencing = random.uniform(0,1)
                if sequencing < sequencing_error_prob:
                    current = copied_clone[pointer]
                    choice_array = ['A', 'C', 'G', 'T']
                    choice_array.remove(current)
                    copied_clone[pointer] = random.choice(choice_array)

                sim_read.append(copied_clone[pointer])

            pointer += 1
        return self.char_to_int(sim_read)

    def generate_data_for_noise_reduction(self, sample_size=1000, image_height=HEIGHT, alignment_error_prob=0.05, sequencing_error_prob=0.05):
        """Function which returns noisy and non-noisy data"""

        noisy_images = []
        clean_images = []
        all_mutation_positions = []

        for _ in range(sample_size):
            self.generate_random_clones()
            noisy = []
            clean = []
            random_array = np.random.randint(len(self.clones), size=image_height)
            for clone_type in random_array:
                noisy.append(self._simulate_read(clone_type, alignment_error_prob, sequencing_error_prob))
                clean.append(self.char_to_int(list(self.clones[clone_type])))
            noisy_images.append(noisy)
            clean_images.append(clean)
            all_mutation_positions.append(self.mutation_positions.copy())

        return np.array(noisy_images), np.array(clean_images), all_mutation_positions

    def generate_random_clones(self, num_of_mutants=3, num_of_mutations_per_clone=3, size=WIDTH):
        clones = []
        bases = ['A', 'C', 'G', 'T']
        random_seq = np.random.choice(bases, size=size).tolist()
        clones.append(random_seq)

        mean_value = (size - 1) / 2  # (0 + 177) / 2
        std_deviation = mean_value / 3.5
        # std_deviation = mean_value / 20
        mutation_pos_1 = np.random.normal(loc=mean_value, scale=std_deviation)
        mutation_pos_1 = np.clip(mutation_pos_1, 0, size - 1)
        mutation_pos_1 = int(round(mutation_pos_1))

        mutation_positions = set()
        for _ in range(num_of_mutants):
            new_variant = random_seq[:]
            for _ in range(num_of_mutations_per_clone):
                mutation_pos_variant = int(np.clip(mutation_pos_1 + np.random.normal(loc=0, scale=0.5)*10, 0, size - 1))
                new_variant[mutation_pos_variant] = np.random.choice(bases)
                mutation_positions.add(mutation_pos_variant)
            clones.append(new_variant)

        self.clones = clones
        self.mutation_positions = list(mutation_positions)
        return clones, mutation_positions

    def generate_data_for_comparator(self, sample_size=1000, image_height=HEIGHT, alignment_error_prob=0.05, sequencing_error_prob=0.05):
        """Function which returns reference, input sequences and labels"""

        ref_seq = []
        inp_seq = []
        labels = []
        for _ in range(sample_size):
            random_array = np.random.randint(len(self.clones), size=image_height)
            for clone_type in random_array:
                for i in range(len(self.clones)):
                    ref_seq.append(self.char_to_int(self.clones[i]))
                    inp_seq.append(self._simulate_read(clone_type, alignment_error_prob, sequencing_error_prob))
                    labels.append(int(clone_type == i))
        return ref_seq, inp_seq, labels
