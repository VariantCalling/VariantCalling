import random
import math
import numpy as np
import matplotlib.pyplot as plt
import sys
from tqdm.notebook import tqdm
import PickleUtil as PU

class VariantCalling:
    def __init__(self, mutation_labels, mutation_types_names, file_name, gen_mode, pkl_sequence, pkl_clones) -> None:
        self.mutation_labels = mutation_labels
        self.mutation_type_names = mutation_types_names
        self.NUCLEOTIDES = "ACGT"
        self.transdict = {"A":0, "C": 1, "G":2, "T":3}
        self.transdict_reverse = {"0":"A", "1": "C", "2":"G", "3":"T"}
        self.pkl_sequence = pkl_sequence
        self.pkl_clones = pkl_clones
        self.pkl_alignment_list = [] # This is the variable that holds the real-world data converted
        self.gen_mode = gen_mode # We define the mode of generation, with 1: Sim, 2: Pickle, we default to 1 to maintain conpatibility
        self.clones = []
        with open(file_name, "r") as f:
            for clone in f:
                alignment = []
                for char in clone.strip():
                    alignment.append(char)
                self.clones.append(alignment)
        self.nb_clones = len(self.clones)
        self.clones_int = self.char_to_int(self.clones)
        if self.gen_mode == 2: # We only load the pickle if necessary
            self._load_pickle()
            self.nb_clones = len(self.pkl_clones)

class VariantCallingData(VariantCalling):
    """Class for simulated data generation"""
    def __init__(self, 
                 mutation_labels={"no_SNP": 0, "heterozygous_SNP": 1, "homozygous_SNP": 2},
                 mutation_types_names={0: "No mutation", 1: "Heterozygous SNP", 2: "Homozygous SNP"},
                 file_name="clones.txt",
                 gen_mode=1,
                 pkl_sequence="CRT",
                 pkl_clones=["3D7"]
                 ) -> None:
        super().__init__(mutation_labels=mutation_labels, mutation_types_names=mutation_types_names,
                         file_name=file_name,gen_mode=gen_mode,pkl_sequence=pkl_sequence,pkl_clones=pkl_clones)
        self.alignments = None
        self.mutation_types = None
                
    def simulate_alignments(self, reference_length=200, 
                        num_alignments = 2000, 
                        coverage = 100, 
                        mutations = None,
                        p_sequencing_error=0.0,
                        p_alignment_error=0.00) -> tuple:
        """Modified provided function to generate simulated data for the project.
        
        Inputs
        ------
        reference_length: Width of the image.
        num_alignments: Size of data samples (1 sample is 1 image).
        coverage: Height of the image.
        mutations: Contains the keys of the mutation types. 
                    If None, use self.mutation_labels.
        p_sequencing_error: Probability of sequencing error.
        p_alignment_error: Probability of alignment error.

        Outputs
        -------
        tuple: (alignments, mutation_types)
        alignments: 3D np.array containing the data in the shape 
                    (num_alignments, coverage + 1, reference_length).
        mutation_types: python list of integers indicating the mutation 
                        type of each sample of data in alignments. 
                        size of num_alignments.
        """
        if mutations is None:
            mutations = self.mutation_labels.keys()

        alignments = []
        mutation_types = []
        
        for i in range(num_alignments):
            snp_index = reference_length // 2 
            if (i % 400 == 0):
                print("Computing alignment ", i)
            reference = [random.choice(self.NUCLEOTIDES) for _ in range(reference_length)]
            reference_base_at_snp = reference[snp_index]
            snp_base = random.choice([i for i in self.NUCLEOTIDES if i != reference_base_at_snp])
            mutation_type=random.choice([self.mutation_labels[m] for m in mutations]) 
            mutation_types.append(mutation_type)
            
            alignment = [reference] #first read is always the reference
            for _ in range(coverage):
                mut_index = snp_index
                new_read = [reference[i] if random.random() > p_sequencing_error 
                            else random.choice(self.NUCLEOTIDES) for i in range(reference_length)]
                if random.random() < p_alignment_error:
                    mut_index = snp_index + random.randint(-1,2)
                if mutation_type == 1 and random.random() > 0.5: # heterozygous SNP
                    new_read[mut_index] = snp_base            
                if mutation_type == 2: #homozygous SNP
                    new_read[mut_index] = snp_base
                if random.random() < p_sequencing_error: #Add errors to SNP region also
                    new_read[mut_index] =  random.choice(self.NUCLEOTIDES)
                alignment.append(new_read)
            alignments.append(alignment)

        self.alignments = alignments
        self.mutation_types = mutation_types
        return np.array(alignments), mutation_types
    

    def char_to_int(self, alignments=None) -> np.ndarray:
        """Maps the char ACGT to the corresponding integers"""
        if alignments is None:
            alignments = self.alignments
        return np.vectorize(self.transdict.get)(alignments)

    def int_to_char(self, alignments=None) -> np.ndarray:
        """Maps the char ACGT to the corresponding integers"""
        if alignments is None:
            alignments = self.alignments
        return np.vectorize(self.transdict_reverse.get)(alignments.astype("str"))
    
    def plot_data(self, alignments_ints=None, mutation_types=None, mutation_index=0) -> None:
        """Use to plot an alignment at a certain mutation type
        
        Note
        ----
        This Function should be refactored, or a new function should be created
        """
        if alignments_ints is None or mutation_types is None:
            alignments_ints = self.char_to_int()
            mutation_types = self.mutation_types
        
        plt.rcParams['figure.dpi'] = 200
        alignment_idx = mutation_types.index(mutation_index)
        plt.title(f"Mutation type: {self.mutation_type_names[mutation_types[alignment_idx]]}")
        plt.imshow(alignments_ints[alignment_idx],cmap='jet')

    def simulate_clones(self, num_alignments = 2000, 
                        coverage = 100, 
                        p_sequencing_error=0.0,
                        p_alignment_error=0.00,
                        verbose=1) -> (np.ndarray, list) :
        """Wrapper to generate n alignments as specified in num_alignments.
        NOTE: Consider to merge this into simulate_alignments in the future.
        
        Parameters
        ----------
        num_alignments : <int>
            Number of alignments to be generated
        coverage : <int>
            Number of read for an alignment
        p_sequencing_error : <double>
            Probability of sequencing error, takes value >= 0, <= 1
        p_alignment_error : <double>
            Probability of alignment error, takes value >= 0, <= 1 
        verbose : <int>
            0 - No progress will be printed
            1 - Progress will be printed for every 400 num_alignments
        Returns
        -------
        np.ndarray :
            Numpy array of alignments generated.
        list :
            List of lists of probability for each of the alignment
        """
        alignments = []
        prob_lists = []
        row_clones = []
        for _ in tqdm(range(num_alignments)):
            alignment, prob_list, row_clone = self.ratio_gen(coverage, p_sequencing_error, p_alignment_error)
            alignments.append(alignment)
            prob_lists.append(prob_list)
            row_clones.append(row_clone)

        self.alignments = alignments
        print(f"Done, Number of alignments: {num_alignments}")
        return np.array(alignments), prob_lists, row_clones

    def ratio_gen(self, coverage, p_sequencing_error=0, p_alignment_error=0) -> (list,list):
        """Wrapper to generate a single alignment based on a randomly generated ratio
        Returns np.ndarray of the alignment and the probability of the distribution
        
        NOTE: Naming of variables shall be improved in the future for easy readability
        
        Parameters
        ----------
        coverage : <int>
            Number of read for an alignment
        p_sequencing_error : <double>, optional
            Probability of sequencing error, takes value >= 0, <= 1
        p_alignment_error : <double>, optinal
            Probability of alignment error, takes value >= 0, <= 1 
        gen_mode : int, optional
            Mode 1: Simulation
            Mode 2: Pickle
        
        Returns
        -------
        list, list
            List of reads for an alignment
        list
            Probability distribution for the alignment read for each of the clone class
        """
        prob_dist = self._gen_prob_list(self.nb_clones, mode=2)

        nb_coverage_list = []
        for prob in prob_dist:
            nb_coverage_list.append(math.floor(prob * coverage))

        for _ in range(coverage - sum(nb_coverage_list)):
            # We randomly increase an element by 1 until we reach the number of coverages specified (n - 1)
            # as the first row is always the reference
            nb_coverage_list[random.randint(0,self.nb_clones - 1)] += 1        

        coverage_list = []
        clone_idx_list = []
        for clone_idx, nb_clone_coverage in enumerate(nb_coverage_list):
            for _ in range(0,nb_clone_coverage):
                match self.gen_mode:
                    case 1: # Simulation
                        coverage_list.append(
                        self._add_errors(self, self.clones[clone_idx],p_sequencing_error,p_alignment_error))
                    case 2: # Load from pickle
                        coverage_list.append(list(self.pkl_alignment_list[clone_idx][self._gen_rand_nb(len(self.pkl_alignment_list[clone_idx])-1)]))
                clone_idx_list.append(clone_idx)

        # This will be the final probability list
        prob_list = [nb_coverage_list[i]/coverage for i in range(0, len(nb_coverage_list))]

        # Here we shuffle the list and concatenate into the final alignment
        choice_indices = np.random.choice(len(coverage_list), coverage, replace=False)
        alignment = [self.clones[0]] # First row is always reference (assumed to be index at 0)
        row_clone = [0]
        alignment += [coverage_list[i] for i in choice_indices] # Concatenate the randomized read to the reference row
        row_clone += [clone_idx_list[i] for i in choice_indices] # Concatenate the randomized read to the reference row

        return alignment, prob_list, row_clone

    @staticmethod
    def _add_errors(self, clone, p_sequencing_error, p_alignment_error) -> list:
        """
        Adds sequencing error and alignment error to a single read, returns clone with error
        Parameters
        ----------
        clone : <List> 
            List of the base-pair encoded in <int>
        p_sequencing_error : <double>
            Probability of sequencing error, takes value >= 0, <= 1
        p_alignment_error : <double>
            Probability of alignment error, takes value >= 0, <= 1
        
        Returns
        -------
        list
            List of bp based on input clone with sequencing and alignment errors added
        """
        # Let's make alignment error applicable to all for now
        new_clone = [clone[i] if random.random() > p_alignment_error 
                     else clone[min(max(0, i + random.randint(-1,2)),len(clone)-1)] for i in range(len(clone))]
        return [new_clone[i] if random.random()> p_sequencing_error
                else random.choice(self.NUCLEOTIDES) for i in range(len(new_clone))]

    def _gen_prob_list(self, nb_class, mode=1) -> list:
        """Generate a list of nb_class elements of probability that sum to 1
        Leaving the nb_class to maintain modularity in case we need to generate 
        probability for other purposes.
        
        Parameters
        ----------
        nb_classes : <int>
            Number of classes of probability to be generated
        
        mode : <int>
            Mode of the randomizer, it seems like the Dirichlet's Distribution outputs overly balanced distribution
            1 : Dirichlet's Distribution - Seems to be quite balanced, often returns within mean of p=0.33
            2 : random.random() - Generates nb_class random numbers that are normalized against the sum
        Returns
        -------
        list
            List of <nb_class> element representing the probability for each of the class.
        """
        match mode:
            case 1: 
                return (np.random.dirichlet(np.ones(nb_class)*1000.,size=1)).flatten().tolist()
            case 2:
                prob_list = [random.random() for _ in range(0,nb_class)]
                return [prob_list[i] / sum(prob_list) for i in range(0, nb_class)]

    def _array_dup(self,arr,coverage) -> np.ndarray: 
        """Produces an array of the reference genome repeated - this is used as a second channel in CNN, Will's Code
        
        Parameters
        ----------
        arr : TYPE
            Input alignment data
        
        coverage : TYPE
            Size of the amplicon
        coverage: 
        
        Returns
        -------
        np.ndarray
            2-dimensional array with second channel as the reference row
        """

        # NOTE: We should be able to get size of coverage from len(arr[0]) but faced np.tile issue (the output shape is different)
        # when passing len(ref_gen) instead of coverage to np.tile, non-blocking but good to be solved in the future.
        ref_gen = arr[0]
        ref_gen_matrix = np.tile(ref_gen,(coverage,1))
        aln_ref_dim = np.array((arr, ref_gen_matrix))
        return aln_ref_dim

    def _array_dup_all(self,arr,coverage) -> np.ndarray:
        """Wrapper to call _array_dup_clone to generate channels for different clones, basically instead of only 1 channel, this wrapper 
        loops and generates additional channels for all the clones considered.
        
        Parameters
        ----------
        arr : TYPE
            Description
        coverage : TYPE
            Description
        """

        # For each clones in the clones.txt file, we create a separate channel for it
        arr_channel = np.array([arr])
        for clone in self.clones:
            arr_channel = self._array_dup_clone(clone,arr_channel,coverage)
        return arr_channel

    def _array_dup_clone(self,clone,arr,coverage) -> np.ndarray: 
        """Produces an array of the reference genome repeated, mc for multiple channel, this is built on top of Will's code and supposed to work
        better as it should work for more than one. The input parameter has been changed to specific alignment row.
        
        Parameters
        ----------
        clone : List
            Row of clone to be duplicate as a separate channel
        arr : List
            Input alignment data
        coverage : Int
            Size of the amplicon
        
        Returns
        -------
        np.ndarray
            2-dimensional array with second channel as the reference row
        """
        clone_gen_matrix = np.array([np.tile(clone,(coverage,1))])
        aln_clone_dim = np.concatenate((arr,clone_gen_matrix))
        return aln_clone_dim

    def _array_binary_mask(self,arr,coverage,ref_mat=False) -> np.ndarray:
        """Generates binary masks based on whether the base pair is the same as the reference sequence
        (0: False   1: True)
        
        Parameters
        ----------
        arr : <List>
            2D input list of the input array with the first row containing the reference sequence
        coverage : <int>
            Size of the coverage for the input array, we should be able to identify the size from arr 
            but seems to be broken right now, may need to be fixed in the future.
        ref_mat : bool, optional
            Boolean option to either add the reference matrix into the return array
        
        Returns
        -------
        np.ndarray
            Np array of multi-channel depends on the nature
        """
        # Here we are transcoding the True to A as it's the highest value
        binary_transdict = {True:"A", False:"T"}
        ref_gen = arr[0]
        #binary_mask_matrix = [np.in1d(arr[i],ref_gen).tolist() for i in range(0,len(arr))]
        binary_mask_matrix = [[arr[j][i]==ref_gen[i] for i in range(0,len(arr[j]))] for j in range(0,len(arr))]
        binary_mask_matrix = np.vectorize(binary_transdict.get)(binary_mask_matrix)

        if ref_mat:
            return np.array((arr, binary_mask_matrix, np.tile(ref_gen,(coverage,1))))
        else:
            return np.array((arr, binary_mask_matrix))

    def _array_constituent_bp(self,arr,coverage, binary_map=False) -> np.ndarray:
        """Generates binary masks based on whether the base pair is the same as the reference sequence
        (0: False   1: True)
        
        Parameters
        ----------
        arr : <List>
            2D input list of the input array with the first row containing the reference sequence
        coverage : <int>
            Size of the coverage for the input array, we should be able to identify the size from arr 
            but seems to be broken right now, may need to be fixed in the future.
        binary_map : bool, optional
            Boolean to either include the binary map
        
        Returns
        -------
        np.ndarray
            Description
        
        Deleted Parameters
        ------------------
        ref_mat : bool, optional
            Boolean option to either add the reference matrix into the return array
        """
        binary_transdict = {True:"T", False:"A"}
        matrix_A = np.vectorize(binary_transdict.get)([[True if arr[j][i] == "A" else False for i in range(0,len(arr[j]))] for j in range(0,len(arr))])
        matrix_T = np.vectorize(binary_transdict.get)([[True if arr[j][i] == "T" else False for i in range(0,len(arr[j]))] for j in range(0,len(arr))])
        matrix_C = np.vectorize(binary_transdict.get)([[True if arr[j][i] == "C" else False for i in range(0,len(arr[j]))] for j in range(0,len(arr))])
        matrix_G = np.vectorize(binary_transdict.get)([[True if arr[j][i] == "G" else False for i in range(0,len(arr[j]))] for j in range(0,len(arr))])

        if binary_map:
            ref_gen = arr[0]
            #binary_mask_matrix = [np.in1d(arr[i],ref_gen).tolist() for i in range(0,len(arr))]
            binary_mask_matrix = [[arr[j][i]==ref_gen[i] for i in range(0,len(arr[j]))] for j in range(0,len(arr))]
            return np.array((matrix_A, matrix_T, matrix_C, matrix_G, np.vectorize(binary_transdict.get)(binary_mask_matrix)))
        else:
            return np.array((matrix_A, matrix_T, matrix_C, matrix_G))

        aln_binary_mask_dim = np.array((arr, binary_mask_matrix))
        return aln_binary_mask_dim

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

    def generate_data_for_noise_reduction(self, sample_size=1000, image_height=100, alignment_error_prob=0.05, sequencing_error_prob=0.05):
        """Function which returns noisy and non-noisy data"""

        noisy_images = []
        clean_images = []

        for _ in tqdm(range(sample_size)):
            noisy = []
            clean = []
            random_array = np.random.randint(len(self.clones), size=image_height)
            for clone_type in random_array:
                noisy.append(self._simulate_read(clone_type, alignment_error_prob, sequencing_error_prob))
                clean.append(self.char_to_int(list(self.clones[clone_type])))
            noisy_images.append(noisy)
            clean_images.append(clean)

        return np.array(noisy_images), np.array(clean_images)

    def _load_pickle(self) -> None:
        """Static method to load pickle files
        """
        pickleloader = PU.PickleLoader(sequence=self.pkl_sequence, clones=self.pkl_clones)
        self.pkl_alignment_list = pickleloader.load_pickle()

    def _gen_rand_nb(self, u_bound, l_bound=0) -> int:
        """Generate a random integer from 0 to limit
        
        Parameters
        ----------
        u_bound : int
            Upper bound of range
        l_bound : int, optional
            Lower bound, default to 0
        
        Returns
        -------
        int
            Description
        """
        return random.randint(0,u_bound)
