import random
import math
import numpy as np
import matplotlib.pyplot as plt

class VariantCalling:
    def __init__(self, mutation_labels, mutation_types_names, file_name) -> None:
        self.mutation_labels = mutation_labels
        self.mutation_type_names = mutation_types_names
        self.NUCLEOTIDES = "ACGT"
        self.transdict = {"A":0, "C": 1, "G":2, "T":3}

        self.clones = []
        with open(file_name, "r") as f:
            for clone in f:
                alignment = []
                for char in clone.strip():
                    alignment.append(char)
                self.clones.append(alignment)
        self.nb_clones = len(self.clones)
        self.clones_int = self.char_to_int(self.clones)

class VariantCallingData(VariantCalling):
    """Class for simulated data generation"""
    def __init__(self, 
                 mutation_labels={"no_SNP": 0, "heterozygous_SNP": 1, "homozygous_SNP": 2},
                 mutation_types_names={0: "No mutation", 1: "Heterozygous SNP", 2: "Homozygous SNP"},
                 file_name="clones.txt"
                 ) -> None:
        super().__init__(mutation_labels=mutation_labels, mutation_types_names=mutation_types_names,
                         file_name=file_name)
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

    def simulate_clones(self, num_alignments=2000, coverage=100, p_sequencing_error_range=(0.0, 0.1), p_alignment_error_range=(0.0, 0.1), verbose=1):
        """
        Wrapper to generate n alignments as specified in num_alignments with random error rates.
        Parameters:
        - num_alignments: Number of alignments to be generated.
        - coverage: Number of reads for an alignment.
        - p_sequencing_error_range: Range of sequencing error rates (tuple of two values).
        - p_alignment_error_range: Range of alignment error rates (tuple of two values).
        - verbose: Verbosity level (0 or 1).
        Returns:
        - alignments: Numpy array of alignments generated.
        - prob_lists: List of lists of probability for each alignment.
        """
        alignments = []
        prob_lists = []

        for i in range(num_alignments):
            if (i % 400 == 0) and (verbose == 1):
                print("Progress: {progress_percentage}%% completed. Computing alignment {current_iter} of {total_iter}".format(
                    progress_percentage=round(i*100/num_alignments, 2), current_iter=i, total_iter=num_alignments
                ))
            
            p_sequencing_error = random.uniform(*p_sequencing_error_range)
            p_alignment_error = random.uniform(*p_alignment_error_range)

            alignment, prob_list = self.ratio_gen(coverage, p_sequencing_error, p_alignment_error)
            alignments.append(alignment)
            prob_lists.append(prob_list)

        self.alignments = alignments
        return np.array(alignments), prob_lists
    

    def ratio_gen(self, coverage, p_sequencing_error, p_alignment_error) -> (list,list):
        """Wrapper to generate a single alignment based on a randomly generated ratio
        Returns np.ndarray of the alignment and the probability of the distribution

        NOTE: Naming of variables shall be improved in the future for easy readability
        
        Parameters
        ----------
        coverage : <int>
            Number of read for an alignment
        p_sequencing_error : <double>
            Probability of sequencing error, takes value >= 0, <= 1
        p_alignment_error : <double>
            Probability of alignment error, takes value >= 0, <= 1 
        Returns
        -------
        list
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
        for clone_idx, nb_clone_coverage in enumerate(nb_coverage_list):
            for _ in range(0,nb_clone_coverage):
                coverage_list.append(
                    self._add_errors(self, self.clones[clone_idx],p_sequencing_error,p_alignment_error))
        
        # This will be the final probability list
        prob_list = [nb_coverage_list[i]/coverage for i in range(0, len(nb_coverage_list))]

        # Here we shuffle the list and concatenate into the final alignment
        choice_indices = np.random.choice(len(coverage_list), coverage, replace=False)
        alignment = [self.clones[0]] # First row is always reference (assumed to be index at 0)
        alignment += [coverage_list[i] for i in choice_indices] # Concatenate the randomized read to the reference row

        return alignment, prob_dist

    @staticmethod
    def _add_errors(clone, error_rate):
        """
        Adds sequencing error and alignment error to a single read based on the given error rate.
        Parameters:
        - clone: List of base pairs.
        - error_rate: Rate of introducing errors (float between 0 and 1).
        Returns:
        - clone with errors.
        """
        error_clone = []
        for base in clone:
            if random.random() < error_rate:
                error_clone.append(random.choice("ACGT"))
            else:
                error_clone.append(base)
        return error_clone
    


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

    def _array_binary_mask(self,arr,coverage) -> np.ndarray:
        """Generates binary masks based on whether the base pair is the same as the reference sequence
        (0: False   1: True)
        
        Parameters
        ----------
        arr : <List>
            2D input list of the input array with the first row containing the reference sequence
        coverage : <int>
            Description
        """
        # Here we are transcoding the True to A as it's the highest value
        binary_transdict = {True:"T", False:"A"}
        ref_gen = arr[0]
        #binary_mask_matrix = [np.in1d(arr[i],ref_gen).tolist() for i in range(0,len(arr))]
        binary_mask_matrix = [[arr[j][i]==ref_gen[i] for i in range(0,len(arr[j]))] for j in range(0,len(arr))]
        print(binary_mask_matrix)
        binary_mask_matrix = np.vectorize(binary_transdict.get)(binary_mask_matrix)
        aln_binary_mask_dim = np.array((arr, binary_mask_matrix))
        return aln_binary_mask_dim
