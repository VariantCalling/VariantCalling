# # -*- coding = utf-8 -*-
# # @time:19/06/2023 09:20
# # Author:Yunbo Long
# # @File:visualisation.py
# # @Software:PyCharm

# import matplotlib.pyplot as plt
# import numpy as np

# import numpy as np
# import matplotlib.pyplot as plt

# import numpy as np
# import matplotlib.pyplot as plt
# # try

# NUCLEOTIDES = ['A', 'C', 'G', 'T']

# def plot_alignment(alignment):
#     if isinstance(alignment, list):
#         alignment = np.array(alignment)

#     transdict = {nucleotide: i for i, nucleotide in enumerate(NUCLEOTIDES)}
#     alignment_ints = np.vectorize(transdict.get)(alignment)

#     nucleotide_colors = ['darkgreen', 'darkblue', 'darkred', 'black']
#     alignment_colors = [nucleotide_colors[n] for n in alignment_ints.flatten()]
#     plt.imshow(alignment_ints, cmap='jet')
#     plt.xticks([])
#     plt.yticks([])
#     plt.grid(visible=False)
#     plt.colorbar(ticks=[0, 1, 2, 3], orientation='vertical')
#     plt.imshow(np.zeros_like(alignment_ints), cmap='jet', alpha=0.4)
#     for i in range(alignment_ints.shape[0]):
#         for j in range(alignment_ints.shape[1]):
#             plt.text(j, i, alignment[i][j], ha='center', va='center', color=alignment_colors[i * alignment_ints.shape[1] + j])
#     plt.show()


# # Example usage
# alignment_example = [
#     ['A', 'C', 'C', 'T'],
#     ['A', 'C', 'G', 'T'],
#     ['A', 'C', 'C', 'T'],
#     ['A', 'C', 'G', 'T']
# ]
# plot_alignment(alignment_example)
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
        for i in range(num_alignments):
            if (i % 400 == 0) and (verbose==1):
                print("Progress:  {progress_percentage}%% completed. \tComputing alignment {current_iter} of {total_iter}".format(progress_percentage=round(i*100/num_alignments,2), current_iter = i, total_iter=num_alignments))
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
        # match mode:
        #     case 1: 
        #         return (np.random.dirichlet(np.ones(nb_class)*1000.,size=1)).flatten().tolist()
        #     case 2:
        #         prob_list = [random.random() for _ in range(0,nb_class)]
        #         return [prob_list[i] / sum(prob_list) for i in range(0, nb_class)]

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


# sfsd

import torch
# import torch.nn as nn
# import torch.nn.functional as F
# import VariantCalling as vc
import numpy as np
import matplotlib.pyplot as plt

"""
This section to test ratio_gen, which is wrapper for alignments
"""
dg = VariantCallingData()
alignments, prob_lists = dg.simulate_clones(10,100,0.01,0.01)
alignments_int = dg.char_to_int(alignments)
alignment_int = dg.char_to_int(alignments[0])
plt.rcParams['figure.dpi'] = 200
#alignment_idx = mutation_types.index(mutation_index)
plt.title("Mutation")
plt.imshow(alignment_int,cmap='jet')
print(prob_lists[0])


"""
This section to test ratio_gen, which is wrapper for alignments
"""
dg = vc.VariantCallingData()
alignments, prob_lists = dg.simulate_clones(10,100,0.01,0.01)

orgn = alignments[0]
print(len(orgn))

print(len(alignments[0]))
alignments = [dg._array_dup(i,101) for i in alignments]
print(len(alignments))

print(alignments[8].shape)
print(alignments[8])
alignments_int = dg.char_to_int(alignments)

rng = np.random.default_rng(seed=42) # use a fixed random generator so runs are consistent
idxs = np.arange(alignments_int.shape[0])

rng.shuffle(idxs)

split_idx = int(alignments_int.shape[0]*0.8)
train_alignments, valid_alignments = alignments_int[idxs[:split_idx]], alignments_int[idxs[split_idx:]]
train_prob_lists, valid_prob_lists = np.array(prob_lists)[idxs[:split_idx]], np.array(prob_lists)[idxs[split_idx:]]
#train_mutation_types, valid_mutation_types = mutation_types[idxs[:split_idx]], mutation_types[idxs[split_idx:]]
print(train_alignments.shape)
print(train_prob_lists.shape)



class CNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = nn.Conv2d(2, 6, 5) # input channels, output channels, kernel size
        self.conv2 = nn.Conv2d(6, 16, 5) # input channels, output channels, kernel size
        self.pool = nn.MaxPool2d(4, 4)
        self.fc1 = nn.Linear(720, 120)
        self.fc2 = nn.Linear(120, 60)
        self.fc3 = nn.Linear(60, 2)  # Update the output size to 2 for two classes


    def forward(self, x):
        #print(x.shape)
        x = self.pool(F.relu(self.conv1(x)))
        #print(x.shape)
        x = self.pool(F.relu(self.conv2(x)))
        #print(x.shape)
        x = torch.flatten(x, 1) # flatten all dimensions except batch
        #print(x.shape)
        x = F.relu(self.fc1(x))
        #print(x.shape)
        x = F.relu(self.fc2(x))
        #print(x.shape)
        x = self.fc3(x)
        #print(x.shape)
        return x


def train(model, train_alignments, train_mutation_types, valid_alignments, valid_mutation_types, epochs=10, lr=0.001):
    crit = torch.nn.CrossEntropyLoss()
    opt = torch.optim.Adam(model.parameters(), lr=lr)

    train_dataset = torch.utils.data.TensorDataset(torch.from_numpy(train_alignments).float(), torch.tensor(train_mutation_types))
    valid_dataset = torch.utils.data.TensorDataset(torch.from_numpy(valid_alignments).float(), torch.tensor(valid_mutation_types))

    train_loader = torch.utils.data.DataLoader(dataset=train_dataset, batch_size=16, shuffle=True)
    valid_loader = torch.utils.data.DataLoader(dataset=valid_dataset, batch_size=16)

    train_losses, valid_losses, valid_accs = [], [], []

    for epoch in range(1, epochs + 1):
        # train for 1 epoch
        model = model.train()
        epoch_loss, total = 0.0, 0
        for i, (batch_alignment, batch_mutation_type) in enumerate(train_loader):
            opt.zero_grad()
            out = model(batch_alignment)
            loss = crit(out, batch_mutation_type)
            loss.backward()
            opt.step()
            epoch_loss += loss.item()
            total += 1
        epoch_loss /= total

        # compute validation loss and accuracy
        model = model.eval()
        valid_loss, n_correct, n, total = 0.0, 0, 0, 0
        for i, (batch_alignment, batch_mutation_type) in enumerate(valid_loader):
            with torch.no_grad():
                out = model(batch_alignment)
                loss = crit(out, batch_mutation_type)

            valid_loss += loss.item()
            total += 1

            predict = torch.nn.functional.softmax(out, dim=1).argmax(dim=1)
            correct = predict == batch_mutation_type
            n += out.shape[0]
            n_correct += correct.sum()
        valid_loss /= total
        accuracy = n_correct / n

        train_losses.append(epoch_loss)
        valid_losses.append(valid_loss)
        valid_accs.append(accuracy)
        print(f"epoch={epoch:2d}, train_loss={epoch_loss:.3f}, valid_loss={valid_loss:.3f}, accuracy={accuracy*100:.2f}%")

    return train_losses, valid_losses, valid_accs


model = CNN()
n_epochs = 15 # number of epochs
lr = 0.001 # learning rate

print(train_alignments.shape)

train_losses, valid_losses, valid_accs = train(model, train_alignments, train_prob_lists, valid_alignments, valid_prob_lists, epochs=n_epochs, lr=lr)
