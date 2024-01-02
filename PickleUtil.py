"""
PickleUtil.py

This file contains the utilities for the processing of the pickles files generated
using Minimap and 
"""

import numpy as np

class PickleLoader:
    def __init__(self, sequence="CRT",clones=["7G8"]) -> None:
        self.supported_sequence = ("CRT","DHPS","DHFR")
        if sequence not in self.supported_sequence:
            raise ValueError("Currently supported are: " + str(self.supported_sequence))
        self.sequence = sequence
        self.clones = clones
        self.pickle_path = r"pickle/" + sequence

    def load_pickle(self):
        """
        This function serves to load the pickle as defined during the initialisation.
        The pickle files are loaded into a 2-D array in the order of [clones]

        For example, if clones=["7G8","DD2"], the clones_sequence will be of size [2,],
        where clones_sequence[0,] will be the silo for 7G8 reads
        while clones_sequence[1,] will be the silo for DD2 reads.

        Returns
        -------
        list
            List of reads loaded from the Pickle files arranged in the order of [clones]
        """
        clones_sequence = []
        for clone in self.clones:
            clones_sequence.append(np.array(np.load(self.pickle_path + "/" + self.map_pickle_file(clone),allow_pickle=True))[:,2].tolist())
        return clones_sequence

    def map_pickle_file(self,clone):
        """
        Map the input clone to the pickle files saved.
        
        Args:
            clone (TYPE): The clone for which the pickle file is to be returned
        """

        # Note: In the future we can modify this part to allow for custom pickle file names
        return "trimmed_data_" + str(clone) + ".pkl"
