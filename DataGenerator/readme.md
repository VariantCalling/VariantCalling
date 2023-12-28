# Variant Calling Data Generator

This folder contains all the functions necessary to generate simulated oxford nanopore reads. There are two main ways to do this
1. Using the VariantCallingMini class
2. Using the cython implementation
Eventhough, the VariantCallingMini Class will be the most up-to-date version incase of updates, but both are stable. The cython implementation is harder to use but should give the user a boost in runtime (15 times faster on the MacBook M1 Pro and google colab!).

## VariantCallingMini

```python

import VariantCallingMini as vc

dg = vc.VariantCallingDataMini()
noisy, clean, mutation_positions = dg.generate_data_for_noise_reduction(1000, 100, 0.01, 0.05)
```

- `noisy` will contain an array of noisy images (ie. a 3D numpy array)
- `clean` will contain an array of clean images with each the same index as noisy
- `mutation_poistions` contains the array which tells you where the mutations are


## Cython Implementation

As of now, there are a few constants you would need to adjust in VCDataGenerator.pyx file which corresponds to 

1. `DEFAULT_WIDTH` : the width of the image
2. `DEFAULT_HEIGHT` : the height of the image
3. `NUMBER_OF_MUTANTS` : the number of mutants to generate for a given unique sequence
4. `NUMBER_OF_MUTATIONS_PER_CLONE` : gives the number of mutations in a clone

As cython is a compiled language, the code in VCDataGenerator will have to first be compiled. To do this run the following.

```bash
python setup.py build_ext --inplace
```

This will generate a binary file to be run. You should run this everytime when you change your machine.

To run the code, import and use it as a normal python function

```python
from VCDataGenerator import generate_data_for_noise_reduction

noisy, clean = generate_data_for_noise_reduction(1000, alignment_error_prob=0.01, sequencing_error_prob=0.05)
```

Note that the API for the two methods are currently a little different. Also, noisy and clean will be a python list not a numpy array.

Now, if you would like to run this on jupyter notebook, then run the following cells in jupyter notebook

```python
!pip install cython
```

```python
%load_ext cython
```

```python
%%cython

<PASTE THE CONTENT of VCDataGenerator.pyx HERE>
```

Then you could just call the function as above.

## Other functions

Eventhough the `generate_data_for_noise_reduction` function should be enough for any usage for the project, both implementations also have the following functionality which you might find useful

1. `generate_random_clones` function to generate a 2D array containing similar DNA sequences with some mutations
2.  `_simulate_read` function to simulate an oxford nanopore read

Again, the API for the two implementations are different but the source code should be simple enough to trace.