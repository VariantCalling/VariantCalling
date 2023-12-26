DEF DEFAULT_WIDTH=178
DEF DEFAULT_HEIGHT=100
DEF NUMBER_OF_MUTANTS=4
DEF NUMBER_OF_MUTATIONS_PER_CLONE=3

from libc.stdlib cimport rand, srand, RAND_MAX
from libc.time cimport time
from libc.math cimport sqrt, log, cos, pi, round
cimport cython

srand(time(NULL))

cdef int[4] NUCLEOTIDE_INDEX = [0, 1, 2, 3]
cdef char[4] NUCLEOTIDE = ['A', 'C', 'G', 'T']  # This also give INDEX_2_NUCLEOTIDE
cdef dict NUCLEOTIDE_2_INDEX = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}

ctypedef int MutantList[NUMBER_OF_MUTANTS][DEFAULT_WIDTH]

cdef int get_random_int(int i):
    """
    Get a random number from 0 to i-1. eg. to draw clones 0 to 3 i will be 4
    - Tested on: Mac and Window
    """
    cdef int random_number = rand() % i
    return random_number

def get_random_int_python(int i):
    return get_random_int(i)


cdef double draw_from_uniform():
    """
    Get a random number in the range 0 and 1 from a uniform distribution
    - Tested on: Mac and Window
    """
    cdef double rand_max = (RAND_MAX + 1)
    cdef double random_draw = rand()
    cdef double random_uniform = abs(random_draw / rand_max)  # getting this weird bug where on mac the numbers are negative
    return random_uniform


cdef double draw_from_normal(double mean, double std_deviation):
    """
    Get a random number in the range 0 and 1 from a uniform distribution
    - Tested on: Mac
    """
    cdef double u1 = draw_from_uniform()
    cdef double u2 = draw_from_uniform()

    z = sqrt(-2 * log(u1)) * cos(2 * pi * u2)
    return std_deviation * z + mean


cdef int clip_to_int(double val, int max_val, int min_val):
    
    cdef int val_int = int(round(val))
    
    if val_int > max_val:
        val_int = max_val
    if val_int < min_val:
        val_int = min_val
        
    return val_int


cdef generate_random_sequence():
    cdef int[DEFAULT_WIDTH] sequence
    for i in range(DEFAULT_WIDTH):
        sequence[i] = get_random_int(4)
    return sequence


def char_to_int():
    return


def int_to_char():
    return


def generate_data_for_noise_reduction():
    return


def generate_random_sequence_api():
    return generate_random_sequence()


def generate_random_clones():
    """Function which returns random sequence with its clones"""
    cdef MutantList clones
    cdef int[DEFAULT_WIDTH] random_seq = generate_random_sequence()
    clones[0] = random_seq

    cdef double mean_val = (DEFAULT_WIDTH - 1) / 2
    cdef double std_deviation = mean_val / 3.5

    cdef double mutation_pos_double = draw_from_normal(mean_val, std_deviation)
    cdef int mutation_pos_int = clip_to_int(mutation_pos_double, DEFAULT_WIDTH - 1, 0)
    print(mutation_pos_int)

    for i in range(NUMBER_OF_MUTANTS - 1):
        new_variant = random_seq.copy()  # should copy the array
        new_variant[0] += 1
        clones[i + 1] = new_variant

    return clones


def simulate_read():
    return


def test_function_output():
    cdef int random_int = get_random_int(4)
    cdef double random_uniform = draw_from_uniform()
    cdef double random_normal = draw_from_normal(0, 1)

    return random_int, random_uniform, random_normal