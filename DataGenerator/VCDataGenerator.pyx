DEF DEFAULT_WIDTH=178
DEF DEFAULT_HEIGHT=100
DEF NUMBER_OF_MUTANTS=4
DEF NUMBER_OF_MUTATIONS_PER_CLONE=3

from libc.stdlib cimport rand, srand, RAND_MAX
from libc.time cimport time
from libc.math cimport sqrt, log, cos, pi, round

srand(time(NULL))

cdef int[4] NUCLEOTIDE_INDEX = [0, 1, 2, 3]
cdef char[4] NUCLEOTIDE = ['A', 'C', 'G', 'T']  # This also give INDEX_2_NUCLEOTIDE
cdef dict NUCLEOTIDE_2_INDEX = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}

ctypedef int MutantList[NUMBER_OF_MUTANTS][DEFAULT_WIDTH]


cdef int _get_random_int(int i):
    """
    Get a random number from 0 to i-1. eg. to draw clones 0 to 3 i will be 4
    - Tested on: Mac and Window
    """
    cdef int random_number = rand() % i
    return random_number


cdef double _draw_from_uniform():
    """
    Get a random number in the range 0 and 1 from a uniform distribution
    - Tested on: Mac and Window
    """
    cdef double rand_max = (RAND_MAX + 1)
    cdef double random_draw = rand()
    cdef double random_uniform = abs(random_draw / rand_max)  # getting this weird bug where on mac the numbers are negative
    return random_uniform


cdef double _draw_from_normal(double mean, double std_deviation):
    """
    Get a random number in the range 0 and 1 from a uniform distribution
    - Tested on: Mac
    """
    cdef double u1 = _draw_from_uniform()
    cdef double u2 = _draw_from_uniform()

    z = sqrt(-2 * log(u1)) * cos(2 * pi * u2)
    return std_deviation * z + mean


cdef int _clip_to_int(double val, int max_val, int min_val):
    """
    Clip a value to a value within max_val and min_val
    - Tested on: Mac
    """
    cdef int val_int = int(round(val))
    
    if val_int > max_val:
        val_int = max_val
    if val_int < min_val:
        val_int = min_val
        
    return val_int


cdef _create_mutant(list sequence, int mutation_pos_main):  # using a list, dont know if this effects runtime.
    """
    Given a sequence and mutation position, return a similar sequence with some mutation
    - Tested on: Mac
    """
    for _ in range(NUMBER_OF_MUTATIONS_PER_CLONE):
        mutation_pos = _clip_to_int(mutation_pos_main + _draw_from_normal(0, 0.5) * 10, DEFAULT_WIDTH - 1, 0)
        sequence[mutation_pos] = _get_random_int(4)
    return sequence


cdef _generate_random_sequence():
    """
    Create a random sequence
    - Tested on: Mac
    """
    cdef int[DEFAULT_WIDTH] sequence
    for i in range(DEFAULT_WIDTH):
        sequence[i] = _get_random_int(4)
    return sequence


cdef _generate_random_clones():
    """Function which returns random sequence with its clones"""
    cdef MutantList clones
    cdef int[DEFAULT_WIDTH] random_seq = _generate_random_sequence()
    clones[0] = random_seq

    cdef double mean_val = (DEFAULT_WIDTH - 1) / 2
    cdef double std_deviation = mean_val / 3.5

    cdef double mutation_pos_double = _draw_from_normal(mean_val, std_deviation)
    cdef int mutation_pos_int = _clip_to_int(mutation_pos_double, DEFAULT_WIDTH - 1, 0)
    print(mutation_pos_int)

    for i in range(NUMBER_OF_MUTANTS - 1):
        new_variant = _create_mutant(random_seq, mutation_pos_int)
        clones[i + 1] = new_variant

    return clones


cdef _simulate_read():
    return


def simulate_read():
    return _simulate_read()


def generate_random_clones():
    return _generate_random_clones()


def generate_data_for_noise_reduction():
    return _generate_random_clones()


def test_function_output():
    cdef int random_int = _get_random_int(4)
    cdef double random_uniform = _draw_from_uniform()
    cdef double random_normal = _draw_from_normal(0, 1)

    return random_int, random_uniform, random_normal