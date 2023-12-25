from libc.stdlib cimport rand, srand, RAND_MAX
from libc.time cimport time

srand(time(NULL))

cdef int[4] NUCLEOTIDE_INDEX = [0, 1, 2, 3]
cdef char[4] NUCLEOTIDE = ['A', 'C', 'G', 'T']  # This also give INDEX_2_NUCLEOTIDE
cdef dict NUCLEOTIDE_2_INDEX = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}


def get_random_int(int i):
    """Get a random number from 0 to i-1. eg. to draw clones 0 to 3 i will be 4"""
    cdef int random_number = rand() % i
    return random_number


def draw_from_uniform():
    """Get a random number in the range 0 and 1 from a uniform distribution"""
    cdef double rand_max = (RAND_MAX + 1)
    cdef double random_draw = rand()
    cdef double random_uniform = abs(random_draw / rand_max)  # getting this weird bug where on mac the numbers are negative
    return random_uniform

