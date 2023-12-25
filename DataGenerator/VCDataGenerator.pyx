# def primes(int nb_primes):
#     cdef int n, i, len_p
#     cdef int[1000] p

#     if nb_primes > 1000:
#         nb_primes = 1000

#     len_p = 0  # The current number of elements in p.
#     n = 2
#     while len_p < nb_primes:
#         # Is n prime?
#         for i in p[:len_p]:
#             if n % i == 0:
#                 break

#         # If no break occurred in the loop, we have a prime.
#         else:
#             p[len_p] = n
#             len_p += 1
#         n += 1

#     # Let's copy the result into a Python list:
#     result_as_list = [prime for prime in p[:len_p]]
#     return result_as_list

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
    cdef double random_uniform = rand() / <double>(RAND_MAX + 1)
    return random_uniform

