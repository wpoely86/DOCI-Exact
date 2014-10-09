#include <assert.h>

#include "Permutation.h"

/**
 * Constructor
 * @param n the number of bit that needs to be set
 */
Permutation::Permutation(unsigned int n)
{
   assert(sizeof(unsigned long)*8 >= n && "Trouble: unsigned long not big enough");

   this->n = n;

   // set n lowest bits to 1
   reset();
}

/**
 * Permutate to the next permutation
 * and return it
 * @return the next permutation
 */
mybitset Permutation::next()
{
   // only efficient way to do this for the moment:
   // convert to unsigned long and convert back in the end
 
   auto v = current.to_ulong(); // current permutation of bits 

   // from https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
   auto t = v | (v - 1); // t gets v's least significant 0 bits set to 1
   // Next set to 1 the most significant bit to change, 
   // set to 0 the least significant ones, and add the necessary 1 bits.
   auto w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1));

   // new/next permutation of bits
   current = mybitset(w);

   return current;
}

/**
 * Get current permutation
 * @return current permutation
 */
mybitset Permutation::get() const
{
   return current;
}

/**
 * Back to the begin position:
 * the lowest n bits are set
 */
void Permutation::reset()
{
   current = mybitset((1L<<n)-1L);
}

/* vim: set ts=3 sw=3 expandtab :*/
