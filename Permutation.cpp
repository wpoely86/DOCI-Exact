#include <iostream>
#include <limits>
#include <stdexcept>
#include <assert.h>

#include "Permutation.h"

/**
 * Constructor
 * @param n the number of bit that needs to be set
 */
Permutation::Permutation(unsigned int n)
{
   if(sizeof(unsigned long long)*8 < n)
      throw std::overflow_error("Cannot store permutations in unsigned long long");

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
 
   auto v = current.to_ullong(); // current permutation of bits 

   // from https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
   auto t = v | (v - 1); // t gets v's least significant 0 bits set to 1
   // Next set to 1 the most significant bit to change, 
   // set to 0 the least significant ones, and add the necessary 1 bits.
   auto w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzll(v) + 1));

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

/**
 * Calculate the number of combinations to choose N out of L
 * From: https://stackoverflow.com/questions/1838368/calculating-the-amount-of-combinations
 * @param L the total number of sites
 * @param N the number of sites to choose
 * @return the number of possible combinations to choose N out of L
 */
unsigned long long Permutation::CalcCombinations(unsigned int L, unsigned int N)
{
   assert(L >= N);

   if(0 == L)
      return 0;
   if(0 == N)
      return 1;
   if(L == N)
      return 1;
   if(1 == N)
      return N;

   unsigned long long result = 1;

   for(unsigned int i=1;i<=N;++i, --L)
   {
      auto g = gcd(result, i);
      result /= g;
      auto t = L / (i / g);

      if(result > std::numeric_limits<unsigned long long>::max() / t)
         throw std::overflow_error("Overflow in CalcCombinations");

      result *= t;
   }

   return result;
}

/**
 * Find greatest common divider
 * @param x the first number
 * @param y the second number
 * @return the greatest common divider from x and y
 */
unsigned long long Permutation::gcd(unsigned long long x, unsigned long long y)
{
   assert(x >= y);

    while (y != 0)
    {
        unsigned long long t = x % y;
        x = y;
        y = t;
    }

    return x;
}

/* vim: set ts=3 sw=3 expandtab :*/
