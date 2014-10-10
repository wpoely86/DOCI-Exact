#ifndef PERMUTATION_H
#define PERMUTATION_H

#include<bitset>

// it not set, set it to 64 bits
#ifndef BITSETSIZE
#define BITSETSIZE 64
#endif

//! we use a std::bitset as underlying representation
typedef std::bitset<BITSETSIZE> mybitset;

/**
 * This class is used to generate all permutations of
 * bitsets with n bits set. These permutations are not
 * stored but generated on the fly.
 * Currently, should work well with up to 32 bits 
 * (64 if you switch to unsigned long long). Beyond that,
 * troubles are waiting.
 *
 * There is also no protect against overflows for the moment.
 */
class Permutation
{
    public:
        Permutation(unsigned int);

        virtual mybitset next();

        virtual mybitset get() const;

        virtual void reset();

        static unsigned long long CalcCombinations(unsigned int, unsigned int);

        static unsigned long long gcd(unsigned long long, unsigned long long);

    private:

        //! the current bitset
        mybitset current;

        //! number of ones needed
        unsigned int n;
};

#endif /* PERMUTATION_H */

/* vim: set ts=3 sw=3 expandtab :*/
