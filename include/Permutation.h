#ifndef PERMUTATION_H
#define PERMUTATION_H

// use unsigned long by default
#if not defined(USELONG) && not defined(USELONGLONG)
#define USELONG
#endif

#if defined(USELONG) && defined(USELONGLONG)
#error "You really have to choose between unsigned long and unsigned long long!"
#endif

#if defined(USELONG)
//! we use a unsigned long as underlying representation
typedef unsigned long mybitset;
#elif defined(USELONGLONG)
//! we use a unsigned long long as underlying representation
typedef unsigned long long mybitset;
#endif



/**
 * We use the namespace doci to put all other stuff in.
 * In this way, the code can be easily used from other code without
 * name slashes or anything
 */
namespace doci {

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

        static unsigned int getMax();

    private:

        //! the current bitset
        mybitset current;

        //! number of ones needed
        unsigned int n;
};

}

#endif /* PERMUTATION_H */

/* vim: set ts=3 sw=3 expandtab :*/
