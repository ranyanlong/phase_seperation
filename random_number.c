
#include<math.h>
#include"random_number.h"


#ifndef TRUE
#define TRUE (!0)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL    /* constant vector a */
#define UMASK 0x80000000UL        /* most significant w-r bits */
#define LMASK 0x7fffffffUL        /* least significant r bits */
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))

/*---------------------------------------------------------------------------*/
static unsigned long state[N];    /* the array for the state vector  */
static int left = 1;
static int initf = 0;
static unsigned long* next;

/*---------------------------------------------------------------------------*/
static void next_state(void);

/*---------------------------------------------------------------------------*/
void mtinit(unsigned long s)
{
    int j;

    /* initializes state[N] with a seed */


    state[0] = s & 0xffffffffUL;
    for (j = 1; j < N; j++)
    {
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array state[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        state[j] = (1812433253UL * (state[j - 1] ^ (state[j - 1] >> 30)) + j);

        /* for >32 bit machines */
        state[j] &= 0xffffffffUL;
    }
    left = 1;
    initf = 1;
}

/*---------------------------------------------------------------------------*/
void mtinit_by_array(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;

    /* initialize by an array with array-length */
    /* init_key is the array for initializing keys */
    /* key_length is its length */


    mtinit(19650218UL);
    i = 1;
    j = 0;
    k = (N > key_length ? N : key_length);
    for (; k; k--)
    {
        state[i] = (state[i] ^ ((state[i - 1] ^ (state[i - 1] >> 30)) * 1664525UL)) + init_key[j] + j;    /* non linear */
            /* for WORDSIZE > 32 machines */
        state[i] &= 0xffffffffUL;
        i++;
        j++;
        if (i >= N)
        {
            state[0] = state[N - 1];
            i = 1;
        }
        if ((unsigned long)j >= key_length)
            j = 0;
    }
    for (k = N - 1; k; k--)
    {
        state[i] = (state[i] ^ ((state[i - 1] ^ (state[i - 1] >> 30)) * 1566083941UL)) - i;    /* non linear */
            /* for WORDSIZE > 32 machines */
        state[i] &= 0xffffffffUL;
        i++;
        if (i >= N)
        {
            state[0] = state[N - 1];
            i = 1;
        }
    }

    /* MSB is 1; assuring non-zero initial array */
    state[0] = 0x80000000UL;
    left = 1;
    initf = 1;
}

/*---------------------------------------------------------------------------*/
void next_state(void)
{
    unsigned long* p = state;
    int j;

    /* if mtinit() has not been called, */
    /* a default initial seed is used         */
    if (initf == 0)
        mtinit(5489UL);

    left = N;
    next = state;

    for (j = N - M + 1; --j; p++)
        *p = p[M] ^ TWIST(p[0], p[1]);

    for (j = M; --j; p++)
        *p = p[M - N] ^ TWIST(p[0], p[1]);

    *p = p[M - N] ^ TWIST(p[0], state[0]);
}

/*---------------------------------------------------------------------------*/
unsigned long mtrand_int32(void)
{
    unsigned long y;

    /* generates a random number on [0,0xffffffff]-interval */


    if (--left == 0)
        next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/*---------------------------------------------------------------------------*/
long mtrand_int31(void)
{
    unsigned long y;

    /* generates a random number on [0,0x7fffffff]-interval */


    if (--left == 0)
        next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (long)(y >> 1);
}

/*---------------------------------------------------------------------------*/
double mtrand_real1(void)
{
    unsigned long y;

    /* generates a random number on [0,1]-real-interval */


    if (--left == 0)
        next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0 / 4294967295.0);

    /* divided by 2^32-1 */
}

/*---------------------------------------------------------------------------*/
double mtrand_real2(void)
{
    unsigned long y;

    /* generates a random number on [0,1)-real-interval */


    if (--left == 0)
        next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0 / 4294967296.0);

    /* divided by 2^32 */
}

/*---------------------------------------------------------------------------*/
double mtrand_real3(void)
{
    unsigned long y;

    /* generates a random number on (0,1)-real-interval */


    if (--left == 0)
        next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return ((double)y + 0.5) * (1.0 / 4294967296.0);
    /* divided by 2^32 */
}

/*---------------------------------------------------------------------------*/
double mtrand_res53(void)
{
    unsigned long a = mtrand_int32() >> 5, b = mtrand_int32() >> 6;

    /* generates a random number on [0,1) with 53-bit resolution*/
    /* These real versions are due to Isaku Wada, 2002/01/09 added */


    return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
}
/*---------------------------------------------------------------------------*/
double nrand(void)
{
    static int number_exists = FALSE;
    static double x, y, r;

    /*  Box-Muller method.  */
    if (number_exists)
    {
        number_exists = FALSE;
        return y * r;
    }
    else
    {
        do
        {
            x = 2.0 * mtrand_real2() - 1.0;
            y = 2.0 * mtrand_real2() - 1.0;
            r = x * x + y * y;
        } while ((r > 1.0) || (r == 0.0));

        r = sqrt(-2.0 * log(r) / r);
        number_exists = TRUE;
        return x * r;
    }
}

/*---------------------------------------------------------------------------*/
