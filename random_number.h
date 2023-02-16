#pragma once
#ifndef _RANDOM_NUMBER_H_
#define _RANDOM_NUMBER_H_



/*---------------------------------------------------------------------------*/
/*!  \brief  initialize random nubmer generator with a seed \a s.
     \param  s  seed for initialization  */
void mtinit(unsigned long s);
/*!  \brief  initialize the random number generator by an array \a init_key with array-length \a key_length.
     \param  init_key  array of seeds for initializing keys
     \param  key_length  length of \a init_key  */
void mtinit_by_array(unsigned long init_key[], unsigned long key_length);
/*!  \brief  generates a random number on [0,0xffffffff]-interval.  */
unsigned long mtrand_int32(void);
/*!  \brief  generates a random number on [0,0x7fffffff]-interval.  */
long mtrand_int31(void);
/*!  \brief  generates a random number on [0,1]-real-interval.  */
double mtrand_real1(void);
/*!  \brief  generates a random number on [0,1)-real-interval.  */
double mtrand_real2(void);
/*!  \brief  generates a random number on (0,1)-real-interval.  */
double mtrand_real3(void);
/*!  \brief  generates a random number on [0,1) with 53-bit resolution.  */
double mtrand_res53(void);
/*---------------------------------------------------------------------------*/
/*!  \brief  returns a normal distribution random number.  */
double nrand(void);

/*---------------------------------------------------------------------------*/

#endif /* _RANDOM_NUMBER_H_ */