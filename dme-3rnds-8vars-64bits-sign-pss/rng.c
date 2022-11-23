#include <stdlib.h>
#include "rng.h"

void randombytes(unsigned char *x, unsigned long long n)
{
    unsigned long long i;
    for (i=0; i<n; i++)
        x[i] = rand() & 0xff;
}
