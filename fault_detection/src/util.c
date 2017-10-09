
#include <math.h>
#include <malloc.h>
#include <stdlib.h>

#include <config.h>
#include <util.h>

// return 2^n en int
int pow2(int n)
{
  int i, p = 1;

  for (i = 0; i<n; i++)
    p *= 2;
  return p;
}

int pow4(int n)
{
  int i, p = 1;

  for (i = 0; i<n; i++)
    p *= 4;
  return p;
}

int pow8(int n)
{
  int i, p = 1;

  for (i = 0; i<n; i++)
    p *= 8;
  return p;
}

long double_equal(double v1, double v2)
{
  if (fabs(v1 - v2) < EPS_EQUAL_DOUBLE) return 1;
  return 0;
}


// TODO
// refaire avec les offset
float **malloc_array_float(int nbre, int *size)
{
  float **ptr = NULL;
  long n;

  ptr = (float**)calloc(nbre, sizeof(float*));
  for (n = 0; n < nbre; n++)
    ptr[n] = (float*)calloc(size[n], sizeof(float));
   return ptr;
}

void **free_array_float(void **ptr, int nbre)
{
  long n;

  if (ptr == NULL) return NULL;
  for (n = 0; n<nbre; n++)
    FREE(ptr[n])
   FREE(ptr)
   return NULL;
}

