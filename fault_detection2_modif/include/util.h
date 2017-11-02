#ifndef __UTIL__
#define __UTIL__

#ifndef MAX
	#define MAX(x,y)		( ( x >= y ) ? x : y )
#endif

#if !defined (PI)
	#define PI				 3.14159265358979
#endif

#if !defined (PI_DOUBLE)
	#define PI_DOUBLE  6.28318530717959
#endif

#ifndef CALLOC
#define CALLOC(size, format) (format*)calloc(size, sizeof(format))
#endif

#ifndef FREE
#define FREE(x){ \
if ( x != NULL ) { \
        free(x); \
        x = NULL;} }
#endif

#ifndef CALLOC_TAB
#define CALLOC_TAB(ptr, nb, size, element)	{ \
												long i; \
												ptr = (element **)calloc(nb, sizeof(element *)); \
												for (i=0; i<nb; i++) (ptr)[i] = (element *)calloc(size, sizeof(element)); \
											}
#endif

#ifndef FREE_TAB
#define FREE_TAB(ptr, nb)   { \
								long i; \
								if ( ptr != NULL ) for (i=0; i<nb; i++) FREE((ptr)[i]) \
								FREE(ptr) \
							}
#endif
    

#ifndef CALLOC_VTAB
#define CALLOC_VTAB(ptr, nb, size, element)	{ \
												long i; \
												element * tmp = NULL; \
												tmp = (element *)calloc(nb*size, sizeof(element)); \
												ptr = (element **)calloc(nb, sizeof(element *));\
												for (i=0; i<nb; i++) (ptr)[i] = tmp+i*size; \
											}
#endif

#ifndef FREE_VTAB
#define FREE_VTAB(ptr)   { \
								if ( ptr != NULL ) FREE(ptr[0]);\
								FREE(ptr) \
							}     
#endif

#endif