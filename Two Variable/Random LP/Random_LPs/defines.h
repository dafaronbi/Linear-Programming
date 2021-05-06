
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

// bailing out

#define ABORTMSG(a)	{ printf("%s", a); \
fprintf(stderr, "\n\nExit at line %d in file %s\n", __LINE__, __FILE__); \
fflush(stderr); fflush(stdout); exit(1); }


// memory allocation

#define MYCALLOC(VAR,TYPE,SIZE) \
    if (((VAR) = (TYPE *) calloc((unsigned)(SIZE), (unsigned)sizeof(TYPE))) == NULL) {\
        fprintf (stderr, "MYCALLOC: Memory allocation failed (%s,%d)\n", __FILE__, __LINE__);\
        exit (1);\
    }

#define MYFREE(VAR) free (VAR);


// file I/O

#define OpenFileRead(X, Y) \
    if ((X = fopen (Y, "r" )) == NULL) { \
        fprintf (stderr, "Cannot open file %s for reading\n", Y); \
        exit (1); \
    }

#define CloseFile(X) \
    fclose (X)

#define TRUE 1
#define FALSE 0