/* definition of major structure for PCx()
 *
 * PCx beta-2.0  10/31/96.
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Steve Wright.
 *
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include "hash.h"
#include <ctype.h>

#ifdef MEX
#define printf          mexPrintf
#endif

#define PINFTY       0     /* normal */
#define NORMAL       0
#define FREE         1
#define UPPER        2
#define LOWER        3
#define UPPERLOWER   4
#define FIX          5
#define MINFTY       6

#define OPTIMAL_SOL      0
/* The value of solution status below are greater than 10 to avoid
   confuse with the main program return code. */
#define SUBOPTIMAL_SOL         11
#define INFEASIBLE_SOL         12
#define UNKNOWN_SOL            13
#define UNKNOWN_SOL_INF_MU     14
#define UNKNOWN_SOL_PHI_SLOW   15

#define ON               1
#define OFF              0

#define YES              1
#define NO               0

/* used to flag error in the input; see source file readmps.c */
#define BAD_INPUT       -10

/* main program returns an integer if error condition is noted. The
   error codes are: */

#define INVOCATION_ERROR     1
#define MEMORY_ERROR         2
#define INPUT_ERROR          3
#define SPECS_ERROR          4
#define PRESOLVE_ERROR       5
#define FACTORIZE_ERROR      6

#define ABS(a)    ((a)>0 ? (a) : -(a))
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <time.h> 

//defines for constants
#define NUMBER_OF_CONTRAINTS 2
#define TOLERANCE_1 1e-10
#define TOLERANCE_2 1

//defines for random problem generator 
#define STARTING_SEED 10000
#define INCREMENTAL_SEED 2000
#define NUMBER_PROBLEMS 10
#define A_MULTIPLE_FACTOR 1000
#define A_SUM_FACTOR 0
#define B_TIGHTNESS_RATIO 0.5
#define C_MULTIPLE_FACTOR 200
#define NEG_NUMBERS 1   // 0 = No | 1 = Yes
#define PROB_NEGATIVES 0.2
#define ZERO_NUMBERS 1   // 0 = No | 1 = Yes
#define PROB_ZERO 0.2
#define SPARSE 0        // 0 = No Sparse Matrix   | 1 = Sparse Matrix
#define SPARSE_PERCENTAGE 0.90

//includesfordyer_algo
#define DPDUBUG 0 //will print values to debug

//define data types
typedef double REAL;

//variables for random problem generator
static double z;


//linked list with two element data
struct Lines{
  double alpha;
  double beta;
  struct Lines* next;
  struct Lines* prev;
};

// Linked list data structure 
struct Node { 
    double data; 
    struct Node* next; 
    struct Node* prev;
}; 

//linked list of pairs of lines
struct Pairs{
  struct Lines* first;
  struct Lines* second;
  double intersection;
  struct Pairs* next;
  struct Pairs* prev;
};

//return variables for two variable function
struct tv_return{
  struct Lines* line;
  double u1;
  double u2;
};

 struct _twoDline{
  double  alpha;
  double  beta;
  double  a1, a2, b;
  int   orig;
  struct _twoDline *next;
};

struct Rows_Double {
    double *cols;
};
 
struct Rows_Integer {
    int *cols;
};
 
struct Instance {
    double *cost_vector;
    double *rhs_vector;
    struct Rows_Double *a_matrix;
    int number_variables;
    int number_constraints;
};
 
typedef struct Instance INSTANCE, *INSTANCEptr;


/****************************************************************/

typedef struct {   /* sparse structure, with pointer information for
                      both the matrix and its transpose */
    int NumRows, NumCols, Nonzeros;
    int *pBeginRow, *pBeginRowT; /* Fortran style indexing is used here! */
    int *pEndRow, *pEndRowT;   /* indices range from 1 to n, NOT 0 to n-1 */
    int *Row, *RowT;           
    double *Value, *ValueT;
} MMTtype;

/****************************************************************/

typedef struct {   /* structure of constraint matrix in Mehrotra format */
    int     NumRows, NumCols, Nonzeros;
    int    *pBeginRow;
    int    *pEndRow;
    int    *Row;
    double *Value;
} sparseMatrix;

/****************************************************************/

typedef struct {
   void     *ptr;               /* solvers can attach data here */
   char     *FactorizationCode;
   int       N, NonzerosL;
   int       SmallDiagonals;
   int      *Perm, *InvPerm;
   MMTtype  *AAT;
   int       Ndense;		/* for dense column handling */
   int      *maskDense;
   double  **W, **Ldense;
} FactorType;
   
/****************************************************************/

typedef struct  {

  int  NumRows;  /* indicates number of entries in data structures */
  int  NumCols;
  int  NumEnts;

  int  RowSize;  /* indicates size of data structures allocated */
  int  ColSize;
  int  EntSize;

  sparseMatrix A;

  double *b, *c;
  double cshift;
  int    *BoundType;
  double *UpBound;
  double *LowBound;
  double *Ranges;
  char   *RowType;

  char   *ProblemName;
  char   *ObjectiveName;
  char   *RHSName;
  char   *RangeName;
  char   *BoundName;

  char  **RowNames;
  char  **ColNames;

  HashTable *RowTable;
  HashTable *ColTable;

} MPStype;

/****************************************************************/

typedef struct  {

  int  Rows;  /* indicates number of entries in data structures */
  int  Cols;

  int  *SignChanges; /* stores 0 or 1, depending on whether the
                        sign of the i-th variable is flipped or not */
  double *VarShifts; /* records shifting of each component. NB sign
                        change is applied AFTER the shift */
  double cshift; /* shift in objective. This value must be ADDED to the
	primal objective the dual objective in the transformed problem 
	to get the original cost */

} MPSchanges;

/****************************************************************/

typedef struct LPtype {

  int  Rows;         /* Number of rows in A matrix */
  int  Cols;         /* Number of columns in A matrix */
  int  Ents;         /* Number of nonzero entries in A matrix */

  sparseMatrix Atranspose;
  sparseMatrix A;

  double  *b;        /* right-hand side vector */
  double  *c;        /* cost coefficient vector */
  double cshift;     /* constant shift for the cost; the objective
                        function is cshift + c.x */

  int    *VarType;   /* For each variable specify type: 
			Upper :  0 <= x <= Upbound,
			Normal:  0 <= x, 
			Free */
  double *UpBound;   /* If upper bound, specify */

  int    NumberBounds;
  int   *BoundIndex; /* List of variables which have upper bounds */

  int      NumScale;
  double  *ColScale; /* Vector of column scalings */
  double  *RowScale; /* Vector of row scalings */

  int    NumberFree;
  int   *FreeIndex;

  int    NumberSplit;
  int   *FreePlus;   /* List of variables of type free */
  int   *FreeMinus;
} LPtype;

/****************************************************************/

typedef struct Parameters {
  
  int     IterationLimit;
  double  OptTol;         /* absolute tolerance */
  double  PriFeasTol, DualFeasTol; /* tolerances for primal and 
                                      dual infeasibilities */

  double  AlphaScale;     /* If AlphaMax is the largest value of Alpha in [0,1]
                             for which (x + Alpha Dx, s + Alpha ds) >= 0,
                             take a step of length AlphaScale * AlphaMax.
                             (Typical values are .99 or .999) */
  int Diagnostics;        /* More detailed messages with higher number */
  int ReportingLevel;     /* More detailed messages with higher number */
  int Refinement;         /* 1 = perform iterative refinement; 0 = don't */
  int Preprocessing;      /* 1 = perform preprocessing; 0 = don't */
  int Scaling;            /* 1 = perform scaling; 0 = don't */
  int HOCorrections;      /* 1 = perform Gondzio corrections; 0 = don't */
  int MaxCorrections;     /* Maximum number of Gondzio corrections; 
                             0 = computed by the code */
  int Minimize;           /* 1 = Minimize the objective; 0= Maximize */
  char *InputDirectory;   /* look for input file in this directory */
  int WriteSolution;      /* 1 = Write solution to a file; 0 = don't */
  int WriteHistory;       /* 1 = Write history to a file; 0 = don't  */

  char *RHSName;          /* If names are given for these values, only  */
  char *ObjectiveName;    /* accept rhs, range, bound, or objectives    */
  char *RangeName;        /* by the given name. Otherwise, accept the   */
  char *BoundName;        /* first name encountered.                    */

  int   CacheSize;        /* CacheSize for Ng-Peyton routines           */
  int   UnrollingLevel;   /* Unrolling level for Ng-Peyton routines     */

  double CenterExponent;  /* Exponent for computing centering parameter */

  int OrderAlg;           /* Ordering algorithm to use. 0 = dont't use any
                             algorithm, 1 = multiple minimum degree (default). */
} Parameters;


/****************************************************************/


typedef struct {
  double PrimalObjective;
  double DualObjective;
  double PriInf;
  double DualInf;
  double logmu;
  int    NumCorrections;
  double phi;
} IterationRecord;

typedef struct {
  int    Nonzeros;
  double Density;
  int    Operations;
  int    NumDenseCols;
} FactorizationRecord;

typedef struct solution {
  
  int Rows;          /* number of rows in A */
  int Columns;       /* number of variables */

  double *x;         /* primal solution */
  double *pi;        /* dual vector corresponding to equality constraints */

  double *DualUpper; /* Dual variables for upper bound constraints */
  double *DualLower; /* Dual variables for lower bound constraints */
                     /* If x_j is fixed, the dual for x_j is stored
                        in both DualUpper[j] and DualLower[j].
                        Components of DualUpper and DualLower are left
                        blank if the variable does not have an upper
                        or lower bound, resp. */

  /* s=DualLower    r=DualUpper       x+w=UpBound      */

  double *Activity;  /* row activity  Ax */

  double PrimalObjective;
  double DualObjective;		
  double Complementarity, RelativeComplementarity;

  double PrimalInfeasibility;
  double DualInfeasibility;
  
  int  Iterations;
  IterationRecord *IterationHistory;

  char *FactorizationCode;
  int   Factorizations;
  FactorizationRecord *FactorizationHistory;
  int   RestoredIteration;

  double  ReadTime;
  double  PreprocessTime;
  double  SolutionTime;
  double  FactorizationTime;
  double  SolveADATTime;
  double  InitTime, LoopTime, PredictorTime, FormADATtime;
  double  CorrectorTime;

  int     Status;      /* Optimal, Suboptimal, Infeasible, Unbounded, etc. */

  /* stuff transferred from ChangeStack data structure */
  int PriorRows;      /* rows in the matrix */
  int PriorColumns;   /* cols in the matrix */
  int Passes;         /* number of passes */
  int ReducedRows;    /* Rows remaining in the matrix after reduction */
  int ReducedColumns; /* Columns remaining in the matrix after reduction */

} solution;

typedef struct {
  int     NumRows, NumCols, NumBounds;

  double *x;       /*                   dimension = Cols */
  double *s;       /* DualLower         dimension = Cols */
  double *pi;      /*                   dimension = Rows */
  double *w;       /* Bound slacks      dimension = Bounds */
  double *r;       /* DualUpper         dimension = Bounds */

  int    *BoundIndex;

  double  PriInf, DualInf;
  double  Mu;
  int     Iteration;
} Iterate;

typedef struct 
{
  FILE  *warningfile;
  FILE  *errorfile;
  char  *libraryused;
  char  *programminglanguage;
  char  *machinetype;
  char  *machinename;
  char  *nativelanguage;
} GeneralInfo;

extern void SetGeneralInfo();

struct Lines* deleteLine(struct Lines*, struct Lines*);

//deleting a node
void deleteNode(struct Node*, struct Node*);

//swap two indexes in array
void swap(double *,double *);

//sort array
void sort(double[],int);
 
// Standard partition process of QuickSort(). 
// It considers the last element as pivot 
// and moves all smaller element to left of 
// it and greater elements to right
//from geeks for geeks
int partition(double[], int, int );

// This function returns k'th smallest  
// element in arr[l..r] using QuickSort  
// based method.  ASSUMPTION: ALL ELEMENTS 
// IN ARR[] ARE DISTINCT 
double kthSmallest(double[], int, int, int);

double quickSelect(double[], int, int,  double (*pivot_fn)(double [],int ));

double nlogn_median(double[],int);


double quickSelect_median(double[],int, double (*pivot_fn)(double [],int ));

double pickPivot(double [],int);

//size of alpha and beta must be equal
double* pairing(struct Lines*,double [],double , double,double []);

struct tv_return* two_variable(double, double ,double [] , double[] , double[] , int);


struct Lines* make_lines(double [], double [],int);

double function_random();

void create_instance(struct Instance *instance);

void create_sparse(struct Instance *instance);

void make_random_problem(int,int, int);
void solve_random_problem(int,int,int);
void solve_problem(char *,int,int );

double getLeftSlope(struct _twoDline*, double , double *);
double getRightSlope(struct _twoDline *, double , double *);
void GetOptPoint(double *, double *, double *, double, double, int, double *, double *);
LPtype* makePCXProblem(double, double ,double [] , double[] , double[] , int);

