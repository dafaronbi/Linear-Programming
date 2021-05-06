
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <time.h> 

//defines for constants
#define NUMBER_OF_CONTRAINTS 2
#define TOLERANCE 1e-6

//defines for random problem generator 
#define STARTING_SEED 10000
#define INCREMENTAL_SEED 2000
#define NUMBER_PROBLEMS 10
#define A_MULTIPLE_FACTOR 1000
#define A_SUM_FACTOR 0
#define B_TIGHTNESS_RATIO 0.5
#define C_MULTIPLE_FACTOR 200
#define NEG_NUMBERS 0   // 0 = No	| 1 = Yes
#define PROB_NEGATIVES 0.1
#define SPARSE 0        // 0 = No Sparse Matrix		| 1 = Sparse Matrix
#define SPARSE_PERCENTAGE 0.90

//includesfordyer_algo
#define DPDUBUG 1	//will print values to debug

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

//return variables for two variable function
struct tv_return{
	struct Lines* line;
	double u1;
	double u2;
};

 struct _twoDline{
	double	alpha;
	double	beta;
	double	a1, a2, b;
	int 	orig;
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

int make_random_problem(int, int);
void solve_random_problem(int,int);

double getLeftSlope(struct _twoDline*, double , double *);
double getRightSlope(struct _twoDline *, double , double *);
void GetOptPoint(double *, double *, double *, double, double, int, double *, double *);
