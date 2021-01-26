#define NUMBER_OF_CONTRAINTS 2
#define TOLERANCE 1e-6

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

//linked list with two element data
struct Lines{
	float alpha;
	float beta;
	struct Lines* next;
	struct Lines* prev;
};

// Linked list data structure 
struct Node { 
    float data; 
    struct Node* next; 
    struct Node* prev;
}; 

//return variables for two variable function
struct tv_return{
	struct Lines* line;
	float u1;
	float u2;
};


struct Lines* deleteLine(struct Lines*, struct Lines*);

//deleting a node
void deleteNode(struct Node*, struct Node*);

//swap two indexes in array
void swap(float *,float *);

//sort array
void sort(float[],int);
 
// Standard partition process of QuickSort(). 
// It considers the last element as pivot 
// and moves all smaller element to left of 
// it and greater elements to right
//from geeks for geeks
int partition(float[], int, int );

// This function returns k'th smallest  
// element in arr[l..r] using QuickSort  
// based method.  ASSUMPTION: ALL ELEMENTS 
// IN ARR[] ARE DISTINCT 
float kthSmallest(float[], int, int, int);

float quickSelect(float[], int, int,  float (*pivot_fn)(float [],int ));

float nlogn_median(float[],int);


float quickSelect_median(float[],int, float (*pivot_fn)(float [],int ));

float pickPivot(float [],int);

//size of alpha and beta must be equal
float* pairing(struct Lines*,float[],float, float);

struct tv_return* two_variable(float, float ,float [] , float[] , float[] , int);


struct Lines* make_lines(float [], float [],int);
