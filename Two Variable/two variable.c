#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#define NUMBER_OF_CONTRAINTS 2
#define TOLERANCE 1e-6

struct Lines{
	float alpha;
	float beta;
	struct Lines* next;
	struct Lines* prev;
};

struct Lines* deleteLine(struct Lines* head_ref, struct Lines* del)
{
    /* base case */
    if (head_ref == NULL || del == NULL)
        return del;
 
    /* If node to be deleted is head node */
    if (head_ref == del){

        head_ref = del->next;
        head_ref->prev = NULL;
    }
 
    /* Change next only if node to be deleted is NOT the last node */
    if (del->next != NULL)
        del->next->prev = del->prev;
 
    /* Change prev only if node to be deleted is NOT the first node */
    if (del->prev != NULL)
        del->prev->next = del->next;


 	//next line to return
	struct Lines* next = del->next;

    /* Finally, free the memory occupied by del*/
    free(del);
    return next;
}

// Linked list data structure 
struct Node { 
    float data; 
    struct Node* next; 
    struct Node* prev;
}; 

//deleting a node
void deleteNode(struct Node** head_ref, struct Node* del)
{
    /* base case */
    if (*head_ref == NULL || del == NULL)
        return;
 
    /* If node to be deleted is head node */
    if (*head_ref == del)
        *head_ref = del->next;
 
    /* Change next only if node to be deleted is NOT the last node */
    if (del->next != NULL)
        del->next->prev = del->prev;
 
    /* Change prev only if node to be deleted is NOT the first node */
    if (del->prev != NULL)
        del->prev->next = del->next;
 
    /* Finally, free the memory occupied by del*/
    free(del);
    return;
}

//swap two indexes in array
void swap(float *p,float *q) {
   float t;
   
   t=*p; 
   *p=*q; 
   *q=t;
}

//sort array
void sort(float a[],int n) { 
   int i,j,temp;

   for(i = 0;i < n-1;i++) {
      for(j = 0;j < n-i-1;j++) {
         if(a[j] > a[j+1])
            swap(&a[j],&a[j+1]);
      }
   }
}
 
// Standard partition process of QuickSort(). 
// It considers the last element as pivot 
// and moves all smaller element to left of 
// it and greater elements to right
//from geeks for geeks
int partition(float arr[], int l, int r) 
{ 
    int x = arr[r], i = l; 
    for (int j = l; j <= r - 1; j++) { 
        if (arr[j] <= x) { 
            swap(&arr[i], &arr[j]); 
            i++; 
        } 
    } 
    swap(&arr[i], &arr[r]); 
    return i; 
} 

// This function returns k'th smallest  
// element in arr[l..r] using QuickSort  
// based method.  ASSUMPTION: ALL ELEMENTS 
// IN ARR[] ARE DISTINCT 
float kthSmallest(float arr[], int l, int r, int k) 
{ 
    // If k is smaller than number of  
    // elements in array 
    if (k > 0 && k <= r - l + 1) { 
  
        // Partition the array around last  
        // element and get position of pivot  
        // element in sorted array 
        int index = partition(arr, l, r); 
  
        // If position is same as k 
        if (index - l == k - 1) 
            return arr[index]; 
  
        // If position is more, recur  
        // for left subarray 
        if (index - l > k - 1)  
            return kthSmallest(arr, l, index - 1, k); 
  
        // Else recur for right subarray 
        return kthSmallest(arr, index + 1, r, k - index + l - 1);
    }
  
    // If k is more than number of  
    // elements in array 
    return (float)INT_MAX; 
} 

float pickPivot(float arr[],int size)
{
	if(size < 5){
		return kthSmallest(arr, 0, size - 1, size/2);
	}

	//make array to slpit value between chunks
	float chunks[size/5][5];

	//fill up chunks array
	for(int i = 0; i < (size/5)*5; i++){

		chunks[i/5][i%5] = arr[i];
	}

	//sort each chunk
	for(int i = 0; i< size/5; i++){
		sort(&chunks[i],5);
	}

	//make an array to store medians of sub arrays
	float medians[size/5];

	for(int i =0; i< size/5; i++){

		medians[i] = chunks[i][2];
	}

	return kthSmallest(medians, 0,size - 1, size/2);

}

//size of alpha and beta must be equal
float* pairing(struct Lines* line,float out[]){

	//coutput of function
	float f_of_x;

	//save starting address of line
	struct Lines* start = line;

	//variables to find min and max of alpha
	int a_min = line->alpha;
	int a_max = line->alpha;

	//step 0
	//find minimum in alpha
	while(line->next != NULL){
		if(line->next->alpha < a_min)
			a_min = line->next->alpha;
		line = line->next;
	}
	printf("Minimum alpha is %d\n", a_min);
	//stop if minimum value in alpha is greater than zero
	if (a_min > 0){
		//x0 = -infinity
		out[0] = -INFINITY;
		//f(x0) = -infinity
		out[1] = -INFINITY;

		return out;	
	}

	//start from beginning
	line = start;

	
	//find maximum in alpha
	while(line->next != NULL){
		if(line->next->alpha > a_max)
			a_max = line->next->alpha;
		line = line->next;
	}

	printf("Maximum alpha is %d\n", a_max);

	//stop if maximum value in alpha is less than zero
	if (a_max < 0){
		//x0 = infinity
		out[0] = INFINITY;
		//f(x0) = infinity MAYBE WRONG CHECK DETAILS
		out[1] = INFINITY;

		return out;	
	}

	//step 1
	//variable to store number of x
	int k;

	while(1){

		//linked list initiallized to zero
		struct Node* x = (struct Node*)malloc(sizeof(struct Node));
		x->next = NULL;
		struct Node* x_start = x;

		do{

			//restart from beginning
			line = start;
/*
*/			
			//step 1
			while(line != NULL){

				if(line->next != NULL){
					if(line->alpha < line->next->alpha + TOLERANCE && line->alpha > line->next->alpha - TOLERANCE){
						//
						//	MAYBE NEED TO WATCH WHICH IS DELTED
						//
						printf("DELETE: ");
						if(line->beta > line->next->beta){
							printf("line: a = %f, b = %f\n", line->next->alpha, line->next->beta);
							//delete value from list
							line->next = deleteLine(start, line->next);	
						}
						else{
							printf("line: a = %f, b = %f\n", line->alpha, line->beta);
							if(line == start){
								//delete next value from list
								line = deleteLine(start, line);
								start = line;
							}
							else
								line = deleteLine(start, line);
						}
						continue;
					}
					else{
						printf("line1: a = %f, b = %f\n", line->alpha, line->beta);
						printf("line2: a = %f, b = %f\n", line->next->alpha, line->next->beta);
						x->next = (struct Node*)malloc(sizeof(struct Node));
						x->next->data = (line->beta - line->next->beta)/(line->next->alpha - line->alpha);
						x->next->next = NULL;
						x = x->next;
						printf("x = %f\n",x->data);
					}
				}

				if(line->next == NULL)
					line = line->next;
				else
					line = line->next->next;

			}
			//Step 2 count k, the number of xi created
			k = 0;
			x = x_start->next;
			while(x != NULL){
				k += 1;
				x = x->next;
				
			}	

			//TODO Check k

			printf("k is %d\n",k);	

			} while(k == 0);

			//make linked list into array
			float x_arr[k];

			int index = 0; 
   			x = x_start->next;

   			// Traverse the Linked List and add the 
			// elements to the array one by one 
			while (x != NULL) { 
			    x_arr[index++] = x->data; 
			    x = x->next; 
			} 

			//display array
			printf("size of x_array is %d\n", sizeof(x_arr)/sizeof(x_arr[0]));

			for(int i = 0; i < sizeof(x_arr)/sizeof(x_arr[0]); i++){
				printf("x%d is %f\n",i,x_arr[i]);
			}

			//calculate median
			float x_median;

			if(k > 1){
				//use linear median algorithm
				x_median = kthSmallest(x_arr, 0, k-1, k/2);
			}
			else{
				//when thereis one value no it is assigned median
				x_median = x_arr[0];
			}

			printf("median is %f\n", x_median);

			//Step 3

			//Sstart from beginning
			line = start;

			float f_max = -INFINITY;
			
			//find maximum f(xmean)
			while(line != NULL){
				float new_value = line->alpha*x_median +line->beta;
				if(new_value > f_max)
					f_max = new_value;
				line = line->next;

			}

			printf("f_max is %f\n",f_max);

			

			//start from beginning
			line = start;


			//calculate left and right gradients
			float lambda =line->alpha;
			float row = line->alpha;

			//find minimum in alpha
			while(line->next != NULL){
				if(line->next->alpha < lambda && ((line->next->alpha*x_median + line->next->beta) < f_max + TOLERANCE) && ((line->next->alpha*x_median + line->next->beta) > f_max - TOLERANCE))
					lambda = line->next->alpha;
				line = line->next;
			}

			printf("lambda = %f\n", lambda);
			//start from beginning
			line = start;


			//find maximum in alpha
			while(line->next != NULL){
				if(line->next->alpha > row && ((line->next->alpha*x_median + line->next->beta) < f_max + TOLERANCE) && ((line->next->alpha*x_median + line->next->beta) > f_max - TOLERANCE))
					row = line->next->alpha;
				line = line->next;
			}

			printf("row = %f\n", row);

			// step 4 exit if optimatl is found
			if( lambda <= TOLERANCE && row >= -TOLERANCE){
				out[0] = x_median;
				out[1] = f_max;
				return out;
			}

			 if(lambda > TOLERANCE){
				//start from beginning
				line = start;
				//delte higher of pair greater than median
				while(line != NULL){
					if(line->next != NULL){
						if((line->beta - line->next->beta)/(line->next->alpha - line->alpha) >= x_median){
								printf("line1: a = %f, b = %f\n", line->alpha, line->beta);
								printf("line2: a = %f, b = %f\n", line->next->alpha, line->next->beta);
								printf("LAMBDA DELETE: line: a = %f, b = %f x = %f\n", line->next->alpha, line->next->beta,(line->beta - line->next->beta)/(line->next->alpha - line->alpha));
								line->next = deleteLine(start, line->next);
								line = line->next;
								continue;
							}
					}
					if(line->next == NULL){
						line = line->next;
					}
					else{
						line = line->next->next;
					}
				}
			}

			if(row < -TOLERANCE){
				//start from beginning
				line = start;
				//delete lower of pair less than median
				while(line != NULL){
					if(line->next != NULL){
						if((line->beta - line->next->beta)/(line->next->alpha - line->alpha) <= x_median){
								//account for when first line is deleted
								if(line == start){
									start = line->next;
								}
								printf("line1: a = %f, b = %f\n", line->alpha, line->beta);
								printf("line2: a = %f, b = %f\n", line->next->alpha, line->next->beta);
								printf("ROW DELETE: line: a = %f, b = %f  x = %f\n", line->alpha, line->beta, (line->beta - line->next->beta)/(line->next->alpha - line->alpha));
								line = deleteLine(start, line);	
								line = line->next;
								continue;
							}
					}
					if(line->next == NULL){
						line = line->next;
					}
					else{
						line = line->next->next;
					}
				}
			}
	}

}

	struct Lines* make_lines(float alphas[], float betas[],int a_size){

		struct Lines* current = (struct Lines*)malloc(sizeof(struct Lines));
		struct Lines* first = current;

		current->alpha = alphas[0];
		current->beta = betas[0];
		current->next = NULL;
		current->prev = NULL;

		for(int i = 1; i < a_size; i++){

			//alocate memory for new linked list element
			struct Lines* newLine = (struct Lines*)malloc(sizeof(struct Lines));

			//link the new element
			current->next = newLine;

			//initiate values
			newLine->alpha = alphas[i];
			newLine->beta = betas[i];
			newLine->next = NULL;
			newLine->prev = current;

			//increment
			current = current->next;

		}

		return first;

	}

	int main(){

		// float alphas[] = {2,-1,.4,-5,.3};
		// float betas[] = {0,5,-1,6,-1};
		// int a_size = 5;

		float alphas[] = {1,1,3,.2,3,3,7,.3,-12,-1,-2,.67,1,-4,45,3.43,-6,-4,-2,4};
		float betas[] = {5,6,15,0,6,1,2,1,3,4,5,6,3,4,2,3,2,4,1,0};
		int a_size = 20;
		
		struct Lines* line = make_lines(alphas,betas,a_size);
		float out[2];
		pairing(line,out);

			printf("x* = %f\n", out[0]);
			printf("f* = %f\n", out[1]);


	}