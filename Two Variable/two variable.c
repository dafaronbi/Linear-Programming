#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUMBER_OF_CONTRAINTS 2


struct Lines{
	float alpha;
	float beta;
	struct Lines* next;
	struct Lines* prev;
};

void deleteLine(struct Lines* head_ref, struct Lines* del)
{
    /* base case */
    if (head_ref == NULL || del == NULL)
        return;
 
    /* If node to be deleted is head node */
    if (head_ref == del)
        head_ref = del->next;
 
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
void swap(int *p,int *q) {
   int t;
   
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

// Function for calculating median 
double findMedian(float a[], int n) 
{ 
    // First we sort the array 
    sort(a, a + n); 
  
    // check for even case 
    if (n % 2 != 0) 
        return (double)a[n / 2]; 
  
    return (double)(a[(n - 1) / 2] + a[n / 2]) / 2.0; 
} 

//size of alpha and beta must be equal
float* pairing(struct Lines* line){
	float out[2];

	//coutput of function
	float f_of_x;

	//save starting address of line
	struct Lines* start = line;

	//linked list initiallized to zero
	struct Node* x = (struct Node*)malloc(sizeof(struct Node));
	struct Node* x_start = x;

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
	//stop if maximum value in alpha is less than zero
	if (a_max < 0){
		//x0 = infinity
		out[0] = INFINITY;
		//f(x0) = infinity
		out[1] = INFINITY;

		return out;	
	}

	//variable to store number of x
	int k;

	while(1){

		do{

			//restart from beginning
			line = start;
			//step 1
			while(line != NULL){
				printf("here\n");
				//if lines are parrallel take the line with higher beta and delete other from list
				if(line->alpha < line->next->alpha + 1e-6 && line->alpha > line->next->alpha - 1e-6){
					if(line->beta > line->next->beta + 1e-6){
						printf("here3\n");
						//delete value from list
						deleteLine(start, line->next);	
					}
					else{
						//delete next value from list
						deleteLine(start, line);
					}
				}
				else{
					printf("here2\n");
					x->next = (struct Node*)malloc(sizeof(struct Node));
					x->next->data = line->beta - line->next->beta/(line->next->alpha - line->alpha);
					x = x->next;

				}
				line = line->next;

			}
			//count k, the number of xi created
			k = 0;

				x = x_start->next;
				while(x != NULL){
					k++;
					x = x->next;
				}
			} while(k == 0);


			//make linked list into array
			float x_arr[k];

			int index = 0; 
   			x = x_start;

   			// Traverse the Linked List and add the 
			// elements to the array one by one 
			while (x != NULL) { 
			    x_arr[index++] = x->data; 
			    x = x->next; 
			} 



			//calculate median
			float x_median = findMedian(x_arr,k);

			//start from beginning
			line = start;

			float f_max = line->alpha*x_median +line->beta;

			//find maximum f(xmean)
			while(line->next != NULL){
				float new_value = line->next->alpha*x_median +line->next->beta;
				if(new_value > f_max)
					f_max = new_value;
				line = line->next;
			}

			

			//start from beginning
			line = start;


			//calculate left and right gradients
			float lambda =line->alpha;
			float row = line->alpha;

			//find minimum in alpha
			while(line->next != NULL){
				if(line->next->alpha < lambda && line->next->alpha*x_median + line->next->beta == f_max)
					lambda = line->next->alpha;
				line = line->next;
			}

			//start from beginning
			line = start;


			//find median in alpha
			while(line->next != NULL){
				if(line->next->alpha > row && line->next->alpha*x_median + line->next->beta == f_max)
					row = line->next->alpha;
				line = line->next;
			}


			//exit if optimatl is found
			if( lambda <= 0 && row >= 0){
				out[0] = x_median;
				out[1] = f_max;
				return out;
			}


			//step 1
			while(line != NULL){
				if((line->beta - line->next->beta)/line->next->alpha - line->alpha > x_median){
		
						deleteLine(start, line);	
					}
				line = line->next;
			}
	}



	}

	struct Lines* make_lines(){
		//save starting address of line
		struct Lines* first = (struct Lines*)malloc(sizeof(struct Lines));
		//save starting address of line
		struct Lines* second = (struct Lines*)malloc(sizeof(struct Lines));

		first->alpha = 2;
		first->beta = 0;
		first->next = second;

		second->alpha = -1;
		second->beta = 5;
		second->next = NULL;

		return first;


	}

	int main(){
		
		struct Lines* line = make_lines();
		float* out = pairing(line);

		for(int i = 0;i < 2;i++) {
			printf("%f", out[i]);
			printf("\n");
		}


	}