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

struct Lines* deleteLine(struct Lines* head_ref, struct Lines* del)
{
    /* base case */
    if (head_ref == NULL || del == NULL)
        return del;
 
    /* If node to be deleted is head node */
    if (head_ref == del)
        head_ref = del->next;
 
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

// Function for calculating median 
double findMedian(float a[], int n) 
{ 
	// int a[11] = {10,45,35,3,24,23,78,89,33,34,11};
 //    int max=0,*ar,co=0,mid=0;
 //    for(int i=0; i<sizeof(a)/4; i++)   // O(n)
 //    {
 //        if(i==0)
 //            max = a[i];
 //        if(a[i]>max)
 //            max = a[i];
 //    }
 //    ar =(int*)malloc((sizeof(int)*max)+1); // O(n)
 //    for(int i=0; i<sizeof(a)/4; i++)
 //        ar[a[i]] = 1;
 //    for(int i=0; i<=max; i++)  // 0(n)
 //        if(ar[i]==1)
 //        {
 //            co++;
 //            if(co == ((sizeof(a)/8)+1))
 //                mid = i;
 //        }
 //    printf("%d",mid);
    // First we sort the array 
    sort(a, n); 
  
    // check for even case 
    if (n % 2 != 0) 
        return (float)a[n / 2]; 
  
    return (float)(a[n/ 2]); 
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
		//f(x0) = infinity
		out[1] = INFINITY;

		return out;	
	}

	//step 1
	//variable to store number of x
	int k;

	while(1){

		//linked list initiallized to zero
		struct Node* x = (struct Node*)malloc(sizeof(struct Node));
		struct Node* x_start = x;

		do{

			//restart from beginning
			line = start;
/*
*/			
			//step 1
			while(line != NULL){
				if(line->next != NULL){
					if(line->alpha < line->next->alpha + 1e-6 && line->alpha > line->next->alpha - 1e-6){
						if(line->beta > line->next->beta + 1e-6){
							//delete value from list
							line->next = deleteLine(start, line->next);	
						}
						else{
							//delete next value from list
							line = deleteLine(start, line);
						}
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
			//count k, the number of xi created
			k = 0;
			x = x_start->next;
			while(x != NULL){
				k += 1;
				x = x->next;
				
			}	

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
			printf("size of x_array is %d\n", sizeof(x_arr)/sizeof(x_arr[0]));
			for(int i = 0; i < sizeof(x_arr)/sizeof(x_arr[0]); i++){
				printf("x%d is %f\n",i,x_arr[i]);
			}

			//calculate median
			float x_median = findMedian(x_arr,k);

			printf("median is %f\n", x_median);

			//start from beginning
			line = start;

			float f_max = line->alpha*x_median +line->beta;
			printf("f_max is %f\n",f_max);

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
				if(line->next->alpha < lambda && (line->next->alpha*x_median + line->next->beta) == f_max)
					lambda = line->next->alpha;
				line = line->next;
			}

			printf("lambda = %f\n", lambda);
			//start from beginning
			line = start;


			//find median in alpha
			while(line->next != NULL){
				if(line->next->alpha > row && (line->next->alpha*x_median + line->next->beta == f_max))
					row = line->next->alpha;
				line = line->next;
			}
			printf("row = %f\n", row);
			// step 4 exit if optimatl is found
			if( lambda <= 0.000001 && row >= -0.000001){
				out[0] = x_median;
				out[1] = f_max;
				return out;
			}

			 if(lambda > 0){
				//start from beginning
				line = start;
				//delte higher of pair greater than median
				while(line != NULL){
					if(line->next != NULL){
						if((line->beta - line->next->beta)/line->next->alpha - line->alpha >= x_median){
								line->next = deleteLine(start, line->next);
							}
					}
					if(line != NULL)
						line = line->next;
				}
			}

			if(row > 0){
				//start from beginning
				line = start;
				//delete lower of pair less than median
				while(line != NULL){
					if(line->next != NULL){
						if((line->beta - line->next->beta)/line->next->alpha - line->alpha <= x_median){
								line = deleteLine(start, line);	
							}
					}
					if(line != NULL)
						line = line->next;
				}
			}
	}

}

	struct Lines* make_lines(){
		//save starting address of line
		struct Lines* first = (struct Lines*)malloc(sizeof(struct Lines));
		//save starting address of line
		struct Lines* second = (struct Lines*)malloc(sizeof(struct Lines));
		//save starting address of line
		struct Lines* third = (struct Lines*)malloc(sizeof(struct Lines));
		//save starting address of line
		struct Lines* fourth = (struct Lines*)malloc(sizeof(struct Lines));
		//save starting address of line
		struct Lines* fifth = (struct Lines*)malloc(sizeof(struct Lines));

		first->alpha = 2;
		first->beta = 0;
		first->next = second;
		first->prev = NULL;

		second->alpha = -1;
		second->beta = 5;
		second->next = third;
		second->prev = first;

		third->alpha = .4;
		third->beta = -1;
		third->next = fourth;
		third->prev = second;

		fourth->alpha = -5;
		fourth->beta = 6;
		fourth->next = fifth;
		fourth->prev = third;

		fifth->alpha = .3;
		fifth->beta = -1;
		fifth->next = NULL;
		fifth->prev = fourth;


		return first;


	}

	int main(){
		
		struct Lines* line = make_lines();
		float out[2];
		pairing(line,out);

			printf("x* = %f\n", out[0]);
			printf("f* = %f\n", out[1]);



	}