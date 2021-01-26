
#include "two variable.h"

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

//deleting a node
void deleteNode(struct Node* head_ref, struct Node* del)
{
    /* base case */
    if (head_ref == NULL || del == NULL)
        return;
 
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

float quickSelect(float arr[], int size, int k,  float (*pivot_fn)(float [],int )){

    // Select the kth element in l (0 based)
    // :param l: List of numerics
    // :param k: Index
    // :param pivot_fn: Function to choose a pivot, defaults to random.choice
    // :return: The kth element of l
    if(size == 1){
        return arr[0];
    }

    float pivot = pivot_fn(arr,size);


    //make linked list of values lower than pivot
	struct Node* lows = (struct Node*)malloc(sizeof(struct Node));
	lows->next = NULL;
	struct Node* lows_start = lows;
	int lows_size = 0;

	//make linked list of values higher than pivot
	struct Node* highs = (struct Node*)malloc(sizeof(struct Node));
	highs->next = NULL;
	struct Node* highs_start = highs;
	int highs_size = 0;

	//make linked list of values equal to pivot
	struct Node* pivots = (struct Node*)malloc(sizeof(struct Node));
	pivots->next = NULL;
	struct Node* pivots_start = pivots;
	int pivots_size = 0;

	for(int i = 0; i < size; i++){
		//add to lows linked list if lower than pivot
		if(arr[i] < pivot){
			lows->next = (struct Node*)malloc(sizeof(struct Node));
			lows->next->data = arr[i];
			lows->next->next = NULL;
			lows = lows->next;
			lows_size++;
		}
		//add to highs linked list if higher than pivot
		if(arr[i] > pivot){
			highs->next = (struct Node*)malloc(sizeof(struct Node));
			highs->next->data = arr[i];
			highs->next->next = NULL;
			highs = highs->next;
			highs_size++;
		}
		//add to pivots linked list if equal to pivot
		if(arr[i] == pivot){
			pivots->next = (struct Node*)malloc(sizeof(struct Node));
			pivots->next->data = arr[i];
			pivots->next->next = NULL;
			pivots = pivots->next;
			pivots_size++;
		}
	}

	//find k in lower array
    if( k < lows_size){
    	//make linked list into array
		float lows_arr[lows_size];

		int index = 0; 
		lows = lows_start->next;

			// Traverse the Linked List and add the 
		// elements to the array one by one 
		while (lows != NULL) { 
		    lows_arr[index++] = lows->data; 
		    lows = lows->next; 
		} 
        return quickSelect(lows_arr, lows_size,k, pivot_fn);
    }

    //We got lucky and guessed the median
    if(k < lows_size + pivots_size){
        return pivots_start->next->data;
    }
    
    //make linked list into array
	float highs_arr[highs_size];

	int index = 0; 
	highs = highs_start->next;

	// Traverse the Linked List and add the 
	// elements to the array one by one 
	while (highs != NULL) { 
	    highs_arr[index++] = highs->data; 
	    highs = highs->next;
	} 
    return quickSelect(highs_arr, highs_size, k - lows_size - pivots_size, pivot_fn);

}

float nlogn_median(float arr[],int size){
	sort(arr,size);
    return arr[size / 2];
}


float quickSelect_median(float arr[],int size, float (*pivot_fn)(float [],int )){
	return quickSelect(arr, size,size / 2, pivot_fn);
}

float pickPivot(float arr[],int size)
{
	if(size < 5){
		return nlogn_median(arr,size);
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

	//fun_ptr is a pointer to function fun()  
    float (*pivot_fn)(float[],int) = &pickPivot; 

	return quickSelect_median(medians, size/5,pivot_fn);

}

//size of alpha and beta must be equal
float* pairing(struct Lines* line,float out[],float u1, float u2){

	//coutput of function
	float f_of_x;

	//save starting address of line
	struct Lines* start = line;



	//step 0
	//find minimum in alpha

	printf("<=============== STEP 0 ===============>\n");

	//out zero when u1 > u2
	if(u1 > u2){
		out[0] = 0;
		out[1] = 0;
		return out;
	}

	float f_of_u1 = -INFINITY;

	//determine f(u1)
	line = start;
	while(line != NULL){
		float new_value = line->alpha*u1 +line->beta;
		if(new_value > f_of_u1 )
			f_of_u1 = new_value;
		line = line->next;
		}

	printf("f(u1) = %f\n", f_of_u1);

	
	//out u1 and F(u1) when u1 = u2
	if(u1 == u2){
		out[0] = u1;
		out[2] = f_of_u1;
		return out;
	}

	//start from beginning of lines
	line = start;

	//calculate left and right gradients
	float u_row = line->alpha;
	float u_lambda =line->alpha;
	

	//find maximum in alpha at u1
	line = start;
	while(line->next != NULL){
		if(line->next->alpha > u_row && ((line->next->alpha*u1 + line->next->beta) < f_of_u1 + TOLERANCE) && ((line->next->alpha*u1 + line->next->beta) > f_of_u1 - TOLERANCE))
			u_row = line->next->alpha;
		line = line->next;
	}

	printf("u_row = %f\n", u_row);

	if(u_row >= -TOLERANCE){
		out[0] = u1;
		out[2] = f_of_u1;
	}


	//determine f(u2)
	float f_of_u2 = -INFINITY;

	line = start;
	while(line != NULL){
		float new_value = line->alpha*u2 +line->beta;
		if(new_value > f_of_u2 )
			f_of_u2 = new_value;
		line = line->next;
		}

	printf("f(u2) = %f\n", f_of_u2);


	//find minum in alpha at u2
	line = start;
	while(line->next != NULL){
		if(line->next->alpha < u_lambda && ((line->next->alpha*u2 + line->next->beta) < f_of_u2 + TOLERANCE) && ((line->next->alpha*u2 + line->next->beta) > f_of_u2 - TOLERANCE))
			u_lambda = line->next->alpha;
		line = line->next;
	}

	printf("u_lambda = %f\n", u_lambda);

	if(u_row <= TOLERANCE){
		out[0] = u2;
		out[2] = f_of_u2;
	}

	//find minimum alpha
	line = start;

	//variables to find min and max of alpha
	float a_min = line->alpha;
	float a_max = line->alpha;

	while(line != NULL){
		if(line->alpha < a_min)
			a_min = line->alpha;
		line = line->next;
	}


	printf("Minimum alpha is %f\n", a_min);
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
	while(line != NULL){
		if(line->alpha > a_max)
			a_max = line->alpha;
		line = line->next;
	}

	printf("Maximum alpha is %f\n", a_max);

	//stop if maximum value in alpha is less than zero
	if (a_max < 0){
		//x0 = infinity
		out[0] = INFINITY;
		//f(x0) = infinity MAYBE WRONG CHECK DETAILS
		out[1] = INFINITY;

		return out;	
	}

	//variable to store number of x
	int k;

	while(1){
		//step 1
		printf("<=============== STEP 1 ===============>\n");
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
						float value = (line->beta - line->next->beta)/(line->next->alpha - line->alpha);

						//delete right half if xi value calculated is greater than u2
						if(value > u2){
							if(line->next->alpha > line->alpha){
								printf("OUT OF U2 DELETE: line: a = %f, b = %f x = %f\n", line->next->alpha, line->next->beta,(line->beta - line->next->beta)/(line->next->alpha - line->alpha));
								line->next = deleteLine(start, line->next);
								line = line->next;
							}
							else{
								printf("OUT OF U2 DELETE: line: a = %f, b = %f x = %f\n", line->alpha, line->beta,(line->beta - line->next->beta)/(line->next->alpha - line->alpha));
								line = deleteLine(start, line);
							}

						}
						//delete left half if xi value calculated is less than u2
						else if(value < u1){
							if(line->next->alpha < line->alpha){
								printf("OUT OF U1 DELETE: line: a = %f, b = %f x = %f\n", line->next->alpha, line->next->beta,(line->beta - line->next->beta)/(line->next->alpha - line->alpha));
								line->next = deleteLine(start, line->next);
								line = line->next;
							}
							else{
								printf("OUT OF U1 DELETE: line: a = %f, b = %f x = %f\n", line->alpha, line->beta,(line->beta - line->next->beta)/(line->next->alpha - line->alpha));
								line = deleteLine(start, line);
							}

						}
						else{
							x->next = (struct Node*)malloc(sizeof(struct Node));
							x->next->data = value;
							x->next->next = NULL;
							x = x->next;
							printf("x = %f\n",x->data);
						}
					}
				}

				if(line->next == NULL)
					line = line->next;
				else
					line = line->next->next;

			}
			//Step 2 count k, the number of xi created
			printf("<=============== STEP 2 ===============>\n");
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
			printf("size of x_array is %lu\n", sizeof(x_arr)/sizeof(x_arr[0]));

			for(int i = 0; i < sizeof(x_arr)/sizeof(x_arr[0]); i++){
				printf("x%d is %f\n",i,x_arr[i]);
			}

			//calculate median
			float x_median;

			if(k > 1){

				if(k%2){
					//use linear median algorithm with k/2 when even number k
					x_median = quickSelect(x_arr, k, k/2, &pickPivot);

				}
				else{
					//use linear median algorithm with k/2 - 1 when odd number k
					x_median = quickSelect(x_arr, k, k/2-1, &pickPivot);

				}
			}
			else{
				//when thereis one value no it is assigned median
				x_median = x_arr[0];
			}

			printf("median is %f\n", x_median);

			//Step 3
			printf("<=============== STEP 3 ===============>\n");

			//Sstart from beginning
			line = start;

			float f_max = -INFINITY;
			
			//find maximum f(xmean)
			while(line != NULL){
				float new_value = line->alpha*x_median +line->beta;
				if(new_value > f_max )
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
			printf("<=============== STEP 4 ===============>\n");

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
								//delete the "right half" or teh value with the greater alpha
								if(line->next->alpha > line->alpha){
									printf("LAMBDA DELETE: line: a = %f, b = %f x = %f\n", line->next->alpha, line->next->beta,(line->beta - line->next->beta)/(line->next->alpha - line->alpha));
									line->next = deleteLine(start, line->next);
									line = line->next;
								}
								else{
									printf("LAMBDA DELETE: line: a = %f, b = %f x = %f\n", line->alpha, line->beta,(line->beta - line->next->beta)/(line->next->alpha - line->alpha));
									line = deleteLine(start, line);
									line = line->next;
								}
								
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
								//delete the "right half" or teh value with the greater alpha
								if(line->next->alpha < line->alpha){
									printf("ROW DELETE: line: a = %f, b = %f x = %f\n", line->next->alpha, line->next->beta,(line->beta - line->next->beta)/(line->next->alpha - line->alpha));
									line->next = deleteLine(start, line->next);
									line = line->next;
								}
								else{
									printf("ROW DELETE: line: a = %f, b = %f x = %f\n", line->alpha, line->beta,(line->beta - line->next->beta)/(line->next->alpha - line->alpha));
									line = deleteLine(start, line);
								}

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

	struct tv_return* two_variable(float c1, float c2,float a1 [], float a2 [], float b [], int num_contraints){

		//set initial starting point
		float x1_start = 0;
		float x2_start = 0;

		//value to return
		struct tv_return* output = (struct tv_return*)malloc(sizeof(struct tv_return));

		//print original problem
		printf("<=============== Original ===============>\n");
		printf("max %fx1 + %fx2\n",c1,c2);
		for(int i =0; i< num_contraints;i++){
			printf("%fx1 + %fx2 <= %f\n",a1[i],a2[i],b[i]);
		}


		printf("<=============== CONVERT ===============>\n");

		//optimal if c1 and c2 coefficients are zero
		if(c1 == 0 && c2 == 0){
			output->line = NULL;
			output->u1 = 0;
			output->u2 = 0;
			return output;
		}

		//allocate memory for contraints
		float A1[num_contraints];
		float A2[num_contraints];
		float B[num_contraints];

		//keep track of I0 and I1 contraints
		int num_I0 = 0;
		int num_I1 = 0;

		printf("NEW CONSTRAINTS\n");

		//caluclue Ai1 and Ai2 from inputted ai1 and ai2
		for(int i =0; i< num_contraints;i++){
			A1[i] = a1[i] - a2[i]*c1/c2;
			A2[i] = a2[i]/c2;
			B[i] = b[i] - a1[i]*x1_start -a2[i]*x2_start;

			printf("%fx + %f <= %f\n",A1[i],A2[i],B[i]);

			if(B[i] == 0){
				num_I0++;
			}
			if(B[i] <0){
				num_I1++;
			}
		}

		printf("THERE ARE %d I0 CONSTRAINTS\n",num_I0);
		printf("THERE ARE %d I1 CONSTRAINTS\n",num_I1);

		//alocated memory for u1, u2, and I1 contraints
		float u1 = -INFINITY;
		float u2 = INFINITY;
		float alphas[num_I1];
		float betas[num_I1];

		//keep track of elemtents in alphas and betas array
		int I1_index = 0;

		printf("I1 CONSTRAINTS\n");

		for(int i =0; i< num_contraints;i++){
			if(B[i] == 0){
				if(A1[i] > 0 && (-A2[i]/A1[i]) > u1)
					u1 = -A2[i]/A1[i];
				if(A1[i] < 0 && (-A2[i]/A1[i]) < u2)
					u2 = -A2[i]/A1[i];
			}
			if(B[i] <0){
				alphas[I1_index] =  A1[i]/B[i];
				betas[I1_index] = A2[i]/B[i];
				printf("y = %fx + %f\n", alphas[I1_index],betas[I1_index]);
				I1_index++;
			}
		}

		printf("u1 = %f\nu2 = %f\n",u1,u2);


		//make lines from generated constraints
		struct Lines* line = make_lines(alphas,betas,num_I1);

		output->line = line;
		output->u1 = u1;
		output->u2 = u2;

		return output;	

	}

	int main(){

		// float alphas[] = {2,-1,.4,-5,.3};
		// float betas[] = {0,5,-1,6,-1};
		// int a_size = 5;

		float a1s [] = {2,-1,.4,-5,.3};
		float a2s [] = {0,5,-1,6,-1};
		float bs [] = {2,3,-1.3,3,-2};

		struct tv_return* pair_parameters = two_variable(1,2,a1s,a2s,bs,5);

		//store output from paring problem
		float out[2];

		if(pair_parameters->line != NULL){
			pairing(pair_parameters->line,out,pair_parameters->u1,pair_parameters->u2);
		}
		else{
			out[0] = 0;
			out[1] = 1;
		}
		

		printf("<=============== RESULTS ===============>\n");
		printf("x* = %f\n", out[0]);
		printf("f* = %f\n", out[1]);

		// float alphas[] = {1,1,3,.2,3,3,7,.3,-12,-1,-2,.67,1,-4,45,3.43,-6,-4,-2,4};
		// float betas[] = {5,6,15,0,6,1,2,1,3,4,5,6,3,4,2,3,2,4,1,0};
		// int a_size = 20;
		
		// struct Lines* line = make_lines(alphas,betas,a_size);
		// float out[2];
		// pairing(line,out);

		// printf("<=============== RESULTS ===============>\n");
		// printf("x* = %f\n", out[0]);
		// printf("f* = %f\n", out[1]);
		// float values[] = {2.3,1.5,3.6,7.6,3.4,2.3,2.1,0.23412,1.234,9.432};
		// int size = 10;

		// printf("median from quickselct is %f\n", quickSelect(values, size, size/2, pickPivot));
		// printf("median from nlogn median is %f\n", nlogn_median(values,size));

	}