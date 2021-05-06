#include "two variable.h"

//swap two indexes in array
void swap(double *p,double *q) {
   double t;
   
   t=*p; 
   *p=*q; 
   *q=t;
}

//sort array
void sort(double a[],int n) { 
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
int partition(double arr[], int l, int r) 
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
double kthSmallest(double arr[], int l, int r, int k) 
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
    return (double)INT_MAX; 
} 

double quickSelect(double arr[], int size, int k,  double (*pivot_fn)(double [],int )){

    // Select the kth element in l (0 based)
    // :param l: List of numerics
    // :param k: Index
    // :param pivot_fn: Function to choose a pivot, defaults to random.choice
    // :return: The kth element of l
    if(size == 1){
        return arr[0];
    }

    double pivot = pivot_fn(arr,size);


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
		double lows_arr[lows_size];

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
	double highs_arr[highs_size];

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

double nlogn_median(double arr[],int size){
	sort(arr,size);
    return arr[size / 2];
}


double quickSelect_median(double arr[],int size, double (*pivot_fn)(double [],int )){
	return quickSelect(arr, size,size / 2, pivot_fn);
}

double pickPivot(double arr[],int size)
{
	if(size < 5){
		return nlogn_median(arr,size);
	}

	//make array to slpit value between chunks
	double chunks[size/5][5];

	//fill up chunks array
	for(int i = 0; i < (size/5)*5; i++){

		chunks[i/5][i%5] = arr[i];
	}

	//sort each chunk
	for(int i = 0; i< size/5; i++){
		sort(&chunks[i],5);
	}

	//make an array to store medians of sub arrays
	double medians[size/5];

	for(int i =0; i< size/5; i++){

		medians[i] = chunks[i][2];
	}

	//fun_ptr is a pointer to function fun()  
    double (*pivot_fn)(double[],int) = &pickPivot; 

	return quickSelect_median(medians, size/5,pivot_fn);

}