
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




	struct Lines* make_lines(double alphas[], double betas[],int a_size){

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

	int main(int argc, char**argv){

		// double alphas[] = {2,-1,.4,-5,.3};
		// double betas[] = {0,5,-1,6,-1};
		// int a_size = 5;

		/* FOR TESTING
		char* fn;
		int min = 1;
		int max = 100;
		int inc = 1;

		if(argc == 4){
			min = atoi(argv[1]);
			max = atoi(argv[2]);
			inc = atoi(argv[3]);
		}

		for(int i = min; i < max + 1; i += inc){
			printf("<====== NEW PROBLEM %d ROWS and 2 COLS ======>\n", i);
			make_random_problem(i,2);
			solve_random_problem(i,2);
		}*/

		//store output from paring problem
        double out[2];

        // if(pair_parameters->line != NULL){
        //     pairing(pair_parameters->line,out,pair_parameters->u1,pair_parameters->u2);
        // }
        // else{
        //     out[0] = 0;
        //     out[1] = 1;
        // }
        double a1s [] = {1,0,3,-1,0};
		double a2s [] = {0,2,2,0,-1};
		double bs [] = {4,12,18,0,0};
		double cost []=  {3,5};

        GetOptPoint(a1s, a2s,bs, cost[0], cost[1], 6, &out[0], &out[1]);

        double x1 = out[0];
        double y1 = out[1];
        

        printf("<=============== RESULTS DYER===============>\n");
        printf("x = %f\n", x1);
        printf("y = %f\n", y1);

        struct tv_return* pair_parameters = two_variable(cost[0],cost[1],a1s,a2s,bs,5);
        pairing(pair_parameters->line,out,pair_parameters->u1,pair_parameters->u2,cost);

        x1 = out[0];
        y1 = out[1];
        

        printf("<=============== RESULTS PERSONAL ===============>\n");
        printf("x = %f\n", x1);
        printf("y = %f\n", y1);
/*
		double a1s [] = {1,0,3};
		double a2s [] = {0,2,2};
		double bs [] = {4,12,18};

		struct tv_return* pair_parameters = two_variable(3,5,a1s,a2s,bs,3);

		//store output from paring problem
		double out[2];

		if(pair_parameters->line != NULL){
			pairing(pair_parameters->line,out,pair_parameters->u1,pair_parameters->u2);
		}
		else{
			out[0] = 0;
			out[1] = 1;
		}

		double x1 = out[0]/out[1];
		double y1 = (1-3*out[0])/5*out[1];
		

		printf("<=============== RESULTS ===============>\n");
		printf("x = %f\n", x1);
		printf("y = %f\n", y1);
		return 0;
*/
		// double alphas[] = {1,1,3,.2,3,3,7,.3,-12,-1,-2,.67,1,-4,45,3.43,-6,-4,-2,4};
		// double betas[] = {5,6,15,0,6,1,2,1,3,4,5,6,3,4,2,3,2,4,1,0};
		// int a_size = 20;
		
		// struct Lines* line = make_lines(alphas,betas,a_size);
		// double out[2];
		// pairing(line,out);

		// printf("<=============== RESULTS ===============>\n");
		// printf("x* = %f\n", out[0]);
		// printf("f* = %f\n", out[1]);
		// double values[] = {2.3,1.5,3.6,7.6,3.4,2.3,2.1,0.23412,1.234,9.432};
		// int size = 10;

		// printf("median from quickselct is %f\n", quickSelect(values, size, size/2, pickPivot));
		// printf("median from nlogn median is %f\n", nlogn_median(values,size));

	}