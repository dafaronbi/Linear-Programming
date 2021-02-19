#include "two variable.h"

//size of alpha and beta must be equal
double* pairing(struct Lines* line,double out[],double u1, double u2,double c[]){

	//save starting address of line
	struct Lines* start = line;



	//step 0
	//find minimum in alpha

	printf("<=============== STEP 0 ===============>\n");

	// out zero when u1 > u2
	if(u1 > u2){
		out[0] = 0;
		out[1] = 0;
		return out;
	}

	double f_of_u1 = -INFINITY;

	//determine f(u1)
	line = start;
	while(line != NULL){
		double new_value = line->alpha*u1 +line->beta;
		if(new_value > f_of_u1 + TOLERANCE)
			f_of_u1 = new_value;
		line = line->next;
		}

	printf("f(u1) = %f\n", f_of_u1);

	
	//out u1 and F(u1) when u1 = u2
	if(u1 == u2){
		out[0] = u1/f_of_u1;
		out[1] = (1-c[0]*u1)/(c[1]*f_of_u1);
		return out;
	}

	//start from beginning of lines
	line = start;

	//calculate left and right gradients
	double u_row = -INFINITY;
	double u_lambda = INFINITY;
	

	//find maximum in alpha at u1
	while(line != NULL){
		//calculate y of current line
		double y = line->alpha*u1 + line->beta;
		if(u_row < line->alpha  && (fabs(f_of_u1 - y) < TOLERANCE )){
			u_row = line->alpha;
			printf("possible u_row = %f\n",u_row);
		}
		line = line->next;
	}
		
	

	printf("u_row = %f\n", u_row);

	//terminate when u_row > 0
	if(u_row > -TOLERANCE){
		out[0] = u1/f_of_u1;
		out[1] = (1-c[0]*u1)/(c[1]*f_of_u1);
		return out;
	}


	//determine f(u2)
	double f_of_u2 = -INFINITY;

	line = start;
	while(line != NULL){
		double new_value = line->alpha*u2 +line->beta;
		if(new_value > f_of_u2 + TOLERANCE)
			f_of_u2 = new_value;
		line = line->next;
		}

	printf("f(u2) = %f\n", f_of_u2);


	//find minum in alpha at u2
	line = start;
	while(line != NULL){
		//calculate y of current line
		double y = line->alpha*u2 + line->beta;
		if(line->alpha < u_lambda && (fabs(f_of_u2 - y) < TOLERANCE)){
			u_lambda = line->alpha;
			printf("possible u_lambda = %f\n",u_lambda);
		}
		line = line->next;
	}

	printf("u_lambda = %f\n", u_lambda);

	if(u_lambda < TOLERANCE){
		out[0] = u2/f_of_u2;
		out[1] = (1-c[0]*u2)/(c[1]*f_of_u2);
		return out;
	}

	//find minimum alpha
	line = start;

	//variables to find min and max of alpha
	double a_min = line->alpha;
	double a_max = line->alpha;

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

		struct Lines* pair = (struct Lines*)malloc(sizeof(struct Lines));
		pair->next = NULL;
		struct Lines* pair_start = pair;

		do{

			//restart from beginning
			line = start;
			
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
						double value = (line->beta - line->next->beta)/(line->next->alpha - line->alpha);

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
							pair->next = (struct Lines*)malloc(sizeof(struct Lines));
							pair = pair->next;
							pair->alpha = line->alpha;
							pair->beta = line->beta;
							pair->next = (struct Lines*)malloc(sizeof(struct Lines));
							pair = pair->next;
							pair->alpha = line->next->alpha;
							pair->beta = line->next->beta;
							pair->next = NULL;

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
			double x_arr[k];

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
			double x_median;

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

			//Start from beginning
			line = start;

			double f_max = -INFINITY;
			
			//find maximum f(xmean)
			while(line != NULL){
				double new_value = line->alpha*x_median +line->beta;
				if(new_value > f_max ){
					f_max = new_value;
					printf("possible f_max is %f from line a=%f b=%f\n",f_max,line->alpha, line->beta);
				}
				line = line->next;

			}

			printf("f_max is %f\n",f_max);

			

			//start from beginning
			line = start;


			//calculate left and right gradients
			double lambda = INFINITY;
			double row = -INFINITY;

			//find minimum in alpha
			while(line != NULL){
				double y = line->alpha*x_median + line->beta;
				if(line->alpha < lambda && (fabs(f_max - y) < TOLERANCE)){
					lambda = line->alpha;
					printf("possible lambda = %f\n", lambda);
				}
				line = line->next;
			}

			printf("lambda = %f\n", lambda);

			//start from beginning
			line = start;

			//find maximum in alpha
			while(line != NULL){
				double y = line->alpha*x_median + line->beta;
				if(line->alpha > row && (fabs(f_max - y) < TOLERANCE)){
					row = line->alpha;
					printf("possible row = %f\n", row);
				}
				line = line->next;
			}

			printf("row = %f\n", row);

			// step 4 exit if optimatl is found
			printf("<=============== STEP 4 ===============>\n");

			if( lambda <= TOLERANCE && row >= -TOLERANCE){
				printf("NO NEED TO DELETE\n");
				out[0] = x_median/f_max;
				out[1] = (1-c[0]*x_median)/(c[1]*f_max);
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
								//account for when first line is deleted
								if(line == start){
									start = line->next;
								}
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

	struct tv_return* two_variable(double c1, double c2,double a1 [], double a2 [], double b [], int num_contraints){

		//set initial starting point
		double x1_start = 0;
		double x2_start = 0;

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
		double A1[num_contraints];
		double A2[num_contraints];
		double B[num_contraints];

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
			if(B[i] >0){
				num_I1++;
			}
		}

		printf("THERE ARE %d I0 CONSTRAINTS\n",num_I0);
		printf("THERE ARE %d I1 CONSTRAINTS\n",num_I1);

		//alocated memory for u1, u2, and I1 contraints
		double u1 = -INFINITY;
		double u2 = INFINITY;
		double alphas[num_I1];
		double betas[num_I1];

		//keep track of elemtents in alphas and betas array
		int I1_index = 0;

		printf("I1 CONSTRAINTS\n");

		for(int i =0; i< num_contraints;i++){
			if(B[i] == 0){
				if(A1[i] > 0 && (-A2[i]/A1[i]) < u2)
					u2 = -A2[i]/A1[i];
				if(A1[i] < 0 && (-A2[i]/A1[i]) > u1)
					u1 = -A2[i]/A1[i];
			}
			if(B[i] > 0){
				alphas[I1_index] =  A1[i]/B[i];
				betas[I1_index] = A2[i]/B[i];
				printf("y => %fx + %f\n", alphas[I1_index],betas[I1_index]);
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