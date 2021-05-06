#include "three variable.h"

three_pairs* deletethree_pairs(three_pairs* head_ref, three_pairs* del)
{
    /* base case */
    if (head_ref == NULL || del == NULL)
        return del;
 
    /* If node to be deleted is head node */
    if (head_ref == del){
        head_ref = del->next;
    }
 
    three_pairs* temp = head_ref;

    while(temp != NULL){

        //remove deleted value from linked list chain
        if(temp->next == del){
            temp->next = del->next;
            break;
        }

        //increment
        temp = temp->next;
    }
    


    //next line to return
    temp = del->next;

    /* Finally, free the memory occupied by del*/
    free(del);
    return temp;
}

q_pairs* deleteq_pairs(q_pairs* head_ref, q_pairs* del)
{
    /* base case */
    if (head_ref == NULL || del == NULL)
        return del;
 
    /* If node to be deleted is head node */
    if (head_ref == del){
        head_ref = del->next;
    }
 
    q_pairs* temp = head_ref;

    while(temp != NULL){

        //remove deleted value from linked list chain
        if(temp->next == del){
            temp->next = del->next;
            break;
        }

        //increment
        temp = temp->next;
    }
    


    //next line to return
    temp = del->next;

    /* Finally, free the memory occupied by del*/
    free(del);
    return temp;
}

Lines_three_I1* deleteLines_three_I1(Lines_three_I1* head_ref, Lines_three_I1* del)
{
    /* base case */
    if (head_ref == NULL || del == NULL)
        return del;
 
    /* If node to be deleted is head node */
    if (head_ref == del){
        head_ref = del->next;
    }
 
    Lines_three_I1* temp = head_ref;

    while(temp != NULL){

        //remove deleted value from linked list chain
        if(temp->next == del){
            temp->next = del->next;
            break;
        }

        //increment
        temp = temp->next;
    }
    


    //next line to return
    temp = del->next;

    /* Finally, free the memory occupied by del*/
    free(del);
    return temp;
}

Lines_three_I0* deleteLines_three_I0(Lines_three_I0* head_ref, Lines_three_I0* del)
{
    /* base case */
    if (head_ref == NULL || del == NULL)
        return del;
 
    /* If node to be deleted is head node */
    if (head_ref == del){
        head_ref = del->next;
    }
 
    Lines_three_I0* temp = head_ref;

    while(temp != NULL){

        //remove deleted value from linked list chain
        if(temp->next == del){
            temp->next = del->next;
            break;
        }

        //increment
        temp = temp->next;
    }
    


    //next line to return
    temp = del->next;

    /* Finally, free the memory occupied by del*/
    free(del);
    return temp;
}

int line_search(Lines_three_I1* F, Lines_three_I0* M, double a1, double a2, double beta, double s0_x){

    if (DPDUBUG) {
      printf("<============= LINE SEARCH ===============>\n");
    }

    Lines_three_I1* F_start = F;
    Lines_three_I0* M_start = M;

    //go to next (first element in lisrt)
    F = F->next;
    M = M->next;

    //make lines linked list
    struct Lines* current = (struct Lines*)malloc(sizeof(struct Lines));
    struct Lines* first = current;

    //variable for output
    double out[2];

    current->alpha = F->f1;
    current->beta =  F->f3;
    current->next = NULL;
    current->prev = NULL;

    double g1 = F->f1;
    double g1m = g1;
    

    while(F->next != NULL){

            //increment to next i1 constraint
            F = F->next;

            double new_g1m = F->f1;

            if(new_g1m*g1 < 0){
                g1m = new_g1m;
            }

            //alocate memory for new linked list element
            struct Lines* newLine = (struct Lines*)malloc(sizeof(struct Lines));

            //link the new element
            current->next = newLine;

            //fill in line values from i1 constraints
            newLine->alpha = F->f1;
            newLine->beta = F->f3;
            newLine->next = NULL;
            newLine->prev = current;

            if (DPDUBUG) {
              printf("F Line: %fx + %f\n", F->f1, F->f3);
            }
    }

    double u1 = -INFINITY;
    double u2 = INFINITY;

    double g2 = M->h1;
    double g2m = g2;


    while(M != NULL){

            double new_g2m = F->f1;

            if(new_g2m*g2 < 0){
                g2m = new_g2m;
            }

            double new_value = -M->h3/M->h1;

            if(new_value > u1 && M->h3 < 0 ){
                u1 = new_value;
            }

            if(new_value < u2 && M->h3 > 0){
                u2 = new_value;
            }

            //increment to next i0 constraint
            M = M->next;
    }


    if (DPDUBUG) {
              printf("u1 = %f;\tu2 = %f\n", u1, u2);
        }

    double cost[2] = {1,1};

    //hihgest f as result of pairing
    double* result = pairing(first, out,u1,u2,cost);

    if (DPDUBUG) {
        printf("<====== Pairing Complete =========>\n");
        printf("Xo = %f, G(Xo) = %f\n",result[0],result[1]);
    }

    if(result[1] == -INFINITY){
        if (DPDUBUG) {
        printf("UNBOUNDED\n");
        }
        return -1;

    }

    int sigma = (g1m*g2 - g1*g2m);

    if(sigma > 0){
        sigma = 1;
    }
    if(sigma < 0){
        sigma = -1;
    }

    if(fabs(sigma) < TOLERANCE){
        sigma = 0;

    }

    if (DPDUBUG) {
        printf("sigma is %d\n", sigma);

    }

    return sigma;

}

void three_variable(double* c, double* a1, double* a2, double* a3, double* b, double* out, int num_constraints){
	/*
	BEFORE STEP 0 PERFORM TRANSFORMATION
	*/
	//starting feasible region point
	double x1_start = 0;
	double x2_start = 0;
	double x3_start = 0;

	//if objective function values are zero, the problem is already optimal
	if(c[0] == 0 && c[1] == 0 && c[2] == 0){
			out[0] = x1_start;
			out[1] = x2_start;
			out[2] = x3_start;
			return;
	}

    //print original problem
    if (DPDUBUG) {
      printf("<=============== Original ===============>\n");
      printf("max %fx1 + %fx2 + %fx3\n",c[0],c[1],c[2]);
      for(int i =0; i< num_constraints;i++){
        printf("%fx1 + %fx2 + %fx3 <= %f\n",a1[i],a2[i],a3[i],b[i]);
      }


      printf("<=============== CONVERT ===============>\n");
    }


	//allocate memory for contraints
	double *A1 = malloc(sizeof(double)* num_constraints);
    double *A2 = malloc(sizeof(double)* num_constraints);
    double *A3 = malloc(sizeof(double)* num_constraints);
    double *B = malloc(sizeof(double)* num_constraints);

	//keep track of I0 and I1 contraints
	int num_I0 = 0;
	int num_I1 = 0;

	if (DPDUBUG) {
      printf("NEW CONSTRAINTS\n");
    }

	//caluclue Ai1 and Ai2 from inputted ai1 and ai2
	for(int i =0; i< num_constraints;i++){
		A1[i] = a1[i] - c[0]*a3[i]/c[2];
		A2[i] = a2[i] - c[1]*a3[i]/c[2];
		A3[i] = a3[i]/c[2];
		B[i] = b[i] - a1[i]*x1_start + a2[i]*x2_start + a3[i]*x3_start;

		printf("%fx + %fy + %f <= %fz\n",A1[i],A2[i],A3[i], B[i]);

		if(B[i] == 0){
			num_I0++;
		}
		if(B[i] >0){
			num_I1++;
		}
	}

	 if (DPDUBUG) {
      printf("THERE ARE %d I0 CONSTRAINTS\n",num_I0);
      printf("THERE ARE %d I1 CONSTRAINTS\n",num_I1);
    }

    //alocated memory for u1, u2, and I1 contraints
    Lines_three_I0* i0_constrains = (Lines_three_I0*)malloc(sizeof(Lines_three_I0));
    i0_constrains->next = NULL;
    Lines_three_I0* I0_start = i0_constrains;
    Lines_three_I1* i1_constrains = (Lines_three_I1*)malloc(sizeof(Lines_three_I1));
    i1_constrains->next = NULL;
    Lines_three_I1* I1_start = i1_constrains;


    if (DPDUBUG) {
      printf("NEW CONSTRAINTS\n");
    }

    for(int i =0; i< num_constraints;i++){
      if(B[i] == 0){

        i0_constrains->next = (Lines_three_I0*)malloc(sizeof(Lines_three_I0));
        i0_constrains = i0_constrains->next;
        i0_constrains->next = NULL;
        i0_constrains->h1 = A1[i];
        i0_constrains->h2 = A2[i];
        i0_constrains->h3 = A3[i];

	        if (DPDUBUG) {
	          printf("%fx + %fy + %f <= 0 - I0 contraint\n", A1[i],A2[i], A3[i]);
	        }	

        }
      if(B[i] > 0){

        i1_constrains->next = (Lines_three_I1*)malloc(sizeof(Lines_three_I1));
        i1_constrains = i1_constrains->next;
        i1_constrains->next = NULL;
        i1_constrains->f1 = A1[i]/B[i];
        i1_constrains->f2 = A2[i]/B[i];
        i1_constrains->f3 = A3[i]/B[i];

        	if (DPDUBUG) {
	          printf("%fx + %fy + %f - I1 contraint\n", i1_constrains->f1, i1_constrains->f2, i1_constrains->f3);
	        }

      }
    }

	//STEP 0
    if(num_I1 == 0){
    	out[0] = x1_start;
		out[1] = x2_start;
		out[2] = x3_start;
		return;
    }

    //STEP 1
    if(num_I1 == 1){
        //run two variable problem with
        double * a1_tv = malloc(sizeof(double)* num_I0);
        double * a2_tv = malloc(sizeof(double)* num_I0);
        double * b_tv = malloc(sizeof(double)* num_I0);

        //fill constraint arrays
        i0_constrains = I0_start->next;
        i1_constrains = I1_start->next;

        //index for arrays
        int i = 0;

        while(i0_constrains != NULL){
            a1_tv[i] = i0_constrains->h1;
            a2_tv[i] = i0_constrains->h2;
            b_tv[i] = -i0_constrains->h3;

            //increment i
            i++;

            //next element in linked list
            i0_constrains = i0_constrains->next;
        }

        //variables for outputing from tv
        double tv_out[2];
        double cost[2];
        cost[0] = i1_constrains->f1;
        cost[1] = i1_constrains->f2;
        //run two variable problem
        struct tv_return* pair_parameters = two_variable(cost[0],cost[1],a1_tv,a2_tv,b_tv,num_I0);
        pairing(pair_parameters->line,tv_out,pair_parameters->u1,pair_parameters->u2,cost);

        printf("From Two Variable: x1 = %f, x2 = %f\n", tv_out[0],tv_out[1]);

        double z = cost[0]*tv_out[0] + cost[1]*tv_out[1] + i1_constrains->f3;

        out[0] = tv_out[0]/z;
        out[1] = tv_out[1]/z;
        out[2] = (-c[0]*tv_out[0] - c[1]*tv_out[1] + 1) / (c[2]*z);
        return;



    }

    //get pointers to index pairs
    Lines_three_I1* i = I1_start->next;
    Lines_three_I1* j = i->next;

    //create linked list to store pairs
    three_pairs* k = (three_pairs*)malloc(sizeof(three_pairs));
    k->next = NULL;
    three_pairs* tp_start = k;

    //keep track of the number of pairs
    int num_pairs = 0;

    //keep track of the number of i0 and i1 constraints
    int num_i0_lines = 0;
    int num_i1_lines = 0;

    while(j != NULL){
        //if parallel only keep dominating
        if(fabs(i->f1 - j->f1) < TOLERANCE && fabs(i->f2 - j->f2) < TOLERANCE){
            //keep first plane
            if(i->f3 > j->f3){
                if (DPDUBUG) {
                  printf("DELETE FIRST LINE: %fx + %fy + %f \n", i->f1, i->f2, i->f3);
                }

            }
            else{
                if (DPDUBUG) {
                  printf("DELETE SECOND LINE: %fx + %fy + %f \n", j->f1, j->f2, j->f3);
                }

            }

        }
        else{
            //calculate  critical line
            k->next = (three_pairs*)malloc(sizeof(three_pairs));
            k = k->next;
            k->next =NULL;
            k->first= i;
            k->second  = j;
            k->cl_alpha = (i->f1 - j->f1) / (j->f2 - i->f2);
            k->cl_beta = (i->f3 - j->f3) / (j->f2 - i->f2);
            num_pairs++;

            if(fabs(i->f2 - j->f2) < TOLERANCE){
                k->type = 1;
                num_i0_lines++;
            }
            else{
                k->type = 0;
                num_i1_lines++;
            }
            

            if (DPDUBUG) {
                  printf("NEW CRITICAL LINE: y = %fx  + %f - TYPE IS I%d\n", k->cl_alpha, k->cl_beta, k->type);
                }

        }

        i = j->next;

        if(i != NULL){
            j = i->next;
        }
        else
            j = NULL;

    }

    if (DPDUBUG) {
                  printf("There are %d pairs \n", num_pairs);
    }

    // go through i0 contraints and add them as lines
    i0_constrains = I0_start->next;

    while(i0_constrains != NULL){

        //calculate  critical line
        k->next = (three_pairs*)malloc(sizeof(three_pairs));
        k = k->next;
        k->next =NULL;
        k->first= NULL;
        k->second  = NULL;
        //when h2 is zero it is i1 constraint so alpha is false
        if(fabs(i0_constrains->h2) < TOLERANCE ){
            k->cl_alpha = NAN;
            k->cl_beta = i0_constrains->h3/(-i0_constrains->h1);
        }
        else{
            k->cl_alpha = i0_constrains->h1/(-i0_constrains->h2);
            k->cl_beta = i0_constrains->h3/(-i0_constrains->h2);
        }

        //set  as i0 constraint if h2 is zero and i1 contraint if h2 isn't
        if(i0_constrains->h2 == 0){
            k->type = 1;
            num_i1_lines++;
        }
        else{
            k->type = 0;
            num_i0_lines++;
        }

        if (DPDUBUG) {
                  printf("NEW LINE FROM I0 CONSTRAINT: y = %fx  + %f - TYPE IS I%d\n", k->cl_alpha, k->cl_beta, k->type);
            }

        //increment to next line
        i0_constrains = i0_constrains->next;

    }

    //make array of i0 contraints alphas
    double * alpha_list = malloc(sizeof(double)* num_i0_lines);

    int index = 0;
    k = tp_start->next;


    if (DPDUBUG) {
                  printf("I0 alphas ");
            }

    //populate list of i0 alpha values
      while (k != NULL) {
        if(k->type == 0){
            alpha_list[index++] = k->cl_alpha;
            if (DPDUBUG) {
                  printf("%f ", k->cl_alpha);
            }
        }
          
          k = k->next;
      }

      if (DPDUBUG) {
                  printf("\n");
            }

    //calculate median
    double median;

    if(num_i0_lines > 1){
        //use linear median algorithm when k > 1
        median = quickSelect(alpha_list, num_i0_lines, num_i0_lines/2, &pickPivot);
    }
    else{
        //when thereis one value no it is assigned median
        median = alpha_list[0];
    }

    if (DPDUBUG) {
        printf("median is %.12f\n", median);
    }

    //storage for i2 and I3 and I4 constraints
    three_pairs* I2_partition = (three_pairs*)malloc(sizeof(three_pairs));
    I2_partition->next = NULL;
    three_pairs* I2_partition_start = I2_partition;

    three_pairs* I3_partition = (three_pairs*)malloc(sizeof(three_pairs));
    I3_partition->next = NULL;
    three_pairs* I3_partition_start = I3_partition;

    three_pairs* I4_partition = (three_pairs*)malloc(sizeof(three_pairs));
    I4_partition->next = NULL;
    three_pairs* I4_partition_start = I4_partition;

    int num_i2_partition = 0;
    int num_i3_partition = 0;
    int num_i4_partition = 0;

    k = tp_start->next;


    if (DPDUBUG) {
        printf("Partioning I0 into I2, I3, and I4\n");
    }

    //populate list of i0 alpha values
      while (k != NULL) {
        if(k->type == 0){   
            //partition based on alpha value
            if(fabs(k->cl_alpha-median) < TOLERANCE){

                I2_partition->next = (three_pairs*)malloc(sizeof(three_pairs));
                I2_partition = I2_partition->next;
                I2_partition->next =NULL;
                I2_partition->first= k->first;
                I2_partition->second  = k->second;
                I2_partition->cl_alpha = k->cl_alpha;
                I2_partition->cl_beta = k->cl_beta;
                I2_partition->type = 2;

                num_i2_partition++;

                if (DPDUBUG) {
                  printf("NEW LINE FROM I2 CONSTRAINT: y = %fx  + %f - TYPE IS I%d\n", I2_partition->cl_alpha, I2_partition->cl_beta, I2_partition->type);
                }


            }
            if(k->cl_alpha < (median - TOLERANCE)){

                I3_partition->next = (three_pairs*)malloc(sizeof(three_pairs));
                I3_partition = I3_partition->next;
                I3_partition->next =NULL;
                I3_partition->first= k->first;
                I3_partition->second  = k->second;
                I3_partition->cl_alpha = k->cl_alpha;
                I3_partition->cl_beta = k->cl_beta;
                I3_partition->type = 3;

                num_i3_partition++;

                if (DPDUBUG) {
                  printf("NEW LINE FROM I3 CONSTRAINT: y = %fx  + %f - TYPE IS I%d\n", I3_partition->cl_alpha, I3_partition->cl_beta, I3_partition->type);
                }

            }

            if(k->cl_alpha > (median + TOLERANCE)){

                I4_partition->next = (three_pairs*)malloc(sizeof(three_pairs));
                I4_partition = I4_partition->next;
                I4_partition->next =NULL;
                I4_partition->first= k->first;
                I4_partition->second  = k->second;
                I4_partition->cl_alpha = k->cl_alpha;
                I4_partition->cl_beta = k->cl_beta;
                I4_partition->type = 4;

                num_i4_partition++;

                if (DPDUBUG) {
                  printf("NEW LINE FROM I4 CONSTRAINT: y = %fx  + %f - TYPE IS I%d\n", I4_partition->cl_alpha, I4_partition->cl_beta, I4_partition->type);
                }

            }
        }
            
          k = k->next;
      }

    //store pairs in q set
    q_pairs* q = (q_pairs*)malloc(sizeof(q_pairs));
    q->next = NULL;
    q_pairs* q_start = q;

    //start from beginning of I3 and I4
    I3_partition = I3_partition_start->next;
    I4_partition = I4_partition_start->next;

    while(I3_partition != NULL && I4_partition != NULL){

        //new data for q pair
        q->next = (q_pairs*)malloc(sizeof(q_pairs));
        q= q->next;
        q->next = NULL;
        q->i3 = I3_partition;
        q->i4 = I4_partition;
        q->u = (I3_partition->cl_beta - I4_partition->cl_beta) / (I4_partition->cl_alpha - I3_partition->cl_alpha);
        q->v = (I3_partition->cl_beta*I4_partition->cl_alpha - I4_partition->cl_beta*I3_partition->cl_alpha) / (I4_partition->cl_alpha - I3_partition->cl_alpha);

        if (DPDUBUG) {
            printf("NEW Q PAIR, I3: y = %fx  + %f I4: y = %fx  + %f point(%f,%f) \n", q->i3->cl_alpha, q->i3->cl_beta, q->i4->cl_alpha, q->i4->cl_beta, q->u, q->v);
        }

        I3_partition = I3_partition->next;
        I4_partition = I4_partition->next;

        num_i1_lines++;

    }

    //storage for I1 u vlaues
    double * i1_us = malloc(sizeof(double)* num_i1_lines);

    //populate i1 lines list
    k = tp_start->next;
    index = 0;

    if (DPDUBUG) {
                  printf("I1 u's ");
            }

    //populate list of i1 u values from i1 contraints
    while (k != NULL) {
        if(k->type == 1){
            i1_us[index++] = k->cl_beta;
            if (DPDUBUG) {
                  printf("%f ", k->cl_beta);
            }
        }
          
        k = k->next;
    }

    q = q_start->next;

    //populate list of i1 u values from q
    while (q != NULL) {
            i1_us[index++] = q->u;
            if (DPDUBUG) {
                  printf("%f ", q->u);
            }
          
        q = q->next;
    }

    if (DPDUBUG){
            printf("\n");
        }

    //calculate median
    double u_median;

    if(num_i1_lines > 1){
        //use linear median algorithm when k > 1
        u_median = quickSelect(i1_us, num_i1_lines, num_i1_lines/2, &pickPivot);
    }
    else{
        //when thereis one value no it is assigned median
        u_median = i1_us[0];
    }

    if (DPDUBUG) {
        printf("u_median is %.12f\n", u_median);
    }


    //check if x0 is found on line
    if(line_search(I1_start, I0_start, 1,0, u_median, 0)){

        q = q_start->next;

        //delete q values lower than u_ median
        while (q != NULL) {

            if(q->u < u_median + TOLERANCE){
                if (DPDUBUG) {
                      printf("Q DELETED: %f\n", q->u);
                  }

                I2_partition->next = (three_pairs*)malloc(sizeof(three_pairs));
                I2_partition = I2_partition->next;
                I2_partition->next =NULL;
                I2_partition->first= NULL;
                I2_partition->second  = NULL;
                I2_partition->cl_alpha = median;
                I2_partition->cl_beta = q->v - (median*q->u);
                I2_partition->type = 2;

                num_i2_partition++;


                if (DPDUBUG) {
                  printf("NEW LINE FROM I2 CONSTRAINT: y = %fx  + %f - TYPE IS I%d\n", I2_partition->cl_alpha, I2_partition->cl_beta, I2_partition->type);
                }
            }
        }
              
        q = q->next;
    }
    else{

    }

    //calculated median of I2 y intercepts
    index = 0;
    I2_partition = I2_partition_start->next;
    double y_intercepts[num_i2_partition];


    if (DPDUBUG) {
                  printf("I2 y intercepts ");
            }

    //populate list of i0 alpha values
      while (k != NULL) {
        if(k->type == 0){
            y_intercepts[index++] = I2_partition->cl_beta;
            if (DPDUBUG) {
                  printf("%f ", I2_partition->cl_beta);
            }
        }
          
          I2_partition = I2_partition->next;
      }

      if (DPDUBUG) {
                  printf("\n");
            }

    //calculate median
    double yi_median;

    if(num_i2_partition > 1){

        //use linear median algorithm when k > 1
        yi_median = quickSelect(y_intercepts, num_i2_partition, num_i2_partition/2, &pickPivot);
    }
    else{

        //when thereis one value no it is assigned median
        median = y_intercepts[0];
    }

    if (DPDUBUG) {
        printf("y intercept median is %.12f\n", yi_median);
    }

    //line search y intercept median
    if(line_search(I1_start, I0_start, median,1, yi_median, 0)){

        q = q_start->next;

        //delete q values lower than u_ median
        while (q != NULL) {

            if(q->u < yi_median + TOLERANCE){
                if (DPDUBUG) {
                      printf("Q DELETED: %f\n", q->u);
                  }
            }
        }
              
        q = q->next;
    }
    else{

    }





}