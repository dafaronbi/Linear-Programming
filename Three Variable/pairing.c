#include "three variable.h"
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




	struct Lines* make_lines(double in_alphas[], double in_betas[],int a_size){

    struct Lines* current = (struct Lines*)malloc(sizeof(struct Lines));
    struct Lines* first = current;

    double *alphas = malloc(sizeof(double)*a_size);
    double *betas = malloc(sizeof(double)*a_size);

    for(int i = 0; i < a_size; i++){
      alphas[i] = in_alphas[a_size-i-1];
      betas[i] = in_betas[a_size-i-1];
    }

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

    free(alphas);
    free(betas);

    return first;

  }
	
double* pairing(struct Lines* line,double out[],double u1, double u2,double c[]){

  //save starting address of line
  struct Lines* start = line;

  double * result  = (double *)malloc(2*sizeof(double));

  //step 0
  //find minimum in alpha
  if (DPDUBUG) {
    printf("<=============== STEP 0 ===============>\n");
  }

  // out zero when u1 > u2
  if(u1 > u2 + TOLERANCE_1){
    if(DPDUBUG){
      printf("u1 > u2\n");
    }
    out[0] = 0;
    out[1] = 0;

    //free memory of lines
    struct Lines* lines_tmp;
    line = start;
    while(line != NULL){
      lines_tmp = line;
      line = line->next;
      free(lines_tmp);
    }
    getchar();

    result[0] = 0;
    result[1] = 0;

    return result;
  }

  double f_of_u1 = -INFINITY;

  //determine f(u1)
  line = start;
  while(line != NULL){
    double new_value = line->alpha*u1 +line->beta;
    if(new_value > f_of_u1)
      f_of_u1 = new_value;
    line = line->next;
    }

  if (DPDUBUG) {
    printf("f(u1) = %.20f\n", f_of_u1);
  }


  //out u1 and F(u1) when u1 = u2
  if(fabs(u1-u2) < TOLERANCE_1){

    if(DPDUBUG){
      printf("u1 = u2\n");
    }
    out[0] = u1/f_of_u1;
    out[1] = (1-c[0]*u1)/(c[1]*f_of_u1);

    //free memory of lines
    struct Lines* lines_tmp;
    line = start;
    while(line != NULL){
      lines_tmp = line;
      line = line->next;
      free(lines_tmp);
    }

    result[0] = u1;
    result[1] = f_of_u1;

    return result;
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
    if(u_row < line->alpha -TOLERANCE_1  && (fabs(f_of_u1 - y) < TOLERANCE_1 )){
      u_row = line->alpha;
      if (DPDUBUG) {
        printf("possible u_row = %f\n",u_row);
      }
    }
    line = line->next;
  }


  if (DPDUBUG) {
    printf("u_row = %f\n", u_row);
  }

  //terminate when u_row > 0
  if(u_row > -TOLERANCE_1){
    out[0] = u1/f_of_u1;
    out[1] = (1-c[0]*u1)/(c[1]*f_of_u1);

    //free memory of lines
    struct Lines* lines_tmp;
    line = start;
    while(line != NULL){
      lines_tmp = line;
      line = line->next;
      free(lines_tmp);
    }

    result[0] = u1;
    result[1] = f_of_u1;

    return result;
  }


  //determine f(u2)
  double f_of_u2 = -INFINITY;

  line = start;
  while(line != NULL){
    double new_value = line->alpha*u2 +line->beta;
    if(new_value > f_of_u2 + TOLERANCE_1)
      f_of_u2 = new_value;
    line = line->next;
    }

  if (DPDUBUG) {
    printf("f(u2) = %f\n", f_of_u2);
  }


  //find minum in alpha at u2
  line = start;
  while(line != NULL){
    //calculate y of current line
    double y = line->alpha*u2 + line->beta;
    if(line->alpha < u_lambda && (fabs(f_of_u2 - y) < TOLERANCE_1)){
      u_lambda = line->alpha;
      if (DPDUBUG) {
        printf("possible u_lambda = %f\n",u_lambda);
      }
    }
    line = line->next;
  }

  if (DPDUBUG) {
    printf("u_lambda = %f\n", u_lambda);
  }

  if(u_lambda < TOLERANCE_1){
    out[0] = u2/f_of_u2;
    out[1] = (1-c[0]*u2)/(c[1]*f_of_u2);

    //free memory of lines
    struct Lines* lines_tmp;
    line = start;
    while(line != NULL){
      lines_tmp  = line;
      line = line->next;
      free(lines_tmp);
    }

    result[0] = u2;
    result[1] = f_of_u2;

    return result;
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

  if (DPDUBUG) {
    printf("Minimum alpha is %f\n", a_min);
  }
  //stop if minimum value in alpha is greater than zero
  if (a_min > 0){

    //x0 = -infinity
    out[0] = -INFINITY;

    //f(x0) = -infinity
    out[1] = -INFINITY;

    //free memory of lines
    struct Lines* lines_tmp;
    line = start;
    while(line != NULL){
      lines_tmp  = line;
      line = line->next;
      free(lines_tmp);
    }

    result[0] = -INFINITY;
    result[1] = -INFINITY;

    return result;
  }

  //start from beginning
  line = start;


  //find maximum in alpha
  while(line != NULL){
    if(line->alpha > a_max)
      a_max = line->alpha;
    line = line->next;
  }

  if (DPDUBUG) {
    printf("Maximum alpha is %f\n", a_max);
  }

  //stop if maximum value in alpha is less than zero
  if (a_max < 0){
    //x0 = infinity
    out[0] = INFINITY;
    //f(x0) = infinity MAYBE WRONG CHECK DETAILS
    out[1] = INFINITY;

    //free memory of lines
    struct Lines* lines_tmp;
    line = start;
    while(line != NULL){
      lines_tmp  = line;
      line = line->next;
      free(lines_tmp);
    }

    result[0] = INFINITY;
    result[1] = INFINITY;

    return result;
  }

  //variable to store number of x
  int k;

  while(1){
    //step 1
    if (DPDUBUG) {
      printf("<=============== STEP 1 ===============>\n");
    }
    //linked list initiallized to zero
    struct Node* x = (struct Node*)malloc(sizeof(struct Node));
    x->next = NULL;
    struct Node* x_start = x;

    //creat linked list for line pair addresses
    struct Pairs* line_pairs = (struct Pairs*)malloc(sizeof(struct Pairs));
    line_pairs->next = NULL;
    struct Pairs* line_pairs_start = line_pairs;

    //restart from beginning
    line = start;
    do{

      //make pointers to acces pairs
      struct Lines* i = start;
      struct Lines* j = start->next;

      while(j != NULL){
        //delete parallel lines
        if(fabs(i->alpha - j->alpha) < TOLERANCE_1){
          //display that a line will be deleted
          if (DPDUBUG) {
              printf("DELETE: ");
            }
          if(i->beta > j->beta){
            //display deleted line
            if (DPDUBUG) {
                printf("line: a = %f, b = %f\n", j->alpha, j->beta);
              }

            //delete line
            j = deleteLine(start,j);

            //go to next pair
            i = j;
            if(i != NULL)
              j = i->next;
            else
              j = NULL;

          }
          else{

            //display deleted line
            if (DPDUBUG) {
                printf("line: a = %f, b = %f\n", i->alpha, i->beta);
              }


            //if deleting first line get a new start
            if(i == start){
              i = deleteLine(start,i);
              start = i;

              if (DPDUBUG) {
                  printf("NEW Start Line: line: address %d, a = %f, b = %f\n", i, i->alpha, i->beta);
                  }
            }
            else{
              i = deleteLine(start, i);
            }

            //go to next pair
            i = j->next;
            if(i != NULL)
              j = i->next;
            else
              j = NULL;

          }
          continue;
        }

        //print current pair
        if (DPDUBUG) {
              printf("line1: a = %f, b = %f\n", i->alpha, i->beta);
              printf("line2: a = %f, b = %f\n", j->alpha, j->beta);
            }


        //get intersection value
        double value = (i->beta - j->beta)/(j->alpha - i->alpha);

        //check if value is greater than u2
        if(value >= u2){
              if(j->alpha > i->alpha){
                //display deleted line
                if (DPDUBUG) {
                  printf("OUT OF U2 DELETE: line: a = %f, b = %f x = %f\n", j->alpha, j->beta,value);
                }

                j = deleteLine(start, j);
                
                //go to next pair
                i = j;
                if(i != NULL)
                  j = i->next;
                else
                  j = NULL;
              }
              else{
                //display deleted line
                if (DPDUBUG) {
                  printf("OUT OF U2 DELETE: line: a = %f, b = %f x = %f\n", i->alpha, i->beta, value);
                }


                if(i == start){
                  //if deleting first line get a new start
                  i = deleteLine(start, i);
                  start = i;

                  if (DPDUBUG) {
                  printf("NEW Start Line: line: address %d, a = %f, b = %f\n", i, i->alpha, i->beta);
                  }

                }
                else{
                  i = deleteLine(start, i);
                }

                //go to next pair
                i = j->next;
                if(i != NULL)
                  j = i->next;
                else
                  j = NULL;
              }
              continue;
        }

        //check if value is less than u1
        if(value <= u1){
          if(j->alpha > i->alpha){
                if (DPDUBUG) {
                  printf("OUT OF U1 DELETE: line: a = %f, b = %f x = %f\n", i->alpha, i->beta,value);
                }

                if(i == start){
                  //if deleting first line get a new start
                  i = deleteLine(start, i);
                  start = i;

                  if (DPDUBUG) {
                  printf("NEW Start Line: line: address %d, a = %f, b = %f\n", i, i->alpha, i->beta);
                  }
                }
                else{
                  i = deleteLine(start, i);
                }

                //go to next pair
                i = j->next;
                if(i != NULL)
                  j = i->next;
                else
                  j = NULL;
                }
              else{
                if (DPDUBUG) {
                  printf("OUT OF U1 DELETE: line: a = %f, b = %f x = %f\n", j->alpha, j->beta,value);
                }

                j = deleteLine(start, j);
                
                //go to next pair
                i = j;
                if(i != NULL)
                  j = i->next;
                else
                  j = NULL;
              }
              continue;
        }

        //assign values to line pairs
        line_pairs->next = (struct Pairs*)malloc(sizeof(struct Pairs));
        line_pairs->next->first = i;
        line_pairs->next->second = j;
        line_pairs->next->intersection = value;
        line_pairs->next->next = NULL;
        line_pairs = line_pairs->next;

        x->next = (struct Node*)malloc(sizeof(struct Node));
        x->next->data = value;
        x->next->next = NULL;
        x = x->next;
        if (DPDUBUG) {
          printf("x = %f\n",x->data);
        }

        //go to next pair
        i = j->next;
        if(i != NULL){
          j = i->next;
        }
        else{
          j = NULL;
        }

      }

      //Step 2 count k, the number of xi created
      if(DPDUBUG){
        line_pairs = line_pairs_start->next;

        while(line_pairs != 0){
          printf("Line 1: a = %f, b = %f, Line 2: a = %f, b = %f, intersection: %f\n", line_pairs->first->alpha, line_pairs->first->beta, line_pairs->second->alpha, line_pairs->second->beta, line_pairs->intersection);
          line_pairs = line_pairs->next;
        }
      }

      if (DPDUBUG) {
        printf("<=============== STEP 2 ===============>\n");
      }
      k = 0;
      x = x_start->next;
      while(x != NULL){
        k += 1;
        x = x->next;

      }

      //TODO Check k
      if (DPDUBUG) {
        printf("k is %d\n",k);
      }

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
      if (DPDUBUG) {
        printf("size of x_array is %lu\n", sizeof(x_arr)/sizeof(x_arr[0]));

        for(int i = 0; i < sizeof(x_arr)/sizeof(x_arr[0]); i++){
          printf("x%d is %.12f\n",i,x_arr[i]);
        }
      }

      //calculate median
      double x_median;

      if(k > 1){

          //use linear median algorithm when k > 1
          x_median = quickSelect(x_arr, k, k/2, &pickPivot);
      }
      else{
        //when thereis one value no it is assigned median
        x_median = x_arr[0];
      }

      if (DPDUBUG) {
        printf("median is %.12f\n", x_median);

        //Step 3
        printf("<=============== STEP 3 ===============>\n");
      }

      //Start from beginning
      line = start;

      double f_max = -INFINITY;

      //find maximum f(xmean)
      while(line != NULL){
        double new_value = line->alpha*x_median +line->beta;
        if(new_value > f_max ){
          f_max = new_value;
          if (DPDUBUG) {
            printf("possible f_max is %f from line a=%f b=%f\n",f_max,line->alpha, line->beta);
          }
        }
        line = line->next;

      }

      if (DPDUBUG) {
        printf("f_max is %.12f\n",f_max);
      }



      //start from beginning
      line = start;
 

      //calculate left and right gradients
      double lambda = INFINITY;
      double row = -INFINITY;

      //find minimum in alpha
      while(line != NULL){
        double y = line->alpha*x_median + line->beta;
        if(line->alpha < lambda && (fabs(f_max - y) < TOLERANCE_1)){
          lambda = line->alpha;
          if (DPDUBUG) {
            printf("possible lambda = %f\n", lambda);
          }
        }
        line = line->next;
      }

      if (DPDUBUG) {
        printf("lambda = %f\n", lambda);
      }

      //start from beginning
      line = start;

      //find maximum in alpha
      while(line != NULL){
        double y = line->alpha*x_median + line->beta;
        if (DPDUBUG) {
          printf("posiible y  = %.12f alpha is %.12f\n", y, line->alpha);
        }
        if(line->alpha > row && (fabs(f_max - y) < TOLERANCE_1)){
          row = line->alpha;
          if (DPDUBUG) {
            printf("possible row = %f\n", row);
          }
        }
        line = line->next;
      }
      if (DPDUBUG) {
        printf("row = %f\n", row);

        // step 4 exit if optimatl is found
        printf("<=============== STEP 4 ===============>\n");
      }

      if( lambda <= TOLERANCE_1 && row >= -TOLERANCE_1){
        if (DPDUBUG) {
          printf("NO NEED TO DELETE\n");
        }
        out[0] = x_median/f_max;
        out[1] = (1-c[0]*x_median)/(c[1]*f_max);

        //free memory of x values
        struct Node* x_tmp;
        x = x_start;
        while(x != NULL){
          x_tmp = x;
          x = x->next;
          free(x_tmp);
        }

        //free memory of pairs
        struct Pairs* pairs_tmp;
        line_pairs = line_pairs_start;
        while(line_pairs != NULL){
          pairs_tmp = line_pairs;
          line_pairs = line_pairs->next;
          free(pairs_tmp);
        }

        //free memory of lines
        struct Lines* lines_tmp;
        line = start;
        while(line != NULL){
          lines_tmp = line;
          line = line->next;
          free(lines_tmp);
        }

	    result[0] = x_median;
	    result[1] = f_max;

	    return result;
      }

       if(lambda > TOLERANCE_1){

        line_pairs = line_pairs_start->next;

        while(line_pairs != NULL){
          if(line_pairs->intersection >= x_median){
            //print current lines
            if (DPDUBUG) {
                printf("line1: a = %f, b = %f\n", line_pairs->first->alpha, line_pairs->first->beta);
                printf("line2: a = %f, b = %f\n", line_pairs->second->alpha, line_pairs->second->beta);
              }

            //delete the "right half" or value with greater alpha
            if(line_pairs->first->alpha > line_pairs->second->alpha){
              if (DPDUBUG) {
                  printf("LAMBDA DELETE: line: a = %f, b = %f x = %.12f\n", line_pairs->first->alpha, line_pairs->first->beta,line_pairs->intersection);
                }

              //if start is deleted make new start
              if(line_pairs->first == start){
                start = line_pairs->first->next;
              }

              //delete line
              deleteLine(start, line_pairs->first);
            }
            else{
              if (DPDUBUG) {
                  printf("LAMBDA DELETE: line: a = %f, b = %f x = %.12f\n", line_pairs->second->alpha, line_pairs->second->beta,line_pairs->intersection);
                }


              deleteLine(start, line_pairs->second);

            }
          }

            line_pairs = line_pairs->next;
        }

      }

      if(row < -TOLERANCE_1){

        line_pairs = line_pairs_start->next;

        while(line_pairs != NULL){
          if(line_pairs->intersection <= x_median){
            //print current lines
            if (DPDUBUG) {
                printf("line1: a = %f, b = %f\n", line_pairs->first->alpha, line_pairs->first->beta);
                printf("line2: a = %f, b = %f\n", line_pairs->second->alpha, line_pairs->second->beta);
              }

            //delete the "right half" or value with greater alpha
            if(line_pairs->first->alpha < line_pairs->second->alpha){
              if (DPDUBUG) {
                  printf("ROW DELETE: line: a = %f, b = %f x = %.12f\n", line_pairs->first->alpha, line_pairs->first->beta,line_pairs->intersection);
                }

              //if start is deleted make new start
              if(line_pairs->first == start){
                start = line_pairs->first->next;
              }

              //delete line
              deleteLine(start, line_pairs->first);
            }
            else{
              if (DPDUBUG) {
                  printf("ROW DELETE: line: a = %f, b = %f x = %.12f\n", line_pairs->second->alpha, line_pairs->second->beta,line_pairs->intersection);
                }


              deleteLine(start, line_pairs->second);

            }
          }

            line_pairs = line_pairs->next;
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
    if (DPDUBUG) {
      printf("<=============== Original ===============>\n");
      printf("max %fx1 + %fx2\n",c1,c2);
      for(int i =0; i< num_contraints;i++){
        printf("%fx1 + %fx2 <= %f\n",a1[i],a2[i],b[i]);
      }


      printf("<=============== CONVERT ===============>\n");
    }


    //optimal if c1 and c2 coefficients are zero
    if(c1 == 0 && c2 == 0){
      output->line = NULL;
      output->u1 = 0;
      output->u2 = 0;
      return output;
    }


    //allocate memory for contraints
    double *A1 = malloc(sizeof(double)* num_contraints);
    double *A2 = malloc(sizeof(double)* num_contraints);
    double *B = malloc(sizeof(double)* num_contraints);


    //keep track of I0 and I1 contraints
    int num_I0 = 0;
    int num_I1 = 0;

    if (DPDUBUG) {
      printf("NEW CONSTRAINTS\n");
    }

    //caluclue Ai1 and Ai2 from inputted ai1 and ai2
    for(int i =0; i< num_contraints;i++){
      A1[i] = a1[i] - a2[i]*c1/c2;
      A2[i] = a2[i]/c2;
      B[i] = b[i] - a1[i]*x1_start -a2[i]*x2_start;


      if (DPDUBUG) {
        printf("%fx + %f <= %f\n",A1[i],A2[i],B[i]);
      }

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
    double u1 = -INFINITY;
    double u2 = INFINITY;
    double *alphas = malloc(sizeof(double)* num_I1);
    double *betas = malloc(sizeof(double)* num_I1);

    //keep track of elemtents in alphas and betas array
    int I1_index = 0;

    if (DPDUBUG) {
      printf("I1 CONSTRAINTS\n");
    }

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

        if (DPDUBUG) {
          printf("y => %fx + %f\n", alphas[I1_index],betas[I1_index]);
        }
        I1_index++;
      }
    }
    if (DPDUBUG) {
      printf("u1 = %.30f\nu2 = %.30f\n",u1,u2);
    }


    //make lines from generated constraints
    struct Lines* line = make_lines(alphas,betas,num_I1);



    output->line = line;
    output->u1 = u1;
    output->u2 = u2;

    free(A1);
    free(A2);
    free(B);
    free(alphas);
    free(betas);

    return output;

  }