                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           /* main program
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 *
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"
#include "memory.h"
#include "pre.h"

char    infile[200];
char    outfile[200];
char    sinfile[200];
char    P_name[100];

void usage(char *argv[]) {
  printf("Usage:\n");
  printf("\t%s mpsfile\n", argv[0]);
  printf("\t%s -s specification_file mpsfile\n", argv[0]);
}

LPtype* makePCXProblem(double c1, double c2,double a1 [], double a2 [], double b [], int num_contraints){

    // make the contraint matrixes
    // double a[num_contraints][2+num_contraints];
    // double Ta[2 + num_contraints][num_contraints];

    // double **a = (double **)calloc(num_contraints * sizeof(double *)); 
    // for (int i=0; i<num_contraints; i++) 
    //      a[i]= (double *)calloc(2+num_contraints * sizeof(double)); 

    // double **Ta = (double **)calloc(2+num_contraints * sizeof(double *)); 
    // for (int i=0; i<2+num_contraints; i++) 
    //      Ta[i] = (double *)calloc(num_contraints * sizeof(double)); 

    




    
    // //fill a matrix
    // for(int i = 0; i < num_contraints; i++){
    //   for(int j = 0; j < 2+num_contraints; j++){
    //     if(j == 0){
    //       a[i][j] = a1[i];
    //     }
    //     else if(j== 1){
    //       a[i][j] = a2[i];
    //     }
    //     else{
    //       if(j-2 == i)
    //         a[i][j] = 1;
    //       else
    //         a[i][j] = 0;
    //     }

    //   }
    // }

    // //fill transpose matrix
    // for(int i = 0; i < 2+num_contraints; i++){
    //   for(int j = 0; j < num_contraints; j++){
    //       // Ta[i][j] = a[j][i];

    //       if(a[j][i] != 0){
    //         num_ents++;
    //       }
    //   }
    // } 

    // for (int i = 0; i < num_contraints; i++)
    //     free(a[i]);
    // free(a);  

    // for (int i = 0; i < 2+num_contraints; i++)
    //     free(Ta[i]);
    // free(Ta); 

    //1 0 1 0 0
    //0 2 0 1 0
    //3 2 0 0 1

    //1 0 3
    //0 2 2
    //1 0 0
    //0 1 1
    //0 0 1
    int num_ents = 0;
    for(int i = 0; i < num_contraints; i++){

      if(a1[i]!= 0)
        num_ents++;

      if(a2[i] != 0)
        num_ents++;
    } 

    num_ents += num_contraints;

    // LPtype* LP = NewLP(num_contraints, 2+num_contraints,num_ents);
    // printf("num rows is %d\n", LP->Rows);

    int Rows = num_contraints;
    int Cols = num_contraints + 2;
    int Ents = num_ents;

    LPtype         *LP;
  

   LP = (LPtype *) Malloc(sizeof(LPtype), "LP");

   LP->Atranspose.pBeginRow = NewInt(Rows, "LP->Atranspose.pBeginRow");
   LP->Atranspose.pEndRow = NewInt(Rows, "LP->Atranspose.pEndRow");
   LP->Atranspose.Row = NewInt(Ents, "LP->Atranspose.Row");
   LP->Atranspose.Value = NewDouble(Ents, "LP->Atranspose.Value");
   
   LP->Rows = Rows;

   LP->A.NumRows = Rows;
   LP->Atranspose.NumRows = Cols;
   LP->Cols = Cols;
   LP->A.NumCols = Cols;
   LP->Atranspose.NumCols = Rows;
   LP->Ents = Ents;
   LP->A.Nonzeros = Ents;
   LP->Atranspose.Nonzeros = Ents;
   LP->cshift = 0.0;
   
   LP->b = NewDouble(Rows, "LP->b");
   LP->c = NewDouble(Cols, "LP->c");

   LP->VarType = NewInt(Cols, "LP->BoundType");
   LP->UpBound = NewDouble(Cols, "LP->UpBound");
   LP->NumberBounds = 0;
   LP->BoundIndex = NewInt(Cols, "LP->BoundIndex");
   LP->NumberFree = 0;
   LP->FreeIndex = NewInt(Cols, "LP->FreeIndex");
   
   LP->A.pBeginRow = NewInt(Cols, "LP->A.pBeginRow");
   LP->A.pEndRow = NewInt(Cols, "LP->A.pEndRow");
   LP->A.Row = NewInt(Ents, "LP->A.Row");
   LP->A.Value = NewDouble(Ents, "LP->A.Value");
   
   LP->FreePlus  = NULL;
   LP->FreeMinus = NULL;
   LP->NumberSplit = 0;
   
   LP->ColScale  = NULL;
   LP->RowScale  = NULL;

    LP->Rows = num_contraints;
    LP->Cols = 2 + num_contraints;
    LP->Ents = num_ents;
    LP->Atranspose.NumRows = LP->Cols;
    LP->Atranspose.NumCols = LP->Rows;


    

    int tEntCount = 0;
    for(int j = 0; j < num_contraints; j++){
      LP->Atranspose.pBeginRow[j] = tEntCount + 1;
      for(int i = 0; i < 2+num_contraints; i++){

          if(i == 0){
            if(a1[j] != 0){
              LP->Atranspose.Value[tEntCount] = a1[j];
              LP->Atranspose.Row[tEntCount] = i+1;
              tEntCount++;

            }
          }
          else if(i== 1){
            if(a2[j] != 0){
              LP->Atranspose.Value[tEntCount] = a2[j];
              LP->Atranspose.Row[tEntCount] = i+1;
              tEntCount++;
            }
          }
          else{
            if(i-2 == j){
              LP->Atranspose.Value[tEntCount] = 1;
              LP->Atranspose.Row[tEntCount] = i+1;
              tEntCount++;

            }
          }

      }
      LP->Atranspose.pEndRow[j] = tEntCount;
 

    }


    LP->Atranspose.Nonzeros = num_ents;

    
    // //1 0 1 0 0
    // //0 2 0 1 0
    // //3 2 0 0 1

    LP->A.NumRows = LP->Rows;
    LP->A.NumCols = LP->Cols;
    LP->A.Nonzeros = num_ents;

    tEntCount = 0;
    for(int j = 0; j < 2+num_contraints; j++){
      LP->A.pBeginRow[j] = tEntCount + 1;
      for(int i = 0; i < num_contraints; i++){
        if(j == 0){
            if(a1[i] != 0){
              LP->A.Value[tEntCount] = a1[i];
              LP->A.Row[tEntCount] = i+1;
              tEntCount++;

            }
          }
          else if(j== 1){
            if(a2[i] != 0){
              LP->A.Value[tEntCount] = a2[i];
              LP->A.Row[tEntCount] = i+1;
              tEntCount++;
            }
          }
          else{
            if(j-2 == i){
              LP->A.Value[tEntCount] = 1;
              LP->A.Row[tEntCount] = i+1;
              tEntCount++;

            }
          }
      }
      LP->A.pEndRow[j] = tEntCount;


    }


    for(int i = 0; i < num_contraints; i++){
      LP->b[i] = b[i];
    }


    //objective function values are negative for to make into maximization problem
    LP->c[0] = c1*-1;
    LP->c[1] = c2*-1;

    for(int i = 2; i < 2 + num_contraints; i++){
      LP->c[i] = 0;
    }



    LP->cshift = 0;

    LP->VarType[0] = 0;
    LP->VarType[1] = 0;
    LP->VarType[2] = 0;
    LP->VarType[3] = 0;
    LP->VarType[4] = 0;

    LP->UpBound[0] = 0;
    LP->UpBound[1] = 0;
    LP->UpBound[2] = 0;
    LP->UpBound[3] = 0;
    LP->UpBound[4] = 0;

    LP->NumberBounds = 0;

    for(int i = 0; i < LP->Cols; i ++){
      LP->BoundIndex[i] = 0;
    }


    for(int i = 0; i < LP->Cols; i ++){
      LP->FreeIndex[i] = 0;
    }

    return LP;
}

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

double getLeftSlope(struct _twoDline* linehead, double x, double *y) {
  double maxy=-INFINITY, fx, slope=INFINITY, diff;
  struct _twoDline *tmp = linehead;
  struct _twoDline *maxl = NULL;

  while (tmp->next != NULL) {
    tmp = tmp->next;
    fx = tmp->alpha*x+tmp->beta;

    if (DPDUBUG) {
      fprintf(stdout, "line %d, alpha=%.12f, beta=%.12f, fx=%.12f, maxy=%.12f; ", tmp->orig, tmp->alpha, tmp->beta, fx, maxy);
    }

    if (maxl != NULL) {
      diff = maxl->a1*tmp->b*x + maxl->a2*tmp->b - tmp->a1*maxl->b*x - tmp->a2*maxl->b;

      if (diff < -1e-10) {
        slope = tmp->alpha;
        maxy = fx;
        maxl = tmp;
        if (DPDUBUG) {
          fprintf(stdout, "maxy < fx, slope = alpha = %.12f\n", slope);
        }

      } else if (fabs(diff) < 1e-10) {
        if (tmp->a1*maxl->b < maxl->a1*tmp->b) {
          slope = tmp->alpha;
          maxl = tmp;
          if (DPDUBUG) {
            fprintf(stdout, "maxy = fx, slope <alpha, slope <- alpha=%.12f\n", slope);
          }
        }
      } else {
        if (DPDUBUG) {
          fprintf(stdout, "\n");
        }
      }
    } else {
      slope = tmp->alpha;
      maxy = fx;
      maxl = tmp;
    }
  }
  (*y) = maxy;

  return slope;
}

double getRightSlope(struct _twoDline *linehead, double x, double *y) {
  double maxy=-INFINITY, fx, slope=-INFINITY, diff;
  struct _twoDline *tmp = linehead;
  struct _twoDline *maxl = NULL;
  while (tmp->next != NULL) {
    tmp = tmp->next;
    fx = tmp->alpha*x+tmp->beta;
    if (DPDUBUG) {
      fprintf(stdout, "line %d, alpha=%.12f, beta=%.12f, fx=%.12f, maxy=%.12f; ", tmp->orig, tmp->alpha, tmp->beta, fx, maxy);
    }
    if (maxl != NULL) {
      diff = maxl->a1*tmp->b*x + maxl->a2*tmp->b - tmp->a1*maxl->b*x - tmp->a2*maxl->b;

      if (diff < -1e-10) {
        slope = tmp->alpha;
        maxy = fx;
        maxl = tmp;
        if (DPDUBUG) {
          fprintf(stdout, "maxy < fx, slope = alpha = %.12f\n", slope);
        }
      } else if (fabs(diff)<1e-10) {
        if (tmp->a1*maxl->b > maxl->a1*tmp->b) {
          slope = tmp->alpha;
          maxl = tmp;
          if (DPDUBUG) {
            fprintf(stdout, "maxy = fx, slope <alpha, slope <- alpha=%.12f\n", slope);
          }
        }
      } else {
        if (DPDUBUG) {
          fprintf(stdout, "\n");
        }
      }
    } else {
      slope = tmp->alpha;
      maxy = fx;
      maxl = tmp;
    }
  }
  (*y) = maxy;

  return slope;
}

void GetOptPoint(double *s1, double *s2, double *b, double c1, double c2, int s2LProws, double *x1, double *x2)
{
  int i, j, nI2=0;
  double a1, a2;
  struct _twoDline *linehead=NULL;
  double u1 = -INFINITY, u2 = INFINITY;
  double saprecision = 1e-10;
  double x, y;
  struct _twoDline *tmp=NULL;

  if ((linehead = malloc(sizeof(struct _twoDline))) == NULL) {
        fprintf(stdout,"unable to allocate memory \n");
        exit(-1);
    }
    linehead->alpha = 0.0;
    linehead->beta = 0.0;
    linehead->orig = 0;
  linehead->next = NULL;

  for (i=0; i<s2LProws; i++) {
    //printf("%d\n",i);
    a1 = s1[i] - s2[i] * c1 / c2;
    a2 = s2[i] / c2;
    if (DPDUBUG) {
      fprintf(stdout, "line %d, s1 = %.8f, s2 = %.8f, b = %.8f, a1 = %.8f, a2 = %.8f; ", i, s1[i], s2[i], b[i], a1, a2);
      fflush(stdout);
    }
    if (b[i] ==0 /*< DEF_EPSD*/) { // b[i] == 0, update u1, u2
      if (a1 > saprecision) { // u2 = min{-a2/a1}
        if (u2 > -a2 / a1) u2 = -a2 / a1;
        if (DPDUBUG) {
          fprintf(stdout, "a1>0; u2 = %.8f; ", u2);
        }
      } else if (a1 < -saprecision) { // u1 = max{-a2/a1}
        if (u1 < -a2 / a1) u1 = -a2 / a1;
        if (DPDUBUG) {
          fprintf(stdout, "a1<0; u1 = %.8f; ", u1);
        }
      } else {
        // what if a1[i] close to 0? -epsilon < a1[i] < epsilon?
        if (a2 > saprecision) { // Y=0, return 0,0
          (*x1) = 0.0;
          (*x2) = 0.0;
          if (DPDUBUG) {
            fprintf(stdout, "Y = 0, return (0,0)\n");
          }
          while (linehead != NULL) {
            tmp = linehead->next;
            free(linehead);
            linehead = tmp;
          }
          return;
        }
        if (DPDUBUG) {
          fprintf(stdout, "a1=0; ");
        }
      }
    } else { // calculate alpha, beta
      if (!(fabs(a1)<saprecision && fabs(a2)<saprecision)) { // a1<>0 && a2<>0
        // struct _twoDline *tmp;

        if (DPDUBUG) {
          printf("<================= i = %d\n ==============>",i);
        }

        tmp = malloc(sizeof(struct _twoDline));
        tmp->alpha = a1/b[i];
        tmp->beta = a2/b[i];
        tmp->a1 = a1;
        tmp->a2 = a2;
        tmp->b = b[i];
        tmp->orig = i;
        tmp->next = linehead->next;
        linehead->next = tmp;
        nI2++;
        if (DPDUBUG) {
          fprintf(stdout, "included; ");
        }
      } else {
        if (DPDUBUG) {
          fprintf(stdout, "a1=a2=0; ");
        }
      }
    }
    if (DPDUBUG) {
      fprintf(stdout, "\n");
    }
  }

  if (DPDUBUG) {
    fprintf(stdout, "finish initialize lines for step 0\n");
    fprintf(stdout, "u1 = %.12f, u2 = %.12f\n", u1, u2);
    fprintf(stdout, "nI2 = %d; following lines:\n", nI2);
    fflush(stdout);
  }

  if (u1 - u2 > saprecision) { // original solution (0,0);
    (*x1) = 0.0;
    (*x2) = 0.0;

    while (linehead != NULL) {
      tmp = linehead->next;
      free(linehead);
      linehead = tmp;
    }
    if (DPDUBUG) {
      fprintf(stdout, "u1>u2, solution (0,0)\n");
    }
    return;
  }

  if (fabs(u1-u2) < saprecision) { // u1 == u2
    x = u1;
    y = -INFINITY;

    while (linehead->next != NULL) {
      tmp = linehead->next;
      if (y < tmp->alpha*x+tmp->beta)
        y = tmp->alpha*x+tmp->beta;
      linehead->next = linehead->next->next;
      free(tmp);
    }
    free(linehead);
    (*x1) = x/y;
    (*x2) = (1-c1*x)/(c2*y);
    if (DPDUBUG) {
      fprintf(stdout, "u1=u2, x=%.12f, y=%.12f, solution (%.12f,%.12f)\n", x, y, (*x1), (*x2));
    }
    linehead = NULL;
    tmp = NULL;
    return;
  }

  // u1 < u2
  double lambda, rho;
  int k;
  double xi[s2LProws/2];
  int deleteList[s2LProws+1];
  int leftIdx[s2LProws/2];
  int rightIdx[s2LProws/2];

  rho = getRightSlope(linehead, u1, (&y));

  if (DPDUBUG) {
    fprintf(stdout, "rho = %.12f\n", rho);
    fprintf(stdout, "y = %.12f\n", y);
    fflush(stdout);
  }
  if (rho > -saprecision) {
    x = u1;
    (*x1) = x/y;
    (*x2) = (1-c1*x)/(c2*y);
    if (DPDUBUG) {
      fprintf(stdout, "opt at u1, x=%.12f, y=%.12f, solution (%.12f,%.12f)\n", x, y, (*x1), (*x2));
      fflush(stdout);
    }
    while (linehead != NULL) {
      tmp = linehead->next;
      free(linehead);
      linehead = tmp;
    }
    linehead = NULL;
    tmp = NULL;
    return;
  }
  lambda = getLeftSlope(linehead, u2, (&y));



  if (DPDUBUG) {
    fprintf(stdout, "lambda = %.12f\n", lambda);
    fprintf(stdout, "y = %.12f\n", y);
    fflush(stdout);
  }
  if (lambda < saprecision) {
    x = u2;
    (*x1) = x/y;
    (*x2) = (1-c1*x)/(c2*y);
    if (DPDUBUG) {
      fprintf(stdout, "opt at u2, x=%.12f, y=%.12f, solution (%.12f,%.12f)\n", x, y, (*x1), (*x2));
      fflush(stdout);
    }
    while (linehead != NULL) {
      tmp = linehead->next;
      free(linehead);
      linehead = tmp;
    }
    linehead = NULL;
    tmp = NULL;
    return;
  }

  while (1) {



    k = 0;
    struct _twoDline * tmph = NULL;
    struct _twoDline * tmpi = NULL;
    struct _twoDline * tmpj = NULL;
    tmph = linehead;
    tmpi = linehead->next;
    if (tmpi != NULL)
      tmpj = tmpi->next;
    else
      tmpj = NULL;
    if (DPDUBUG) {
      fprintf(stdout, "iteration start\n");
      fprintf(stdout, "u1 = %.12f, u2 = %.12f\n", u1, u2);
      fprintf(stdout, "nI2 = %d; following lines:\n", nI2);
      fflush(stdout);
    }

    if (nI2 == 2) {
      // get solution
      struct _twoDline *tmp = linehead;
      tmp = linehead;
      while (tmp->next != NULL) {
        tmp = tmp->next;
        fprintf(stdout, "line %d, alpha = %.12f, beta = %.12f\n", tmp->orig, tmp->alpha, tmp->beta);
      }
      fprintf(stdout, "\n");
      x = (tmpi->beta-tmpj->beta)/(tmpj->alpha-tmpi->alpha);
      y = tmpi->alpha*x+tmpi->beta;
      if (y > saprecision) {
        (*x1) = x/y;
        (*x2) = (1-c1*x)/(c2*y);
      } else {
        (*x1) = 0.0;
        (*x2) = 0.0;
      }
      // struct _twoDline *tmp;
      while (linehead != NULL) {
        tmp = linehead->next;
        free(linehead);
        linehead = tmp;

      }
      linehead = NULL;
      tmp = NULL;
      if (DPDUBUG) {
        fprintf(stdout, "opt, x=%.12f, y=%.12f, solution (%.12f,%.12f)\n", x, y, (*x1), (*x2));
        fflush(stdout);
      }
      return;
    }
    for (i=0; i<nI2; i++)
      deleteList[i] = 0;

    // step1: form pair
    if (DPDUBUG) {
      fprintf(stdout, "step1: form pair\n");
      fflush(stdout);
    }
    i = 0; j = 1;
    while (tmpj != NULL) {
      if (fabs(tmpi->alpha - tmpj->alpha)<saprecision) { // parallel
        if (tmpi->beta > tmpj->beta) {// remove tmpj
          deleteList[j] = 1;
          if (DPDUBUG) {
            fprintf(stdout, "parallel; remove tmpj; tmpi:alpha = %.12f, beta = %.12f, orig = %d; tmpj:alpha = %.12f, beta=%.12f, orig = %d; \n", tmpi->alpha, tmpi->beta,tmpi->orig, tmpj->alpha, tmpj->beta, tmpj->orig);
          }
        } else {  // remove tmpi
          deleteList[i] = 1;
          if (DPDUBUG) {
            fprintf(stdout, "parallel; remove tmpi; tmpi:alpha = %.12f, beta = %.12f, orig = %d; tmpj:alpha = %.12f, beta=%.12f, orig = %d; \n", tmpi->alpha, tmpi->beta, tmpi->orig, tmpj->alpha, tmpj->beta, tmpj->orig);
          }
        }
      } else { // calculate interseting point
        double tmpx = (tmpi->beta-tmpj->beta)/(tmpj->alpha-tmpi->alpha);
        if (DPDUBUG) {
          fprintf(stdout, "intersecting; x[%d]=%.12f; tmpi:alpha = %.12f, beta = %.12f, orig = %d; tmpj:alpha = %.12f, beta=%.12f, orig = %d; ", k+1, tmpx, tmpi->alpha, tmpi->beta, tmpi->orig, tmpj->alpha, tmpj->beta, tmpj->orig);
        }
        if (tmpx <= u1) {
          if (DPDUBUG) {
            fprintf(stdout, "x<= u1, ");
          }
          if (tmpi->alpha<tmpj->alpha) {
            deleteList[i] = 1;
            if (DPDUBUG) {
              fprintf(stdout, "remove tmpi\n");
            }
          } else {
            deleteList[j] = 1;
            if (DPDUBUG) {
              fprintf(stdout, "remove tmpj\n");
            }
          }
        } else if (tmpx >= u2) {
          if (DPDUBUG) {
            fprintf(stdout, "x>= u2, ");
          }
          if (tmpi->alpha>tmpj->alpha) {
            deleteList[i] = 1;
            if (DPDUBUG) {
              fprintf(stdout, "remove tmpi\n");
            }
          } else {
            deleteList[j] = 1;
            if (DPDUBUG) {
              fprintf(stdout, "remove tmpj\n");
            }
          }
        } else {
          if (DPDUBUG) {
            fprintf(stdout, "u1<x<u2, included in xi[]\n");
          }
          xi[k] = tmpx;
          if (tmpi->alpha < tmpj->alpha) {
            leftIdx[k] = i;
            rightIdx[k] = j;
          } else {
            leftIdx[k] = j;
            rightIdx[k] = i;
          }
          k++;
        }
      }
      i += 2;
      j += 2;
      tmpi = tmpj->next;
      if (tmpi != NULL)
        tmpj = tmpi->next;
      else
        tmpj = NULL;
    }

    if (k>0) {
      if (DPDUBUG) {
        fprintf(stdout, "k>0\n");
        fprintf(stdout, "k = %d\n", k);
        for (i=0; i<k; i++) {
          fprintf(stdout, "x%d = %.12f; ", i, xi[i]);
        }
        fflush(stdout);
      }
      //step2: find median
      double xtmp[k+1];
      for (i=0; i<k; i++)
        xtmp[i] = xi[i];
      double xstar = quickSelect(xtmp, k, k/2, &pickPivot);//kthSmallest(xtmp, 0, k-1, k/2+1);
      if (DPDUBUG) {
        fprintf(stdout, "k>0\n");
        fprintf(stdout, "k = %d\n", k);
        for (i=0; i<k; i++) {
          fprintf(stdout, "x%d = %.12f; ", i, xi[i]);
        }
        fflush(stdout);
      }
      if (DPDUBUG) {

        fprintf(stdout, "\n x* = %.12f\n", xstar);
        fflush(stdout);
      }

      // step3: evaluate lambda, rho at xstar
      rho = getRightSlope(linehead, xstar, (&y));
      if (DPDUBUG) {
        fprintf(stdout, "y from rightslope=%.12f\n", y);
      }
      lambda = getLeftSlope(linehead, xstar, (&y));
      if (DPDUBUG) {
        fprintf(stdout, "y from leftslope=%.12f\n", y);
      }
      if (DPDUBUG) {
        fprintf(stdout, "rho = %.12f, lambda %.12f\n", rho, lambda);
      }

      // step4:
      if (rho>-saprecision && lambda < saprecision) { // optimal
        (*x1) = xstar/y;
        (*x2) = (1-c1*xstar)/(c2*y);
        if (DPDUBUG) {
          fprintf(stdout, "opt, xstar=%.12f, y=%.12f, solution (%.12f,%.12f)\n", xstar, y, (*x1), (*x2));
        }

        while (linehead != NULL) {
          tmp = linehead->next;
          free(linehead);
          linehead = tmp;
        }
        linehead = NULL;
        tmp = NULL;
        return;
      } else if (rho < -saprecision) {
        if (DPDUBUG) {
          fprintf(stdout, "rho < 0, removing left part\n");
        }
        for (i=0; i<k; i++)
          if (xi[i] < xstar + saprecision) {
            deleteList[leftIdx[i]] = 1;
            if (DPDUBUG) {
              fprintf(stdout, "remove %d left part %d; xi=%.12f, xstar = %.12f\n", i, leftIdx[i], xi[i], xstar);
            }
          }
      } else if (lambda > saprecision) {
        if (DPDUBUG) {
          fprintf(stdout, "lambda > 0, removing right part\n");
        }
        for(i=0; i<k; i++){
          if (xi[i] >= xstar - saprecision ) {
            deleteList[rightIdx[i]] = 1;
            if (DPDUBUG) {
              fprintf(stdout, "remove %d right part %d; xi=%.12f, xstar = %.12f\n", i, rightIdx[i], xi[i], xstar);
            }
          }
        }
      }
    }

    // delete lines
    if (DPDUBUG) {
      fprintf(stdout, "deleteList:\n");
      for (i=0; i<nI2; i++)
        fprintf(stdout, "%d %d; ", i, deleteList[i]);
      fprintf(stdout, "\n");
      fflush(stdout);
    }
    struct _twoDline *del = linehead;
    i=0;
    while (del->next != NULL) {
      if(DPDUBUG){
        printf("i = %d line = %d linked value = %d\n",i,deleteList[i],del->next->orig);
      }
      if (deleteList[i] == 1) {
        tmp = del->next;
        del->next = del->next->next;
        free(tmp);
        nI2--;
      } else {
        del = del->next;
      }
      i++;
    }


  }

  // struct _twoDline *tmp;
  while (linehead != NULL) {
    tmp = linehead->next;
    free(linehead);
    linehead = tmp;
  }
  linehead = NULL;
  tmp = NULL;
  return;
}

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

//size of alpha and beta must be equal
double* pairing(struct Lines* line,double out[],double u1, double u2,double c[]){

  //save starting address of line
  struct Lines* start = line;



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

    return out;
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

    return out;
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

    return out;
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

        return out;
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




void make_random_problem(int seed, int ROWS, int COLUMNS){
    srand (seed);
    int i, j, k, inc_seed;
  char fn[256];

    FILE *fp = NULL;
    struct Instance *instance;
  instance = (struct Instance*)calloc(1, sizeof(struct Instance));
  instance->cost_vector = (double*)calloc(COLUMNS, sizeof(double));
    instance->rhs_vector = (double*)calloc(ROWS, sizeof(double));
    instance->a_matrix = (struct Rows_Double*)calloc(ROWS, sizeof(struct Rows_Double));
  for (i = 0; i < ROWS; i++) {
    instance->a_matrix[i].cols = (double*)calloc(COLUMNS, sizeof(double));
    }
    instance->number_variables = COLUMNS;
    instance->number_constraints = ROWS;

  inc_seed = 0;
  for (i = 0; i < NUMBER_PROBLEMS; i++) {
    // seed_value = STARTING_SEED + inc_seed;
  //       inc_seed += INCREMENTAL_SEED;
    sprintf(fn, "RP_TXT/%d_%d.txt", instance->number_constraints, instance->number_variables);
    if (!(fp = fopen(fn, "w"))) {
      fprintf(stderr, "There was a problem opening the output file!\n");
      goto TERMINATE;
    }
    if (SPARSE) {
      create_sparse(instance);
    }
    else {
      create_instance(instance);
    }
    for (j = 0; j < instance->number_variables; j++) {
      fprintf(fp, "%f ", instance->cost_vector[j]);
    }
    fprintf(fp, "\n\n");
    for (j = 0; j < instance->number_constraints; j++) {
      for (k = 0; k < instance->number_variables; k++) {
        fprintf(fp, "%f ", instance->a_matrix[j].cols[k]);
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (j = 0; j < instance->number_constraints; j++) {
      fprintf(fp, "%f\n", instance->rhs_vector[j]);
    }
    fclose(fp);
    }
  free(instance);

TERMINATE:
  return;

}

double function_random() {
    // const int a = 100801;
    // const int c = 103319;
    // const int m = 193723;
    // seed_value = (a*seed_value*c) % m;
    // random_number = (double)(seed_value) / (double)(m);
    // return random_number;
    return (double)rand()/(double)(RAND_MAX);
}

void create_instance(struct Instance *instance) {
    int i, j;
    double auxiliar, sum, store;
    for (i = 0; i < instance->number_constraints; i++) {
        sum = 0;
        for (j = 0; j < instance->number_variables; j++) {
            auxiliar = function_random();
            store = ceil(auxiliar * A_MULTIPLE_FACTOR + A_SUM_FACTOR);
            if (store == 0) {
                store = 1;
            }
            if (NEG_NUMBERS == 1) {
                auxiliar = function_random();
                if (auxiliar <= PROB_NEGATIVES) {
                    store = -store;
                }
                else{
                  if (ZERO_NUMBERS == 1){
                    auxiliar = function_random(); 
                    if (auxiliar <= PROB_ZERO) {
                      store = 0.0;
                    }
                  }
                }
            }
            else{
              if (ZERO_NUMBERS == 1){
                auxiliar = function_random();
                if (auxiliar <= PROB_ZERO) {
                  store = 0.0;
                }
              }
            }
            instance->a_matrix[i].cols[j] = store;
            sum += instance->a_matrix[i].cols[j];
        }
    instance->rhs_vector[i] = fabs(ceil(sum * B_TIGHTNESS_RATIO));
    }
    for (j = 0; j < instance->number_variables; j++) {
        sum = 0;
        for (i = 0; i < instance->number_constraints; i++) {
            sum += instance->a_matrix[i].cols[j];
        }
        auxiliar = function_random();
        store = fabs(ceil(auxiliar * C_MULTIPLE_FACTOR));
        instance->cost_vector[j] = sum + store;
    }
}

void create_sparse(struct Instance *instance) {
    int i, j;
    double auxiliar, sum, store;

    for (i = 0; i < instance->number_constraints; i++) {
        for (j = 0; j < (instance->number_variables + instance->number_constraints); j++) {
            instance->a_matrix[i].cols[j] = 0;
            instance->cost_vector[i] = 0;
        }
        instance->rhs_vector[i] = 0;
    }
    for (i = 0; i < instance->number_constraints; i++) {
        sum = 0;
        for (j = 0; j < instance->number_variables; j++) {
            auxiliar = function_random();
            if (auxiliar > SPARSE_PERCENTAGE) {
                if (NEG_NUMBERS == 1) {
                    auxiliar = function_random();
                    if (auxiliar <= PROB_NEGATIVES) {
                        instance->a_matrix[i].cols[j] = -1;
                    }
                    else {
                        instance->a_matrix[i].cols[j] = 1;
                    }
                }
                else {
                    instance->a_matrix[i].cols[j] = 1;
                }
            }
            sum += instance->a_matrix[i].cols[j];
        }
        instance->rhs_vector[i] = ceil(sum * B_TIGHTNESS_RATIO);
    }
    for (j = 0; j < instance->number_variables; j++) {
        sum = 0;
        for (i = 0; i < instance->number_constraints; i++) {
            sum += instance->a_matrix[i].cols[j];
        }
        auxiliar = function_random();
        store = ceil(auxiliar * C_MULTIPLE_FACTOR);
        instance->cost_vector[j] = sum + store;
    }
}

void solve_random_problem(int seed, int in_ROWS,int COLLUMNS){

  FILE           *fp, *OpenInputFile();
  int             Preprocess(), Postprocess(), passes, PCx();
  int             CheckParameters();
  LPtype         *LP, *ReducedLP, *Convert_MPS_LP();
  MPStype        *MPS, *ReadMPS();
  Parameters     *Inputs, *NewParameters();
  ChangeStack    *Record;
  MPSchanges     *Changes;
  solution       *Solution, *NewSolution();
  void            SplitFreeVars(), UnSplitFreeVars(), PrintSolution();
  int             filesize;
  int             status;
  double          readtime = 0.0, pretime = 0.0, SysTime;
  double          UserTime, OldSysTime, OldUserTime;

  extern        char            infile[200];
  extern        char            outfile[200];
  extern        char            sinfile[200];

  double time_start;
  double time_end;

  int ROWS = in_ROWS + 2;

  //cost function vales
  double cost[2];

  // //read in from file
  char fn[256];
  if(DPDUBUG){
    printf("RP_TXT/%d_%d.txt\n", in_ROWS, COLLUMNS);
  }
  sprintf(fn, "RP_TXT/%d_%d.txt", in_ROWS, COLLUMNS);
  FILE *fp2 = fopen(fn, "r");
  double *a1s = malloc(sizeof(double)*ROWS);
  double *a2s = malloc(sizeof(double)*ROWS);
  double *bs = malloc(sizeof(double)*ROWS);



        if(fp2){
            //read from file until end of fiel EOF
            fscanf(fp2, "%lf %lf\n", &cost[0],&cost[1]);

            //fill up linked list with file valuesxfgb
            for(int i = 0; i < in_ROWS; i++){

                fscanf(fp2, "%lf %lf\n", &a1s[i], &a2s[i]);

            }

            for(int i = 0; i < in_ROWS; i++){

                fscanf(fp2, "%lf\n", &bs[i]);

            }

            //fill last row with x1 > 0 and x2 > 0 constraints
            a1s[ROWS-2] = -1;
            a2s[ROWS-2] = 0;
            bs[ROWS-2]  = 0;

            a1s[ROWS-1] = 0;
            a2s[ROWS-1] = -1;
            bs[ROWS-1]  = 0;

        }
        fclose(fp2);
        if(DPDUBUG){
          printf("%fx1 %fx2\n",cost[0],cost[1]);
          for(int i = 0; i < ROWS; i++){
              printf("%fx1 + %fx2 <= %f\n", a1s[i],a2s[i],bs[i]);
          }
        }


   /* load the default parameters */
   Inputs = NewParameters();

   /* read modified parameters (if any) from specs file */
   ParseSpecsFile(Inputs, infile, sinfile);

   /* check for errors */
   if (CheckParameters(Inputs))
      {
      if (DPDUBUG) {
     printf("Error return from CheckParameters\n");
     }
    }

    double out[2];



    // double a1s [] = {1,0,3,-1,0};
    // double a2s [] = {0,2,2,0,-1};
    // double bs [] = {4,12,18,0,0};
    // double cost []=  {3,5};
    // ROWS = 5;

     LPtype *newLP = makePCXProblem(cost[0], cost[1],a1s, a2s, bs, ROWS);

     /* run the preprocessor, which converts LP to ReducedLP, and split the free
      * variables into positive and negative parts.  Double-check that there are
      * no free variables in ReducedLP.  */

     GetTime(&OldUserTime, &OldSysTime);

     if (Inputs->Preprocessing)
        {
     passes = Preprocess(newLP, &ReducedLP, &Record, Inputs);
     if (passes < 0)
        {
              if (DPDUBUG) {
              printf("Error return from Preprocessor:");
              }
           exit(PRESOLVE_ERROR);
        }
        }
     else
        ReducedLP = newLP;

     /* PrintLP(ReducedLP); */

     if (Inputs->Scaling)
        ScaleLP(ReducedLP, Inputs);

     SplitFreeVars(ReducedLP);

     GetTime(&UserTime, &SysTime);
     fflush(stdout);
     pretime = UserTime - OldUserTime + SysTime - OldSysTime;

     /* set up a solution data structure, and put the timing information that
      * we've gathered so far into it.  */

     Solution = NewSolution(ReducedLP->Rows, ReducedLP->Cols,
          Inputs->IterationLimit);
     Solution->ReadTime = readtime;
     Solution->PreprocessTime = pretime;

     /* solve the problem, keeping track of CPU times.  */

     time_start = clock();
     status = PCx(ReducedLP, Solution, Inputs);
     time_end = clock(); 


     double pcx_time = (time_end - time_start) / CLOCKS_PER_SEC;

     double pcx_results;
     double dyer_results;
     double personal_results;


     pcx_results = Solution->PrimalObjective*-1;

     DeleteLP(newLP);
     DeleteLP(ReducedLP);
     FreeSolution(Solution);
     FreeParameters(Inputs);

    time_start = clock();
    GetOptPoint(a1s, a2s,bs, cost[0], cost[1], ROWS, &out[0], &out[1]);
    time_end = clock(); 

    double dyer_time = (time_end - time_start) / CLOCKS_PER_SEC;

    double x1 = out[0];
    double y1 = out[1];
    dyer_results = cost[0]*x1 + cost[1]*y1;

    time_start = clock();
    struct tv_return* pair_parameters = two_variable(cost[0],cost[1],a1s,a2s,bs,ROWS);
    pairing(pair_parameters->line,out,pair_parameters->u1,pair_parameters->u2,cost);
    time_end = clock(); 

    double personal_time = (time_end - time_start) / CLOCKS_PER_SEC;

    x1 = out[0];
    y1 = out[1];
    personal_results = cost[0]*x1 + cost[1]*y1;

    freopen("/dev/tty","w",stdout);  

    printf("SEED              PCX       DYER       PERSONAL       D_CORRECT?       P_CORRECT       D<->P_CORRECT\n");
    printf("%d       %f %f  %f      ", seed, pcx_results, dyer_results, personal_results);

    if(fabs(pcx_results - dyer_results) < TOLERANCE_2){
      printf("YES              ");
    }
    else{
      printf("NO               ");
    }
    if(fabs(pcx_results - personal_results) < TOLERANCE_2){
      printf("YES              ");
    }
    else{
      printf("NO               ");
    }
    if(fabs(dyer_results - personal_results) < TOLERANCE_2){
      printf("YES");
    }
    else{
      printf("NO");
    }
    printf("\n");

    printf("TIMES: %fs %fs %fs\n",pcx_time,dyer_time,personal_time);

    free(a1s);
    free(a2s);
    free(bs);


}

void solve_problem(char * file_name, int in_ROWS, int COLLUMNS){
    FILE           *fp, *OpenInputFile();
  int             Preprocess(), Postprocess(), passes, PCx();
  int             CheckParameters();
  LPtype         *LP, *ReducedLP, *Convert_MPS_LP();
  MPStype        *MPS, *ReadMPS();
  Parameters     *Inputs, *NewParameters();
  ChangeStack    *Record;
  MPSchanges     *Changes;
  solution       *Solution, *NewSolution();
  void            SplitFreeVars(), UnSplitFreeVars(), PrintSolution();
  int             filesize;
  int             status;
  double          readtime = 0.0, pretime = 0.0, SysTime;
  double          UserTime, OldSysTime, OldUserTime;

  double time_start;
  double time_end;

  int ROWS = in_ROWS + 2;

  //cost function vales
  double cost[2];

  // //read in from file
  char fn[256];
  if(DPDUBUG){
    printf("RP_TXT/%s", file_name);
  }
  sprintf(fn, "RP_TXT/%s", file_name);
  FILE *fp2 = fopen(fn, "r");
  double *a1s = malloc(sizeof(double)*ROWS);
  double *a2s = malloc(sizeof(double)*ROWS);
  double *bs = malloc(sizeof(double)*ROWS);



        if(fp2){
            //read from file until end of fiel EOF
            fscanf(fp2, "%lf %lf\n", &cost[0],&cost[1]);

            //fill up linked list with file valuesxfgb
            for(int i = 0; i < in_ROWS; i++){

                fscanf(fp2, "%lf %lf\n", &a1s[i], &a2s[i]);

            }

            for(int i = 0; i < in_ROWS; i++){

                fscanf(fp2, "%lf\n", &bs[i]);

            }

            //fill last row with x1 > 0 and x2 > 0 constraints
            a1s[ROWS-2] = -1;
            a2s[ROWS-2] = 0;
            bs[ROWS-2]  = 0;

            a1s[ROWS-1] = 0;
            a2s[ROWS-1] = -1;
            bs[ROWS-1]  = 0;

        }
        fclose(fp2);
        if(DPDUBUG){
          printf("%fx1 %fx2\n",cost[0],cost[1]);
          for(int i = 0; i < ROWS; i++){
              printf("%fx1 + %fx2 <= %f\n", a1s[i],a2s[i],bs[i]);
          }
        }

   /* load the default parameters */
   Inputs = NewParameters();

   /* read modified parameters (if any) from specs file */
   ParseSpecsFile(Inputs, infile, sinfile);

   /* check for errors */
   if (CheckParameters(Inputs))
      {
      if (DPDUBUG) {
     printf("Error return from CheckParameters\n");
     }
    }

    double out[2];

     LPtype *newLP = makePCXProblem(cost[0], cost[1],a1s, a2s, bs, ROWS);

     /* run the preprocessor, which converts LP to ReducedLP, and split the free
      * variables into positive and negative parts.  Double-check that there are
      * no free variables in ReducedLP.  */

     GetTime(&OldUserTime, &OldSysTime);

     if (Inputs->Preprocessing)
        {
     passes = Preprocess(newLP, &ReducedLP, &Record, Inputs);
     if (passes < 0)
        {
              if (DPDUBUG) {
              printf("Error return from Preprocessor:");
              }
           exit(PRESOLVE_ERROR);
        }
        }
     else
        ReducedLP = newLP;

     /* PrintLP(ReducedLP); */

     if (Inputs->Scaling)
        ScaleLP(ReducedLP, Inputs);

     SplitFreeVars(ReducedLP);

     GetTime(&UserTime, &SysTime);
     fflush(stdout);
     pretime = UserTime - OldUserTime + SysTime - OldSysTime;

     /* set up a solution data structure, and put the timing information that
      * we've gathered so far into it.  */

     Solution = NewSolution(ReducedLP->Rows, ReducedLP->Cols,
          Inputs->IterationLimit);
     Solution->ReadTime = readtime;
     Solution->PreprocessTime = pretime;

     /* solve the problem, keeping track of CPU times.  */

     time_start = clock();
     status = PCx(ReducedLP, Solution, Inputs);
     time_end = clock(); 


     double pcx_time = (time_end - time_start) / CLOCKS_PER_SEC;

     double pcx_results;
     double dyer_results;
     double personal_results;


     pcx_results = Solution->PrimalObjective*-1;

     DeleteLP(newLP);
     DeleteLP(ReducedLP);
     FreeSolution(Solution);
     FreeParameters(Inputs);


    

    time_start = clock();
    GetOptPoint(a1s, a2s,bs, cost[0], cost[1], ROWS, &out[0], &out[1]);
    time_end = clock(); 

    double dyer_time = (time_end - time_start) / CLOCKS_PER_SEC;

    double x1 = out[0];
    double y1 = out[1];
    dyer_results = cost[0]*x1 + cost[1]*y1;

    time_start = clock();
    struct tv_return* pair_parameters = two_variable(cost[0],cost[1],a1s,a2s,bs,ROWS);
    pairing(pair_parameters->line,out,pair_parameters->u1,pair_parameters->u2,cost);
    time_end = clock(); 

    double personal_time = (time_end - time_start) / CLOCKS_PER_SEC;

    x1 = out[0];
    y1 = out[1];
    personal_results = cost[0]*x1 + cost[1]*y1;

    freopen("/dev/tty","w",stdout); 

    printf("SEED              PCX       DYER       PERSONAL       D_CORRECT?       P_CORRECT       D<->P_CORRECT\n");
    printf("       %f %f  %f      ", pcx_results, dyer_results, personal_results);

    if(fabs(pcx_results - dyer_results) < TOLERANCE_2){
      printf("YES              ");
    }
    else{
      printf("NO               ");
    }
    if(fabs(pcx_results - personal_results) < TOLERANCE_2){
      printf("YES              ");
    }
    else{
      printf("NO               ");
    }
    if(fabs(dyer_results - personal_results) < TOLERANCE_2){
      printf("YES");
    }
    else{
      printf("NO");
    }
    printf("\n");

    printf("TIMES: %fs %fs %fs\n",pcx_time,dyer_time,personal_time);

    free(a1s);
    free(a2s);
    free(bs);

}

int main(int argc, char**argv){

   
 
    //FOR TESTING
    int min = 1;
    int max = 100;
    int inc = 1;

    if(argc == 4){
      if(!strcmp(argv[1],"-p")){
        if(!DPDUBUG){
            fclose(stdout);
          }
        solve_problem(argv[2],atoi(argv[3]),2);
      }
      else{
        min = atoi(argv[1]);
        max = atoi(argv[2]);
        inc = atoi(argv[3]);
      

        for(int i = min; i < max + 1; i += inc){
          printf("\n<====== NEW PROBLEM %d ROWS and 2 COLS ======>\n", i);
          if(!DPDUBUG){
            fclose(stdout);
          }
          int seed = time(NULL);
          make_random_problem(seed, i,2);
          solve_random_problem(seed, i,2);
        }
      }
    }

}


//SEED Number   PCx Obj. Value    Your Obj. Value   Their Obj. Value    YES/NO    YES/NO



