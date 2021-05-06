/***********************************************************

   Linear Programming Solver (LPSOL) 1995 -- 2020.

   Copyright by Author, Eva K. Lee, All rights reserved

   Purpose:

   Last Update : Zhuonan (Leon) Li

***********************************************************/
#include "two variable.h"
  
// It searches for x in arr[l..r], and partition_news the array
// around x. 

/*Comment out median finding for now
int partition_new(double arr[], int l, int r,  double x)
{ 
    // Search for x in arr[l..r] and move it to end 
    int i, j; 
    for (i=l; i<r; i++) 
        if (arr[i] == x) 
           break; 
    // if (DPDUBUG) {
    // 	fprintf(stdout, "partition_new: switch i = %d, l=%d, r=%d\n", i, l, r);
    // 	fflush(stdout);


double getShortMedian(double *xi, int k) { // return position ceil(k/2)
	if (k <= 5) {
		qsort(xi, k, sizeof(double), cmpfunc);
		return xi[k/2];
	}
}

double kthSmallest(double arr[], int l, int r, int k)
{ 
    // If k is smaller than number of elements in array 
    if (k > 0 && k <= r - l + 1) 
    { 
        int n = r-l+1; // Number of elements in arr[l..r] 
  
        // Divide arr[] in groups of size 5, calculate median 
        // of every group and store it in median[] array. 
        int i;
        double median[(n+4)/5]; // There will be floor((n+4)/5) groups;
        for (i=0; i<n/5; i++) 
            median[i] = getShortMedian(arr+l+i*5, 5); 
        if (i*5 < n) //For last group with less than 5 elements 
        { 
            median[i] = getShortMedian(arr+l+i*5, n%5);  
            i++; 
        }     
  
        // Find median of all medians using recursive call. 
        // If median[] has only one element, then no need 
        // of recursive call 
        double medOfMed = (i == 1)? median[i-1]:kthSmallest(median, 0, i-1, i/2);
  
        // partition_new the array around a random element and
        // get position of pivot element in sorted array 
        int pos = partition_new(arr, l, r, medOfMed);
  
        // If position is same as k 
        if (pos-l == k-1) 
            return arr[pos]; 
        if (pos-l > k-1)  // If position is more, recur for left 
            return kthSmallest(arr, l, pos-1, k); 
  
        // Else recur for right subarray 
        return kthSmallest(arr, pos+1, r, k-pos+l-1); 
    } 
  
    // If k is more than number of elements in array 
    return INFINITY; 
} 
 */ 


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

				printf("<================= i = %d\n ==============>",i);
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
				} else {	// remove tmpi
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
				fprintf(stdout, "rho < 0, removing left part\n");
				for (i=0; i<k; i++)
					if (xi[i] < xstar + saprecision) {
						deleteList[leftIdx[i]] = 1;
						if (DPDUBUG) {
							fprintf(stdout, "remove %d left part %d; xi=%.12f, xstar = %.12f\n", i, leftIdx[i], xi[i], xstar);
						}
					}
			} else if (lambda > saprecision) {
				fprintf(stdout, "lambda > 0, removing right part\n");
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
			printf("i = %d line = %d linked value = %d\n",i,deleteList[i],del->next->orig);
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
