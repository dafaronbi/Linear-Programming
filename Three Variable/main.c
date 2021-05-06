
#include "two variable.h"

	int main(int argc, char**argv){

        double c[] = {2,3,4};
        double a1[] = {-1,0,0,1,3.4,4.5,2.34};
        double a2[] = {0,-1,0,1,2.3,4.3,2.34};
        double a3[] = {0,0,-1,1,0.34,1.2,3.4};
        double b[] = {0,0,0,1,2.3,1,4};
        double out[3];

	three_variable(c, a1, a2, a3, b, out, 7);

        printf("x1 = %f, x2 = %f, x3 = %f\n", out[0],out[1],out[2]);
        printf("OBJECTIVE FUNCTION VALUE IS: %f\n", c[0]*out[0] + c[1]*out[1] +c[2]*out[2]);

	}