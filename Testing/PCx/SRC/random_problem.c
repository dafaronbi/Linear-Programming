

#include "two variable.h"
 
int make_random_problem(double seed, int ROWS, int COLUMNS){
    srand ( seed );
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
		fprintf(fp, "\n");		
		for (j = 0; j < instance->number_constraints; j++) {
			for (k = 0; k < instance->number_variables; k++) {
				fprintf(fp, "%f ", instance->a_matrix[j].cols[k]);
			}
			fprintf(fp, "\n");
		}
		for (j = 0; j < instance->number_constraints; j++) {
			fprintf(fp, "%f\n", instance->rhs_vector[j]);
		}
		fclose(fp);
    }
	free(instance);
 
TERMINATE:

    system("PAUSE");
    return 0;
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
 
void create_instance(double seed, struct Instance *instance) {
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
            }
            instance->a_matrix[i].cols[j] = store;
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

void solve_random_problem(int ROWS,int COLLUMNS){

        //cost function vales
        double cost[2];

        // //read in from file
        char fn[256];
        sprintf(fn, "RP_TXT/%d_%d.txt", ROWS, COLLUMNS);
        FILE * fp = fopen(fn, "r");
        double a1s[ROWS];
        double a2s[ROWS];
        double bs[ROWS];



        if(fp){
            //read from file until end of fiel EOF
            int r;
            r =fscanf(fp, "%f %f\n", &cost[0],&cost[1]);
            
            //fill up linked list with file values
            for(int i = 0; i < ROWS; i++){

                fscanf(fp, "%f %f\n", &a1s[i], &a2s[i]);

            }

            for(int i = 0; i < ROWS; i++){

                fscanf(fp, "%f\n", &bs[i]);

            }

        }
        fclose(fp);

        // printf("%fx1 %fx2\n",cost[0],cost[1]);
        // for(int i = 0; i < ROWS; i++){
        //     printf("%fx1 + %fx2 <= %f\n", a1s[i],a2s[i],bs[i]);
        // }

        // struct tv_return* pair_parameters = two_variable(cost[0],cost[1],a1s,a2s,bs,ROWS);

        

        //store output from paring problem
        double out[2];

        // if(pair_parameters->line != NULL){
        //     pairing(pair_parameters->line,out,pair_parameters->u1,pair_parameters->u2);
        // }
        // else{
        //     out[0] = 0;
        //     out[1] = 1;
        // }
        GetOptPoint(a1s, a2s,bs, cost[0], cost[1], ROWS, &out[0], &out[1]);

        double x1 = out[0]/out[1];
        double y1 = (1-3*out[0])/5*out[1];
        

        printf("<=============== RESULTS ===============>\n");
        printf("x = %f\n", x1);
        printf("y = %f\n", y1);

}