/* Copyright 2012-2014 The Pennsylvania State University
 *
 * This software was written by David Hadka and others.
 * 
 * The use, modification and distribution of this software is governed by the
 * The Pennsylvania State University Research and Educational Use License.
 * You should have received a copy of this license along with this program.
 * If not, contact <dmh309@psu.edu>.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <errno.h>
#include "borg.h"
#include "mt19937ar.h"

/* Fix for Microsoft visual studio versions before MSVS 2013 */
#ifdef _MSC_VER 
	#if (_MSC_VER < 1800)
		#define isnan _isnan
	#endif	
#endif


#ifdef DBL_DECIMAL_DIG
  #define BORG_DIGITS (DBL_DECIMAL_DIG)
#else  
  #ifdef DECIMAL_DIG
    #define BORG_DIGITS (DECIMAL_DIG)
  #else  
    #define BORG_DIGITS (DBL_DIG + 3)
  #endif
#endif

static int BORG_Debug_enabled = 0;
static char* BORG_Debug_name = NULL;

struct BORG_Problem_t {
	int numberOfVariables;
	int numberOfObjectives;
	int numberOfConstraints;
	double* lowerBounds;
	double* upperBounds;
	double* epsilons;
	const char** names;
	void (*function)(double*, double*, double*);
};

struct BORG_Solution_t {
	BORG_Problem problem;
	double* variables;
	double* objectives;
	double* constraints;
	int operatorIndex;
};

struct BORG_Operator_t {
	const char* name;
	int numberOfParents;
	int numberOfOffspring;
	int numberOfParameters;
	double selectionProbability;
	double* parameters;
	void (*function)(BORG_Operator, BORG_Solution*, BORG_Solution*);
	BORG_Operator mutation;
};

struct BORG_Population_t {
	int size;
	int capacity;
	BORG_Solution* members;
};

typedef struct BORG_Entry_t *BORG_Entry;

struct BORG_Entry_t {
	BORG_Solution solution;
	BORG_Entry next;
	BORG_Entry prev;
	BORG_Archive archive;
};

struct BORG_Archive_t {
	int size;
	BORG_Entry head;
	BORG_Entry tail;
	int numberOfImprovements;
	int* recencyList;
	int recencyListSize;
	int recencyListPosition;
};

struct BORG_Algorithm_t {
	BORG_Problem problem;
	int numberOfEvaluations;
	int numberOfOperators;
	int tournamentSize;
	int windowSize;
	int maximumWindowSize;
	int initialPopulationSize;
	int minimumPopulationSize;
	int maximumPopulationSize;
	double populationRatio;
	double selectionRatio;
	int evaluationsAtLastCheck;
	int evaluationsAtLastRestart;
	int updateInterval;
	int operationsSinceLastUpdate;
	int improvementsAtLastCheck;
	BORG_Operator* operators;
	BORG_Population population;
	BORG_Archive archive;
	BORG_Restart restartMode;
	int restartedLastCheck;
	int baseMutationIndex;
	int maxMutationIndex;	
	BORG_Probabilities probabilityMode;
	int numberOfRestarts;
};

void BORG_Copyright(FILE* fp) {
	BORG_Validate_pointer(fp);

	fprintf(fp, "Copyright 2012-2014 The Pennsylvania State University\n");
	fprintf(fp, "This software was written by David Hadka and others.\n");
	fprintf(fp, "\n");
	fprintf(fp, "The use, modification and distribution of this software is governed by the\n");
	fprintf(fp, "The Pennsylvania State University Research and Educational Use License.\n");
	fprintf(fp, "You should have received a copy of this license along with this program.\n");
	fprintf(fp, "If not, contact <dmh309@psu.edu>.\n");
}

void BORG_Validate_file(FILE* fp) {
	if (!fp) {
		BORG_Error("unable to open file: %s\n", strerror(errno));
	}
}

void BORG_Validate_pointer(const void* ptr) {
	if (!ptr) {
		BORG_Error("pointer to NULL object\n");
	}
}

void BORG_Validate_index(int index, int size) {
	if ((index < 0) || (index >= size)) {
		BORG_Error("index out of bounds\n");
	}
}

void BORG_Validate_malloc(void* ptr) {
	if (!ptr) {
		BORG_Error("memory allocation failed\n");
	}
}

void BORG_Validate_positive(double value) {
	if (value < 0) {
		BORG_Error("value must be positive\n");
	}
}

void BORG_Debug_on() {
	BORG_Debug_enabled = 1;
}

void BORG_Debug_off() {
	BORG_Debug_enabled = 0;
}

void BORG_Debug_set_name(char* name) {
	BORG_Debug_name = name;
}

void BORG_Debug(const char* format, ...) {
	if (BORG_Debug_enabled) {
		va_list arguments;

		va_start(arguments, format);
		if (BORG_Debug_name) {
			fprintf(stderr, "%s: ", BORG_Debug_name);
		}
		vfprintf(stderr, format, arguments);
		fflush(stderr);
		va_end(arguments);
	}
}

void BORG_Error(const char* format, ...) {
	va_list arguments;

	va_start(arguments, format);
	if (BORG_Debug_name) {
		fprintf(stderr, "%s: ", BORG_Debug_name);
	}
	vfprintf(stderr, format, arguments);
	va_end(arguments);

	exit(EXIT_FAILURE);
}

BORG_Problem BORG_Problem_create(
		int numberOfVariables,
		int numberOfObjectives,
		int numberOfConstraints,
		void (*function)(double*, double*, double*)) {
	BORG_Validate_positive(numberOfVariables);
	BORG_Validate_positive(numberOfObjectives);
	BORG_Validate_positive(numberOfConstraints);
	BORG_Validate_pointer((void*)function);

	BORG_Problem problem = (BORG_Problem)malloc(sizeof(struct BORG_Problem_t));

	BORG_Validate_malloc(problem);

	problem->numberOfVariables = numberOfVariables;
	problem->numberOfObjectives = numberOfObjectives;
	problem->numberOfConstraints = numberOfConstraints;
	problem->lowerBounds = (double*)calloc(problem->numberOfVariables, sizeof(double));
	problem->upperBounds = (double*)calloc(problem->numberOfVariables, sizeof(double));
	problem->epsilons = (double*)calloc(problem->numberOfObjectives, sizeof(double));
	problem->names = (const char**)calloc(problem->numberOfObjectives, sizeof(const char*));
	problem->function = function;

	BORG_Validate_malloc(problem->lowerBounds);
	BORG_Validate_malloc(problem->upperBounds);
	BORG_Validate_malloc(problem->epsilons);
	BORG_Validate_malloc(problem->names);

	memset(problem->lowerBounds, 0, numberOfVariables*sizeof(double));
	memset(problem->upperBounds, 0, numberOfVariables*sizeof(double));
	memset(problem->epsilons, 0, numberOfObjectives*sizeof(double));
	memset(problem->names, 0, numberOfObjectives*sizeof(char*));

	return problem;
}

void BORG_Problem_destroy(BORG_Problem problem) {
	if (problem) {
		free(problem->lowerBounds);
		free(problem->upperBounds);
		free(problem->epsilons);
		free(problem->names);
	}

	free(problem);
}

void BORG_Problem_validate(BORG_Problem problem) {
	BORG_Validate_pointer(problem);

	int i;

	for (i = 0; i < problem->numberOfVariables; i++) {
		if (problem->lowerBounds[i] == 0.0 && problem->upperBounds[i] == 0.0) {
			BORG_Error("lower and upper bounds not set or both are set to 0.0\n");
		}
	}

	for (i = 0; i < problem->numberOfObjectives; i++) {
		if (problem->epsilons[i] == 0.0) {
			BORG_Error("epsilon values not set or assigned the value of 0.0\n");
		}
	}
}

void BORG_Problem_set_epsilon(BORG_Problem problem, int index, double epsilon) {
	BORG_Validate_pointer(problem);
	BORG_Validate_index(index, problem->numberOfObjectives);
	BORG_Validate_positive(epsilon);

	problem->epsilons[index] = epsilon;
}

void BORG_Problem_set_name(BORG_Problem problem, int index, const char* name) {
	BORG_Validate_pointer(problem);
	BORG_Validate_index(index, problem->numberOfObjectives);
	BORG_Validate_pointer((void*)name);

	problem->names[index] = name;
}

void BORG_Problem_set_epsilons(BORG_Problem problem, double* epsilons) {
	BORG_Validate_pointer(problem);
	BORG_Validate_pointer(epsilons);

	int i;

	for (i=0; i<problem->numberOfObjectives; i++) {
		BORG_Problem_set_epsilon(problem, i, epsilons[i]);
	}
}

void BORG_Problem_set_bounds(BORG_Problem problem, int index, double lowerBound, double upperBound) {
	BORG_Validate_pointer(problem);
	BORG_Validate_index(index, problem->numberOfVariables);

	problem->lowerBounds[index] = lowerBound;
	problem->upperBounds[index] = upperBound;
}

int BORG_Problem_number_of_variables(BORG_Problem problem) {
	BORG_Validate_pointer(problem);
	return problem->numberOfVariables;
}

int BORG_Problem_number_of_objectives(BORG_Problem problem) {
	BORG_Validate_pointer(problem);
	return problem->numberOfObjectives;
}

int BORG_Problem_number_of_constraints(BORG_Problem problem) {
	BORG_Validate_pointer(problem);
	return problem->numberOfConstraints;
}

BORG_Solution BORG_Solution_create(BORG_Problem problem) {
	BORG_Validate_pointer(problem);

	BORG_Solution solution = (BORG_Solution)malloc(sizeof(struct BORG_Solution_t));

	BORG_Validate_malloc(solution);

	solution->problem = problem;

	solution->variables = (double*)calloc(problem->numberOfVariables, sizeof(double));
	BORG_Validate_malloc(solution->variables);

	solution->objectives = (double*)calloc(problem->numberOfObjectives, sizeof(double));
	BORG_Validate_malloc(solution->objectives);

	if (problem->numberOfConstraints > 0) {
		solution->constraints = (double*)calloc(problem->numberOfConstraints, sizeof(double));
		BORG_Validate_malloc(solution->constraints);
	} else {
		solution->constraints = NULL;
	}

	solution->operatorIndex = -1;

	return solution;
}

void BORG_Solution_destroy(BORG_Solution solution) {
	if (solution) {
		free(solution->variables);
		free(solution->objectives);
		free(solution->constraints);
	}

	free(solution);
}

BORG_Solution BORG_Solution_clone(BORG_Solution original) {
	BORG_Validate_pointer(original);

	BORG_Problem problem = original->problem;
	BORG_Solution clone = BORG_Solution_create(problem);

	memcpy(clone->variables, original->variables, problem->numberOfVariables * sizeof(double));
	memcpy(clone->objectives, original->objectives, problem->numberOfObjectives * sizeof(double));
	memcpy(clone->constraints, original->constraints, problem->numberOfConstraints * sizeof(double));
	clone->operatorIndex = original->operatorIndex;

	return clone;
}

double BORG_Solution_get_variable(BORG_Solution solution, int index) {
	BORG_Validate_pointer(solution);
	BORG_Validate_index(index, solution->problem->numberOfVariables);

	return solution->variables[index];
}

double BORG_Solution_get_objective(BORG_Solution solution, int index) {
	BORG_Validate_pointer(solution);
	BORG_Validate_index(index, solution->problem->numberOfObjectives);

	return solution->objectives[index];
}

double BORG_Solution_get_constraint(BORG_Solution solution, int index) {
	BORG_Validate_pointer(solution);
	BORG_Validate_index(index, solution->problem->numberOfConstraints);

	return solution->constraints[index];
}

void BORG_Solution_set_variable(BORG_Solution solution, int index, double value) {
	BORG_Validate_pointer(solution);
	BORG_Problem problem = solution->problem;
	BORG_Validate_index(index, problem->numberOfVariables);

	if ((value < problem->lowerBounds[index]) || (value > problem->upperBounds[index])) {
		BORG_Error("value outside variable bounds\n");
	}

	solution->variables[index] = value;
}

void BORG_Solution_set_variables(BORG_Solution solution, double* variables) {
	BORG_Validate_pointer(solution);
	BORG_Validate_pointer(variables);

	int i;
	BORG_Problem problem = solution->problem;

	for (i=0; i<problem->numberOfVariables; i++) {
		BORG_Solution_set_variable(solution, i, variables[i]);
	}
}

void BORG_Solution_set_objective(BORG_Solution solution, int index, double value) {
	BORG_Validate_pointer(solution);
	BORG_Validate_index(index, solution->problem->numberOfObjectives);

	solution->objectives[index] = value;
}

void BORG_Solution_set_constraint(BORG_Solution solution, int index, double value) {
	BORG_Validate_pointer(solution);
	BORG_Validate_index(index, solution->problem->numberOfConstraints);

	solution->constraints[index] = value;
}

BORG_Problem BORG_Solution_get_problem(BORG_Solution solution) {
	BORG_Validate_pointer(solution);

	return solution->problem;
}

void BORG_Solution_evaluate(BORG_Solution solution) {
	BORG_Validate_pointer(solution);

	int i;
	BORG_Problem problem = solution->problem;
	int numberOfObjectives = problem->numberOfObjectives;
	int numberOfConstraints = problem->numberOfConstraints;

	solution->problem->function(solution->variables, solution->objectives, solution->constraints);

	/* guard against NaN values, which will always be non-dominated in the archive */
	for (i = 0; i < numberOfObjectives; i++) {
		if (isnan(solution->objectives[i])) {
			BORG_Error("Encountered a NaN objective value\n");
		}
	}

	for (i = 0; i < numberOfConstraints; i++) {
		if (isnan(solution->constraints[i])) {
			BORG_Error("Encountered a NaN constraint value\n");
		}
	}
}

void BORG_Solution_print(BORG_Solution solution, FILE* fp) {
	BORG_Validate_pointer(solution);
	BORG_Validate_pointer(fp);

	int i;

	if (solution->problem->numberOfVariables > 0) {
		fprintf(fp, "%.*g", BORG_DIGITS, solution->variables[0]);
	}

	for (i=1; i<solution->problem->numberOfVariables; i++) {
		fprintf(fp, " %.*g", BORG_DIGITS, solution->variables[i]);
	}

	for (i=0; i<solution->problem->numberOfObjectives; i++) {
		fprintf(fp, " %.*g", BORG_DIGITS, solution->objectives[i]);
	}

	for (i=0; i<solution->problem->numberOfConstraints; i++) {
		fprintf(fp, " %.*g", BORG_DIGITS, solution->constraints[i]);
	}

	fprintf(fp, "\n");
}

void BORG_Solution_initialize(BORG_Solution solution) {
	BORG_Validate_pointer(solution);

	int i;
	BORG_Problem problem = solution->problem;
	int numberOfVariables = problem->numberOfVariables;
	double* lowerBounds = problem->lowerBounds;
	double* upperBounds = problem->upperBounds;
	
	for (i=0; i<numberOfVariables; i++) {
		solution->variables[i] = BORG_Random_uniform(lowerBounds[i], upperBounds[i]);
	}
}

int BORG_Solution_violates_constraints(BORG_Solution solution) {
	BORG_Validate_pointer(solution);

	int i;
	int numberOfConstraints = solution->problem->numberOfConstraints;

	for (i=0; i<numberOfConstraints; i++) {
		if (solution->constraints[i] != 0.0) {
			return 1;
		}
	}

	return 0;
}

BORG_Dominance BORG_Dominance_pareto(BORG_Solution solution1, BORG_Solution solution2) {
	BORG_Validate_pointer(solution1);
	BORG_Validate_pointer(solution2);

	int i;
	int dominate1 = 0;
	int dominate2 = 0;
	int numberOfObjectives = solution1->problem->numberOfObjectives;
	
	for (i=0; i<numberOfObjectives; i++) {
		if (solution1->objectives[i] < solution2->objectives[i]) {
			dominate1 = 1;

			if (dominate2) {
				return NONDOMINATED;
			}
		} else if (solution1->objectives[i] > solution2->objectives[i]) {
			dominate2 = 1;

			if (dominate1) {
				return NONDOMINATED;
			}
		}
	}

	if (dominate1 == dominate2) {
		return NONDOMINATED;
	} else if (dominate1) {
		return DOMINATES;
	} else {
		return DOMINATED;
	}
}

BORG_Dominance BORG_Dominance_epsilon(BORG_Solution solution1, BORG_Solution solution2) {
	BORG_Validate_pointer(solution1);
	BORG_Validate_pointer(solution2);

	int i;
	int dominate1 = 0;
	int dominate2 = 0;
	int numberOfObjectives = solution1->problem->numberOfObjectives;
	double* epsilons = solution1->problem->epsilons;

	for (i=0; i<numberOfObjectives; i++) {
		double epsilon = epsilons[i];
		double index1 = floor(solution1->objectives[i] / epsilon);
		double index2 = floor(solution2->objectives[i] / epsilon);

		if (index1 < index2) {
			dominate1 = 1;

			if (dominate2) {
				return NONDOMINATED;
			}
		} else if (index1 > index2) {
			dominate2 = 1;

			if (dominate1) {
				return NONDOMINATED;
			}
		}
	}

	if (!dominate1 && !dominate2) {
		double dist1 = 0;
		double dist2 = 0;

		for (i=0; i<numberOfObjectives; i++) {
			double epsilon = epsilons[i];
			double index1 = floor(solution1->objectives[i] / epsilon);
			double index2 = floor(solution2->objectives[i] / epsilon);

			dist1 += pow(solution1->objectives[i] - index1*epsilon, 2.0);
			dist2 += pow(solution2->objectives[i] - index2*epsilon, 2.0);
		}

		// Calculating the sqrt is not necessary, the end result is the same
		// dist1 = sqrt(dist1);
		// dist2 = sqrt(dist2);

		if (dist1 < dist2) {
			return DOMINATES_SAME_BOX;
		} else {
			return DOMINATED_SAME_BOX;
		}
	} else if (dominate1) {
		return DOMINATES;
	} else {
		return DOMINATED;
	}
}

BORG_Dominance BORG_Dominance_constraints(BORG_Solution solution1, BORG_Solution solution2) {
	BORG_Validate_pointer(solution1);
	BORG_Validate_pointer(solution2);

	int i;
	double sum1 = 0.0;
	double sum2 = 0.0;
	int numberOfConstraints = solution1->problem->numberOfConstraints;

	for (i=0; i<numberOfConstraints; i++) {
		sum1 += fabs(solution1->constraints[i]);
		sum2 += fabs(solution2->constraints[i]);
	}

	if (sum1 == sum2) {
		return NONDOMINATED;
	} else if (sum1 < sum2) {
		return DOMINATES;
	} else {
		return DOMINATED;
	}
}

BORG_Dominance BORG_Dominance_compound(
		BORG_Solution solution1,
		BORG_Solution solution2,
		BORG_Dominance (*comparator1)(BORG_Solution, BORG_Solution),
		BORG_Dominance (*comparator2)(BORG_Solution, BORG_Solution)) {
	BORG_Validate_pointer((void*)comparator1);
	BORG_Validate_pointer((void*)comparator2);

	BORG_Dominance dominance = comparator1(solution1, solution2);

	if (dominance == NONDOMINATED) {
		dominance = comparator2(solution1, solution2);
	}

	return dominance;
}

void BORG_Random_seed(unsigned long seed) {
	init_genrand(seed);
}

double BORG_Random_uniform(double lowerBound, double upperBound) {
	return genrand_real1()*(upperBound - lowerBound) + lowerBound;
}

int BORG_Random_int(int n) {
	BORG_Validate_positive(n);

	int value = (int)BORG_Random_uniform(0.0, n);

	if (value > n-1) {
		value = n-1;
	}

	return value;
}

static double nextNextGaussian;
static int haveNextNextGaussian = 0;

double BORG_Random_gaussian(double mean, double stdev) {
	BORG_Validate_positive(stdev);

	double r;

	if (haveNextNextGaussian) {
		haveNextNextGaussian = 0;
		r = nextNextGaussian;
	} else {
		double v1, v2, s, m;

		do {
			v1 = BORG_Random_uniform(-1.0, 1.0);
			v2 = BORG_Random_uniform(-1.0, 1.0);
			s = v1*v1 + v2*v2;
		} while (s >= 1 || s == 0);

		m = sqrt(-2 * log(s)/s);
		nextNextGaussian = v2 * m;
		haveNextNextGaussian = 1;
		r = v1 * m;
	}

	return stdev*r + mean;
}

BORG_Operator BORG_Operator_create(
		const char* name,
		int numberOfParents,
		int numberOfOffspring,
		int numberOfParameters,
		void (*function)(BORG_Operator, BORG_Solution*, BORG_Solution*)) {
	BORG_Validate_pointer((void*)name);
	BORG_Validate_positive(numberOfParents);
	BORG_Validate_positive(numberOfOffspring);
	BORG_Validate_positive(numberOfParameters);
	BORG_Validate_pointer((void*)function);

	BORG_Operator variation = (BORG_Operator)malloc(sizeof(struct BORG_Operator_t));

	BORG_Validate_malloc(variation);

	variation->name = name;
	variation->numberOfParents = numberOfParents;
	variation->numberOfOffspring = numberOfOffspring;
	variation->numberOfParameters = numberOfParameters;
	variation->selectionProbability = 1.0;
	variation->parameters = (double*)calloc(numberOfParameters, sizeof(double));
	variation->function = function;
	variation->mutation = NULL;

	BORG_Validate_malloc(variation->parameters);

	memset(variation->parameters, 0, numberOfParameters*sizeof(double));

	return variation;
}

void BORG_Operator_destroy(BORG_Operator variation) {
	if (variation) {
		free(variation->parameters);
	}

	free(variation);
}

double BORG_Operator_get_probability(BORG_Operator variation) {
	BORG_Validate_pointer(variation);

	return variation->selectionProbability;
}

void BORG_Operator_set_parameter(BORG_Operator variation, int index, double parameter) {
	BORG_Validate_pointer(variation);
	BORG_Validate_index(index, variation->numberOfParameters);

	variation->parameters[index] = parameter;
}

void BORG_Operator_set_parameters(BORG_Operator variation, double* parameters) {
	BORG_Validate_pointer(variation);
	BORG_Validate_pointer(parameters);

	memcpy(variation->parameters, parameters, variation->numberOfParameters * sizeof(double));
}

int BORG_Operator_get_number_of_offspring(BORG_Operator variation) {
	BORG_Validate_pointer(variation);

	return variation->numberOfOffspring;
}

int BORG_Operator_get_number_of_parents(BORG_Operator variation) {
	BORG_Validate_pointer(variation);

	return variation->numberOfParents;
}

void BORG_Operator_set_mutation(BORG_Operator variation, BORG_Operator mutation) {
	BORG_Validate_pointer(variation);
	BORG_Validate_pointer(mutation);

	if (mutation->numberOfParents != 1 || mutation->numberOfOffspring != 1) {
		BORG_Error("mutation operator may only have 1 parent and 1 offspring\n");
	}

	variation->mutation = mutation;
}

void BORG_Operator_apply(BORG_Operator variation, BORG_Solution* parents, BORG_Solution* offspring) {
	BORG_Validate_pointer(variation);
	BORG_Validate_pointer(parents);
	BORG_Validate_pointer(offspring);

	variation->function(variation, parents, offspring);

	if (variation->mutation) {
		int i;
		BORG_Solution oldOffspring;

		for (i=0; i<variation->numberOfOffspring; i++) {
			oldOffspring = offspring[i];
			BORG_Operator_apply(variation->mutation, &offspring[i], &offspring[i]);
			BORG_Solution_destroy(oldOffspring);
		}
	}
}

void BORG_Operator_UM(BORG_Operator um, BORG_Solution* parents, BORG_Solution* offspring) {
	BORG_Validate_pointer(um);
	BORG_Validate_pointer(parents);
	BORG_Validate_pointer(offspring);

	int i;
	double probability = um->parameters[0];
	BORG_Problem problem = parents[0]->problem;
	int numberOfVariables = problem->numberOfVariables;

	offspring[0] = BORG_Solution_clone(parents[0]);

	for (i=0; i<numberOfVariables; i++) {
		if (BORG_Random_uniform(0.0, 1.0) <= probability) {
			offspring[0]->variables[i] = BORG_Random_uniform(problem->lowerBounds[i], problem->upperBounds[i]);
		}
	}
}

void BORG_Operator_SBX(BORG_Operator sbx, BORG_Solution* parents, BORG_Solution* offspring) {
	BORG_Validate_pointer(sbx);
	BORG_Validate_pointer(parents);
	BORG_Validate_pointer(offspring);

	int i;
	double probability = sbx->parameters[0];
	double distributionIndex = sbx->parameters[1];
	BORG_Problem problem = parents[0]->problem;
	int numberOfVariables = problem->numberOfVariables;

	offspring[0] = BORG_Solution_clone(parents[0]);
	offspring[1] = BORG_Solution_clone(parents[1]);

	if (BORG_Random_uniform(0.0, 1.0) <= probability) {
		for (i=0; i<numberOfVariables; i++) {
			if (BORG_Random_uniform(0.0, 1.0) <= 0.5) {
				double x0 = parents[0]->variables[i];
				double x1 = parents[1]->variables[i];
				double dx = fabs(x1 - x0);

				if (dx > EPSILON) {
					double lb = problem->lowerBounds[i];
					double ub = problem->upperBounds[i];
					double bl;
					double bu;

					if (x0 < x1) {
						bl = 1 + 2 * (x0 - lb) / dx;
						bu = 1 + 2 * (ub - x1) / dx;
					} else {
						bl = 1 + 2 * (x1 - lb) / dx;
						bu = 1 + 2 * (ub - x0) / dx;
					}

					if (bl < bu) {
						bu = bl;
					} else {
						bl = bu;
					}

					double p_bl = 1 - 1 / (2 * pow(bl, distributionIndex + 1));
					double p_bu = 1 - 1 / (2 * pow(bu, distributionIndex + 1));
					double u = genrand_real2(); /* need random number between [0,1) */
					double u0 = u * p_bl;
					double u1 = u * p_bu;
					double b0;
					double b1;

					if (u0 <= 0.5) {
						b0 = pow(2 * u0, 1 / (distributionIndex + 1));
					} else {
						b0 = pow(0.5 / (1 - u0), 1 / (distributionIndex + 1));
					}

					if (u1 <= 0.5) {
						b1 = pow(2 * u1, 1 / (distributionIndex + 1));
					} else {
						b1 = pow(0.5 / (1 - u1), 1 / (distributionIndex + 1));
					}

					double v1;
					double v2;

					if (x0 < x1) {
						v1 = 0.5 * (x0 + x1 + b0 * (x0 - x1));
						v2 = 0.5 * (x0 + x1 + b1 * (x1 - x0));
					} else {
						v1 = 0.5 * (x0 + x1 + b1 * (x0 - x1));
						v2 = 0.5 * (x0 + x1 + b0 * (x1 - x0));
					}

					if (BORG_Random_uniform(0.0, 1.0) <= 0.5) {
						double temp = v1;
						v1 = v2;
						v2 = temp;
					}

					if (v1 < lb) {
						v1 = lb;
					} else if (v1 > ub) {
						v1 = ub;
					}

					if (v2 < lb) {
						v2 = lb;
					} else if (v2 > ub) {
						v2 = ub;
					}

					offspring[0]->variables[i] = v1;
					offspring[1]->variables[i] = v2;
				}
			}
		}
	}
}

void BORG_Operator_PM(BORG_Operator pm, BORG_Solution* parents, BORG_Solution* offspring) {
	BORG_Validate_pointer(pm);
	BORG_Validate_pointer(parents);
	BORG_Validate_pointer(offspring);

	int i;
	double probability = pm->parameters[0];
	double distributionIndex = pm->parameters[1];
	BORG_Problem problem = parents[0]->problem;
	int numberOfVariables = problem->numberOfVariables;

	offspring[0] = BORG_Solution_clone(parents[0]);

	for (i=0; i<numberOfVariables; i++) {
		if (BORG_Random_uniform(0.0, 1.0) <= probability) {
			double u = BORG_Random_uniform(0.0, 1.0);
			double x = offspring[0]->variables[i];
			double lb = problem->lowerBounds[i];
			double ub = problem->upperBounds[i];
			double dx = ub - lb;
			double delta;

			if (u < 0.5) {
				double bl = (x - lb) / dx;
				double b = 2 * u + (1 - 2 * u)
					* pow(1 - bl, (distributionIndex + 1));
				delta = pow(b, (1.0 / (distributionIndex + 1))) - 1.0;
			} else {
				double bu = (ub - x) / dx;
				double b = 2 * (1 - u) + 2 * (u - 0.5)
					* pow(1 - bu, (distributionIndex + 1));
				delta = 1.0 - pow(b, (1.0 / (distributionIndex + 1)));
			}

			x = x + delta * dx;

			if (x < lb) {
				x = lb;
			} else if (x > ub) {
				x = ub;
			}

			offspring[0]->variables[i] = x;
		}
	}
}

void BORG_Operator_DE(BORG_Operator de, BORG_Solution* parents, BORG_Solution* offspring) {
	BORG_Validate_pointer(de);
	BORG_Validate_pointer(parents);
	BORG_Validate_pointer(offspring);

	int i;
	double CR = de->parameters[0];
	double F = de->parameters[1];
	BORG_Problem problem = parents[0]->problem;
	int numberOfVariables = problem->numberOfVariables;
	int irand = BORG_Random_int(numberOfVariables);

	offspring[0] = BORG_Solution_clone(parents[0]);

	for (i=0; i<numberOfVariables; i++) {
		if ((BORG_Random_uniform(0.0, 1.0) <= CR) || (i == irand)) {
			double v = parents[3]->variables[i] + F * (parents[1]->variables[i] - parents[2]->variables[i]);

			if (v < problem->lowerBounds[i]) {
				v = problem->lowerBounds[i];
			} else if (v > problem->upperBounds[i]) {
				v = problem->upperBounds[i];
			}

			offspring[0]->variables[i] = v;
		}
	}
}

void BORG_Operator_SPX(BORG_Operator spx, BORG_Solution* parents, BORG_Solution* offspring) {
	BORG_Validate_pointer(spx);
	BORG_Validate_pointer(parents);
	BORG_Validate_pointer(offspring);

	int i;
	int j;
	int k;
	int numberOfParents = spx->numberOfParents;
	int numberOfOffspring = spx->numberOfOffspring;
	double epsilon = spx->parameters[0];
	BORG_Problem problem = parents[0]->problem;
	int numberOfVariables = problem->numberOfVariables;
	double* G = (double*)calloc(numberOfVariables, sizeof(double));
	double* r = (double*)calloc(numberOfParents-1, sizeof(double));
	double** x = (double**)calloc(numberOfParents, sizeof(double*));
	double** C = (double**)calloc(numberOfParents, sizeof(double*));

	BORG_Validate_malloc(G);
	BORG_Validate_malloc(r);
	BORG_Validate_malloc(x);
	BORG_Validate_malloc(C);

	for (i=0; i<numberOfParents; i++) {
		x[i] = (double*)calloc(numberOfVariables, sizeof(double));
		C[i] = (double*)calloc(numberOfVariables, sizeof(double));

		BORG_Validate_malloc(x[i]);
		BORG_Validate_malloc(C[i]);
	}

	/* compute center of mass */
	memset(G, 0, numberOfVariables*sizeof(double));

	for (i=0; i<numberOfParents; i++) {
		for (j=0; j<numberOfVariables; j++) {
			G[j] += parents[i]->variables[j];
		}
	}

	for (j=0; j<numberOfVariables; j++) {
		G[j] /= numberOfParents;
	}

	/* compute simplex vertices expanded by epsilon */
	for (i=0; i<numberOfParents; i++) {
		for (j=0; j<numberOfVariables; j++) {
			x[i][j] = G[j] + epsilon * (parents[i]->variables[j] - G[j]);
		}
	}

	/* generate offspring */
	for (k=0; k<numberOfOffspring; k++) {
		offspring[k] = BORG_Solution_clone(parents[numberOfParents-1]);

		for (i=0; i<numberOfParents-1; i++) {
			r[i] = pow(BORG_Random_uniform(0.0, 1.0), 1.0 / (i + 1.0));
		}

		for (i=0; i<numberOfParents; i++) {
			for (j=0; j<numberOfVariables; j++) {
				if (i == 0) {
					C[i][j] = 0;
				} else {
					C[i][j] = r[i-1] * (x[i-1][j] - x[i][j] + C[i-1][j]);
				}
			}
		}

		for (j=0; j<numberOfVariables; j++) {
			double value = x[numberOfParents-1][j] + C[numberOfParents-1][j];

			if (value < problem->lowerBounds[j]) {
				value = problem->lowerBounds[j];
			} else if (value > problem->upperBounds[j]) {
				value = problem->upperBounds[j];
			}

			offspring[k]->variables[j] = value;
		}
	}
	
	/* release resources */
	for (i=0; i<numberOfParents; i++) {
		free(x[i]);
		free(C[i]);
	}

	free(G);
	free(r);
	free(x);
	free(C);
}

double* BORG_Vector_clone(int length, double* original) {
	int i;
	double* clone = (double*)calloc(length, sizeof(double));

	BORG_Validate_malloc(clone);

	for (i=0; i<length; i++) {
		clone[i] = original[i];
	}

	return clone;
}

void BORG_Vector_destroy(double* v) {
	free(v);
}

double* BORG_Vector_add(int length, double* u, double* v) {
	int i;

	for (i=0; i<length; i++) {
		u[i] += v[i];
	}

	return u;
}

double* BORG_Vector_subtract(int length, double* u, double* v) {
	int i;

	for (i=0; i<length; i++) {
		u[i] -= v[i];
	}

	return u;
}

double* BORG_Vector_multiply(int length, double* v, double c) {
	int i;

	for (i=0; i<length; i++) {
		v[i] *= c;
	}

	return v;
}

int BORG_Vector_iszero(int length, double* v) {
	int i;

	for (i=0; i<length; i++) {
		if (fabs(v[i]) > EPSILON) {
			return 0;
		}
	}

	return 1;
}

double BORG_Vector_dot(int length, double* u, double* v) {
	int i;
	double dot = 0.0;

	for (i=0; i<length; i++) {
		dot += u[i] * v[i];
	}

	return dot;
}

double BORG_Vector_magnitude(int length, double* u) {
	return sqrt(BORG_Vector_dot(length, u, u));
}

double* BORG_Vector_project(int length, double* u, double* v) {
	return BORG_Vector_multiply(length, v, BORG_Vector_dot(length, u, v)/BORG_Vector_dot(length, v, v));
}

double* BORG_Vector_normalize(int length, double* u) {
	return BORG_Vector_multiply(length, u, 1.0/BORG_Vector_magnitude(length, u));
}

double* BORG_Vector_orthogonalize(int length, double* v, int size, double** basis) {
	int i;

	for (i=0; i<size; i++) {
		double* u = BORG_Vector_clone(length, basis[i]);
		v = BORG_Vector_subtract(length, v, BORG_Vector_project(length, v, u));
		BORG_Vector_destroy(u);
	}

	return v;
}

BORG_Solution BORG_Operator_PCX_internal(BORG_Operator pcx, BORG_Solution* parents) {
	int i;
	int j;
	int numberOfParents = pcx->numberOfParents;
	double eta = pcx->parameters[0];
	double zeta = pcx->parameters[1];
	BORG_Problem problem = parents[0]->problem;
	int numberOfVariables = problem->numberOfVariables;

	if (numberOfParents < 2) {
		BORG_Error("PCX requires at least 2 parents");
	}

	double* g = (double*)calloc(numberOfVariables, sizeof(double));
	double** x = (double**)calloc(numberOfParents, sizeof(double*));

	BORG_Validate_malloc(g);
	BORG_Validate_malloc(x);

	memset(g, 0, numberOfVariables*sizeof(double));

	for (i=0; i<numberOfParents; i++) {
		x[i] = (double*)calloc(numberOfVariables, sizeof(double));

		BORG_Validate_malloc(x[i]);

		for (j=0; j<numberOfVariables; j++) {
			x[i][j] = parents[i]->variables[j];
			g[j] += x[i][j];
		}
	}

	for (j=0; j<numberOfVariables; j++) {
		g[j] /= numberOfParents;
	}

	double** e_eta = (double**)calloc(numberOfParents, sizeof(double*));
	int size_eta = 0;
	double D = 0.0;

	BORG_Validate_malloc(e_eta);

	e_eta[size_eta++] = BORG_Vector_subtract(numberOfVariables, 
		BORG_Vector_clone(numberOfVariables, x[numberOfParents-1]), g);

	/* basis vectors defined by parents */
	for (i=0; i<numberOfParents-1; i++) {
		double* d = BORG_Vector_subtract(numberOfVariables, 
			BORG_Vector_clone(numberOfVariables, x[i]), g);

		if (!BORG_Vector_iszero(numberOfVariables, d)) {
			d = BORG_Vector_orthogonalize(numberOfVariables, d, size_eta, e_eta);

			if (!BORG_Vector_iszero(numberOfVariables, d)) {
				D += BORG_Vector_magnitude(numberOfVariables, d);
				e_eta[size_eta++] = BORG_Vector_normalize(numberOfVariables, d);
			} else {
				BORG_Vector_destroy(d);
			}
		} else {
			BORG_Vector_destroy(d);
		}
	}

	D /= numberOfParents-1;

	/* construct the offspring */
	double rzeta = BORG_Random_gaussian(0.0, zeta);
	double reta = BORG_Random_gaussian(0.0, eta);
	double* variables = BORG_Vector_add(numberOfVariables, x[numberOfParents-1], 
		BORG_Vector_multiply(numberOfVariables, e_eta[0], rzeta));


	for (i=1; i<size_eta; i++) {
		variables = BORG_Vector_add(numberOfVariables, variables, 
			BORG_Vector_multiply(numberOfVariables, e_eta[i], reta*D));
	}

	BORG_Solution result = BORG_Solution_clone(parents[numberOfParents-1]);

	for (j=0; j<numberOfVariables; j++) {
		double value = variables[j];

		if (value < problem->lowerBounds[j]) {
			value = problem->lowerBounds[j];
		} else if (value > problem->upperBounds[j]) {
			value = problem->upperBounds[j];
		}

		result->variables[j] = value;
	}

	/* release resources */
	for (i=0; i<numberOfParents; i++) {
		free(x[i]);
	}

	for (i=0; i<size_eta; i++) {
		free(e_eta[i]);
	}

	free(e_eta);
	free(x);
	free(g);

	return result;
}

void BORG_Operator_PCX(BORG_Operator pcx, BORG_Solution* parents, BORG_Solution* offspring) {
	BORG_Validate_pointer(pcx);
	BORG_Validate_pointer(parents);
	BORG_Validate_pointer(offspring);

	int i;
	int numberOfParents = pcx->numberOfParents;
	int numberOfOffspring = pcx->numberOfOffspring;

	for (i=0; i<numberOfOffspring; i++) {
		int index = BORG_Random_int(numberOfParents);
		BORG_Solution temp = parents[index];
		parents[index] = parents[numberOfParents-1];
		parents[numberOfParents-1] = temp;

		offspring[i] = BORG_Operator_PCX_internal(pcx, parents);
	}
}

BORG_Solution BORG_Operator_UNDX_internal(BORG_Operator undx, BORG_Solution* parents) {
	int i;
	int j;
	int numberOfParents = undx->numberOfParents;
	double zeta = undx->parameters[0];
	double eta = undx->parameters[1];
	BORG_Problem problem = parents[0]->problem;
	int numberOfVariables = problem->numberOfVariables;

	if (numberOfParents < 2) {
		BORG_Error("PCX requires at least 2 parents");
	}

	double* g = (double*)calloc(numberOfVariables, sizeof(double));
	double** x = (double**)calloc(numberOfParents, sizeof(double*));

	BORG_Validate_malloc(g);
	BORG_Validate_malloc(x);

	memset(g, 0, numberOfVariables*sizeof(double));

	for (i=0; i<numberOfParents; i++) {
		x[i] = (double*)calloc(numberOfVariables, sizeof(double));

		BORG_Validate_malloc(x[i]);

		for (j=0; j<numberOfVariables; j++) {
			x[i][j] = parents[i]->variables[j];
			g[j] += x[i][j];
		}
	}

	for (j=0; j<numberOfVariables; j++) {
		g[j] /= numberOfParents;
	}

	double** e_zeta = (double**)calloc(numberOfParents, sizeof(double*));
	double** e_eta = (double**)calloc(numberOfVariables, sizeof(double*));
	int size_zeta = 0;
	int size_eta = 0;

	BORG_Validate_malloc(e_zeta);
	BORG_Validate_malloc(e_eta);

	double* d = BORG_Vector_subtract(numberOfVariables, 
		BORG_Vector_clone(numberOfVariables, x[numberOfParents-1]), g);
	double D = BORG_Vector_magnitude(numberOfVariables, d);
	BORG_Vector_destroy(d);

	/* basis vectors defined by parents */
	for (i=0; i<numberOfParents-1; i++) {
		d = BORG_Vector_subtract(numberOfVariables, 
			BORG_Vector_clone(numberOfVariables, x[i]), g);

		if (!BORG_Vector_iszero(numberOfVariables, d)) {
			double dbar = BORG_Vector_magnitude(numberOfVariables, d);
			d = BORG_Vector_orthogonalize(numberOfVariables, d, size_zeta, e_zeta);

			if (!BORG_Vector_iszero(numberOfVariables, d)) {
				e_zeta[size_zeta++] = BORG_Vector_multiply(numberOfVariables, 
					BORG_Vector_normalize(numberOfVariables, d), dbar);
			} else {
				BORG_Vector_destroy(d);
			}
		} else {
			BORG_Vector_destroy(d);
		}
	}

	/* create the remaining basis vectors */
	for (i=0; i<numberOfVariables-size_zeta; i++) {
		d = (double*)calloc(numberOfVariables, sizeof(double));

		BORG_Validate_malloc(d);

		for (j=0; j<numberOfVariables; j++) {
			d[j] = BORG_Random_gaussian(0.0, 1.0);
		}

		if (!BORG_Vector_iszero(numberOfVariables, d)) {
			d = BORG_Vector_orthogonalize(numberOfVariables, d, size_eta, e_eta);

			if (!BORG_Vector_iszero(numberOfVariables, d)) {
				e_eta[size_eta++] = BORG_Vector_multiply(numberOfVariables, 
					BORG_Vector_normalize(numberOfVariables, d), D);
			} else {
				BORG_Vector_destroy(d);
			}
		} else {
			BORG_Vector_destroy(d);
		}
	}

	/* construct the offspring */
	double* variables = g;

	for (i=0; i<size_zeta; i++) {
		variables = BORG_Vector_add(numberOfVariables, variables, 
			BORG_Vector_multiply(numberOfVariables, e_zeta[i], BORG_Random_gaussian(0.0, zeta)));
	}

	for (i=0; i<size_eta; i++) {
		variables = BORG_Vector_add(numberOfVariables, variables, 
			BORG_Vector_multiply(numberOfVariables, e_eta[i], BORG_Random_gaussian(0.0, eta / sqrt((double)numberOfVariables))));
	}

	BORG_Solution result = BORG_Solution_clone(parents[numberOfParents-1]);

	for (j=0; j<numberOfVariables; j++) {
		double value = variables[j];

		if (value < problem->lowerBounds[j]) {
			value = problem->lowerBounds[j];
		} else if (value > problem->upperBounds[j]) {
			value = problem->upperBounds[j];
		}

		result->variables[j] = value;
	}

	/* release resources */
	for (i=0; i<numberOfParents; i++) {
		free(x[i]);
	}

	for (i=0; i<size_zeta; i++) {
		free(e_zeta[i]);
	}

	for (i=0; i<size_eta; i++) {
		free(e_eta[i]);
	}

	free(e_zeta);
	free(e_eta);
	free(x);
	free(g);

	return result;
}

void BORG_Operator_UNDX(BORG_Operator undx, BORG_Solution* parents, BORG_Solution* offspring) {
	BORG_Validate_pointer(undx);
	BORG_Validate_pointer(parents);
	BORG_Validate_pointer(offspring);

	int i;
	int numberOfOffspring = undx->numberOfOffspring;

	for (i=0; i<numberOfOffspring; i++) {
		offspring[i] = BORG_Operator_UNDX_internal(undx, parents);
	}
}

BORG_Population BORG_Population_create(int capacity) {
	BORG_Validate_positive(capacity);

	BORG_Population population = (BORG_Population)malloc(sizeof(struct BORG_Population_t));

	BORG_Validate_malloc(population);

	population->size = 0;
	population->capacity = capacity;
	population->members = (BORG_Solution*)calloc(capacity, sizeof(struct BORG_Solution_t*));

	BORG_Validate_malloc(population->members);

	return population;
}

void BORG_Population_destroy(BORG_Population population) {
	int i;

	if (population) {
		for (i=0; i<population->size; i++) {
			BORG_Solution_destroy(population->members[i]);
		}

		free(population->members);
	}

	free(population);
}

void BORG_Population_reset(BORG_Population population, int newCapacity) {
	BORG_Validate_pointer(population);
	BORG_Validate_positive(newCapacity);

	int i;

	for (i=0; i<population->size; i++) {
		BORG_Solution_destroy(population->members[i]);
	}

	population->members = (BORG_Solution*)realloc(population->members, newCapacity * sizeof(struct BORG_Solution_t*));

	BORG_Validate_malloc(population->members);

	population->size = 0;
	population->capacity = newCapacity;
}

/**
 * Private method for adding a solution to the population when the population
 * is at capacity.
 */
void BORG_Population_replace(BORG_Population population, BORG_Solution solution) {
	BORG_Validate_pointer(population);
	BORG_Validate_pointer(solution);

	int i;
	int dominated = 0;
	int size = 0;
	int* dominates = (int*)calloc(population->size, sizeof(int));

	BORG_Validate_malloc(dominates);

	for (i=0; i<population->size; i++) {
		BORG_Dominance flag = BORG_Dominance_compound(solution, population->members[i], 
			BORG_Dominance_constraints, BORG_Dominance_pareto);

		if (flag < 0) {
			dominates[size++] = i;
		} else if (flag > 0) {
			dominated = 1;
		}
	}

	if (size > 0) {
		int index = dominates[BORG_Random_int(size)];
			
		BORG_Solution_destroy(population->members[index]);
		population->members[index] = BORG_Solution_clone(solution);
	} else if (!dominated) {
		int index = BORG_Random_int(population->size);

		BORG_Solution_destroy(population->members[index]);
		population->members[index] = BORG_Solution_clone(solution);
	}

	free(dominates);
}

void BORG_Population_add(BORG_Population population, BORG_Solution solution) {
	BORG_Validate_pointer(population);
	BORG_Validate_pointer(solution);

	if (population->size < population->capacity) {
		population->members[population->size] = BORG_Solution_clone(solution);
		population->size++;
	} else {
		BORG_Population_replace(population, solution);
	}
}

BORG_Solution BORG_Population_tournament(BORG_Population population, int tournamentSize) {
	BORG_Validate_pointer(population);
	
	int i;
	BORG_Solution winner = population->members[BORG_Random_int(population->size)];

	for (i=1; i<tournamentSize; i++) {
		BORG_Solution candidate = population->members[BORG_Random_int(population->size)];
		BORG_Dominance flag = BORG_Dominance_compound(winner, candidate,
			BORG_Dominance_constraints, BORG_Dominance_pareto);

		if (flag > 0) {
			winner = candidate;
		}
	}

	return winner;
}

void BORG_Population_select(BORG_Population population, int arity, BORG_Solution* parents, int tournamentSize) {
	BORG_Validate_pointer(population);
	BORG_Validate_positive(arity);
	BORG_Validate_pointer(parents);
	int i;

	for (i=0; i<arity; i++) {
		parents[i] = BORG_Population_tournament(population, tournamentSize);
	}
}

BORG_Entry BORG_Entry_create(BORG_Archive archive, BORG_Solution solution) {
	BORG_Validate_pointer(archive);
	BORG_Validate_pointer(solution);

	BORG_Entry entry = (BORG_Entry)malloc(sizeof(struct BORG_Entry_t));

	BORG_Validate_malloc(entry);

	entry->solution = BORG_Solution_clone(solution);
	entry->next = NULL;
	entry->prev = NULL;
	entry->archive = archive;

	return entry;
}

void BORG_Entry_destroy(BORG_Entry entry) {
	if (entry) {
		BORG_Solution_destroy(entry->solution);
	}

	free(entry);
}

BORG_Archive BORG_Archive_create() {
	int i;
	BORG_Archive archive = (BORG_Archive)malloc(sizeof(struct BORG_Archive_t));

	BORG_Validate_malloc(archive);

	archive->size = 0;
	archive->head = NULL;
	archive->tail = NULL;
	archive->numberOfImprovements = 0;
	archive->recencyListSize = 50;
	archive->recencyList = (int*)calloc(archive->recencyListSize, sizeof(int));
	archive->recencyListPosition = 0;

	BORG_Validate_malloc(archive->recencyList);

	for (i=0; i<archive->recencyListSize; i++) {
		archive->recencyList[i] = -1;
	}

	return archive;
}

BORG_Archive BORG_Archive_clone(BORG_Archive original) {
	BORG_Validate_pointer(original);

	BORG_Archive clone = BORG_Archive_create();
	BORG_Entry entry = original->head;

	while (entry) {
		BORG_Archive_add(clone, entry->solution);
		entry = entry->next;
	}

	clone->numberOfImprovements = original->numberOfImprovements;

	return clone;
}

BORG_Entry BORG_Archive_remove(BORG_Entry entry) {
	BORG_Validate_pointer(entry);

	BORG_Archive archive = entry->archive;
	BORG_Entry next = entry->next;

	if (archive->head == entry) {
		archive->head = entry->next;
	}

	if (archive->tail == entry) {
		archive->tail = entry->prev;
	}

	if (entry->next) {
		entry->next->prev = entry->prev;
	}

	if (entry->prev) {
		entry->prev->next = entry->next;
	}

	BORG_Entry_destroy(entry);
	archive->size--;

	return next;
}

void BORG_Archive_destroy(BORG_Archive archive) {
	if (archive) {
		BORG_Entry entry = archive->head;

		while (entry) {
			entry = BORG_Archive_remove(entry);
		}

		free(archive->recencyList);
	}

	free(archive);
}

void BORG_Archive_add(BORG_Archive archive, BORG_Solution solution) {
	BORG_Validate_pointer(archive);
	BORG_Validate_pointer(solution);

	int isSameBox = 0;
	BORG_Entry entry = archive->head;

	while (entry) {
		BORG_Dominance flag = BORG_Dominance_compound(solution, entry->solution, 
			BORG_Dominance_constraints, BORG_Dominance_epsilon);

		if (flag < 0) {
			entry = BORG_Archive_remove(entry);

			if (flag == DOMINATES_SAME_BOX) {
				isSameBox = 1;
			}
		} else if (flag > 0) {
			return;
		} else {
			entry = entry->next;
		}
	}

	entry = BORG_Entry_create(archive, solution);

	if (archive->tail) {
		archive->tail->next = entry;
		entry->prev = archive->tail;
		archive->tail = entry;
	} else {
		archive->head = entry;
		archive->tail = entry;
	}

	archive->size++;

	if (!isSameBox) {
		archive->numberOfImprovements++;
	}

	archive->recencyList[archive->recencyListPosition] = solution->operatorIndex;
	archive->recencyListPosition = (archive->recencyListPosition + 1) % archive->recencyListSize;
}

int BORG_Archive_get_size(BORG_Archive archive) {
	BORG_Validate_pointer(archive);

	return archive->size;
}

BORG_Solution BORG_Archive_get(BORG_Archive archive, int index) {
	BORG_Validate_pointer(archive);
	BORG_Validate_index(index, archive->size);

	BORG_Entry entry = archive->head;

	while (index > 0) {
		entry = entry->next;
		index--;
	}

	return entry->solution;
}

int BORG_Archive_get_improvements(BORG_Archive archive) {
	BORG_Validate_pointer(archive);

	return archive->numberOfImprovements;
}

BORG_Solution BORG_Archive_select(BORG_Archive archive) {
	BORG_Validate_pointer(archive);

	return BORG_Archive_get(archive, BORG_Random_int(archive->size));
}

void BORG_Archive_print(BORG_Archive archive, FILE* fp) {
	BORG_Validate_pointer(archive);
	BORG_Validate_pointer(fp);

	BORG_Entry entry = archive->head;

	while (entry) {
		BORG_Solution_print(entry->solution, fp);
		entry = entry->next;
	}
}

void BORG_Archive_append(BORG_Archive archive, FILE* file) {
	int i;
	int j;
	int size = 0;

	for (i=0; i<BORG_Archive_get_size(archive); i++) {
		BORG_Solution solution = BORG_Archive_get(archive, i);

		if (BORG_Solution_violates_constraints(solution)) {
			continue;
		} else {
			size++;
		}

		for (j=0; j<solution->problem->numberOfVariables; j++) {
			if (j > 0) {
				fprintf(file, " ");
			}

			fprintf(file, "%.*g", BORG_DIGITS, solution->variables[j]);
		}

		for (j=0; j<solution->problem->numberOfObjectives; j++) {
			if (j > 0 || solution->problem->numberOfVariables > 0) {
				fprintf(file, " ");
			}

			fprintf(file, "%.*g", BORG_DIGITS, solution->objectives[j]);
		}

		fprintf(file, "\n");
	}

	if (size <= 0) {
		fprintf(file, "//Empty\n");
	}

	fprintf(file, "#\n");
	
	fflush(file);
}

BORG_Algorithm BORG_Algorithm_create(BORG_Problem problem, int numberOfOperators) {
	BORG_Validate_pointer(problem);
	BORG_Validate_positive(numberOfOperators);


	BORG_Algorithm algorithm = (BORG_Algorithm)malloc(sizeof(struct BORG_Algorithm_t));

	BORG_Validate_malloc(algorithm);

	algorithm->problem = problem;
	algorithm->numberOfEvaluations = 0;
	algorithm->numberOfOperators = numberOfOperators;
	algorithm->tournamentSize = 2;
	algorithm->windowSize = 200;
	algorithm->maximumWindowSize = 20000;
	algorithm->initialPopulationSize = 100;
	algorithm->minimumPopulationSize = 100;
	algorithm->maximumPopulationSize = 10000;
	algorithm->populationRatio = 4.0;
	algorithm->selectionRatio = 0.02;
	algorithm->evaluationsAtLastCheck = 0;
	algorithm->evaluationsAtLastRestart = 0;
	algorithm->updateInterval = 100;
	algorithm->operationsSinceLastUpdate = 0;
	algorithm->improvementsAtLastCheck = 0;
	algorithm->operators = (BORG_Operator*)calloc(numberOfOperators, sizeof(struct BORG_Operator_t*));
	algorithm->population = BORG_Population_create(100);
	algorithm->archive = BORG_Archive_create();
	algorithm->restartMode = RESTART_DEFAULT;
	algorithm->restartedLastCheck = 0;
	algorithm->baseMutationIndex = 0;
	algorithm->maxMutationIndex = 10;
	algorithm->probabilityMode = PROBABILITIES_DEFAULT;
	algorithm->numberOfRestarts = 0;

	BORG_Validate_malloc(algorithm->operators);
	memset(algorithm->operators, 0, numberOfOperators*sizeof(struct BORG_Operator_t*));

	return algorithm;
}

void BORG_Algorithm_destroy(BORG_Algorithm algorithm) {
	if (algorithm) {
		free(algorithm->operators);
		BORG_Population_destroy(algorithm->population);
		BORG_Archive_destroy(algorithm->archive);
	}

	free(algorithm);
}

void BORG_Algorithm_validate(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	int i;

	if (algorithm->windowSize > algorithm->maximumWindowSize) {
		BORG_Error("window size is greater than max window size\n");
	}

	if (algorithm->minimumPopulationSize > algorithm->maximumPopulationSize) {
		BORG_Error("minimum population size is greater than maximum population size\n");
	}

	for (i=0; i<algorithm->numberOfOperators; i++) {
		if (!algorithm->operators[i]) {
			BORG_Error("operator %d is undefined\n", i);
		}
	}
}

void BORG_Algorithm_update_tournament_size(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	BORG_Population population = algorithm->population;
	algorithm->tournamentSize = (int)(algorithm->selectionRatio * population->size);

	if (algorithm->tournamentSize < 2) {
		algorithm->tournamentSize = 2;
	}
}

int BORG_Algorithm_get_number_improvements(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	return algorithm->archive->numberOfImprovements;
}

int BORG_Algorithm_get_number_restarts(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	return algorithm->numberOfRestarts;
}

void BORG_Algorithm_set_window_size(BORG_Algorithm algorithm, int windowSize) {
	BORG_Validate_pointer(algorithm);
	BORG_Validate_positive(windowSize);

	algorithm->windowSize = windowSize;
}

void BORG_Algorithm_set_maximum_window_size(BORG_Algorithm algorithm, int maximumWindowSize) {
	BORG_Validate_pointer(algorithm);
	BORG_Validate_positive(maximumWindowSize);

	algorithm->maximumWindowSize = maximumWindowSize;
}

void BORG_Algorithm_set_initial_population_size(BORG_Algorithm algorithm, int initialPopulationSize) {
	BORG_Validate_pointer(algorithm);
	BORG_Validate_positive(initialPopulationSize);

	algorithm->initialPopulationSize = initialPopulationSize;
}

void BORG_Algorithm_set_minimum_population_size(BORG_Algorithm algorithm, int minimumPopulationSize) {
	BORG_Validate_pointer(algorithm);
	BORG_Validate_positive(minimumPopulationSize);

	algorithm->minimumPopulationSize = minimumPopulationSize;
}

void BORG_Algorithm_set_maximum_population_size(BORG_Algorithm algorithm, int maximumPopulationSize) {
	BORG_Validate_pointer(algorithm);
	BORG_Validate_positive(maximumPopulationSize);

	algorithm->maximumPopulationSize = maximumPopulationSize;
}

void BORG_Algorithm_set_population_ratio(BORG_Algorithm algorithm, double populationRatio) {
	BORG_Validate_pointer(algorithm);
	BORG_Validate_positive(populationRatio);

	algorithm->populationRatio = populationRatio;
}

void BORG_Algorithm_set_selection_ratio(BORG_Algorithm algorithm, double selectionRatio) {
	BORG_Validate_pointer(algorithm);
	BORG_Validate_positive(selectionRatio);

	algorithm->selectionRatio = selectionRatio;
}

void BORG_Algorithm_set_restart_mode(BORG_Algorithm algorithm, BORG_Restart restartMode) {
	BORG_Validate_pointer(algorithm);

	algorithm->restartMode = restartMode;
}

void BORG_Algorithm_set_probability_mode(BORG_Algorithm algorithm, BORG_Probabilities probabilityMode) {
	BORG_Validate_pointer(algorithm);

	algorithm->probabilityMode = probabilityMode;
}

void BORG_Algorithm_set_operator(BORG_Algorithm algorithm, int index, BORG_Operator variation) {
	BORG_Validate_pointer(algorithm);
	BORG_Validate_index(index, algorithm->numberOfOperators);
	BORG_Validate_pointer(variation);

	algorithm->operators[index] = variation;
}

void BORG_Algorithm_update(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	int i;
	int index;
	double sum = 0.0;
	BORG_Archive archive = algorithm->archive;

	for (i=0; i<algorithm->numberOfOperators; i++) {
		algorithm->operators[i]->selectionProbability = 1.0;
	}

	if ((algorithm->probabilityMode == PROBABILITIES_DEFAULT) ||
		(algorithm->probabilityMode == PROBABILITIES_BOTH) ||
		(algorithm->probabilityMode == PROBABILITIES_ADAPTIVE))  {
		BORG_Entry entry = archive->head;

		while (entry) {
			if (entry->solution->operatorIndex >= 0) {
				algorithm->operators[entry->solution->operatorIndex]->selectionProbability++;
			}

			entry = entry->next;
		}
	}

	if ((algorithm->probabilityMode == PROBABILITIES_RECENCY) ||
		(algorithm->probabilityMode == PROBABILITIES_BOTH)) {
		for (i=0; i<archive->recencyListSize; i++) {
			if (archive->recencyList[i] >= 0) {
				algorithm->operators[archive->recencyList[i]]->selectionProbability++;
			}
		}
	} else if (algorithm->probabilityMode == PROBABILITIES_ADAPTIVE) {
		for (i=0; i<archive->recencyListSize-archive->size; i++) {
			index = archive->recencyListPosition-i;

			if (index < 0) {
				index += archive->recencyListSize;
			}

			if (archive->recencyList[index] >= 0) {
				algorithm->operators[archive->recencyList[index]]->selectionProbability++;
			}
		}
	}

	for (i=0; i<algorithm->numberOfOperators; i++) {
		sum += algorithm->operators[i]->selectionProbability;
	}

	BORG_Debug("Update (%d NFE) - Probabilities:", algorithm->numberOfEvaluations);

	for (i=0; i<algorithm->numberOfOperators; i++) {
		algorithm->operators[i]->selectionProbability /= sum;
		BORG_Debug(" %lg", algorithm->operators[i]->selectionProbability);
	}

	BORG_Debug("; Count: %lg\n", sum);
}

int BORG_Algorithm_select(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	int i;
	double rand = BORG_Random_uniform(0.0, 1.0);
	double sum = 0.0;
	algorithm->operationsSinceLastUpdate++;
	
	if (algorithm->operationsSinceLastUpdate >= algorithm->updateInterval) {
		algorithm->operationsSinceLastUpdate = 0;
		BORG_Algorithm_update(algorithm);
	}

	for (i=0; i<algorithm->numberOfOperators; i++) {
		sum += algorithm->operators[i]->selectionProbability;

		if (sum > rand) {
			return i;
		}
	}

	return algorithm->numberOfOperators-1;
}

void BORG_Algorithm_shuffle(int numberOfParents, BORG_Solution* parents) {
	BORG_Validate_positive(numberOfParents);
	BORG_Validate_pointer(parents);

	int i;
	int j;
	BORG_Solution temp;

	for (i=numberOfParents-1; i>=1; i--) {
		j = BORG_Random_int(i+1);

		if (i != j) {
			temp = parents[i];
			parents[i] = parents[j];
			parents[j] = temp;
		}
	}
}

int BORG_Algorithm_get_mutation_index(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	return algorithm->baseMutationIndex;
}

void BORG_Algorithm_set_max_mutation_index(BORG_Algorithm algorithm, int maxMutationIndex) {
	BORG_Validate_pointer(algorithm);

	algorithm->maxMutationIndex = maxMutationIndex;
}

int BORG_Algorithm_get_population_size(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	return algorithm->population->size;
}

int BORG_Algorithm_get_archive_size(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	return algorithm->archive->size;
}

void BORG_Algorithm_injection(BORG_Algorithm algorithm, int size) {
	BORG_Validate_pointer(algorithm);

	int i;
	int remaining = size;
	BORG_Population population = algorithm->population;
	BORG_Archive archive = algorithm->archive;
	BORG_Operator injectionOperator = BORG_Operator_create("UM", 1, 1, 1, BORG_Operator_UM);
	int numberOfParents = injectionOperator->numberOfParents;
	int numberOfOffspring = injectionOperator->numberOfOffspring;
	BORG_Solution* parents = (BORG_Solution*)calloc(numberOfParents, sizeof(struct BORG_Solution_t*));
	BORG_Solution* offspring = (BORG_Solution*)calloc(numberOfOffspring, sizeof(struct BORG_Solution_t*));

	BORG_Validate_malloc(parents);
	BORG_Validate_malloc(offspring);

	while (remaining > 0) {
		if (algorithm->restartMode == RESTART_DEFAULT) {
			BORG_Operator_set_parameter(injectionOperator, 0, 1.0 / algorithm->problem->numberOfVariables);
		} else if (algorithm->restartMode == RESTART_RANDOM) {
			BORG_Operator_set_parameter(injectionOperator, 0, 1.0);
		} else if (algorithm->restartMode == RESTART_RAMPED) {
			BORG_Operator_set_parameter(injectionOperator, 0, remaining / (double)size);
		} else if (algorithm->restartMode == RESTART_ADAPTIVE) {
			double d = 0.0;

			if (algorithm->maxMutationIndex > 0) {
				d = algorithm->baseMutationIndex/(double)algorithm->maxMutationIndex;
			}

			double rate = d + (1.0 - d)/algorithm->problem->numberOfVariables;

			BORG_Operator_set_parameter(injectionOperator, 0, rate);
		} else if (algorithm->restartMode == RESTART_INVERTED) {
			double d = 0.0;

			if (algorithm->maxMutationIndex > 0) {
				d = (algorithm->maxMutationIndex-algorithm->baseMutationIndex)/(double)algorithm->maxMutationIndex;
			}

			double rate = d + (1.0 - d)/algorithm->problem->numberOfVariables;

			BORG_Operator_set_parameter(injectionOperator, 0, rate);
		} else {
			BORG_Error("unknown restart mode\n");
		}

		for (i=0; i<numberOfParents; i++) {
			parents[i] = BORG_Archive_select(archive);
		}

		BORG_Operator_apply(injectionOperator, parents, offspring);

		for (i=0; i<numberOfOffspring; i++) {
			offspring[i]->operatorIndex = -1;
			BORG_Solution_evaluate(offspring[i]);
			BORG_Population_add(population, offspring[i]);
			BORG_Archive_add(archive, offspring[i]);
			BORG_Solution_destroy(offspring[i]);

			algorithm->numberOfEvaluations++;
			remaining--;
		}
	}

	BORG_Operator_destroy(injectionOperator);
	free(parents);
	free(offspring);
}

void BORG_Algorithm_restart(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);


	BORG_Population population = algorithm->population;
	BORG_Archive archive = algorithm->archive;
	int newPopulationSize = (int)(algorithm->populationRatio * archive->size);

	if (newPopulationSize < algorithm->minimumPopulationSize) {
		newPopulationSize = algorithm->minimumPopulationSize;
	} else if (newPopulationSize > algorithm->maximumPopulationSize) {
		newPopulationSize = algorithm->maximumPopulationSize;
	}

	BORG_Debug("Restart (%d NFE) - New Population Size: %d\n",
		algorithm->numberOfEvaluations, newPopulationSize);

	/* clear and resize population */
	BORG_Population_reset(population, newPopulationSize);

	/* fill population with all individuals from the archive */
	BORG_Entry entry = archive->head;

	while (entry) {
		BORG_Population_add(population, entry->solution);
		entry = entry->next;
		newPopulationSize--;
	}

	/* inject remaining individuals */
	BORG_Algorithm_injection(algorithm, newPopulationSize);

	/* adjust selection pressure */
	BORG_Algorithm_update_tournament_size(algorithm);

	algorithm->evaluationsAtLastRestart = algorithm->numberOfEvaluations;
	algorithm->numberOfRestarts++;
}

int BORG_Algorithm_check(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	int isRestart = 0;
	BORG_Population population = algorithm->population;
	BORG_Archive archive = algorithm->archive;
	double targetSize = algorithm->populationRatio * archive->size;

	if (algorithm->numberOfEvaluations-algorithm->evaluationsAtLastRestart > algorithm->maximumWindowSize) {
		BORG_Debug("Check (%d NFE) - Window Limit\n", algorithm->numberOfEvaluations);
		isRestart = 1;
	} else if ((targetSize >= algorithm->minimumPopulationSize) &&
			(targetSize <= algorithm->maximumPopulationSize) && 
			(fabs(population->size - targetSize) > (0.25 * targetSize))) {
		BORG_Debug("Check (%d NFE) - Population to Archive Ratio\n", algorithm->numberOfEvaluations);
		isRestart = 1;
	} else if (archive->numberOfImprovements <= algorithm->improvementsAtLastCheck) {
		BORG_Debug("Check (%d NFE) - Number of Improvements\n", algorithm->numberOfEvaluations);
		isRestart = 1;
	} else {
		BORG_Debug("Check (%d NFE) - Ok\n", algorithm->numberOfEvaluations);
	}

	return isRestart;
}

void BORG_Algorithm_step(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	int i;
	BORG_Problem problem = algorithm->problem;
	BORG_Population population = algorithm->population;
	BORG_Archive archive = algorithm->archive;

	/* check if a restart is necessary */
	if ((algorithm->numberOfEvaluations-algorithm->evaluationsAtLastCheck) >= algorithm->windowSize) {
		if (BORG_Algorithm_check(algorithm)) {
			if (algorithm->restartedLastCheck) {
				algorithm->baseMutationIndex++;

				if (algorithm->baseMutationIndex > algorithm->maxMutationIndex) {
					algorithm->baseMutationIndex = algorithm->maxMutationIndex;
				}
			}

			BORG_Algorithm_restart(algorithm);

			algorithm->restartedLastCheck = 1;
		} else {
			if (algorithm->restartedLastCheck) {
				algorithm->baseMutationIndex--;

				if (algorithm->baseMutationIndex < 0) {
					algorithm->baseMutationIndex = 0;
				}
			}

			algorithm->restartedLastCheck = 0;	
		}

		algorithm->improvementsAtLastCheck = archive->numberOfImprovements;
		algorithm->evaluationsAtLastCheck = algorithm->numberOfEvaluations;
	}

	if (algorithm->numberOfEvaluations == 0) {
		/* initialization */
		BORG_Problem_validate(algorithm->problem);
		BORG_Algorithm_validate(algorithm);

		if (population->capacity != algorithm->initialPopulationSize) {
			BORG_Population_reset(population, algorithm->initialPopulationSize);
		}

		for (i=0; i<population->capacity; i++) {
			BORG_Solution solution = BORG_Solution_create(problem);
			BORG_Solution_initialize(solution);
			BORG_Solution_evaluate(solution);
			BORG_Population_add(population, solution);
			BORG_Archive_add(archive, solution);
			BORG_Solution_destroy(solution);

			algorithm->numberOfEvaluations++;
		}

		BORG_Algorithm_update(algorithm);
		BORG_Algorithm_update_tournament_size(algorithm);
	} else {
		/* iteration */
		int operatorIndex = BORG_Algorithm_select(algorithm);
		BORG_Operator variation = algorithm->operators[operatorIndex];
		int numberOfParents = variation->numberOfParents;
		int numberOfOffspring = variation->numberOfOffspring;
		BORG_Solution* parents = (BORG_Solution*)calloc(numberOfParents, sizeof(struct BORG_Solution_t*));
		BORG_Solution* offspring = (BORG_Solution*)calloc(numberOfOffspring, sizeof(struct BORG_Solution_t*));

		BORG_Validate_malloc(parents);
		BORG_Validate_malloc(offspring);

		if (archive->size <= 1) {
			BORG_Population_select(population, numberOfParents, parents, 
				algorithm->tournamentSize);
		} else {
			BORG_Population_select(population, numberOfParents-1, parents,
				algorithm->tournamentSize);
			parents[numberOfParents-1] = BORG_Archive_select(archive);
		}

		BORG_Algorithm_shuffle(numberOfParents, parents);
		BORG_Operator_apply(variation, parents, offspring);

		for (i=0; i<numberOfOffspring; i++) {
			offspring[i]->operatorIndex = operatorIndex;
			BORG_Solution_evaluate(offspring[i]);
			BORG_Population_add(population, offspring[i]);
			BORG_Archive_add(archive, offspring[i]);
			BORG_Solution_destroy(offspring[i]);

			algorithm->numberOfEvaluations++;
		}

		free(parents);
		free(offspring);
	}
}

int BORG_Algorithm_get_nfe(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	return algorithm->numberOfEvaluations;
}

BORG_Archive BORG_Algorithm_get_result(BORG_Algorithm algorithm) {
	BORG_Validate_pointer(algorithm);

	return BORG_Archive_clone(algorithm->archive);
}

BORG_Archive BORG_Algorithm_run(BORG_Problem problem, int maxEvaluations) {
	BORG_Validate_pointer(problem);
	BORG_Validate_positive(maxEvaluations);

	BORG_Operator pm = BORG_Operator_create("PM", 1, 1, 2, BORG_Operator_PM);
	BORG_Operator_set_parameter(pm, 0, 1.0 / problem->numberOfVariables);
	BORG_Operator_set_parameter(pm, 1, 20.0);

	BORG_Operator sbx = BORG_Operator_create("SBX", 2, 2, 2, BORG_Operator_SBX);
	BORG_Operator_set_parameter(sbx, 0, 1.0);
	BORG_Operator_set_parameter(sbx, 1, 15.0);
	BORG_Operator_set_mutation(sbx, pm);

	BORG_Operator de = BORG_Operator_create("DE", 4, 1, 2, BORG_Operator_DE);
	BORG_Operator_set_parameter(de, 0, 0.1);
	BORG_Operator_set_parameter(de, 1, 0.5);
	BORG_Operator_set_mutation(de, pm);

	BORG_Operator um = BORG_Operator_create("UM", 1, 1, 1, BORG_Operator_UM);
	BORG_Operator_set_parameter(um, 0, 1.0 / problem->numberOfVariables);

	BORG_Operator spx = BORG_Operator_create("SPX", 10, 2, 1, BORG_Operator_SPX);
	BORG_Operator_set_parameter(spx, 0, 3.0);

	BORG_Operator pcx = BORG_Operator_create("PCX", 10, 2, 2, BORG_Operator_PCX);
	BORG_Operator_set_parameter(pcx, 0, 0.1);
	BORG_Operator_set_parameter(pcx, 1, 0.1);

	BORG_Operator undx = BORG_Operator_create("UNDX", 10, 2, 2, BORG_Operator_UNDX);
	BORG_Operator_set_parameter(undx, 0, 0.5);
	BORG_Operator_set_parameter(undx, 1, 0.35);

	BORG_Algorithm algorithm = BORG_Algorithm_create(problem, 6);
	BORG_Algorithm_set_operator(algorithm, 0, sbx);
	BORG_Algorithm_set_operator(algorithm, 1, de);
	BORG_Algorithm_set_operator(algorithm, 2, pcx);
	BORG_Algorithm_set_operator(algorithm, 3, spx);
	BORG_Algorithm_set_operator(algorithm, 4, undx);
	BORG_Algorithm_set_operator(algorithm, 5, um);

	while (BORG_Algorithm_get_nfe(algorithm) < maxEvaluations) {
		BORG_Algorithm_step(algorithm);
	}

	BORG_Archive result = BORG_Algorithm_get_result(algorithm);

	BORG_Operator_destroy(sbx);
	BORG_Operator_destroy(de);
	BORG_Operator_destroy(pm);
	BORG_Operator_destroy(um);
	BORG_Operator_destroy(spx);
	BORG_Operator_destroy(pcx);
	BORG_Operator_destroy(undx);
	BORG_Algorithm_destroy(algorithm);

	return result;
}

void BORG_Algorithm_checkpoint(BORG_Algorithm algorithm, FILE *file) {
	int i;
	int j;

	BORG_Validate_pointer(algorithm);
	BORG_Validate_file(file);

	BORG_Problem problem = algorithm->problem;
	BORG_Population population = algorithm->population;
	BORG_Archive archive = algorithm->archive;

	/* save problem settings */
	fprintf(file, "Version: %d\n", BORG_VERSION_NUMBER);
	fprintf(file, "Number of Variables: %d\n", problem->numberOfVariables);
	fprintf(file, "Number of Objectives: %d\n", problem->numberOfObjectives);
	fprintf(file, "Number of Constraints: %d\n", problem->numberOfConstraints);
	fprintf(file, "Lower Bounds:");
	
	for (i = 0; i < problem->numberOfVariables; i++) {
		fprintf(file, " %.*g", BORG_DIGITS, problem->lowerBounds[i]);
	}

	fprintf(file, "\n");
	fprintf(file, "Upper Bounds:");

	for (i = 0; i < problem->numberOfVariables; i++) {
		fprintf(file, " %.*g", BORG_DIGITS, problem->upperBounds[i]);
	}

	fprintf(file, "\n");
	fprintf(file, "Epsilons:");

	for (i = 0; i < problem->numberOfObjectives; i++) {
		fprintf(file, " %.*g", BORG_DIGITS, problem->epsilons[i]);
	}

	fprintf(file, "\n");
	fprintf(file, "Number of Operators: %d\n", algorithm->numberOfOperators);
	fprintf(file, "Operator Selection Probabilities:");

	for (i = 0; i < algorithm->numberOfOperators; i++) {
		fprintf(file, " %.*g", BORG_DIGITS, algorithm->operators[i]->selectionProbability);
	}

	fprintf(file, "\n");
	fprintf(file, "Number of Evaluations: %d\n", algorithm->numberOfEvaluations);
	fprintf(file, "Tournament Size: %d\n", algorithm->tournamentSize);
	fprintf(file, "Window Size: %d\n", algorithm->windowSize);
	fprintf(file, "Maximum Window Size: %d\n", algorithm->maximumWindowSize);
	fprintf(file, "Initial Population Size: %d\n", algorithm->initialPopulationSize);
	fprintf(file, "Minimum Population Size: %d\n", algorithm->minimumPopulationSize);
	fprintf(file, "Maximum Population Size: %d\n", algorithm->maximumPopulationSize);
	fprintf(file, "Population Ratio: %.*g\n", BORG_DIGITS, algorithm->populationRatio);
	fprintf(file, "Selection Ratio: %.*g\n", BORG_DIGITS, algorithm->selectionRatio);
	fprintf(file, "Evaluations At Last Check: %d\n", algorithm->evaluationsAtLastCheck);
	fprintf(file, "Evaluations At Last Restart: %d\n", algorithm->evaluationsAtLastRestart);
	fprintf(file, "Update Interval: %d\n", algorithm->updateInterval);
	fprintf(file, "Operations Since Last Update: %d\n", algorithm->operationsSinceLastUpdate);
	fprintf(file, "Improvements At Last Check: %d\n", algorithm->improvementsAtLastCheck);
	fprintf(file, "Restart Mode: %d\n", algorithm->restartMode);
	fprintf(file, "Restarted Last Check: %d\n", algorithm->restartedLastCheck);
	fprintf(file, "Base Mutation Index: %d\n", algorithm->baseMutationIndex);
	fprintf(file, "Max Mutation Index: %d\n", algorithm->maxMutationIndex);
	fprintf(file, "Probability Mode: %d\n", algorithm->probabilityMode);
	fprintf(file, "Number of Restarts: %d\n", algorithm->numberOfRestarts);

	/* save population members */
	fprintf(file, "Population Size: %d\n", population->size);
	fprintf(file, "Population Capacity: %d\n", population->capacity);
	fprintf(file, "Population:\n");

	for (i = 0; i < population->size; i++) {
		BORG_Solution solution = population->members[i];

		for (j = 0; j < problem->numberOfVariables; j++) {
			fprintf(file, " %.*g", BORG_DIGITS, solution->variables[j]);
		}

		for (j = 0; j < problem->numberOfObjectives; j++) {
			fprintf(file, " %.*g", BORG_DIGITS, solution->objectives[j]);
		}

		for (j = 0; j < problem->numberOfConstraints; j++) {
			fprintf(file, " %.*g", BORG_DIGITS, solution->constraints[j]);
		}

		fprintf(file, " %d\n", solution->operatorIndex);
	}

	/* save archive members */
	fprintf(file, "Archive Size: %d\n", archive->size);
	fprintf(file, "Archive:\n");

	BORG_Entry entry = archive->head;

	while (entry) {
		BORG_Solution solution = entry->solution;

		for (j = 0; j < problem->numberOfVariables; j++) {
			fprintf(file, " %.*g", BORG_DIGITS, solution->variables[j]);
		}

		for (j = 0; j < problem->numberOfObjectives; j++) {
			fprintf(file, " %.*g", BORG_DIGITS, solution->objectives[j]);
		}

		for (j = 0; j < problem->numberOfConstraints; j++) {
			fprintf(file, " %.*g", BORG_DIGITS, solution->constraints[j]);
		}

		fprintf(file, " %d\n", solution->operatorIndex);
		entry = entry->next;
	}

	fprintf(file, "Number of Improvements: %d\n", archive->numberOfImprovements);
	fprintf(file, "Recency List Size: %d\n", archive->recencyListSize);
	fprintf(file, "Recency List Position: %d\n", archive->recencyListPosition);
	fprintf(file, "Recency List:");

	for (i = 0; i < archive->recencyListSize; i++) {
		fprintf(file, " %d", archive->recencyList[i]);
	}

	fprintf(file, "\n");

	/* save random number generator state */
	int state_length;
	unsigned long* state;
	int* state_index;

	get_state(&state_length, &state, &state_index);

	fprintf(file, "RNG State Length: %d\n", state_length);
	fprintf(file, "RNG State:");

	for (i = 0; i < state_length; i++) {
		fprintf(file, " %lu", state[i]);
	}

	fprintf(file, "\n");
	fprintf(file, "RNG Index: %d\n", *state_index);
	fprintf(file, "Next Guassian: %.*g\n", BORG_DIGITS, nextNextGaussian);
	fprintf(file, "Have Next Guassian: %d\n", haveNextNextGaussian);
}

void BORG_Check_scan(int expected, int result) {
	if (result == EOF) {
		BORG_Error("End of file reached before reading entire checkpoint file\n");
	} else if (result != expected) {
		BORG_Error("Unable to parse value from checkpoint file\n");
	}
}

void BORG_Check_scan1(int result) {
	BORG_Check_scan(1, result);
}

void BORG_Check_scan0(int result) {
	BORG_Check_scan(0, result);
}

void BORG_Algorithm_restore(BORG_Algorithm algorithm, FILE *file) {
	int i;
	int j;
	int readInt;
	double readDouble;
	int populationSize;
	int populationCapacity;
	int archiveSize;
	BORG_Entry entry;

	BORG_Validate_pointer(algorithm);
	BORG_Validate_file(file);

	BORG_Problem problem = algorithm->problem;
	BORG_Population population = algorithm->population;
	BORG_Archive archive = algorithm->archive;

	/* if the problem settings do not match, fail */
	BORG_Check_scan1(fscanf(file, "Version: %d\n", &readInt));

	if (readInt != BORG_VERSION_NUMBER) {
		BORG_Error("Checkpoint file from different version\n");
	}

	BORG_Check_scan1(fscanf(file, "Number of Variables: %d\n", &readInt));

	if (readInt != problem->numberOfVariables) {
		BORG_Error("Checkpoint file does not match number of variabes\n");
	}

	BORG_Check_scan1(fscanf(file, "Number of Objectives: %d\n", &readInt));

	if (readInt != problem->numberOfObjectives) {
		BORG_Error("Checkpoint file does not match number of objectives\n");
	}

	BORG_Check_scan1(fscanf(file, "Number of Constraints: %d\n", &readInt));

	if (readInt != problem->numberOfConstraints) {
		BORG_Error("Checkpoint file does not match number of constraints\n");
	}

	BORG_Check_scan0(fscanf(file, "Lower Bounds:"));

	for (i = 0; i < problem->numberOfVariables; i++) {
		BORG_Check_scan1(fscanf(file, " %lg", &readDouble));

		if (readDouble != problem->lowerBounds[i]) {
			BORG_Error("Checkpoint file does not match lower bounds\n");
		}
	}

	BORG_Check_scan0(fscanf(file, "\n"));
	BORG_Check_scan0(fscanf(file, "Upper Bounds:"));

	for (i = 0; i < problem->numberOfVariables; i++) {
		BORG_Check_scan1(fscanf(file, " %lg", &readDouble));

		if (readDouble != problem->upperBounds[i]) {
			BORG_Error("Checkpoint file does not match upper bounds\n");
		}
	}

	BORG_Check_scan0(fscanf(file, "\n"));
	BORG_Check_scan0(fscanf(file, "Epsilons:"));

	for (i = 0; i < problem->numberOfObjectives; i++) {
		BORG_Check_scan1(fscanf(file, " %lg", &readDouble));

		if (readDouble != problem->epsilons[i]) {
			BORG_Error("Checkpoint file does not match epsilon values\n");
		}
	}

	BORG_Check_scan0(fscanf(file, "\n"));
	BORG_Check_scan1(fscanf(file, "Number of Operators: %d\n", &readInt));

	if (readInt != algorithm->numberOfOperators) {
		BORG_Error("Checkpoint file does not match number of operators\n");
	}

	BORG_Check_scan0(fscanf(file, "Operator Selection Probabilities:"));

	for (i = 0; i < algorithm->numberOfOperators; i++) {
		BORG_Check_scan1(fscanf(file, " %lg", &algorithm->operators[i]->selectionProbability));
	}

	BORG_Check_scan0(fscanf(file, "\n"));

	/* load remaining settings, some will be ignored */
	BORG_Check_scan1(fscanf(file, "Number of Evaluations: %d\n", &algorithm->numberOfEvaluations));
	BORG_Check_scan1(fscanf(file, "Tournament Size: %d\n", &algorithm->tournamentSize));
	BORG_Check_scan1(fscanf(file, "Window Size: %d\n", &algorithm->windowSize));
	BORG_Check_scan1(fscanf(file, "Maximum Window Size: %d\n", &readInt));
	BORG_Check_scan1(fscanf(file, "Initial Population Size: %d\n", &readInt));
	BORG_Check_scan1(fscanf(file, "Minimum Population Size: %d\n", &readInt));
	BORG_Check_scan1(fscanf(file, "Maximum Population Size: %d\n", &readInt));
	BORG_Check_scan1(fscanf(file, "Population Ratio: %lg\n", &readDouble));
	BORG_Check_scan1(fscanf(file, "Selection Ratio: %lg\n", &readDouble));
	BORG_Check_scan1(fscanf(file, "Evaluations At Last Check: %d\n", &algorithm->evaluationsAtLastCheck));
	BORG_Check_scan1(fscanf(file, "Evaluations At Last Restart: %d\n", &algorithm->evaluationsAtLastRestart));
	BORG_Check_scan1(fscanf(file, "Update Interval: %d\n", &readInt));
	BORG_Check_scan1(fscanf(file, "Operations Since Last Update: %d\n", &algorithm->operationsSinceLastUpdate));
	BORG_Check_scan1(fscanf(file, "Improvements At Last Check: %d\n", &algorithm->improvementsAtLastCheck));
	BORG_Check_scan1(fscanf(file, "Restart Mode: %d\n", &readInt));
	BORG_Check_scan1(fscanf(file, "Restarted Last Check: %d\n", &algorithm->restartedLastCheck));
	BORG_Check_scan1(fscanf(file, "Base Mutation Index: %d\n", &algorithm->baseMutationIndex));
	BORG_Check_scan1(fscanf(file, "Max Mutation Index: %d\n", &readInt));
	BORG_Check_scan1(fscanf(file, "Probability Mode: %d\n", &readInt));
	BORG_Check_scan1(fscanf(file, "Number of Restarts: %d\n", &algorithm->numberOfRestarts));

	/* load the saved population */
	BORG_Check_scan1(fscanf(file, "Population Size: %d\n", &populationSize));
	BORG_Check_scan1(fscanf(file, "Population Capacity: %d\n", &populationCapacity));
	BORG_Check_scan0(fscanf(file, "Population:\n"));

	BORG_Population_reset(population, populationCapacity);

	for (i = 0; i < populationSize; i++) {
		BORG_Solution solution = BORG_Solution_create(problem);

		for (j = 0; j < problem->numberOfVariables; j++) {
			BORG_Check_scan1(fscanf(file, " %lg", &solution->variables[j]));
		}

		for (j = 0; j < problem->numberOfObjectives; j++) {
			BORG_Check_scan1(fscanf(file, " %lg", &solution->objectives[j]));
		}

		for (j = 0; j < problem->numberOfConstraints; j++) {
			BORG_Check_scan1(fscanf(file, " %lg", &solution->constraints[j]));
		}

		BORG_Check_scan1(fscanf(file, " %d\n", &solution->operatorIndex));
		BORG_Population_add(population, solution);
	}

	/* clear any solutions in the existing archive */
	entry = archive->head;

	while (entry) {
		entry = BORG_Archive_remove(entry);
	}

	/* load the saved archive */
	BORG_Check_scan1(fscanf(file, "Archive Size: %d\n", &archiveSize));
	BORG_Check_scan0(fscanf(file, "Archive:\n"));

	for (i = 0; i < archiveSize; i++) {
		BORG_Solution solution = BORG_Solution_create(problem);

		for (j = 0; j < problem->numberOfVariables; j++) {
			BORG_Check_scan1(fscanf(file, " %lg", &solution->variables[j]));
		}

		for (j = 0; j < problem->numberOfObjectives; j++) {
			BORG_Check_scan1(fscanf(file, " %lg", &solution->objectives[j]));
		}

		for (j = 0; j < problem->numberOfConstraints; j++) {
			BORG_Check_scan1(fscanf(file, " %lg", &solution->constraints[j]));
		}

		BORG_Check_scan1(fscanf(file, " %d\n", &solution->operatorIndex));
		BORG_Archive_add(archive, solution);
	}

	BORG_Check_scan1(fscanf(file, "Number of Improvements: %d\n", &archive->numberOfImprovements));
	BORG_Check_scan1(fscanf(file, "Recency List Size: %d\n", &readInt));
	BORG_Check_scan1(fscanf(file, "Recency List Position: %d\n", &archive->recencyListPosition));
	BORG_Check_scan0(fscanf(file, "Recency List:"));

	for (i = 0; i < archive->recencyListSize; i++) {
		BORG_Check_scan1(fscanf(file, " %d", &archive->recencyList[i]));
	}

	BORG_Check_scan0(fscanf(file, "\n"));

	/* load random number generator state */
	int state_length;
	unsigned long* state;
	int* state_index;

	get_state(&state_length, &state, &state_index);

	BORG_Check_scan1(fscanf(file, "RNG State Length: %d\n", &readInt));

	if (readInt != state_length) {
		BORG_Error("Checkpoint file has invalid random number generator state\n");
	}

	BORG_Check_scan0(fscanf(file, "RNG State:"));

	for (i = 0; i < state_length; i++) {
		BORG_Check_scan1(fscanf(file, " %lu", &state[i]));
	}
	
	BORG_Check_scan0(fscanf(file, "\n"));
	BORG_Check_scan1(fscanf(file, "RNG Index: %d\n", state_index));
	BORG_Check_scan1(fscanf(file, "Next Guassian: %lg\n", &nextNextGaussian));
	BORG_Check_scan1(fscanf(file, "Have Next Guassian: %d\n", &haveNextNextGaussian));
}

