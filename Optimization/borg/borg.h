/* Copyright 2012-2014 The Pennsylvania State University
 *
 * This software was written by David Hadka and others.
 * 
 * The use, modification and distribution of this software is governed by the
 * The Pennsylvania State University Research and Educational Use License.
 * You should have received a copy of this license along with this program.
 * If not, contact <dmh309@psu.edu>.
 */
#ifndef _BORG_H_
#define _BORG_H_

/**
 * Flags the methods comprising the public API.
 */
#ifdef _MSC_VER
	#ifdef BORG_EXPORTS
		#define BORG_API __declspec(dllexport)
	#else
		#define BORG_API __declspec(dllimport)
	#endif
#else
	#define BORG_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Human-readable version number.
 */
#define BORG_VERSION "1.8"

/**
 * Numeric version number that is compatible with comparison operators (e.g., <, =, >).
 */
#define BORG_VERSION_NUMBER 108000

/**
 * The minimum precision for numeric values (i.e., machine precision).
 */
#define EPSILON 1e-9

/**
 * Opaque structure representing a problem.
 */
typedef struct BORG_Problem_t *BORG_Problem;

/**
 * Opaque structure representing a solution.
 */
typedef struct BORG_Solution_t *BORG_Solution;

/**
 * Opaque structure representing a crossover, mutation or composite operator.
 */
typedef struct BORG_Operator_t *BORG_Operator;

/**
 * Opaque structure representing a fixed-size population of solutions.
 */
typedef struct BORG_Population_t *BORG_Population;

/**
 * Opaque structure representing an epsilon-dominance archive.
 */
typedef struct BORG_Archive_t *BORG_Archive;

/**
 * Opaque structure representing the Borg algorithm.
 */
typedef struct BORG_Algorithm_t *BORG_Algorithm;

/**
 * The dominance relation flags.  Note that the values distinguish between
 * dominance that occurs within and outside the same epsilon-box.
 */
typedef enum BORG_Dominance {
	DOMINATES = -2,
	DOMINATES_SAME_BOX = -1,
	NONDOMINATED = 0,
	DOMINATED_SAME_BOX = 1,
	DOMINATED = 2
} BORG_Dominance;

/**
 * The available restart options.
 */
typedef enum BORG_Restart {
	RESTART_DEFAULT = 0,
	RESTART_RANDOM = 1,
	RESTART_RAMPED = 2,
	RESTART_ADAPTIVE = 3,
	RESTART_INVERTED = 4
} BORG_Restart;

/**
 * The available operator probability options.
 */
typedef enum BORG_Probabilities {
	PROBABILITIES_DEFAULT = 0,
	PROBABILITIES_RECENCY = 1,
	PROBABILITIES_BOTH = 2,
	PROBABILITIES_ADAPTIVE = 3
} BORG_Probabilities;

/**
 * Prints the copyright and other legal notices to the specified file.
 */
BORG_API void BORG_Copyright(FILE* fp);

/**
 * Aborts with an error if the file pointer is invalid.
 */
BORG_API void BORG_Validate_file(FILE* fp);

/**
 * Aborts with an error if the pointer is invalid.
 */
BORG_API void BORG_Validate_pointer(const void* ptr);

/**
 * Aborts with an error if the index is out of bounds.
 */
BORG_API void BORG_Validate_index(int index, int size);

/**
 * Aborts with an error if the memory allocation failed.
 */
BORG_API void BORG_Validate_malloc(void* ptr);

/**
 * Aborts with an error if the argument is negative.
 */
BORG_API void BORG_Validate_positive(double value);

/**
 * Enables the logging of debug statements to stderr.
 */
BORG_API void BORG_Debug_on();

/**
 * Disables the logging of debug statements to stderr.
 */
BORG_API void BORG_Debug_off();

/**
 * Sets the name of this node displayed in the debug output, or
 * NULL if no name is displayed.
 */
BORG_API void BORG_Debug_set_name(char* name);

/**
 * Displays debugging information.
 */
BORG_API void BORG_Debug(const char* format, ...);

/**
 * Displays an error message and aborts the program.
 */
#ifdef __GNUC__
__attribute__((noreturn))
#endif
void BORG_Error(const char* format, ...);

/**
 * Creates a new problem with the specified number of decision variables, objectives and
 * constraints.  The evaluation function for the problem is also provided, which inputs
 * as arguments the arrays holding the decision variables, objectives and constraints.
 */
BORG_API BORG_Problem BORG_Problem_create(
		int numberOfVariables,
		int numberOfObjectives,
		int numberOfConstraints,
		void (*function)(double*, double*, double*));

/**
 * Destroys a problem.  No methods may be invoked on the problem after it is destroyed.
 */
BORG_API void BORG_Problem_destroy(BORG_Problem problem);

/**
 * Sets the epsilon value for the specified objective index.
 */
BORG_API void BORG_Problem_set_epsilon(BORG_Problem problem, int index, double epsilon);

/**
 * Sets the objective name for the specified objective index.
 */
BORG_API void BORG_Problem_set_name(BORG_Problem problem, int index, const char* name);

/**
 * Sets all epsilon values.
 */
BORG_API void BORG_Problem_set_epsilons(BORG_Problem problem, double* epsilons);

/**
 * Sets the lower and upper bounds for the decision variable specified by the index.
 */
BORG_API void BORG_Problem_set_bounds(BORG_Problem problem, int index, double lowerBound, double upperBound);

/**
 * Returns the number of decision variables defined by the problem.
 */
BORG_API int BORG_Problem_number_of_variables(BORG_Problem problem);

/**
 * Returns the number of objectives defined by the problem.
 */
BORG_API int BORG_Problem_number_of_objectives(BORG_Problem problem);

/**
 * Returns the number of constraints defined by the problem.
 */
BORG_API int BORG_Problem_number_of_constraints(BORG_Problem problem);

/**
 * Creates a new, uninitialized solution for the specified problem.
 */
BORG_API BORG_Solution BORG_Solution_create(BORG_Problem problem);

/**
 * Destroys a solution.  No methods may be invoked on the solution after it is destroyed.
 */
BORG_API void BORG_Solution_destroy(BORG_Solution solution);

/**
 * Returns a clone of the specified solution.  The cloned solution represents the same value
 * as the original solution, but occupies a new memory location.
 */
BORG_API BORG_Solution BORG_Solution_clone(BORG_Solution original);

/**
 * Returns the decision variable value at the specified index for the given solution.
 */
BORG_API double BORG_Solution_get_variable(BORG_Solution solution, int index);

/**
 * Returns the objective value at the specified index for the given solution.
 */
BORG_API double BORG_Solution_get_objective(BORG_Solution solution, int index);

/**
 * Returns the constraint value at the specified index for the given solution.
 */
BORG_API double BORG_Solution_get_constraint(BORG_Solution solution, int index);

/**
 * Sets the value of the decision variable specified by the index.
 */
BORG_API void BORG_Solution_set_variable(BORG_Solution solution, int index, double value);

/**
 * Sets the value of all decision variables.
 */
BORG_API void BORG_Solution_set_variables(BORG_Solution solution, double* variables);

/**
 * Sets the objective value at the specified index for the given solution.  Note that
 * all objectives are assumed to be minimized.
 */
BORG_API void BORG_Solution_set_objective(BORG_Solution solution, int index, double value);

/**
 * Sets the value of the constraint at the specified index for the given solution.
 */
BORG_API void BORG_Solution_set_constraint(BORG_Solution solution, int index, double value);

/**
 * Returns the optimization problem associated with this solution.
 */
BORG_API BORG_Problem BORG_Solution_get_problem(BORG_Solution solution);

/**
 * Evaluates the solution, assigning all objectives and constraints.
 */
BORG_API void BORG_Solution_evaluate(BORG_Solution solution);

/**
 * Prints the specified solution to a file stream.
 */
BORG_API void BORG_Solution_print(BORG_Solution solution, FILE* fp);

/**
 * Randomly initializes the decision variables of the specified solution.
 */
BORG_API void BORG_Solution_initialize(BORG_Solution solution);

/**
 * Returns true if the solution violates any constraints; false otherwise.
 */
BORG_API int BORG_Solution_violates_constraints(BORG_Solution solution);

/**
 * Returns DOMINATES if solution1 Pareto dominates solution2; DOMINATED if solution2
 * Pareto dominates solution1; or NONDOMINATED otherwise.
 */
BORG_API BORG_Dominance BORG_Dominance_pareto(BORG_Solution solution1, BORG_Solution solution2);

/**
 * Returns DOMINATES if the epsilon-box vector of solution1 Pareto dominates solution2;
 * DOMINATED if the epsilon-box vector of solution2 Pareto dominates solution1;
 * DOMINATES_SAME_BOX if solution1 and solution2 reside in the same epsilon-box but solution1
 * is nearer to the optimal corner; DOMINATED_SAME_BOX if solution1 and solution2 reside in
 * the same epsilon-box but solution2 is nearer to the optimal corner; or NONDOMINATED
 * otherwise.
 */
BORG_API BORG_Dominance BORG_Dominance_epsilon(BORG_Solution solution1, BORG_Solution solution2);

/**
 * Returns DOMINATES if the aggregate constraint violation of solution1 is less than solution2;
 * DOMINATED if the aggregate constraint violation of solution2 is less than solution2; or
 * NONDOMINATED if solution1 and solution2 share equal aggregate constraint violations.
 */
BORG_API BORG_Dominance BORG_Dominance_constraints(BORG_Solution solution1, BORG_Solution solution2);

/**
 * Applies two dominance checks consecuatively, returning the result of the first if is not
 * NONDOMINATED; otherwise the result of the second dominance check is returned.
 */
BORG_API BORG_Dominance BORG_Dominance_compound(
		BORG_Solution solution1,
		BORG_Solution solution2,
		BORG_Dominance (*comparator1)(BORG_Solution, BORG_Solution),
		BORG_Dominance (*comparator2)(BORG_Solution, BORG_Solution));

/**
 * Seeds the internal pseudo-random number generator with the specified seed.
 */
BORG_API void BORG_Random_seed(unsigned long seed);

/**
 * Returns a random number sampled uniformly from within the specified lower and upper
 * bounds, inclusively.
 */
BORG_API double BORG_Random_uniform(double lowerBound, double upperBound);

/**
 * Returns a random integer sampled uniformly from the range [0, n).
 */
BORG_API int BORG_Random_int(int n);

/**
 * Returns a random number sampled from a Gaussian (normal) distribution with the specified
 * mean and standard deviation.
 */
BORG_API double BORG_Random_gaussian(double mean, double stdev);

/**
 * Creates a new operator with the specified number of parents, offspring and parameters.
 * Also specified is the function for performing the operation, which inputs as arguments
 * a reference to this operator, an array containing the parent solutions,
 * and an array to store the resulting offspring.
 */
BORG_API BORG_Operator BORG_Operator_create(
		const char* name,
		int numberOfParents,
		int numberOfOffspring,
		int numberOfParameters,
		void (*function)(BORG_Operator, BORG_Solution*, BORG_Solution*));

/**
 * Destroys an operator.  No methods may be invoked on an operator after it is destroyed.
 */
BORG_API void BORG_Operator_destroy(BORG_Operator variation);

/**
 * Returns the selection probability of the given operator.
 */
BORG_API double BORG_Operator_get_probability(BORG_Operator variation);

/**
 * Sets the value of the parameter specified by the index.
 */
BORG_API void BORG_Operator_set_parameter(BORG_Operator variation, int index, double parameter);

/**
 * Sets all parameter values for the specified operator.
 */
BORG_API void BORG_Operator_set_parameters(BORG_Operator variation, double* parameters);

/**
 * Returns the number of offspring produced by the specified operator.
 */
BORG_API int BORG_Operator_get_number_of_offspring(BORG_Operator variation);

/**
 * Returns the number of parents required by the specified operator.
 */
BORG_API int BORG_Operator_get_number_of_parents(BORG_Operator variation);

/**
 * Sets an additional mutation operator which is applied to each individual offspring.
 */
BORG_API void BORG_Operator_set_mutation(BORG_Operator variation, BORG_Operator mutation);

/**
 * Applies the specified variation operator and its optional mutation to the parents, filling the
 * offspring array with the result of the operator.
 */
BORG_API void BORG_Operator_apply(BORG_Operator variation, BORG_Solution* parents, BORG_Solution* offspring);

/**
 * The uniform mutation operator.
 */
BORG_API void BORG_Operator_UM(BORG_Operator um, BORG_Solution* parents, BORG_Solution* offspring);

/**
 * The simulated binary crossover operator.
 */
BORG_API void BORG_Operator_SBX(BORG_Operator sbx, BORG_Solution* parents, BORG_Solution* offspring);

/**
 * The polynomial mutation operator.
 */
BORG_API void BORG_Operator_PM(BORG_Operator pm, BORG_Solution* parents, BORG_Solution* offspring);

/**
 * The differential evolution operator.
 */
BORG_API void BORG_Operator_DE(BORG_Operator de, BORG_Solution* parents, BORG_Solution* offspring);

/**
 * The simplex crossover operator.
 */
BORG_API void BORG_Operator_SPX(BORG_Operator spx, BORG_Solution* parents, BORG_Solution* offspring);

/**
 * The parent centric crossover operator.
 */
BORG_API void BORG_Operator_PCX(BORG_Operator pcx, BORG_Solution* parents, BORG_Solution* offspring);

/**
 * The unimodal normal distribution crossover operator.
 */
BORG_API void BORG_Operator_UNDX(BORG_Operator undx, BORG_Solution* parents, BORG_Solution* offspring);

/**
 * Creates a new population of the specified capacity.
 */
BORG_API BORG_Population BORG_Population_create(int capacity);

/**
 * Destroys a population and all solutions it contains.  No methods may be invoked on a population
 * after it is destroyed.
 */
BORG_API void BORG_Population_destroy(BORG_Population population);

/**
 * Adds the solution to the population, replacing a random dominated member; or if no members are
 * dominated, a random nondominated member; otherwise the solution is discarded.
 */
BORG_API void BORG_Population_add(BORG_Population population, BORG_Solution solution);

/**
 * Clears and resizes the population.
 */
BORG_API void BORG_Population_reset(BORG_Population population, int newCapacity);

/**
 * Performs tournament selection on the population.
 */
BORG_API BORG_Solution BORG_Population_tournament(BORG_Population population, int tournamentSize);

/**
 * Selects one or more parents from the population.
 */
BORG_API void BORG_Population_select(BORG_Population population, int arity, BORG_Solution* parents, int tournamentSize);

/**
 * Creates a new epsilon-dominance archive.
 */
BORG_API BORG_Archive BORG_Archive_create();

/**
 * Returns a clone of the specified epsilon-dominance archive.  The solutions contained in
 * the archive are also cloned.
 */
BORG_API BORG_Archive BORG_Archive_clone(BORG_Archive original);

/**
 * Destroys an archive and all solutions it contains.  No methods may be invoked on the
 * archive after it is destroyed.
 */
BORG_API void BORG_Archive_destroy(BORG_Archive archive);

/**
 * Adds the specified solution to the epsilon-dominance archive.  The solution will not
 * be added if the solution is dominated by existing members of the archive.
 */
BORG_API void BORG_Archive_add(BORG_Archive archive, BORG_Solution solution);

/**
 * Returns the size of the epsilon-dominance archive.
 */
BORG_API int BORG_Archive_get_size(BORG_Archive archive);

/**
 * Returns the solution at the specified index in the epsilon-dominance archive.
 */
BORG_API BORG_Solution BORG_Archive_get(BORG_Archive archive, int index);

/**
 * Returns the number of epsilon-progress improvements recorded by the epsilon-dominance
 * archive.
 */
BORG_API int BORG_Archive_get_improvements(BORG_Archive archive);

/**
 * Selects one solution uniformly at random from the archive.
 */
BORG_API BORG_Solution BORG_Archive_select(BORG_Archive archive);

/**
 * Prints all solutions contained in the specified archive to a file stream.
 */
BORG_API void BORG_Archive_print(BORG_Archive archive, FILE* fp);

/**
 * Prints all feasible solutions in the specified archive to the specified result file.
 */
BORG_API void BORG_Archive_append(BORG_Archive archive, FILE* file);

/**
 * Creates a new instance of the Borg algorithm.
 */
BORG_API BORG_Algorithm BORG_Algorithm_create(BORG_Problem problem, int numberOfOperators);

/**
 * Destroys a Borg algorithm instance.  No methods may be invoked on the algorithm after it is
 * destroyed.
 */
BORG_API void BORG_Algorithm_destroy(BORG_Algorithm algorithm);

/**
 * Validates the algorithm settings, aborting the program if invalid settings are detected.
 */
BORG_API void BORG_Algorithm_validate(BORG_Algorithm algorithm);

/**
 * Returns the number of restarts that have occurred.
 */
BORG_API int BORG_Algorithm_get_number_restarts(BORG_Algorithm algorithm);

/**
 * Returns the number of epsilon-progress improvements that have occurred.
 */
BORG_API int BORG_Algorithm_get_number_improvements(BORG_Algorithm algorithm);

/**
 * Sets the minimum number of evaluations between epsilon-progress checks.
 */
BORG_API void BORG_Algorithm_set_window_size(BORG_Algorithm algorithm, int windowSize);

/**
 * Sets the maximum number of evaluations between epsilon-progress checks.
 */
BORG_API void BORG_Algorithm_set_maximum_window_size(BORG_Algorithm algorithm, int maxWindowSize);

/**
 * Sets the initial population size for adaptive population sizing.
 */
BORG_API void BORG_Algorithm_set_initial_population_size(BORG_Algorithm algorithm, int initialPopulationSize);

/**
 * Sets the minimum population size for adaptive population sizing.
 */
BORG_API void BORG_Algorithm_set_minimum_population_size(BORG_Algorithm algorithm, int minimumPopulationSize);

/**
 * Sets the maximum population size for adaptive population sizing.
 */
BORG_API void BORG_Algorithm_set_maximum_population_size(BORG_Algorithm algorithm, int maximumPopulationSize);

/**
 * Sets the population-to-archive ratio for adaptive population sizing.
 */
BORG_API void BORG_Algorithm_set_population_ratio(BORG_Algorithm algorithm, double populationRatio);

/**
 * Sets the ratio of the tournament selection size to the population size.
 */
BORG_API void BORG_Algorithm_set_selection_ratio(BORG_Algorithm algorithm, double selectionRatio);

/**
 * Sets the restart mode.
 */
BORG_API void BORG_Algorithm_set_restart_mode(BORG_Algorithm algorithm, BORG_Restart restartMode);

/**
 * Set the operator probability update mode.
 */
BORG_API void BORG_Algorithm_set_probability_mode(BORG_Algorithm algorithm, BORG_Probabilities probilityMode);

/**
 * Returns the current mutation index when adaptive mutation is used.
 */
BORG_API int BORG_Algorithm_get_mutation_index(BORG_Algorithm algorithm);

/**
 * Sets the max mutation index when adaptive mutation is used.  This determine the number of
 * intermediate mutation probabilities between 1/L and 100%.
 */
BORG_API void BORG_Algorithm_set_max_mutation_index(BORG_Algorithm algorithm, int maxMutationIndex);

/**
 * Assigns one of the variation operators used by the Borg algorithm.
 */
BORG_API void BORG_Algorithm_set_operator(BORG_Algorithm algorithm, int index, BORG_Operator variation);

/**
 * Updates the operator selection probabilities.
 */
BORG_API void BORG_Algorithm_update(BORG_Algorithm algorithm);

/**
 * Returns the index of the operator selected for auto-adaptive multioperator search.
 */
BORG_API int BORG_Algorithm_select(BORG_Algorithm algorithm);

/**
 * Randomly shuffles the array of solutions.
 */
BORG_API void BORG_Algorithm_shuffle(int numberOfParents, BORG_Solution* parents);

/**
 * Performs a restart with injection.
 */
BORG_API void BORG_Algorithm_restart(BORG_Algorithm algorithm);

/**
 * Checks if the conditions for epsilon-progress are satisfied, returning true
 * if a restart should occur; false otherwise.
 */
BORG_API int BORG_Algorithm_check(BORG_Algorithm algorithm);

/**
 * Performs one logical step of the Borg algorithm.
 */
BORG_API void BORG_Algorithm_step(BORG_Algorithm algorithm);

/**
 * Returns the number of objective function evaluations expended by the Borg algorithm.
 */
BORG_API int BORG_Algorithm_get_nfe(BORG_Algorithm algorithm);

/**
 * Returns the size of the population.
 */
BORG_API int BORG_Algorithm_get_population_size(BORG_Algorithm algorithm);

/**
 * Returns the size of the archive.
 */
BORG_API int BORG_Algorithm_get_archive_size(BORG_Algorithm algorithm);

/**
 * Returns the Pareto approximation set discovered by the Borg algorithm.
 */
BORG_API BORG_Archive BORG_Algorithm_get_result(BORG_Algorithm algorithm);

/**
 * Performs a complete run of the Borg algorithm.
 */
BORG_API BORG_Archive BORG_Algorithm_run(BORG_Problem problem, int maxEvaluations);

/**
 * Save a checkpoint file.
 */
BORG_API void BORG_Algorithm_checkpoint(BORG_Algorithm algorithm, FILE* file);

/**
 * Loads a checkpoint file.
 */
BORG_API void BORG_Algorithm_restore(BORG_Algorithm algorithm, FILE* file);


#ifdef __cplusplus
}
#endif

#endif
