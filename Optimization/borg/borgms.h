/* Copyright 2012-2014 The Pennsylvania State University
 *
 * This software was written by David Hadka and others.
 * 
 * The use, modification and distribution of this software is governed by the
 * The Pennsylvania State University Research and Educational Use License.
 * You should have received a copy of this license along with this program.
 * If not, contact <dmh309@psu.edu>.
 */
#ifndef _BORG_MSDRIVER_H_
#define _BORG_MSDRIVER_H_

#ifdef _BORG_H_
#error do not include borg.h if you are including borgms.h
#endif

#include "borg.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Enumeration of initialization strategies.
 */
typedef enum BORG_Initialization {
	INITIALIZATION_UNIFORM,
	INITIALIZATION_LATIN,
} BORG_Initialization;

/**
 * Returns the difference between two high-resolution time values.  The
 * differences is reported in microseconds.
 */
long BORG_Timer_diff(struct timeval* start, struct timeval* end);

/**
 * Sets the referenced timeval structure to the current high-resolution
 * time.
 */
void BORG_Timer_now(struct timeval* time);

/**
 * Output timing data to files with the specified prefix.  The filename
 * should contain one %%d to be replaced by the world rank.
 */
void BORG_Algorithm_output_timing(char* filename);

/**
 * Output runtime dynamics to the specified file.
 */
void BORG_Algorithm_output_runtime(char* filename);

/**
 * Output runtime dynamics with the specified frequency.
 */
void BORG_Algorithm_output_frequency(int frequency);

/**
 * Output all evaluated solutions to the specified file.
 *
 * CAUTION: This feature can quickly generate large data files.
 */
void BORG_Algorithm_output_evaluations(char* filename);

/**
 * Output runtime dynamics in the format of Aerovis.
 */
void BORG_Algorithm_output_aerovis();

/**
 * Sets the initial population of the algorithm.  If set, MS Borg will
 * not perform any random initialization of the population.
 */
void BORG_Algorithm_initial_population(BORG_Population population);

/**
 * Returns the next solution to be evaluated.
 */
BORG_Solution BORG_Master_next(BORG_Algorithm algorithm);

/**
 * Sends a shutdown signal to the specified worker node.
 */
void BORG_Master_shutdown(int rank);

/**
 * Sends a solution for evaluation to the specified worker node.
 */
void BORG_Master_send(int rank, BORG_Solution solution);

/**
 * Receives an evaluated solution from any worker node, returning the rank of
 * the worker node which provided the solution.
 */
int BORG_Master_receive();

/**
 * Injects the specified number of individuals, adding them onto the solution
 * queue to be evaluated and added into the population.
 */
void BORG_Master_injection(BORG_Algorithm algorithm, int size);

/**
 * Performs a restart on the MS Borg algorithm; the only difference is any new,
 * randomly-generated solutions are added to the queue.
 */
void BORG_Master_restart(BORG_Algorithm algorithm);

/**
 * Initializes the MS Borg algorithm, sending the initial population for evaluation.
 */
void BORG_Master_initialize(BORG_Algorithm algorithm);

/**
 * Runs the initialized MS Borg algorithm until all allocated NFE are consumed.
 */
void BORG_Master_loop(BORG_Algorithm algorithm);

/**
 * Main thread for the master node, returning the end-of-run approximation set.
 */
BORG_Archive BORG_Master_run();

/**
 * Main thread for the worker node, evaluating solutions until a shutdown
 * signal is received.
 */
void BORG_Worker_run();

/**
 * Starts the MS Borg algorithm, assigning each node a role (master or worker).  
 */
void BORG_Algorithm_ms_startup(int* argc, char*** argv);

/**
 * Allocates the specified number of objective function evaluations when running the MS Borg
 * algorithm.  The algorithm will terminate automatically if the NFE exceeds this threshold.
 */
void BORG_Algorithm_ms_max_evaluations(int maxEvaluations);

/**
 * Allocates the specified number of hours when running the MS Borg algorithm.  The algorithm
 * will terminate automatically if the wallclock time exceeds this threshold.
 */
void BORG_Algorithm_ms_max_time(double maxHours);

/**
 * Sets the initialization strategy.
 */
void BORG_Algorithm_ms_initialization(BORG_Initialization initialization);

/**
 * Runs the MS Borg algorithm on the specified problem given the allocated
 * resources.  This method must be invoked on all nodes, but only the master node will
 * return the result; all other nodes return with NULL.
 */
BORG_Archive BORG_Algorithm_ms_run(BORG_Problem problem);

/**
 * Shuts down the MS Borg algorithm.  No additional MS_Algorithm_* methods should be invoked.
 */
void BORG_Algorithm_ms_shutdown();

#ifdef __cplusplus
}
#endif

#endif
