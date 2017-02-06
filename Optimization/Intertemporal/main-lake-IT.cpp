/* main-lake_DPS.cpp
   
  Riddhi Singh, May, 2014
  The Pennsylvania State University
  rus197@psu.edu

  Adapted by Tori Ward, July 2014 
  Cornell University
  vlw27@cornell.edu

  Adapted by Jonathan Herman and David Hadka, Sept-Dec 2014
  Cornell University and The Pennsylvania State University

  Adapted by Julianne Quinn, July 2015
  Cornell University
  jdq8@cornell.edu

  A multi-objective represention of the lake model from Carpenter et al., 1999
  This simulation is designed for optimization with Multi-Master Borg.

  Stochasticity is introduced by natural phosphorous inflows. 
  These follow a lognormal distribution with specified mu and sigma.

  Decision variable
    vars : vector of 100 years of P emissions

  Objectives
  1: max avg P concentration
  2: mean economic benefits
  3: mean inertia
  4: mean reliability

  Constraints
  Reliability must be > 85%

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sstream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/tools/roots.hpp>
#include "./../borg/moeaframework.h"
#include "./../boostutil.h"
#include "./../borg/borgms.h"

#define nYears 100
#define nSamples 100
#define inertia_threshold -0.02
#define reliability_threshold 0.85

double b, q, alpha, delta, pCrit;

int nvars = nYears;
int nobjs = 4;
int nconstrs = 1;
double nat_flowmat [10000][nYears]; // create a matrix of [ 10000 x nYears ]

namespace ublas = boost::numeric::ublas;
namespace tools = boost::math::tools;
using namespace std;

ublas::vector<double> average_annual_P(nYears);
ublas::vector<double> discounted_benefit(nSamples);
ublas::vector<double> yrs_inertia_met(nSamples);
ublas::vector<double> yrs_pCrit_met(nSamples);

ublas::vector<double> lake_state(nYears+1);

void lake_problem(double* vars, double* objs, double* constrs) 
{
  // initialize variables
  zero(average_annual_P);
  zero(discounted_benefit);
  zero(yrs_inertia_met);
  zero(yrs_pCrit_met);
  zero(lake_state);

  int linesToUse [nSamples];
  srand (time(NULL)); //gives PRNG a seed based on time
  for (int s=0; s < nSamples; s++) {
    //pick a random number based on time
    //choose 100 of 10,000 available inflow value lines
    linesToUse[s] = rand() % 10000;
  }

  // run lake model simulation
  for (int s = 0; s < nSamples; s++){
    // randomly generated natural phosphorous inflows
    double *nat_flow = new double [nYears];
    int index = linesToUse[s];
    // get the random natural flow from the States of the world file
    //each line of SOW file covers 100 days of inflow
    for (int i=0; i < nYears; i++){
      nat_flow[i] = nat_flowmat[index][i]; 
    }

    // initialize lake_state
    lake_state(0) = 0;
    // find initial policy-derived release

    //implement the lake model from Carpenter et al. 1999
    for (int i = 0; i < nYears; i++)
    {
      // new state: previous state - decay + recycling + pollution
      lake_state(i+1) = lake_state(i)*(1-b) + pow(lake_state(i),q)/(1+pow(lake_state(i),q)) + vars[i] + nat_flow[i];
      average_annual_P(i) += lake_state(i+1)/nSamples;
      discounted_benefit(s) += alpha*vars[i]*pow(delta,i);

      if (i>=1 && (vars[i] - vars[i-1]) > inertia_threshold)
        yrs_inertia_met(s) += 1;

      if(lake_state(i+1) < pCrit)
        yrs_pCrit_met(s) += 1;
    }
  }
  
  objs[0] = -1*vsum(discounted_benefit)/nSamples; // average economic benefit
  objs[1] = vmax(average_annual_P);; // max average annual P concentration
  objs[2] = -1*vsum(yrs_inertia_met)/((nYears-1)*nSamples); // average inertia
  objs[3] = -1*vsum(yrs_pCrit_met)/(nYears*nSamples); // average reliability

  constrs[0] = max(0.0, reliability_threshold - (-1*objs[3]));
  
  average_annual_P.clear();
  discounted_benefit.clear();
  yrs_inertia_met.clear();
  yrs_pCrit_met.clear();
  lake_state.clear();

}

double root_function(double x) {
  return pow(x,q)/(1+pow(x,q)) - b*x;
}

bool root_termination(double min, double max) {
  return abs(max - min) <= 0.000001;
}

int main(int argc, char* argv[]) 
{  
  // initialize defaults
  b = 0.42;
  q = 2;
  alpha = 0.4;
  delta = 0.98;

  std::pair<double, double> root = tools::bisect(root_function, 0.01, 1.0, root_termination);
  pCrit = (root.first + root.second) / 2;

  for (int i=0;i<10000;i++){   //this is 10,000 to match nat_flowmat's size
    for (int j=0;j<nYears;j++){
      nat_flowmat[i][j] = 0.0; 
    }
  }
  
  FILE * myfile;
  myfile = fopen("./../SOWs_Type6.txt","r");
  
  int linenum = 0;
  int maxSize = 5000;
  
  if (myfile==NULL){
    perror("Error opening file");
  } else {
    char buffer [maxSize];
    while (fgets(buffer, maxSize, myfile)!=NULL){
      linenum++;
      if (buffer[0]!='#'){
        char *pEnd;
        char *testbuffer = new char [maxSize];
        for (int i=0; i <maxSize; i++){
          testbuffer[i] = buffer[i];
        }
  
        for (int cols=0; cols < nYears; cols++){
          nat_flowmat[linenum-1][cols] = strtod(testbuffer, &pEnd);
          testbuffer  = pEnd; 
        }       
      }
    }
  }
  
  fclose(myfile);

  // setting random seed
  unsigned int seed = atoi(argv[1]);
  srand(seed);
  int NFE = atoi(argv[2]);

  // interface with Borg-MS
  BORG_Algorithm_ms_startup(&argc, &argv);
  BORG_Algorithm_ms_max_evaluations(NFE);
  BORG_Algorithm_output_frequency(NFE/200);

  // Define the problem with decisions, objectives, constraints and the evaluation function
  BORG_Problem problem = BORG_Problem_create(nvars, nobjs, nconstrs, lake_problem);

  // Set all the parameter bounds and epsilons
  for (int j=0; j < nvars; j++){
    BORG_Problem_set_bounds(problem, j, 0.01, 0.1);
  }

  BORG_Problem_set_epsilon(problem, 0, 0.01); // average economic benefit
  BORG_Problem_set_epsilon(problem, 1, 0.01); // max average annual P concentration
  BORG_Problem_set_epsilon(problem, 2, 0.0001); // average inertia
  BORG_Problem_set_epsilon(problem, 3, 0.0001); // average reliability

  //This is set up to run only one seed at a time
  char outputFilename[256];
  char runtime[256];
  FILE* outputFile = NULL;
  sprintf(outputFilename, "./sets/LakeIT_S%d.set", seed);
  sprintf(runtime, "./runtime/LakeIT_S%d.runtime", seed);

  BORG_Algorithm_output_runtime(runtime);

  BORG_Random_seed(seed);
  BORG_Archive result = BORG_Algorithm_ms_run(problem); // this actually runs the optimization

  //If this is the master node, print out the final archive
  if (result != NULL){
    outputFile = fopen(outputFilename, "w");
    if(!outputFile){
      BORG_Debug("Unable to open final output file\n");
    }
    BORG_Archive_print(result, outputFile);
    BORG_Archive_destroy(result);
    fclose(outputFile);
  }

  BORG_Algorithm_ms_shutdown();
  BORG_Problem_destroy(problem);

  return EXIT_SUCCESS;

}
