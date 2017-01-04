/***************************************************************************
 *   copyright (C) 2013 by Marco Caserta                                   *
 *   marco dot caserta at ie dot edu                                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/**! \file chart.cpp
   \brief Algorithm for the Seating Chart Assignment Problem

   Author: Marco Caserta (marco dot caserta at ie dot edu)\n
   Started: 13.08.13\n
   Ended:   26.08.13\n 

   Compile with: make
   Execute with: ./bin/chart -f "input_data_file" 
   Some useful options are:
   # -w nr.workgroups : This option defines the value of variable "nWG". It defines the number 
                        of workgroups assignment already known. Depending on when the software
			is used, the workgroup assigments for terms 1, 2, or 3 could be known. 
			If this is the case, one of the criteria of the algorithm is to avoid 
			having students from the same term seating next to each other. 

   # -l neighborhood type : This option defines the value of variable "largeNeighborhood". It
                        define the type of neighborhood used for the most important 
                        constraints. Two types of neighborhoods are currently defined, i.e.:
			SMALL ->      o------X------o

			LARGE ->             o
			                     |
				      o------X------o
				             |
					     o

			The diffence lies in the size of the neighborhood itself. Of course, the
			larger the neighborhood, the harder the problem (and, consequently,
			we experience an increase in running time). 
   # -t cplex time    : This option defines the value of variable IloMaxTime. It establishes the 
                        maximum running time given to cplex to solve each phase. Phase 1 (the
                        cluster assignment problem) is a network problem and takes no time. Phase 2,
			however, is a much harder problem and often requires around 500 seconds for
			the typical IMBA instance. A way to reduce the running time is to use a 
			SMALL neighborhood.
			
   # -p nPeriodsGiven : This option allows to read previouosly defined solutions. For example,
                        imagine that, for a given section, the seating chart of the first term
			has already been defined. We need to set "-p 1", to ensure that the previous
			solution is read and that the corresponding variables in the models of phase 1
			and phase 2 are fixed to the appropriate values (0 or 1). The previous
			solution(s) should be added as last column(s) in the data file.

   \param A text file selected via command line with the flag "-f input_data_file". The format
          of the file is specified below and must be strictly respected.
	  
   \param A list of disk files, required to characterize and properly define the layout of the
          room along with other parameters:
       a) _ADJACENCY = "data/networkLarge.csv". The adjacency matrix for the LARGE version of 
          the neighborhood. Note that the matrix is symmetric.
       b) _ADJACENCYREDUCED = "data/networkReduced.csv". The adjacency matrix for the SMALL
          version of the neighborhood. Note that the matrix is symmetric.
       c) _FORBIDDEN = "data/forbidden.csv". A binary matrix that defines whether students "i"
          and "j" can be seated in the same neighborhood (where the definition of the neighborhood
	  depends on parameters "largeNeighborhood"). Note that this is a hard constraint, i.e.,
	  in the formulation, this set of constraints is never relaxed and, therefore, depending
	  on the structure of the matrix, it could lead to an infeasible problem.
       d) _CLUSTERS = "data/clusters.csv". For each cluster (currently E, B, F), it provides the
          list of seats included in the cluster itself. The seat number corresponds to the numbers
	  given in the file _MATRIX -- see point e) below.
       e) _MATRIX = "data/matrixChart.csv". The MAP that defines the layout of the room itself.
          It also provides a unique identifier for each seat. The seat numbering used in the code
	  follows the one provided in this file.
       f) _EMPTY = "data/emptySeats.csv". A list of spared seat. This list is used whenever the
          number of seats exceeds the number of students. In decreasing order of priorities, we
	  establish here which seats should be left empty.

	  
   \return Two disk files:
       i. _SOLFILE = "sol/solution.csv". A file containing the best solution found for the three
          terms. For each student (each row) the seats assigned to each term (1, 2, and 3) are 
	  provided. The seat numbering is again given by the file _MATRIX. 
	  (This solution can be manually copied into the initial Excel file.)
   
      ii. _OUTPUTFILE = "sol/output.csv". A summary of most relevant statistics that allow
          to get an idea of the solution quality. This is what should be used to generate a report 
	  for the user.

   ** OVERALL ALGORITHM

   The algorithm is divided in two parts:
   - stage 1 : Assignment of students to clustes (E, B, F). We defined two different cost functions
     that lead to different pattens of assignments. The type of solution chosen here affects the next
     step.
    
   - state 2 : Assignment of students to seats (within the cluster determined in stage 1). Three 
     types of constraints have been currently defined:

     i. Hard Contraints : Constraints that must be necessarily satisfied. A violation of such 
     constraints leads to an infeasible solution. These constraints are:
        - Forbidden matrix: Some students cannot be neighbors (currently due to same lastname).

    ii. Important Constraints : Constraints whose violation is expensive, thus implying a high 
     priority. We use a binary variable to "relax" each of these constraints. We set a cost
     in the objective function and attemps to minimize the total cost of violation. Such constraints
     are:
        - Neighborhood : No repetition of students in the neighborhood over the three time periods
	- Nationality  : Students with same nationality should not be in the same neighborhood

   iii. Soft Constraints : Constaints whose violation is cheap. We do not add these constraints
     in the fist run of cplex. Thus, the first MIP model defined only includes constraints of type
     i. and ii. We then count then number of violations of soft constraints and, if any exists, 
     we add these soft constraints to the formulation and re-run cplex. The method used is the same
     of ii., i.e., we define a binary variable to "relax" each soft constraint and minimize the cost
     of violation. Such constraints are:
        - Gender : Women should not be in the same neighborhood
	- Workgroup : In each period, students belonging to the same workgroup should be in the same
	  neighborhood.
	  
     NOTE : The type of neighborhood used for constraints i. and ii. can be chosen using option -l
            from the command line (which, in turn, determines the value of variable 
	    "largeNeighborhood" - 0 or 1). However, constraints iii. only use the small neighborhood,
	    in an attempt to have, e.g., women not seating next to each other.

  ** MODIFICATIONS

  Add the option of reading previvously defined solutions. Via command line (option -p), we can define
  the number of seating assigments already defined. Thus, given that the total number of time periods
  is nPeriods, we only need to solve the problem for the remaining semesters. These solutions, if
  available, are passed via the same data file, as last column(s). For example, if for a given section
  we alredy have the seating assigment for the first term, we just need to add a new column to the 
  data file (at the end of the row for each student) that contains the seat number for that student.
  We also need to set "-p 1" to ensure that this new column is read.

  The strategy to deal with previous assigments is to fix the corresponding variables to 1, both in
  the Phase 1 and Phase 2 model. The other variables are thus fixed to zero.
*/
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <cstdlib>
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <stdexcept>

#include "options.h"
#include "timer.h"

#include <RInside.h>                    // for the embedded R via RInside

using namespace std;

// disk files needed to define data structure (see description above)
char* _FILENAME;                //!< name of the instance file
char* _FORBIDDEN;                //!< name of the instance file
const char * _ADJACENCY        = "data/networkLarge.csv";
const char * _ADJACENCYREDUCED = "data/networkReduced.csv";
//const char * _FORBIDDEN        = "data/forbidden.csv";
const char * _CLUSTERS         = "data/clusters.csv";
const char * _MATRIX           = "data/matrixChart.csv";
const char * _EMPTY            = "data/emptySeats.csv";
const char * _SOLFILE          = "sol/solution.csv";
const char * _OUTPUTFILE       = "sol/output.csv";

typedef IloArray < IloNumVarArray > TwoD;
typedef IloArray < TwoD > ThreeD;
typedef IloArray < ThreeD > FourD;

double INFTY = std::numeric_limits<double>::infinity();
const long _MAXRANDOM = 2147483647;
const double ZERO     = 0.0e0;
const double EPSI      = 0.00000001;

int nHash = 13;
string clusterLabels[] = {"E", "B", "F"};
string clusterNames[] = {"NA", "NA", "NA", "0-0-3", "0-1-2", "0-2-1", "1-0-2", "1-1-1", "1-2-0", "2-0-1", "2-1-0", "NA","3-0-0"};

// cplex parameters
int IloMaxTime;			//<! max time given to each run of cplex (option -t to change it)
const int IloMaxSol  = 150;
const int IloOutput  = 2;

// input file structure
const int genderPos  = 0;
const int natPos     = 9;

// general parameters
const int nClusters  = 3;	//!< E, B, F
const int nPeriods   = 3;	//!< thee terms of MBA program
const int nSeats     = 57;	//!< total number of seats in the room
int nStudents;
int nAttributes;

int nrowsChart;			//!< seating chart info : nr rows
int ncolsChart;			//!< seating chart info : nr cols

double time_limit;	//!< cpu-clock time limit
int wLS;		//!< with (1) or without (0) local search
int nWG;                //!< number of worgroups already defined (0,...,3)
int largeNeighborhood;	//!< type of neighborhood to be used by main constraints (1->large, 0->small)
int nPeriodsGiven;	//!< number of time periods already solved (some variables must be fixed)
timer tTime;            //!< Ojbect clock to measure REAL and VIRTUAL (cpu) time
int withForbidden;	//!< Defines whether a forbidden matrix should be read

int  ** y;			//!< seat assigned in period t to student i
int  ** s;			//!< student assigned to seat in each period
int  ** x;			//!< cluster assignment of student i in each period t
int  ** A;	                //!< adjacency matrix
int  ** F;			//!< incompatibilities matrix
int   * ub;	                //!< number of seats of each type (E, B, F)
int  ** clusters;               //!< list of seat of each type (E, B, F)
int  ** bin_clusters;		//!< binary indicating whether seat k belongs to cluster c
int  ** matrixChart;		//!< seating chart using a matrix structure
int  ** attr;			//!< attributes of students (read from disk file)
int   * freqHash;		//!< patterns frequencies (phase 1)
int  ** givenSol;		//!< solution (seat assignment) for previous semesters (if given)

std::vector < std::vector < int> > N_set;
std::vector < std::vector < int> > N_setReduced;
std::vector < std::vector < int> > nationalities;
std::vector < std::vector < std::vector < int > > > workgroups;
std::map < int, std::vector <int> > genderMap;

int infeasNeigh;  //!< count number of infeasibilities w.r.t. neighborhoods
int infeasNat;    //!< count number of infeasibilities w.r.t. nationalities
int infeasGender; //!< count number of infeasibilities w.r.t. gender
int infeasWG;	  //!< count number of infeasibilities w.r.t. workgroups
/************************ FUNCTIONS ******************************/
void read_problem_data();
void create_list_nationalities();
void create_list_workgroups();
void create_list_gender();
void define_network_problem(int nPeriods, int nStudents, int nClusters, int ** x, int * freqHash);
void compute_summary(IloModel model, IloCplex cplex, FourD x_ilo, int nPeriods, int nStudents, 
		     int nClusters, int ** x, int * freqHash);
int solve_cplex(IloModel model, IloCplex cplex, int solLimit, int displayLimit, int timeLimit);
int hash_pattern(int * stats, int nClusters);
void read_adjacency_matrix();
void read_clusters();
void create_neighborhoods(std::vector < std::vector < int > > & N_set);
void create_reduced_neighborhoods(std::vector < std::vector < int> > & N_setReduced);
void define_chart_model(int ** x);
void read_incompatibilities();
void get_chart(IloModel model, IloCplex cplex, ThreeD x_ilo);
void print_options();
void initialize_sol_data_structure();
void write_solution(int ** y);
void write_statistics(int ** s, int ** y, int * freqHash);
void get_cplex_solution(IloModel model, IloCplex cplex, ThreeD x_ilo, int ** y, int ** s);
void print_chart(int ** s);
void add_previous_charts(IloModel & model, FourD & x_ilo);
int get_student_cluster(int seat_number);
void add_previous_seats(IloModel & model, ThreeD & x_ilo);
/************************ FUNCTIONS ******************************/


/************************ MAIN PROGRAM ******************************/
/// Main program
/************************ MAIN PROGRAM ******************************/
int main(int argc, char *argv[])
{
   //freopen("debug.txt", "w", stdout); //!< redirect output to a file
   int err = parseOptions(argc, argv);
   if ( err != 0)
   { if (err != -1) cout << "Error argument " << err+1 << endl; exit(1); }

   string aa = "NA";
   if (aa.compare(_FORBIDDEN) == 0)
   {
      withForbidden = 0;
   }
   else
      withForbidden = 1;

   cout << "Forbidden is " << _FORBIDDEN << endl;
   IloMaxTime = time_limit;	// max time assigned to cplex (-t option of command line)
   
   read_problem_data();
   create_list_nationalities();
   create_list_workgroups();
   create_list_gender();
   initialize_sol_data_structure();

   read_adjacency_matrix();   
   read_clusters();
   if (withForbidden)
      read_incompatibilities();
   print_options();

   tTime.resetTime();		      // start clock

   // phase 1. assign students to each cluster (E, B, F) in each period
   define_network_problem(nPeriods, nStudents, nClusters, x, freqHash);

   // phase 2. assign students to each seat, respecting the cluster assignment
   create_neighborhoods(N_set);
   create_reduced_neighborhoods(N_setReduced);

   define_chart_model(x);

   cout << "Solution found in " << tTime.elapsedTime(timer::VIRTUAL) << " seconds." << endl;
   
   write_solution(y);
   write_statistics(s, y, freqHash);

   exit(0);
}


/************************ FUNCTIONS ******************************/
/// Functions
/************************ FUNCTIONS ******************************/
/// Read the file _FILENAME containing the instance to be solved
/** The format of the input file must be the following:
    - row 1 : number_of_students (nStudents)  number_of_attributes
    (nAttributes)
    - rows 2 .. nStudents+1 : 
    student_ID   'nAttributes' columns with the value of each
    attribute. The meaning of the column is defined as:
    - ...
    
    Note: The columns beyond 10, if present, are workgroup assignments.
*/
void read_problem_data()
{
   int temp;
   ifstream fdata(_FILENAME, ios::in);
   if (!fdata)
   { cerr << "Cannot open file " << _FILENAME << endl; exit(1); }

   fdata >> nStudents >> nAttributes;

   // read attributes for each student
   attr = new (int*[nStudents]);
   for (int i = 0; i < nStudents; i++)
      attr[i] = new int[nAttributes];

   // read previous periods solution
   if (nPeriodsGiven > 0)
   {
      givenSol = new (int*[nStudents]);
      for (int i = 0; i < nStudents; i++)
	 givenSol[i] = new int[nPeriodsGiven];
   }

   for (int i = 0; i < nStudents; i++)
   {
      fdata >> temp;		// read student number
      for (int j = 0; j < nAttributes; j++)
	 fdata >> attr[i][j];	// read attributes
      
      // if previous semestrers charts have been solved already, we read the solutions
      for (int t = 0; t < nPeriodsGiven; t++)
	 fdata >> givenSol[i][t];
      
#ifdef M_DEBUG
      cout << "s(" << i << ") :: ";
      for (int j = 0; j < nAttributes; j++)
	 cout << setw(4) << attr[i][j];
      cout << endl;
      int aaa;
      cin >> aaa;
#endif
   }
   fdata.close();
}

/// Print instance info and algorithmic parameters.
void print_options()
{
   cout << "-------------------------------------" << endl;
   cout << "- OPTIONS : " << endl;
   cout << "-------------------------------------" << endl;
   cout << " Chart MAP\t :: " << _MATRIX << endl;
   cout << " Clusters MAP\t :: " << _CLUSTERS << endl;
   cout << " Nr. Clusters\t :: " << nClusters << endl;
   cout << " Nr. Seats\t :: " << nSeats << endl;
   for (int c = 0; c < nClusters; c++)
      cout << "\t " << clusterLabels[c] << " : " << ub[c] << endl;
   cout << " Nr. Periods\t :: " << nPeriods << endl;
   cout << " Nr. Students\t :: " << nStudents << endl;
   cout << " Nr. Attributes\t :: " << nAttributes << endl;
   cout << " Nr. Workgroups\t :: " << nWG << endl;
   cout << " Nr. Charts Given\t :: " << nPeriodsGiven << endl;
   cout << endl;
   cout << "-------------------------------------" <<  endl;   
   cout << "-------------------------------------" <<  endl;   
   cout << "* Algorithmic Parameters :: " << endl;
   cout << "    Neighborhood Type\t :: ";
   if (largeNeighborhood == 1)
      cout << "Large" << endl;
   else
      cout << "Small" << endl;
   cout << "    Cplex Max Time\t :: " << IloMaxTime << endl;;
   cout << "    Incompatibilities\t :: ";
   if (withForbidden)
      cout << "Yes" << endl;
   else
      cout << "No" << endl;
   cout << "-------------------------------------" <<  endl << endl;   
}


/// For each student, create list of students with same nationality (i.e., incompatible)
void create_list_nationalities()
{
   for (int i = 0; i < nStudents-1; i++)
   {
      int nation = attr[i][natPos];
      std::vector <int> aux;
      for (int j = i+1; j < nStudents; j++)
	 if (nation == attr[j][natPos])
	    aux.push_back(j);
      nationalities.push_back(aux);      
      aux.clear();
   }

#ifdef M_DEBUG
   for (int i = 0; i < nStudents-1; i++)
   {
      cout << "Student " << i << " has same nationality of ";
      for (unsigned j = 0; j < nationalities[i].size(); j++)
	 cout << setw(4) << nationalities[i][j];
      cout << endl;
   }
   int aaa;
   cin >> aaa;
#endif
}

/// For each student, create list of students belonging to same workgroup
/** The list is three-dimensional, since we need to specify:
    - the number of time periods for which we have a workgroup information available (0, ..., 3)
    - the student "i"
    - the list of students sharing the workgroup with student i in period t
 */ 
void create_list_workgroups()
{
   int wgPos;
   std::vector < std::vector < int > > wgAux;
   for (int t = 0; t < nWG; t++)
   {
      wgPos = nAttributes - nWG + t;
      for (int i = 0; i < nStudents-1; i++)
      {
	 int wg = attr[i][wgPos];
	 std::vector <int> aux;
	 for (int j = i+1; j < nStudents; j++)
	    if (wg == attr[j][wgPos])
	       aux.push_back(j);

	 wgAux.push_back(aux);      
	 aux.clear();
      }
      workgroups.push_back(wgAux);
      wgAux.clear();
   }

#ifdef M_DEBUG
   for (unsigned t = 0; t < workgroups.size(); t++)
   {
      cout << "WG(" << t << ") ################################ " << endl;
      for (int i = 0; i < nStudents-1; i++)
      {
	 cout << "Student " << i << " is in same WG with ";
	 for (unsigned j = 0; j < workgroups[t][i].size(); j++)
	    cout << setw(4) << workgroups[t][i][j];
	 cout << endl;
	 int aaa;
	 cin >> aaa;
      }
   }
#endif
}

/// Create a list of forbidden pairs of women
void create_list_gender()
{

   for (int i = 0; i < nStudents-1; i++)
   {
      if (attr[i][genderPos] == 0) continue;

      std::vector <int> aux;      
      for (int j = i+1; j < nStudents; j++)
      {
	 if (attr[j][genderPos] == 0) continue;

	 aux.push_back(j);	// i and j are women
      }
      genderMap.insert(std::pair <int, std::vector <int> > (i, aux));
   }

#ifdef M_DEBUG   
   std::map <int, std::vector <int> >::iterator it;
   for (it = genderMap.begin(); it != genderMap.end(); it++)
   {
      cout << "Student " << it->first << " incompatible with :: ";
      for (unsigned j = 0; j < it->second.size(); j++)
	 cout << setw(4) << it->second[j];
      cout << endl;
   }
#endif
}

/// Solve network problem using cplex
int solve_cplex(IloModel model, IloCplex cplex, int solLimit, int displayLimit, int timeLimit)
{
   try
   {
      IloEnv env = model.getEnv();
      //cplex.setOut(env.getNullStream());
      cplex.setParam(IloCplex::MIPInterval,5000);
      cplex.setParam(IloCplex::MIPDisplay, displayLimit);
      //cplex.setParam(IloCplex::RINSHeur, 1);
      cplex.setParam(IloCplex::TiLim, timeLimit);

      // set number of solutions to be obtained before stopping
      cplex.setParam(IloCplex::IntSolLim, solLimit);

      // Optimize the problem and obtain solution.
      if ( !cplex.solve() ) 
      {
	 //env.error() << "Failed to optimize MIP" << endl;
	 throw(-1);
      }
      return 1;
   }
   catch (...) 
   {
      cout << "Exception caught :: Problem is infeasible " << endl;
      cout << ".";
      return -1;
   }
}

/// Phase 1 : Define network problem to allocate students to clusters
/** In this phase, we allocate each student to a cluster for each period. However,
    a specific seat is not assigned.

    Note: Different definitions of the costs of the objective function lead to different
    solutions. We could re-run the problem more than once to give the user the chance to
    choose among different solutions.
 */
void define_network_problem(int nPeriods, int nStudents, int nClusters, int ** x, int * freqHash)
{
   IloEnv env;
   IloModel model(env);

   FourD x_ilo(env, nPeriods);
   for (int t = 0; t < nPeriods; t++)
   {
      x_ilo[t] = ThreeD(env, nStudents);
      for (int i = 0; i < nStudents; i++)
      {
	 x_ilo[t][i] = TwoD(env, nClusters);
	 for (int k = 0; k < nClusters; k++)
	    x_ilo[t][i][k] = IloNumVarArray(env, nClusters, 0, 1, ILOINT);
      }
   }

   // max inflow capacity at each node (including node 1)
   for (int t = 0; t < nPeriods; t++)
      for (int k = 0; k < nClusters; k++)
      { 
	 IloExpr sum(env);
	 for (int l = 0; l < nClusters; l++)
	    for (int i = 0; i < nStudents; i++)
	       sum += x_ilo[t][i][l][k];
	 model.add(sum <= ub[k]);
      }

   // in each period, each student to only one level
   for (int t = 0; t < nPeriods; t++)
      for (int i = 0; i < nStudents; i++)
      {
	 IloExpr sum(env);
	 for (int k = 0; k < nClusters; k++)
	    for (int l = 0; l < nClusters; l++)
	       sum += x_ilo[t][i][k][l];
	 model.add(sum == 1.0);
      }

   // level 0 balance constraint
   for (int i = 0; i < nStudents; i++)
   {
      for (int k = 0; k < nClusters; k++)
      {
	 IloExpr sum(env);
	 sum += x_ilo[0][i][0][k];
	 for (int l = 0; l < nClusters; l++)
	    sum -= x_ilo[1][i][k][l];
	 model.add(sum == 0.0);
      }
   }

   // level 1 balance constraint
   for (int i = 0; i < nStudents; i++)
   {
      for (int k = 0; k < nClusters; k++)
      {
	 IloExpr lhs(env);
	 IloExpr rhs(env);
	 for (int l = 0; l < nClusters; l++)
	    lhs += x_ilo[1][i][l][k];
	 for (int l = 0; l < nClusters; l++)
	    rhs += x_ilo[2][i][k][l];
	 model.add(lhs == rhs);
      }
   }

   // all students out of node 0
   IloExpr sum(env);
   for (int i = 0; i < nStudents; i++)
      for (int k = 0; k < nClusters; k++)
	 sum += x_ilo[0][i][0][k];
   model.add(sum == nStudents);
   

   // cost definition
   int ** c = new (int*[nClusters]);
   for (int k = 0; k < nClusters; k++)
   {
      c[k] = new int[nClusters];
      for (int l = 0; l < nClusters; l++)
	 c[k][l] = 10000;
   }
   
   c[0][0] = 10000;
   c[0][1] = 10;
   c[0][2] = 1;

   c[1][0] = 100;
   c[1][1] = 10000;
   c[1][2] = 1;

   c[2][0] = 100;
   c[2][1] = 10;
   c[2][2] = 200;

//#ifdef OLD
   IloExpr totCost(env);
   for (int t = 0; t < nPeriods; t++)
      for (int k = 0; k < nClusters; k++)
	 for (int l = 0; l < nClusters; l++)
	    for (int i = 0; i < nStudents; i++)
	       totCost += c[k][l]*x_ilo[t][i][k][l];
//#endif

#ifdef OLD
   int * w = new int[nClusters];
   w[0] = 100;
   w[1] = 10;
   w[2] = 1;

   IloExpr totCost(env);
   for (int t = 0; t < nPeriods; t++)
      for (int k = 0; k < nClusters; k++)
	 for (int l = 0; l < nClusters; l++)
	    for (int i = 0; i < nStudents; i++)
	       totCost += w[l]*x_ilo[t][i][k][l];
#endif

   model.add(IloMinimize(env, totCost));

   // no more than once in lowest level
   for (int i = 0; i < nStudents; i++)
   {
      IloExpr sum(env);
      sum += x_ilo[0][i][0][0];
      for (int j = 0; j < nClusters; j++)
	 sum += x_ilo[1][i][j][0];
      for (int j = 0; j < nClusters; j++)
	 sum += x_ilo[2][i][j][0];

      model.add(sum <= 1.0);
   }

   if (nPeriodsGiven > 0)
      add_previous_charts(model, x_ilo);

   IloCplex cplex(model);
   solve_cplex(model, cplex, IloMaxSol, 0, IloMaxTime);

   // summary for each period
   compute_summary(model, cplex, x_ilo, nPeriods, nStudents, nClusters, x, freqHash);
}

void compute_summary(IloModel model, IloCplex cplex, FourD x_ilo, int nPeriods, int nStudents, 
		     int nClusters, int ** x, int * freqHash)
{
   IloEnv env = model.getEnv();

   // count how many time L, M, and H appears
   int ** statsStudent = new (int*[nStudents]);
   for (int i = 0; i < nStudents; i++)
   {
      statsStudent[i] = new int[nClusters];
      for (int k = 0; k < nClusters; k++)
	 statsStudent[i][k] = 0;	     
   }

   // read solution
   for (int i = 0; i < nStudents; i++)
      for (int t = 0; t < nPeriods; t++)
      {
	 for (int k = 0; k < nClusters; k++)
	    for (int l = 0; l < nClusters; l++)
	       if (cplex.getValue(x_ilo[t][i][k][l]) >= 1.0 - EPSI)
	       {
		  x[i][t] = l;
		  statsStudent[i][l]++;
	       }
      }

#ifdef M_DEBUG
   for (int i = 0; i < nStudents; i++)
   {
      cout << "s("<< setw(2) << i<<") = \t";
      for (int k = 0; k < nClusters; k++)
      {
	 for (int m = 0; m < statsStudent[i][k]; m++)
	    cout << clusterLabels[k];
	 cout << "\t";
      }
      cout << endl;
   }      
   cout << endl << endl;


   for (int i = 0; i < nStudents; i++)
      for (int t = 0; t < nPeriodsGiven; t++)
	 cout << "k(" << setw(2) << i << ") = " << x[i][t] << endl;
#endif

   cout << endl << "Number of students in each cluster" << endl;
   cout << "==================================" << endl;
   cout << endl;

   int * totInCluster = new int[nClusters];
   cout << setw(28) << "E" << setw(5) << "B" << setw(5) << "F" << endl;
   for (int t = 0; t < nPeriods; t++)
   {
      cout << "Period t = " << t+1 << " :: \t";
      for (int k = 0; k < nClusters; k++)
      {
	 totInCluster[k] = 0;
	 for (int i = 0; i < nStudents; i++)
	    for (int l = 0; l < nClusters; l++)
	       if (cplex.getValue(x_ilo[t][i][l][k]) >= 1.0 - EPSI)
		  totInCluster[k]++;
      }	 
      for (int k = 0; k < nClusters; k++)
	 cout << setw(5) << totInCluster[k];
      cout << endl;
   }
   
   // compute patterns frequency using hash function
   for (int i = 0; i < nHash; i++)
      freqHash[i] = 0;
   for (int i = 0; i < nStudents; i++)
   {
      int hash = hash_pattern(statsStudent[i], nClusters);
      freqHash[hash]++;
   }
   
   cout << endl << "Patterns Frequency" << endl;
   cout << "==================" << endl;
   cout << endl;
   int cum = 0;
   cout << setw(20) << "----------" << endl;
   cout << setw(15) << "E-B-F" << setw(5) << "#"  << endl;
   cout << setw(15) << "-----" << setw(5) << "--" << endl;
   for (int i = 0; i < nHash; i++)
      if (freqHash[i] > 0)
      {
	 cout << setw(15) << clusterNames[i] << setw(5) << freqHash[i] << endl;
	 cum += freqHash[i];
      }
   cout << setw(20) << "----------" << endl << endl;
   assert(cum == nStudents);
}

/// Compute hash value for each possible pattern (E-B-F)
int hash_pattern(int * stats, int nClusters)
{

   int hash = 0;
   for (int k = 0; k < nClusters; k++)
      hash += pow(2, nClusters-k-1)*stats[k];
   return hash;
}

/// Read the adjacency matrix that defines the main Neighborhood       
/** Note: The variable "largeNeighborhood" determines what type of neighborhood is used
    in the definition of the main constraints (forbidden pairs, and previous neighbors)

    - largeNeighborhood = 1 => for each seat, the neighborhood is defined as a crux, i.e,
      four seats (east, west, north, south) are considered neighbors

    - largeNeighborhood = 0 => for each seat, the neighborhood is defined as east, west only.
 */
void read_adjacency_matrix()
{
   int aux;
   ifstream fadjacency;

   if (largeNeighborhood)
      fadjacency.open(_ADJACENCY, ios::in);
   else
      fadjacency.open(_ADJACENCYREDUCED, ios::in);

   if (!fadjacency)
   {
      cout << "Could not open file " << _ADJACENCY << endl;
      exit(100);
   }

   fadjacency >> aux;
   assert(aux == nSeats);
   
   A = new (int*[nSeats]);
   for (int j = 0; j < nSeats; j++)
      A[j] = new int[nSeats];

   for (int j = 0; j < nSeats; j++)
      for (int k = 0; k < nSeats; k++)
	 fadjacency >> A[j][k];
   fadjacency.close();
}

/// A matrix of forbidden, or incompatible, students is read from a disk file
/** The matrix F is read from a disk file and is symmetric.
 */ 
void read_incompatibilities()
{
   int aux;
   ifstream fforbidden(_FORBIDDEN, ios::in);
   if (!fforbidden)
   { cout << "Could not open file " << _FORBIDDEN << endl; exit(300); }

   fforbidden >> aux;
   assert(aux == nStudents);

   F = new (int*[nStudents]);
   for (int i = 0; i < nStudents; i++)
   {
      F[i] = new int[nStudents];
      for (int j = 0; j < nStudents; j++)
	 fforbidden >> F[i][j];
   }
   fforbidden.close();
}

/// Read file containing the clusters of seats
/** The format of the file is the following:
    row 1 : list of seats of type E
    row 2 : list of seats of type B
    row 3 : list of seats of type F
*/ 
void read_clusters()
{
   ifstream fcluster(_CLUSTERS, ios::in);
   if (!fcluster)
   { cout << "Could not open file " << _CLUSTERS << endl; exit(400); }

   clusters = new (int*[nClusters]);
   ub       = new int[nClusters];

   for (int k = 0; k < nClusters; k++)
   {
      fcluster >> ub[k];
      clusters[k] = new int[ub[k]];
      for (int j = 0; j < ub[k]; j++)
	 fcluster >> clusters[k][j];
   }
   fcluster.close();

   bin_clusters  = new (int*[nClusters]);
   for (int c = 0; c < nClusters; c++)
   {
      bin_clusters[c] = new int[nSeats];
      for (int k = 0; k < nSeats; k++)
	 bin_clusters[c][k] = 0;
   }
   
   for (int c = 0; c < nClusters; c++)
      for (int k = 0; k < ub[c]; k++)
	 bin_clusters[c][clusters[c][k]] = 1;
}

/// Use adjacency matrix to define N(k), the list of neighbors of seat k
void create_neighborhoods(std::vector < std::vector < int> > & N_set)
{
   for (int i = 0; i < nSeats; i++)
   {
      std::vector <int> aux;
      for (int j = 0; j < nSeats; j++)
	 if (A[i][j] == 1) 
	    aux.push_back(j);
      N_set.push_back(aux);
   }
#ifdef M_DEBUG
   for (int i = 0; i < nSeats; i++)
   {
      cout << "N("<< setw(2) << i << ") \t :: ";
      for (int j = 0; j < N_set[i].size(); j++)
	 cout << setw(4) << N_set[i][j]+1;
      cout << endl;
   }
#endif
}

/// Use reduced adjacency matrix to define smaller neighborhoods (see introduction)
void create_reduced_neighborhoods(std::vector < std::vector < int> > & N_setReduced)
{
   int temp;
   ifstream fadjacencyReduced(_ADJACENCYREDUCED, ios::in);
   if (!fadjacencyReduced)
   { cout << "Could not open file " << _ADJACENCYREDUCED << endl; exit(200); }

   fadjacencyReduced >> temp;
   assert(temp == nSeats);

   for (int j = 0; j < nSeats; j++)
   {
      std::vector <int> aux;
      for (int k = 0; k < nSeats; k++)
      {
	 fadjacencyReduced >> temp;
	 if (temp == 1)
	    aux.push_back(k);
      }
      N_setReduced.push_back(aux);
   }

   fadjacencyReduced.close();

#ifdef M_DEBUG
   for (int i = 0; i < nSeats; i++)
   {
      cout << "NRed("<< setw(2) << i << ") \t :: ";
      for (int j = 0; j < N_setReduced[i].size(); j++)
	 cout << setw(4) << N_setReduced[i][j]+1;
      cout << endl;
   }
#endif
}

/// Phase 2: Define model to solve the seating chart problem
/** A note on the objective function:
    The most difficult constraint to satisfy is ensuring that students that were neighbors in any
    of the previous terms are no longer seating nearby (according to the definition of neighborhood).
    Since this set of constraints might lead to an infeasible model, we setup an objective function
    aimed at minimizing the infeasibility, via:
    1. creation of a (binary) gamma variable for each "no-repeating-neigh" constraint
    2. such variable takes value 1 to ensure satisfaction of a constraint, thus allowing two
    students to repeat neighborh
    3. the objective function of the model is now minimization of the sum of the gamma variables,
    i.e., the number of violated constraints.
*/
void define_chart_model(int ** x)
{
   IloEnv env;
   IloModel model(env);

   IloNumVarArray gamma(env);	// infeasibility variables (for neighborhood constraints)
   IloNumVarArray delta(env);	// infeasibility variables (for nationalities constraints)
   IloNumVarArray epsilon(env);	// infeasibility variables (for gender constraints)
   IloNumVarArray eta(env);	// infeasibility variables (for workgroup constraints)
   IloObjective zValue   = IloAdd(model, IloMinimize(env));

   ThreeD x_ilo(env, nPeriods);
   for (int t = 0; t < nPeriods; t++)
   {
      x_ilo[t] = TwoD(env, nStudents);
      for (int i = 0; i < nStudents; i++)
	 x_ilo[t][i] = IloNumVarArray(env, nSeats, 0, 1, ILOINT);
   }

   int seat, student_i, student_j, cluster_of_i, cluster_of_j, neigh_of_seat;
   int c_i_t1, c_j_t1, seat1, neigh_of_seat1;

   // set here empty seats (if number of seats exceeds number of students)
   int diff = nSeats - nStudents;
   int availableEmpty;
   int emptySeat;
   ifstream fempty(_EMPTY, ios::in);
   if (!fempty)
   { cout << "Could not open file " << _EMPTY << endl; exit(600); }

   fempty >> availableEmpty;
   assert(availableEmpty >= diff);
   for (int j = 0; j < diff; j++)
   {
      fempty >> emptySeat;
      for (int t = 0; t < nPeriods; t++)
	 for (int i = 0; i < nStudents; i++)
	    model.add(x_ilo[t][i][emptySeat] == 0.0); // keep these seats empty
   }
   fempty.close();

   // in each period, have each student seated in the right cluster (and only one place)
   for (int t = 0; t < nPeriods; t++)
      for (int i = 0; i < nStudents; i++)
      {
	 IloExpr sum(env);
	 int c = x[i][t];
	 for (int k = 0; k < nSeats; k++)
	    if (bin_clusters[c][k] == 1)
	       sum += x_ilo[t][i][k];
	    else
	       model.add(x_ilo[t][i][k] == ZERO); // this seat cannot be used
	 model.add(sum == 1.0);			  // only one seat can be taken
      }
   
   // in each period, each seat taken at most once (but only by students in the corresponding cluster)
   for (int t = 0; t < nPeriods; t++)
      for (int k = 0; k < nSeats; k++)
      {
	 IloExpr sum(env);
	 for (int i = 0; i < nStudents; i++)
	 {
	    cluster_of_i = x[i][t];
	    if (bin_clusters[cluster_of_i][k] == 1)
	       sum += x_ilo[t][i][k];
	 }
	 model.add(sum <= 1.0);
      }


   // incompatibilities among students (depends on the "forbidden" matrix)
   // Note: This constraint is hard, i.e., it cannot be violated (thus, it might lead to infeasible
   // problems).
   if (withForbidden)		// only if a forbidden matrix is defined
   {
      for (int t = 0; t < nPeriods; t++)
	 for (int i = 0; i < nStudents; i++)
	 {
	    int c = x[i][t];
	    for (int k = 0; k < ub[c]; k++)
	    {
	       seat = clusters[c][k]; // assume student i is assigned to seat "seat"
	       for (int j = i+1; j < nStudents; j++)
	       {
		  if (F[i][j] == 1) // assume student j is incompatible
		  {
		     cluster_of_j = x[j][t];
		     for (unsigned p = 0; p < N_set[seat].size(); p++)
		     {
			neigh_of_seat = N_set[seat][p];
			if (bin_clusters[cluster_of_j][neigh_of_seat] == 1)
			{
			   IloExpr sum(env);
			   sum += x_ilo[t][i][seat];
			   // student j could be assigned to "neigh_of_seat"
			   sum += x_ilo[t][j][neigh_of_seat];
			
			   model.add(sum <= 1.0);
			}
		     }
		  }
	       }
	    }
	 }
   }

   // not repeating neighbors from previous periods
   // we need to check and avoid repetitions between periods 0 and 1, 0 and 2, as well as 1 and 2
   int count_gamma     = 0;
   double weight_gamma = 0;	// weight of this constraint (value defined below)
   for (int t = 0; t < nPeriods-1; t++)
      for (int tp = t+1; tp < nPeriods; tp++)
	 for (int i = 0; i < nStudents; i++)
	 {
	    int c_i_t = x[i][t];	// cluster of student i in period t
	    for (int j = i + 1; j < nStudents; j++)
	    {
	       if (withForbidden && F[i][j] == 1) continue; // if they are incompatible, they cannot be neighbors 
	       
	       int c_j_t = x[j][t];	// cluster of student j
	       
	       for (int k = 0; k < ub[c_i_t]; k++)
	       {
		  int foundOne = 0;
		  IloExpr periodT(env);
		  seat = clusters[c_i_t][k]; // assume student i is assigned to "seat" in period t
		  periodT += x_ilo[t][i][seat];
		  for (unsigned p = 0; p < N_set[seat].size(); p++)
		  {
		     neigh_of_seat = N_set[seat][p];
		     if (bin_clusters[c_j_t][neigh_of_seat] == 1) // stud j could seat here
		     {
			periodT += x_ilo[t][j][neigh_of_seat];
			foundOne++;
		     }
		  }

		  if (foundOne == 0)
		     continue;

		  c_i_t1 = x[i][tp]; // cluster of student i in period t+1
		  c_j_t1 = x[j][tp]; // cluster of student j in period t+1
		  
		  foundOne = 0;
		  for (int w = 0; w < ub[c_i_t1]; w++)
		  {
		    
		     IloExpr periodT1(env);
		     seat1 = clusters[c_i_t1][w]; // assume i is assigned to w in period t+1
		     periodT1 += x_ilo[tp][i][seat1];
		     
		     for (unsigned p = 0; p < N_set[seat1].size(); p++)
		     {
			neigh_of_seat1 = N_set[seat1][p];
			if (bin_clusters[c_j_t1][neigh_of_seat1] == 1)
			{
			   periodT1 += x_ilo[tp][j][neigh_of_seat1];
			   foundOne++;
			}
		     }

		     if (foundOne > 0)
		     {
		      	IloRange rng(env, 0.0, periodT + periodT1, 3.0);
		      	model.add(rng);
		      	// add gamma variable to account for infeasibility
		      	// we give double weight to violations in two consecutive periods
		      	weight_gamma = 2.0*((double)nPeriods - (double)(tp - t)); 
		      	gamma.add(IloNumVar(zValue(weight_gamma) + rng(-1.0)));
		     }
		  }
	       }
	    }
	 }


   // nationalities constraints (students with same nationality cannot be seated in same neigh)
   double weight_delta = 0.1;
   for (int t = 0; t < nPeriods; t++)
   {
      for (int i = 0; i < nStudents-1; i++)
      {
	 int c = x[i][t];
	 for (int k = 0; k < ub[c]; k++)
	 {
	    seat = clusters[c][k]; // assume student i seats here
	    for (unsigned j = 0; j < nationalities[i].size(); j++)
	    {
	       student_j    = nationalities[i][j];
	       cluster_of_j = x[student_j][t];
	       for (unsigned p = 0; p < N_set[seat].size(); p++)
	       {

		  neigh_of_seat = N_set[seat][p];
		  if (bin_clusters[cluster_of_j][neigh_of_seat] == 1)
		  {
		     IloExpr sum(env);
		     sum += x_ilo[t][i][seat];

		     // student j could be assigned to "neigh_of_seat"
		     sum += x_ilo[t][student_j][neigh_of_seat];

		     IloRange rng(env, 0.0, sum, 1.0);
		     model.add(rng);
		     // add delta variables to account for nationality infeasibility
		     delta.add(IloNumVar(zValue(weight_delta) + rng(-1.0)));
		  }
	       }
	    }
	 }
      }
   }

   if (nPeriodsGiven > 0)
      add_previous_seats(model, x_ilo);

   IloCplex cplex(model);


   // run cplex without gender and workgroup constraints
   solve_cplex(model, cplex, IloMaxSol, IloOutput, IloMaxTime);
   get_cplex_solution(model, cplex, x_ilo, y, s);
   //print_chart(s);
   get_chart(model, cplex, x_ilo);
   
   // infeasGender = 1;
   // infeasWG     = 1;
   if (infeasGender > 0)     // ADD gender constraint and re-run cplex
   {
      double weight_epsilon = 0.05; // weight in the objective function
      std::map <int, std::vector <int> >::iterator it;
      for (int t = 0; t < nPeriods; t++)
      {
	 for (it = genderMap.begin(); it != genderMap.end(); it++)
	 {
	    student_i = it->first;
	    int c = x[student_i][t];
	    for (int k = 0; k < ub[c]; k++)
	    {
	       seat = clusters[c][k]; // assume student i seats here

	       for (unsigned j = 0; j < it->second.size(); j++)
	       {
		  student_j    = it->second[j];
		  cluster_of_j = x[student_j][t];

		  for (unsigned p = 0; p < N_setReduced[seat].size(); p++)
		  {
		     neigh_of_seat = N_setReduced[seat][p];
		     if (bin_clusters[cluster_of_j][neigh_of_seat] == 1)
		     {
			IloExpr sum(env);
			sum += x_ilo[t][student_i][seat];

			// student j could be assigned to "neigh_of_seat"
			sum += x_ilo[t][student_j][neigh_of_seat];

			IloRange rng(env, 0.0, sum, 1.0);
			model.add(rng);
			// add epsilon variables to account for gender infeasibility
			epsilon.add(IloNumVar(zValue(weight_epsilon) + rng(-1.0)));
		     }
		  }
	       }
	    }
	 }
      }
   }

   if (infeasWG > 0)
   {
      double weight_eta = 0.01;	// weight in the objective function
      for (int t = 0; t < workgroups.size(); t++)
      {
	 for (int i = 0; i < nStudents-1; i++)
	 {
	    int c = x[i][t];
	    for (int k = 0; k < ub[c]; k++)
	    {
	       seat = clusters[c][k]; // assume student i seats here

	       for (unsigned j = 0; j < workgroups[t][i].size(); j++)
	       {
		  student_j    = workgroups[t][i][j];
		  cluster_of_j = x[student_j][t];
		  for (unsigned p = 0; p < N_setReduced[seat].size(); p++)
		  {
		     neigh_of_seat = N_setReduced[seat][p];
		     if (bin_clusters[cluster_of_j][neigh_of_seat] == 1)
		     {
			IloExpr sum(env);
			sum += x_ilo[t][i][seat];
			
			// student j could be assigned to "neigh_of_seat"
			sum += x_ilo[t][student_j][neigh_of_seat];

			IloRange rng(env, 0.0, sum, 1.0);
			model.add(rng);
			// add eta variables to account for workgroups infeasibility
			eta.add(IloNumVar(zValue(weight_eta) + rng(-1.0)));
		     }
		  }
	       }
	    }
	 }
      }
   }

   solve_cplex(model, cplex, IloMaxSol, IloOutput, 2*IloMaxTime);
   get_chart(model, cplex, x_ilo);
}

/// Create output of the final solution and print some basic statistics about infeasibility
void get_chart(IloModel model, IloCplex cplex, ThreeD x_ilo)
{

   get_cplex_solution(model, cplex, x_ilo, y, s);
   print_chart(s);

   infeasNeigh  = 0;
   infeasNat    = 0;
   infeasGender = 0;
   infeasWG     = 0;
   
   // verify that there is no repetition of neighborhs
   int ** counting = new (int*[nStudents]);
   for (int i = 0; i < nStudents; i++)
   {
      counting[i] = new int[nStudents];
      for (int j = 0; j < nStudents; j++)
	 counting[i][j] = 0;
   }

   for (int t = 0; t < nPeriods; t++)
      for (int i = 0; i < nStudents; i++)
      {
	 int seat_i = y[t][i];
	 for (unsigned k = 0; k < N_set[seat_i].size(); k++)
	 {
	    int seat_neigh = N_set[seat_i][k];
	    counting[i][s[t][seat_neigh]]++;
	 }
      }
   
   // count number of infeasibilities w.r.t. neighborhoods
   for (int i = 0; i < nStudents-1; i++)
      for (int j = i+1; j < nStudents; j++)
	 if (counting[i][j] > 1)
	 {
	    infeasNeigh++;
	    cout << "WARNING :: Students " << i << " and " << j << " repeat neighborhood " << endl;
	 }

   // count number of infeasibilities due to nationalities
   for (int i = 0; i < nStudents-1; i++)
      for (unsigned j = 0; j < nationalities[i].size(); j++)
      {
	 int student_j = nationalities[i][j];
	 if (counting[i][student_j] > 0)
	 {
	    infeasNat++;
	    cout << "WARNING :: Students " << i << " and " << student_j << " are neighbors and have same nationality (" << attr[i][natPos] << ") " << endl; 
	    assert(attr[i][natPos] == attr[student_j][natPos]);
	 }
      }

   int ** countingReduced = new (int*[nStudents]);
   for (int i = 0; i < nStudents; i++)
   {
      countingReduced[i] = new int[nStudents];
      for (int j = 0; j < nStudents; j++)
	 countingReduced[i][j] = 0;
   }

   for (int t = 0; t < nPeriods; t++)
      for (int i = 0; i < nStudents; i++)
      {
	 int seat_i = y[t][i];
	 for (unsigned k = 0; k < N_setReduced[seat_i].size(); k++)
	 {
	    int seat_neigh = N_setReduced[seat_i][k];
	    countingReduced[i][s[t][seat_neigh]]++;
	 }
      }


   // count number of infeasibilities due to gender
   std::map <int, std::vector <int> >::iterator it;
   for (it = genderMap.begin(); it != genderMap.end(); it++)
   {
      int student_i = it->first;
      for (unsigned j = 0; j < it->second.size(); j++)
      {
	 int student_j = it->second[j];
	 if (countingReduced[student_i][student_j] > 0)
	 {
	    infeasGender++;
	    cout << "WARNING :: Students " << student_i << " and " << student_j << " are neighbors and women (" << attr[student_i][genderPos] << ") " << endl; 
	    //assert(attr[student_i][genderPos] == attr[student_j][genderPos]);
	 }
      }
   }


   // count number of infeasibilities due to workgroups
   for (int t = 0; t < workgroups.size(); t++)
   {
      int wgPos = nAttributes - nWG + t;
      for (int i = 0; i < nStudents-1; i++)
      {
	 int seat_i = y[t][i];
	 for (unsigned j = 0; j < workgroups[t][i].size(); j++)
	 {
	    int student_j = workgroups[t][i][j];
	    int seat_j = y[t][student_j];
	    for (unsigned p = 0; p < N_setReduced[seat_i].size(); p++)
	    {
	       int neigh_of_seat = N_setReduced[seat_i][p];
	       if (seat_j == neigh_of_seat)
	       {
		  infeasWG++;
		  cout << "WARNING :: Students " << i << " and " << student_j << " are neighbors and same WG (" << attr[i][wgPos] << ") in period " << t << endl; 
		  assert(attr[i][wgPos] == attr[student_j][wgPos]);
		  break;
	       }
	    }
	 }
      }
   }
   cout << " ====================" << endl;
   cout << "| Summary Statistics |" << endl;
   cout << " ====================" << endl;
   cout << " *  Nr. Neighbors Infeasibilities \t :: " << infeasNeigh << endl;
   cout << " *  Nr. Nationality Infeasibilities \t :: " << infeasNat << endl;
   cout << " *  Nr. Gender Infeasibilities \t\t :: " << infeasGender << endl;
   cout << " *  Nr. WG Infeasibilities \t\t :: " << infeasWG << endl;
}

/// Write summary statistics to a disk file
void write_statistics(int ** s, int ** y, int * freqHash)
{
   int val;

   ofstream fout(_OUTPUTFILE, ios::out);

   fout << setw(30) << "SUMMARY STATISTICS " << endl << endl << endl;

   fout << "==================" << endl;
   fout << endl << "Patterns Frequency" << endl;
   fout << "==================" << endl;
   fout << endl;
   fout << setw(20) << "----------" << endl;
   fout << setw(15) << "E-B-F" << setw(5) << "#"  << endl;
   fout << setw(15) << "-----" << setw(5) << "--" << endl;

   for (int i = 0; i < nHash; i++)
      if (freqHash[i] > 0)
	 fout << setw(15) << clusterNames[i] << setw(5) << freqHash[i] << endl;
   fout << setw(20) << "----------" << endl << endl;
   
   
   fout << endl << endl << "=========================" << endl;
   fout << "** ** Seating chart ** **" << endl;
   fout << "=========================" << endl;
   fout << endl;
   for (int t = 0; t < nPeriods; t++)
   {
      fout << endl << endl << " -------------------------------PERIOD " << t+1 <<"-------------------------------" << endl;
      for (int i = 0; i < nrowsChart; i++)
      {
	 for (int j = 0; j < ncolsChart; j++)
	 {
	    val = matrixChart[i][j]; 
	    if (val == -1)
	       fout << setw(4) << " ";
	    else
	       if (val == -2)
		  fout << setw(3) << "  | ";
	       else
		  if (s[t][val] == -1)
		     fout << setw(4) << " ";
		  else
		     fout << setw(4) << s[t][val];
	 }
	 fout << endl;
      }
   }

   fout << endl << endl;
   fout << " ====================" << endl;
   fout << "| Summary Statistics |" << endl;
   fout << " ====================" << endl;
   fout << " *  Nr. Neighbors Infeasibilities \t :: " << infeasNeigh << endl;
   fout << " *  Nr. Nationality Infeasibilities \t :: " << infeasNat << endl;
   fout << " *  Nr. Gender Infeasibilities \t\t :: " << infeasGender << endl;
   fout << " *  Nr. WG Infeasibilities \t\t :: " << infeasWG << endl;

   fout.close();
}

/// Data structure initialization, for phase 1 and 2 of the problem
void initialize_sol_data_structure()
{
   // hash function to compute frequencies of cluster patterns
   freqHash = new int[nHash];

   // cluster assigned to each student in each period (phase 1 of the problem)
   x = new (int*[nStudents]);
   for (int i = 0; i < nStudents; i++)
      x[i] = new int[nPeriods];

   // assignment of each seat in each period (phase 2 of the problem)
   s = new (int*[nPeriods]);
   for (int t = 0; t < nPeriods; t++)
   {
      s[t] = new int[nSeats];
      for (int k = 0; k < nSeats; k++)
	 s[t][k] = -1;
   }

   // position of each student in each period (phase 2 of the problem)
   y = new (int*[nPeriods]);
   for (int t = 0; t < nPeriods; t++)
      y[t] = new int[nStudents];


   // read matrix chart (the physical MAP of the classroom)
   ifstream fmatrix(_MATRIX, ios::in);
   if (!fmatrix)
   { cout << "Could not open file " << _MATRIX << endl; exit(500); }

   fmatrix >> nrowsChart >> ncolsChart;
   matrixChart = new (int*[nrowsChart]);
   for (int i = 0; i < nrowsChart; i++)
   {
      matrixChart[i] = new int[ncolsChart];
      for (int j = 0; j < ncolsChart; j++)
	 fmatrix >> matrixChart[i][j];
   }
   fmatrix.close();
}

/// Write final solution to disk
void write_solution(int ** y)
{
   ofstream fsol(_SOLFILE, ios::out);
   for (int i = 0; i < nStudents; i++)
   {
      for (int t = 0; t < nPeriods; t++)
	 fsol << y[t][i] << "\t";
      fsol << endl;
   }
   fsol.close();
}

/// Get cplex solution
void get_cplex_solution(IloModel model, IloCplex cplex, ThreeD x_ilo, int ** y, int ** s)
{
   IloEnv env = model.getEnv();

   // assignment of each seat in each period
   for (int t = 0; t < nPeriods; t++)
      for (int k = 0; k < nSeats; k++)
	 s[t][k] = -1;

   // position of each student in each period
   for (int t = 0; t < nPeriods; t++)
      for (int i = 0; i < nStudents; i++)
	 for (int k = 0; k < nSeats; k++)
	    if (cplex.getValue(x_ilo[t][i][k] >= 1.0 - EPSI))
	    {
	       y[t][i] = k;
	       s[t][k] = i;
	    }
}

void print_chart(int ** s)
{
   cout << endl << endl;
   cout << "=========================" << endl;
   cout << "** ** Seating chart ** **" << endl;
   cout << "=========================" << endl;
   cout << endl;
   for (int t = 0; t < nPeriods; t++)
   {
      cout << endl << endl << " ------------------------------PERIOD " << t+1 <<"-------------------------------" << endl;
      for (int i = 0; i < nrowsChart; i++)
      {
	 for (int j = 0; j < ncolsChart; j++)
	 {
	    int val = matrixChart[i][j]; 
	    if (val == -1)
	       cout << setw(4) << " ";
	    else
	       if (val == -2)
		  cout << setw(3) << "  | ";
	       else
		  if (s[t][val] == -1)
		     cout << setw(4) << " ";
		  else
		     cout << setw(4) << s[t][val];
	 }
	 cout << endl << endl;
      }
   }
}

/// Add solution (seat assigments) for previous semesters
void add_previous_charts(IloModel & model, FourD & x_ilo)
{
   IloEnv env = model.getEnv();

   for (int i = 0; i < nStudents; i++)
   {
      int previous_k_i = 0;
      for (int t = 0; t < nPeriodsGiven; t++)
      {
	 int k_i = get_student_cluster(givenSol[i][t]);

	 model.add(x_ilo[t][i][previous_k_i][k_i] == 1.0); // fix cluster

	 for (int k = 0; k < nClusters; k++)
	 {
	    if (k == k_i) continue;
	    model.add(x_ilo[t][i][previous_k_i][k] == 0.0); // fix all the others to zero
	 }
	 
	 previous_k_i = k_i;
      }
   }
}

/// Given the seat number, identify the cluster type
int get_student_cluster(int seat_number)
{
   for (int k = 0; k < nClusters; k++)
      if (bin_clusters[k][seat_number] == 1)
	 return k;
}

/// Add seats defined in previous semesters to Phase 2 model
void add_previous_seats(IloModel & model, ThreeD & x_ilo)
{
   IloEnv env = model.getEnv();

   for (int i = 0; i < nStudents; i++)
   {
      for (int t = 0; t < nPeriodsGiven; t++)
      {
	 int seat_i = givenSol[i][t];
	 model.add(x_ilo[t][i][seat_i] == 1.0);

	 for (int k = 0; k < nSeats; k++)
	 {
	    if (k == seat_i) continue;
	    model.add(x_ilo[t][i][k] == 0.0);
	 }
      }
   }
}
