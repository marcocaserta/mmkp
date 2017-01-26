/***************************************************************************
 *   copyright (C) 2015 by Marco Caserta                                   *
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
  \brief Algorithm for the Multiple-choice Multidimensional Knapsack Problem

Author: Marco Caserta (marco dot caserta at ie dot edu)\n
Started: 05.02.15\n
Ended:   19.03.15\n 

Compile with: make
Execute with: ./bin/mmkp -f "input_data_file" 

To have a look at the flags that can be defined via command line, type:
> ./bin/mmkp -h

The options for the command line are defined in the file "options.cpp".

PROJECT STRUCTURE
=================

The project is composed of the following files:
- mmkp.cpp: The main file, which reads the instance data, calls the 
lagrangean phase, and implements the refinement procedure.
- lagrange.cpp: The implementation of the lagrangean relaxation with
subgradient optimization. 
- cplex.cpp : It manages the call to cplex (model creation+call to solver)
- options.cpp : Command-line parameters management.
*/

//#define M_DEBUG  // activate this to see the output for debugging

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cassert>


#include "timer.h"
#include "options.h"

#include "SampleDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"

#include "lagrange.h"
#include "cplex.h"

using namespace std;
char* _FILENAME;          // !< name of the instance file
double time_limit;        // !< wall-clock time limit
double corridorWidthBase; // base value for corridor width
int nSolBase;             // nr of sols to be found by cplex within corridor
double propFixed0;        // % of vars to be fixed to zero
double propFixed1;        // % of vars to be fixed to one
int add_oldSol_cut;       // boolean: cut out old feasible solution
int max_iter;             // max number lagrangean iterations
double Omega;             // parameter of robust formulation
double sigmaSq;           // parameter of individual item

ofstream flagr("lagrange.txt", ios::out);
ofstream fsol("solution.txt", ios::out);

/// Implementation of the basic data structure to hold the instance data
struct INSTANCE {
		int nR;       // number of resources (constraints)
		int nC;       // number of classes
		int * ri;     // number of items in each class i

		double  ** c; // obj function coefficients
		int *** w;    // knapsack constraints coefficients
		int   * R;    // r.h.s. values of the knapsack constraints
};

INSTANCE inp;         // INSTANCE DATA

double * sigma2;     // sigma^2 for robust formulation
int  * F;             // F[i] = 3 ==> var 3 in class "i" is set to 1
int **F0;             // F0[i][j] = 1 ==> var x_{ij} is set to 0

double best_time;     // time to best solution
int bestIter;         // lagr iteration in which best was found
double zBest;         // best feasible solution (lower bound)
double ubStar;        // best upper bound (lagr solution)


// definition of constant values
double INFTY = std::numeric_limits<double>::infinity();
const long _MAXRANDOM  = 2147483647;
const double ZERO      = 0.0e0;
const double EPSI      = 0.00001;


timer tTime;            //!< Ojbect clock to measure REAL and VIRTUAL (cpu) time


/************************ FUNCTIONS ******************************/
void read_problem_data(INSTANCE & inp);
void print_options(INSTANCE inp);
void lagrangean_phase(INSTANCE & inp, IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
				IloObjective & obj, IloNumVar & Q_ilo, double Omega, double sigma);
double get_first_lb(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj, 
				int nSol, IloNumVar & Q_ilo, double Omega, double sigma);
void update_best(int * xLBest, int * xL, double & ubStar, double zL, double * lambda, double * lambdaBest);
double corridor_method(IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
				IloObjective & obj, int * xLBest, double & zBest, double ** rc, int * xIlo);
double refine_solution(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj, int * xL, 
				double & cWidth, int & nSol, double & zBest, int * xIlo, int & iterImproved, 
				int iter, double ** rc, bool fix2zero, bool fix2one, IloRangeArray & cutSol,
				int lagrIter);
void find_best(int * F, double ** rc);
void find_worst(int ** F0, double ** rc, double propFixed0);
void cut_out_solution(IloModel & model, IloCplex & cplex, TwoD & x_ilo, int * xIlo, IloRangeArray & cutSol);
void add_z_cut(IloModel & model, IloCplex & cplex, TwoD & x_ilo, double zBest);

double compute_fitness_value(double corridorWidthBase, int nSolBase, double propFixed0, 
				double propFixed1, int add_z_cut);
int stopping_criteria(int iter, int max_iter, double time_limit);

double solve_robust_problem(INSTANCE & inp, IloModel & model, IloCplex & cplex, 
				TwoD & x_ilo, IloObjective & obj, IloNumVar & Q_ilo, 
				double Omega, double * sigma2, int * xNominal, double zNominal);

double solve_nominal_problem(INSTANCE & inp, IloModel & model, IloCplex & cplex, 
    TwoD & x_ilo, IloObjective & obj, int * xNominal, double & zNominal);

void computeSigmaSq(INSTANCE inp, double * sigma2);
/************************ FUNCTIONS ******************************/




/************************ MAIN PROGRAM ******************************/
/// Main program
/************************ MAIN PROGRAM ******************************/
int main(int argc, char *argv[])
{
 
    srand(time(0));

		//freopen("debug.txt", "w", stdout); //!< redirect output to a file
		int err = parseOptions(argc, argv);
		if ( err != 0)
		{ if (err != -1) cout << "Error argument " << err+1 << endl; exit(1); }


		tTime.resetTime();		      // start clock

		read_problem_data(inp); // read instance data
		print_options(inp);

        sigma2 = new double[inp.nC];
        computeSigmaSq(inp, sigma2);
		tTime.resetTime();

		// vectors used to defined variables fixed to 0 and 1
		F = new int[inp.nC];
		F0 = new int*[inp.nC];
		for (int i = 0; i < inp.nC; i++)
				F0[i] = new int[inp.ri[i]];
		zBest     = ZERO;    
		ubStar    = INFTY;

		/// cplex stuff
		IloEnv env;
		IloModel model(env);
		IloCplex cplex(model);
		TwoD x_ilo(env, inp.nC);
		for (int i = 0; i < inp.nC; i++)
				x_ilo[i] = IloNumVarArray(env, inp.ri[i], 0, 1, ILOINT);
		IloObjective obj = IloAdd(model, IloMaximize(env, 0));

		IloNumVar Q_ilo(env, 0.0, IloInfinity, IloNumVar::Float);

        double zNominal = -1.0; // obj function value nominal pbr
        int * xNominal = new int[inp.nC];
        solve_nominal_problem(inp, model, cplex, x_ilo, obj, xNominal, zNominal); 

		solve_robust_problem(inp, model, cplex, x_ilo, obj, Q_ilo, Omega, 
                sigma2, xNominal, zNominal);


		// call lagrangean phase (subgradient optimization)
		lagrangean_phase(inp, model, cplex, x_ilo, obj, Q_ilo, Omega, sigmaSq);

        robust_lagrangean_phase();

		fsol << _FILENAME << "\t" << zBest << "\t" << best_time << "\t" 
				<< ubStar << "\t" << bestIter << "\t" << corridorWidthBase 
				<< "\t" << nSolBase << "\t" << propFixed0 << "\t" 
				<< propFixed1 << "\t" << add_oldSol_cut << endl;

		flagr.close();
		fsol.close();

		env.end();

		return zBest;		// return value to brkGA for 

}
/************************ FUNCTIONS ******************************/
/// Functions
/************************ FUNCTIONS ******************************/
void read_problem_data(INSTANCE & inp)
{
   


		int temp;

		ifstream fdata(_FILENAME, ios::in);
		if (!fdata)
		{
				cerr << "Cannot open file " << _FILENAME << endl; exit(1);
		}

		fdata >> inp.nC >> temp >> inp.nR;

		// all the classes have the same nr of elements
		inp.ri = new int[inp.nC];
		inp.c  = new double*[inp.nC];
		inp.w  = new int**[inp.nC];
		for (int i = 0; i < inp.nC; i++)
		{
				inp.ri[i] = temp;
				inp.c[i] = new double[inp.ri[i]];
				inp.w[i] = new int*[inp.ri[i]];
				for (int j = 0; j < inp.ri[i]; j++)
						inp.w[i][j] = new int[inp.nR];
		}

		// read rhs of knapsack constraints
		inp.R = new int[inp.nR];
		for (int k = 0; k < inp.nR; k++)
				fdata >> inp.R[k];

		// read data for each class
		for (int i = 0; i < inp.nC; i++)
		{
				fdata >> temp;
				assert(temp == (i+1));

				for (int j = 0; j < inp.ri[i]; j++)
				{
						fdata >> inp.c[i][j];
						for (int k = 0; k < inp.nR; k++)
								fdata >> inp.w[i][j][k];
				}
		}
        

}

void computeSigmaSq(INSTANCE inp, double * sigma2)
{
    // some statistics
    double * avg = new double[inp.nC];
    double tot = 0.0;

    for (int i = 0; i < inp.nC; i++)
    {
        // avg profit per group
        avg[i]= 0.0;
        for (int j = 0; j < inp.ri[i]; j++)
        {
            avg[i] += inp.c[i][j];
        }
        avg[i] /= (double)inp.ri[i];
        tot += avg[i];

        cout << "Avg Profit Class " << i << " = " << avg[i] << endl;
    }

    int abc;
    cin >> abc;

    for (int i = 0; i < inp.nC; i++)
    {
       sigma2[i] =  0.8 + avg[i]/tot;
       sigma2[i] = 1.0;
        cout << "SigmaSq[" << i << "] = " << sigma2[i] << endl;
    }
}

/// Print instance info and algorithmic parameters.
void print_options(INSTANCE inp)
{
		cout << "-------------------------------------" << endl;
		cout << "- INSTANCE DATA : " << endl;
		cout << "-------------------------------------" << endl;
		cout << " Nr. Classes\t \t :: " << inp.nC << endl;
		cout << " Nr. Resources\t \t :: " << inp.nR << endl;
		cout << " Nr. Items per class\t :: " << inp.ri[0] << endl;
		cout << "-------------------------------------" <<  endl << endl;   
		cout << "-------------------------------------" << endl;
		cout << "- ALGORITHMIC PARAMETERS : " << endl;
		cout << "-------------------------------------" << endl;
		cout << " Time Limit \t \t :: " << time_limit << endl;
		cout << " Iter Limit \t \t :: " << max_iter << endl;
		cout << " Corridor Base\t \t :: " << corridorWidthBase << endl;
		cout << " Base Sol. Nr.\t \t :: " << nSolBase << endl;
		cout << " % Fix to Zero \t \t :: " << propFixed0 << endl;
		cout << " % Fix to One \t \t :: " << propFixed1 << endl;
		cout << " With Cut Sol \t \t :: " << add_oldSol_cut << endl;
        cout << " Omega \t \t \t :: " << Omega << endl;
		cout << " Sigma^2 \t \t :: " << sigmaSq << endl;
		cout << "-------------------------------------" << endl;
		cout << "-------------------------------------" << endl << endl;;


}

/// Define stopping criteria of the algorithm
int stopping_criteria(int iter, int max_iter, double time_limit)
{
		return ( (tTime.elapsedTime(timer::REAL) >= time_limit) ||
						(iter >= max_iter) );
}

/// Lagrangean phase
/*  The implementation of the lagrangean relaxation coupled with subgradient
	optimization is defined in the file "lagrange.cpp." 

	The Lagrangean relaxation method produces both UPPER and LOWER bounds:
	- at each lagrange iteration, i.e., for each set of lagrange multipliers, we
	solve to optimality the corresponding relaxed problem (we obtain a valid
	UPPER BOUND zL)
	- starting from the infeasible lagrangean solution, we apply the corridor
	method to look for feasible solutions in the neighborhood of the lagrangean
	solution. Any feasible solution obtained during this refining procedure is
	a valid LOWER BOUND.

	The key idea is:
	- we repeat the lagrangean cycle 3 times (each time with a different
	starting set of lagrange multipliers)
	- each cycle is divided into two main parts:
	i.  the first, e.g., 100 iterations are used to quickly get near-optimal
	lagrange multipliers
	ii. once near-optimal multipliers are available, we periodically call a 
	refining procedure (based on the corridor method) to achieve better 
	lower bounds.
 **/
void lagrangean_phase(INSTANCE & inp, IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
				IloObjective & obj, IloNumVar & Q_ilo, double Omega, double sigma)
{

		IloEnv env = model.getEnv();
		IloRangeArray cutSol(env);

		map<int,int> zList;
		// initialize data structure for Lagrangean phase
		double * lambda     = new double[inp.nR];
		double * lambdaBest = new double[inp.nR];
		int * xL            = new int[inp.nC]; // current lagrangean solution
		int * xLBest        = new int[inp.nC]; // best lagrangean solution (infeasible)
		int * xIlo          = new int[inp.nC]; // best feasible solution from cplex
		double lb           = ZERO;

		zBest = get_first_lb(model, cplex, x_ilo, obj, 1, Q_ilo, Omega, sigma);// nr feasible solutions to be found by cplex
		cout << "First feasible solution = " << zBest << endl;

		lambda_initialization(lambda);

#ifdef M_DEBUG
		for (int k = 0; k < inp.nR; k++)
				cout << "l(" << k << ") = " << lambda[k] << endl;
#endif

		double ** rc = new double*[inp.nC]; // lagrangean reduced costs
		for (int i = 0; i < inp.nC; i++)
				rc[i] = new double[inp.ri[i]];


		int lagrIter = 0;
		int iter     = 0;
		while (lagrIter < 3 && !stopping_criteria(iter, max_iter, time_limit))
		{
				double cWidth     = corridorWidthBase;   // <---------- used to be 0.8
				int iterImproved  = 200;
				int nSol          = nSolBase;		// nr of solutions for cplex
				int Freq          = 40;

				double zL, bestLagr, worstLagr;

				double start300  = INFTY;
				double best300   = INFTY;
				double delta     = 0.1;
				bool stopping    = false;
				bool fix2zero    = false;
				bool fix2one     = false;

				while (!stopping_criteria(iter, max_iter, time_limit) && !stopping)
				{
						// note: zL in an "upper bound" of the optimal value
						zL = lagrange_step(lambda, xL, iter, rc);
						if (zL < ubStar)
								update_best(xLBest, xL, ubStar, zL, lambda, lambdaBest);

						stopping = lambda_update(lambda, delta, xL, bestLagr, worstLagr, zL, zBest, iter, 
										best300, start300);

						if (iter > 199) 	// start fixing schemes
						{
								fix2one = true;
								fix2zero = true;
						}

						// refine lagrange solution to obtain a valid lower bound
						if (iter > 100 && (iter % Freq) == 0)
								lb = refine_solution(model, cplex, x_ilo, obj, xL, cWidth, nSol, zBest, 
												xIlo, iterImproved, iter, rc, fix2zero, fix2one, 
												cutSol, lagrIter);

						if ((iter > 250) && (iter - iterImproved) > 100)
						{
								nSol++;
								cWidth       *= 0.9;
								iterImproved  = iter;
								Freq         /= 2;
								if (Freq < 5) Freq = 5;

								//cout << "Enhancing corridor to : " << cWidth << " and nsol to " << nSol << endl;
						}

						iter++;
				} // end lagrangean cycle

				cout << "** ** ** SUMMARY ** ** ** " << endl;
				cout << "** ** ** Lagrangean Iteration Nr. " << lagrIter << " :: Best LB = " << zBest << endl;
				cout << setw(49) << " :: Best UB = " << ubStar << endl;
#ifdef M_DEBUG
				cout << "Before removing cuts " << cplex.getNrows() << " constraints " << endl;
#endif

				model.remove(cutSol);

#ifdef M_DEBUG
				cout << "After removing cuts " << cplex.getNrows() << " constraints " << endl;
#endif

				lagrIter++;
				iter = 0;		// reset counter

				if (lagrIter < 3)
						lambda_random_perturbation(lambda, lambdaBest);

		} // END OF LAGRANGEAN CYCLE -- repeat 3 times

		//fix_to_zero(rc, ubStar, zBest);
		//corridor_method(model, cplex, x_ilo, obj, xIlo, zBest, rc, xIlo);
		//corridor_method(model, cplex, x_ilo, obj, xLBest, zBest, (double **)inp.c);
}

/// Update best solution
void update_best(int * xLBest, int * xL, double & ubStar, double zL, double * lambda, double * lambdaBest)
{
		ubStar = zL;
		for (int i = 0; i < inp.nC; i++)
				xLBest[i] = xL[i];

		for (int k = 0; k < inp.nR; k++)
				lambdaBest[k] = lambda[k];
}

/// Get first feasible solution
double get_first_lb(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj, 
				int nSol, IloNumVar & Q_ilo, double Omega, double sigma)
{

		//defineModel(model, cplex, x_ilo, obj); // get the first cplex solution

		//defineRobustModel(model, cplex, x_ilo, obj, Q_ilo, Omega, sigma);

		// cycle for all the possible values of u in W	
		for (int cc = 0; cc < inp.nC; cc++)
		{


				double u = (double)(cc+1)*sigma;	
				cout << "** ** ** U[" << cc+1 <<"] = " << u << " ** ** ** " << endl;
				defineRobustDet(model, cplex, x_ilo, obj, Omega, sigma, u);
				double statusBin = solve_KNAP(model, cplex, 99999, 4, 10000); // call cplex
				int abc;
				cin >> abc;
		}


		//double statusBin = solve_KNAP(model, cplex, nSol, 0, 10); // call cplex

		double statusBin = solve_KNAP(model, cplex, 99999, 4, 10000); // call cplex
		cout << "Status :: " << statusBin << endl;
		cout << "CPLEX = " << cplex.getStatus() << endl;

		// get solution
		int * xIlo = new int[inp.nC];

		for (int i = 0; i < inp.nC; i++)
				for (int j = 0; j < inp.ri[i]; j++)
						if (cplex.getValue(x_ilo[i][j]) >= 1.0 - EPSI)
								xIlo[i] = j;		
		for (int i = 0; i < inp.nC; i++)
				cout << "x[" << i <<  "] = " << xIlo[i] << endl;

		exit(145);
		return statusBin;
}

/// NO LONGER USED
double corridor_method(IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
				IloObjective & obj, int * xLBest, double & zBest, double ** rc, int * xIlo)
{
		IloEnv env = model.getEnv();

#ifdef W_LAGR_COSTS
		// solve using lagrangean costs + feasibility (i.e., knapsack) constraints
		for (int i = 0; i < inp.nC; i++)
				for (int j = 0; j < inp.ri[i]; j++)
						obj.setLinearCoef(x_ilo[i][j], rc[i][j]);
		model.add(obj);
#endif


		// add corridor around best feasible solution
		IloExpr lhsCorridor(env);
		IloRangeArray neighborhood(env);
		double rhs = 0.98*(double)(inp.nC);

		for (int i = 0; i < inp.nC; i++)
				lhsCorridor += x_ilo[i][xIlo[i]];

		neighborhood.add(lhsCorridor >= rhs);
		model.add(neighborhood);


		int iter = 0;
		bool stopping = false;
		while(!stopping)
		{
				double statusBin = solve_KNAP(model, cplex, 10000, 3, 300);	// call cplex

				cout << "changing objective coefficients to original costs and restart ... " << endl;
				// evaluate solution (this is needed in the first iter)
				cout << "z[" << iter << "]  = " << get_cplex_sol(model, cplex, x_ilo, xIlo) << endl;

				cout << "CPLEX Status :: " << cplex.getStatus() << endl;

#ifdef W_LAGR_COSTS
				for (int i = 0; i < inp.nC; i++)
						for (int j = 0; j < inp.ri[i]; j++)
								obj.setLinearCoef(x_ilo[i][j], inp.c[i][j]);
				model.add(obj);
#endif

				// remove neighborhood constraint
				model.remove(neighborhood);
				lhsCorridor.end();
				cout << " .... corridor removed ... " << endl;

				iter++;

				if (iter > 1) stopping = true;

		}

		return  -1.0;
}

/// Refine lagrange solution to obtain a valid lower bound (i.e., a feasible solution)
double refine_solution(IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
				IloObjective & obj, int * xL, double & cWidth, int & nSol, 
				double & zBest, int * xIlo, int & iterImproved, int iter, 
				double ** rc, bool fix2zero, bool fix2one, IloRangeArray & cutSol, 
				int lagrIter)
{

		//int  * F;
		//int ** F0;
		//int nFixed1       = 0;
		//int propFix1      = 0;
		//double propFixed0 = ZERO;

		int nFixed1;
		double width1;

		if (fix2zero)
				propFixed0 = 0.75;	// used to be 0.75

		if (fix2one)
		{
				nFixed1  = ceil(propFixed1*(double)inp.nC); // used to be 0.25
				width1 = ceil(cWidth*(double)nFixed1);
		}


		IloEnv env = model.getEnv();


		// add corridor constraint
		IloExpr lhsCorridor(env);
		IloRangeArray neighborhood(env);
		double rhs = cWidth*(double)(inp.nC);

		for (int i = 0; i < inp.nC; i++)
				lhsCorridor += x_ilo[i][xL[i]];

		neighborhood.add(lhsCorridor >= rhs);
		model.add(neighborhood);


		IloExpr lhs0(env);
		IloRangeArray fix_to_zero(env);
		IloExpr lhs(env);
		IloRangeArray cut(env);
		int count0 = 0;


		if (fix2zero)
		{

				for (int i = 0; i < inp.nC; i++)
						for (int j = 0; j < inp.ri[i]; j++)
								F0[i][j] = 0;

				find_worst(F0, rc, propFixed0);    

				for (int i = 0; i < inp.nC; i++)
						for (int j = 0; j < inp.ri[i]; j++)
								if (F0[i][j] == 1)
								{
										lhs0 += x_ilo[i][j];
										fix_to_zero.add(x_ilo[i][j] == ZERO);
										count0++;
								}
				//fix_to_zero.add(lhs0 == ZERO);
				//fix_to_zero.add(lhs0 <= (1.0-cWidth)*count0); //hard
				//fix_to_zero.add(lhs0 <= cWidth*count0);  // soft
				model.add(fix_to_zero);
		}


		if (fix2one)
		{
				for (int i = 0; i < inp.nC; i++)
						F[i] = -1;

				for (int k = 0; k < nFixed1; k++)
						find_best(F, rc);

				for (int i = 0; i < inp.nC; i++)
						if (F[i] != -1)
								lhs += x_ilo[i][F[i]];

				//cout << " ---- FIxED = " << ceil(cWidth*(double)nFixed) << " over " << nFixed << endl << endl;

				cut.add(lhs >= width1);
				model.add(cut);
		}

		//cout << "Nr constraints after adding ALL CUTS is " << cplex.getNrows() << endl;

		double statusBin = solve_KNAP(model, cplex, nSol, 0, 10);	// call cplex

		if (statusBin >= -1.0+EPSI) 
		{
#ifdef M_DEBUG
				cout << " ... repaired solution is " << statusBin << "(z* = " <<  zBest << ") with cplex status :: " 
						<< cplex.getStatus() << ". Params: [ zeros = " << count0 << ", ones = " << width1 << "/" 
						<< nFixed1 << "]" << endl;

				cout << "Nr. Constraints before cutting out solutions is : " << cplex.getNrows() << endl;
#endif

				int * xRepaired = new int[inp.nC];
				get_cplex_sol(model, cplex, x_ilo, xRepaired);
				//cut_out_solution(model, cplex, x_ilo, xRepaired, cutSol);

#ifdef M_DEBUG
				cout << "Nr. Constraints after cutting out solutions is : " << cplex.getNrows() << endl;
#endif

				if (statusBin > (zBest+EPSI))
				{
						zBest        = statusBin;
						best_time    = tTime.elapsedTime(timer::REAL); // measure wall-clock time
						bestIter     = lagrIter;
						iterImproved = iter; // save iteration of last improvement
						cWidth       = corridorWidthBase;
						nSol         = nSolBase;
						for (int i = 0; i < inp.nC; i++) 
								xIlo[i] = xRepaired[i];

						if(add_oldSol_cut)
								add_z_cut(model, cplex, x_ilo, statusBin);

						cout << "... improved best lb(" << iter << ") = " << zBest << " found after " 
								<< best_time << " seconds." << endl;

#ifdef M_DEBUG
						cout << "... corridor width restored to " << cWidth << endl;
						cout << "... nSol restored to " << nSol << endl;
#endif


#ifdef M_DEBUG
						if (fix2one)
						{
								// study cplex vs lagrangean solution
								int cc = 0;
								for (int i = 0; i < inp.nC; i++)
								{
										cout << "here .... " << xL[i] << " vs " << xIlo[i] << endl;
										if (xL[i] != xIlo[i])
												cc++;

										cout << "COMPONENT " << i << " fixed element is " << F[i] << endl;
										if (F[i] != -1 && F[i] != xIlo[i])
										{
												cout << " ++++++++++++++++++++++++++++++++++" << endl;
												int aka;
												cin >> aka;
										}

										if (F0[i][xIlo[i]] == 1)
										{
												cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
												int aka;
												cin >> aka;
										}


										for (int j = 0; j < inp.ri[i]; j++)
										{
												if (xIlo[i] == j)
														cout << "x(" << i << "," << j << ") = 1  vs rc = " << rc[i][j] << " vs " << F0[i][j] << endl;
												else
														cout << "x(" << i << "," << j << ") = 0  vs rc = " << rc[i][j] << " vs " << F0[i][j] << endl;
										}
								}

								cout << "Nr. different between lagr and cplex = " << cc << endl;
								int abc;
								cin >> abc;
						}
#endif
				}

		}

		if (fix2one)
		{
				// remove fix_to_one constraint
				model.remove(cut);
				cut.end();
				lhs.end();
		}

		if (fix2zero)
		{
				// remove fix to zero
				model.remove(fix_to_zero);
				lhs0.end();
				fix_to_zero.end();
		}


		// remove neighborhood constraint
		model.remove(neighborhood);
		lhsCorridor.end();

		//cout << "Nr constraints after removing ALL CUTS is " << cplex.getNrows() << endl;

		return statusBin;
}

/// Find variables with best lagrangean cost
void find_best(int * F, double ** rc)
{
		double bestRC = -INFTY;
		int group = -1;
		int item  = -1;
		for (int i = 0; i < inp.nC; i++)
		{
				if (F[i] != -1) continue;
				for (int j = 0; j < inp.ri[i]; j++)
						if (rc[i][j] > bestRC)
						{
								bestRC = rc[i][j];
								group  = i;
								item   = j;
						}
		}
		if (group != -1)
				F[group] = item;
		else
		{
				cout << "Unable to fix variable " << endl;
				exit(144);
		}
}

/// Find variables with worst lagrangean costs
void find_worst(int ** F0, double ** rc, double propFixed0)
{

		for (int i = 0; i < inp.nC; i++)
		{
				int nFixed = ceil(propFixed0*(double)inp.ri[i]);	
				//cout << "for module " << i << " we fix " << nFixed << " elements " << endl;

				for (int counter = 0; counter < nFixed; counter++)
				{
						double worstRC = INFTY;
						int item       = -1;
						for (int j = 0; j < inp.ri[i]; j++)
						{
								if (F0[i][j] == 1) continue;

								if (rc[i][j] < worstRC)
								{
										worstRC = rc[i][j];
										item    = j;
								}
						}
						if (item != -1)
								F0[i][item] = 1;
						else
						{
								cout << "something wrong here ... " << endl;
								exit(123);
						}
				}
		}
}


/// Add cut to exclude current solution from feasible space
void cut_out_solution(IloModel & model, IloCplex & cplex, TwoD & x_ilo, int * xIlo, IloRangeArray & cutSol)
{
		IloEnv env = model.getEnv();

		IloExpr lhs(env);
		//IloRangeArray cutSol1(env);
		for (int i = 0; i < inp.nC; i++)
				lhs += x_ilo[i][xIlo[i]];	 

		// for (int j = 0; j < inp.ri[i]; j++)
		//     if (cplex.getValue(x_ilo[i][j]) >= 1.0 - EPSI)
		//  	lhs += x_ilo[i][j];



		cutSol.add(lhs <= inp.nC-1);
		model.add(cutSol);

}

/// Add cut with obj function value of best found solution
void add_z_cut(IloModel & model, IloCplex & cplex, TwoD & x_ilo, double zBest)
{
		IloEnv env = model.getEnv();

		IloExpr lhs(env);
		//IloRangeArray cutSol1(env);
		for (int i = 0; i < inp.nC; i++)
				for (int j = 0; j < inp.ri[i]; j++)
						lhs += (double)inp.c[i][j]*x_ilo[i][j];

		model.add(lhs >= zBest);    
}
