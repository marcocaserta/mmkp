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

/*! \mainpage Multiple-choice Multidimensional Knapsack Problem
 
  \author Marco Caserta 2015 (c) 
  \version v. 1.0.0
  \date Start Date : 30.01.15
  \date Last Update: 
**/


#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cassert>

#include "cplex.h"


const double EPSI      = 0.00000001;

struct INSTANCE {
    int nR;			// number of resources (constraints)
    int nC;			// number of classes
    int * ri;			// number of items in each class i
  
    double  ** c;
    int *** w;
    int   * R;
};

extern INSTANCE inp;

void add_cut_corridor(IloModel & model, IloCplex & cplex, TwoD & x_ilo, int * xCorridor, int rhs, IloRangeArray & corridor, IloExpr lhs)
{
    IloEnv env = model.getEnv();

    for (int i = 0; i < inp.nC; i++)
	lhs += x_ilo[i][xCorridor[i]];
    
    corridor.add(lhs >= rhs);

    model.add(corridor);
}


double defineModel(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj)
{
    IloEnv env = model.getEnv();

    // knapsack constraints
    for (int k = 0; k < inp.nR; k++)
    {
	IloExpr sum(env);
	for (int i = 0; i < inp.nC; i++)
	    for (int j = 0; j < inp.ri[i]; j++)
		sum += (inp.w[i][j][k]*x_ilo[i][j]);
	
	model.add(sum <= inp.R[k]);
	sum.end();
    }

    // multi-choice constraint
    for (int i = 0; i < inp.nC; i++)
    {
	IloExpr sum(env);
	for (int j = 0; j < inp.ri[i]; j++)
	    sum += x_ilo[i][j];
	
	model.add(sum == 1.0);
    }

    // add objective function
    for (int i = 0; i < inp.nC; i++)
     	for (int j = 0; j < inp.ri[i]; j++)
	    obj.setLinearCoef(x_ilo[i][j], inp.c[i][j]);
    model.add(obj);
     
}

/// Solve Multiple-choice Multidimensional Knapsack Problem
double solve_KNAP(IloModel model, IloCplex cplex, int solLimit, int displayLimit, int timeLim)
{
    try
    {

	IloEnv env = model.getEnv();

	//cplex.setOut(env.getNullStream());

	//cplex.setParam(IloCplex::ItLim, 1000); // Iterations limit
	//cplex.setParam(IloCplex::TiLim, c_time); // Time limit
	//cplex.setParam(IloCplex::ClockType, 2); // Clock type (Wall clock)
	//cplex.setParam(IloCplex::AdvInd, 2); // MIP Start
	//cplex.setParam(IloCplex::RINSHeur, 1); // how often to apply the RINS heuristic

        cplex.setParam(IloCplex::MIPInterval,10000); // after how many iterations to print
	cplex.setParam(IloCplex::MIPDisplay, displayLimit); 
	
	// set number of solutions to be obtained before stopping
	cplex.setParam(IloCplex::IntSolLim, solLimit);
	cplex.setParam(IloCplex::TiLim, timeLim); // Time limit

	// Optimize the problem and obtain solution.
	if ( !cplex.solve() ) 
	{
	    //env.error() << "Failed to optimize MIP." << endl;
	    throw(-1.0);
	}
	
	//cout << "Best OBj Value is " << cplex.getBestObjValue() << endl;
	//cout << "Status :: " << cplex.getStatus() << endl;

	return cplex.getObjValue();
    }
    catch (IloException & e) 
    {
	//cout << "Cplex exception caught: " << e << "." << endl;
	return -1.0;
    }
    catch (...) 
    {
	//cout << "Cplex unknown exception caught. " << endl;
	return -1.0;
    }
}

double get_cplex_sol(IloModel & model, IloCplex & cplex, TwoD & x_ilo, int * xIlo)
{

    for (int i = 0; i < inp.nC; i++)
	for (int j = 0; j < inp.ri[i]; j++)
	    if (cplex.getValue(x_ilo[i][j]) >= 1.0 - EPSI)
		xIlo[i] = j;		
    
    // verify obj function value and return solution
    int z = 0;
    for (int i = 0; i < inp.nC; i++)
	z += inp.c[i][xIlo[i]];

    assert( (cplex.getObjValue() - (double)z) <= EPSI);

#ifdef W_SOL	    
    cout << "CPLEX solution has z = " << z << " and components are :: " << endl;
    for (int i = 0; i < inp.nC; i++)
     	cout << setw(4) << xIlo[i];
    cout << endl;
#endif
    
    return cplex.getObjValue();
}
