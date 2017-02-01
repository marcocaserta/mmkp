/********************************en*******************************************
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

double solve_nominal_problem(INSTANCE & inp, IloModel & model, IloCplex & cplex, 
        TwoD & x_ilo, IloObjective & obj, int * xNominal, double & zNominal)
{
    double statusBin = -1;
    defineModel(model, cplex, x_ilo, obj); // get the first cplex solution
    statusBin = solve_KNAP(model, cplex, 99999, 0, 10000); // call cplex
    zNominal = cplex.getObjValue();
    cout << "Nominal Problem : ====== " << endl;
    cout << "Status :: " << statusBin << endl;
    cout << "CPLEX = " << cplex.getStatus() << endl;
    cout << "z*    = " << cplex.getObjValue() << endl; 
    // get solution

    for (int i = 0; i < inp.nC; i++)
        for (int j = 0; j < inp.ri[i]; j++)
            if (cplex.getValue(x_ilo[i][j]) >= 1.0 - EPSI)
                xNominal[i] = j;		
    for (int i = 0; i < inp.nC; i++)
        cout << "x[" << i <<  "] = " << xNominal[i] << endl;

}



double solve_robust_problem(INSTANCE & inp, IloModel & model, IloCplex & cplex, 
        TwoD & x_ilo, IloObjective & obj, IloNumVar & Q_ilo, 
        double Omega, double * sigma2, int * xNominal, double zNominal)
{
    //cout << "Defining problem instance DET " << endl;
    //defineModel(model, cplex, x_ilo, obj); // get the first cplex solution
    //defineRobustModel(model, cplex, x_ilo, obj, Q_ilo, Omega, sigmaSq);
    //double statusBin = solve_KNAP(model, cplex, 99999, 4, 10000); // call cplex

    // compute sum of squared sigma
    double SSsigma  = 0.0;
    for (int i = 0; i < inp.nC; i++)
        SSsigma += sigma2[i];


    // cycle for all the possible values of u in W	
    int * xIlo = new int[inp.nC];
    double statusBin = -1;
    for (int cc = inp.nC; cc <= inp.nC; cc++)
    {
        double u = (double)(cc);	
        cout << "** ** ** U[" << cc <<"] = " << u << " ** ** ** " << endl;
        try
        {

            defineRobustDet(model, cplex, x_ilo, obj, Omega, SSsigma, u);
            //defineRobustBertsimas(model, cplex, x_ilo, obj, Omega, sigmaSq, u);

            statusBin = solve_KNAP(model, cplex, 99999, 0, 10000); // call cplex
        }
        catch(IloException& e)
        {
            cerr << "Exception CPLEX " << endl;
            e.end();
            throw;
        }
        catch (...)
        {
            cerr << "Unknown exception " << endl;
            throw;
        }


        cout << "Status :: " << statusBin << endl;
        cout << "CPLEX = " << cplex.getStatus() << endl;
        int totInfeasible = 0;
        int totInfeasibleNominal = 0;
        int nRuns = 0;
        if (statusBin == -1)
        {
            cout << "No feasible solution found " << endl; 
        }

        else // cplex has found a solution
        {

            cout << "z*    = " << cplex.getObjValue() << endl; 
            for (int i = 0; i < inp.nC; i++)
                for (int j = 0; j < inp.ri[i]; j++)
                    if (cplex.getValue(x_ilo[i][j]) >= 1.0 - EPSI)
                        xIlo[i] = j;

            // cout << "z*[" << u << "," << sigmaSq*u <<"] = " 
            // << cplex.getObjValue();
            cout << "z*[" << Omega << "," <<  SSsigma/(double)inp.nC <<"] = " 
                << cplex.getObjValue();

            for (int i = 0; i < inp.nC; i++)
                cout << " " << xIlo[i];
            cout << endl;

            // get solution
            //int * xIlo = new int[inp.nC];

            for (int i = 0; i < inp.nC; i++)
                for (int j = 0; j < inp.ri[i]; j++)
                    if (cplex.getValue(x_ilo[i][j]) >= 1.0 - EPSI)
                        xIlo[i] = j;
            for (int i = 0; i < inp.nC; i++)
                cout << "x[" << i <<  "] = " << xIlo[i] << endl;


            // now evaluate solution with random generated 
            int ** wN = new int*[inp.nC];
            for (int i = 0; i < inp.nC; i++)
                wN[i] = new int[inp.nR];

            int ** wR = new int*[inp.nC];
            for (int i = 0; i < inp.nC; i++)
                wR[i] = new int[inp.nR];

            // =============== EVAL ======================
            double * rVector = new double[inp.nC*nRuns];
            for (int cc = 0; cc < nRuns; cc++)
            {
                cout << "Iter " << cc << endl;

                for (int i = 0; i < inp.nC; i++)
                {
                    for (int k = 0; k < inp.nR; k++)
                    {
                        double rr = (double)(rand()+1)/((double)(RAND_MAX)+1.0);
                        rVector[cc*inp.nR + k] = rr;
                        if (rr <= 0.5)
                        {
                            wR[i][k] = inp.w[i][xIlo[i]][k] - sigma2[i];
                            wN[i][k] = inp.w[i][xNominal[i]][k] - sigma2[i];
                        }
                        else
                        {
                            wN[i][k] = inp.w[i][xNominal[i]][k] + sigma2[i];
                            wR[i][k] = inp.w[i][xIlo[i]][k] + sigma2[i];
                        }
                    }
                }

                // this is what we have done
                cout << "Comparison of nominal vs real :: " << endl;
                for (int i = 0; i < inp.nC; i++)
                {
                    cout << "Item " << xIlo[i] << " ... " << endl;
                    for (int k = 0; k < inp.nR; k++)
                    {
                        cout << inp.w[i][xIlo[i]][k] << " vs " << wR[i][k] << endl;
                    }
                }

                // is it feasible?
                for (int k = 0; k < inp.nR; k++)
                {
                    double tot = 0.0;
                    for (int i = 0; i < inp.nC; i++) 
                        tot += wR[i][k];
                    // cout << "lhs = " << tot << " vs rhs = " << inp.R[k] << endl;
                    if (tot > inp.R[k])
                    {
                        // cout << "Infeasible in constraint " << k << endl;
                        totInfeasible++;

                        break;
                    }
                }

                // is nominal solution feasible?
                for (int k = 0; k < inp.nR; k++)
                {
                    double tot = 0.0;
                    for (int i = 0; i < inp.nC; i++) 
                        tot += wN[i][k];
                    // cout << "lhs = " << tot << " vs rhs = " << inp.R[k] << endl;
                    if (tot > inp.R[k])
                    {
                        // cout << "Infeasible in constraint " << k << endl;
                        totInfeasibleNominal++;

                        break;
                    }
                }
            }
        }
        cout << " [" << totInfeasible << "/" << nRuns << "]" 
            << endl;

        cout << " [" << totInfeasibleNominal << "/" << nRuns << "]" 
            << endl;

        string sRatio;
        if (statusBin != -1)
        {
            double ratio = (zNominal - statusBin)/zNominal;
            stringstream ss;
            ss << ratio;
            sRatio = ss.str();
        }
        else
            sRatio = "-";



        // write info to diskfile
        ofstream fRandom("robust.txt", ios::out);
        fRandom << zNominal << "\t" << totInfeasibleNominal << "\t" 
            << statusBin << "\t" << sRatio << "\t" << Omega << "\t" << 
            SSsigma/(double)inp.nC << "\t" 
            << totInfeasible << "\t" << nRuns << endl; 
        fRandom.close();
    }

}



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
// Quadratic Model SOCP
double defineRobustModel(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj, IloNumVar & Q_ilo, double Omega, double sigma)
{
    IloEnv env = model.getEnv();

    // robust constraint
    cout << "Omega and sigma are ;" << Omega << " " << sigma << endl;

    IloExpr lhsRobust(env);
    double rhs = 0.0;
    lhsRobust = -Q_ilo*Q_ilo;
    for (int i = 0; i < inp.nC; i++)
        for (int j = 0; j < inp.ri[i]; j++)
            lhsRobust += x_ilo[i][j]*x_ilo[i][j]*sigma;

    model.add(lhsRobust <= rhs);

    for (int k = 0; k < inp.nR; k++)
    {
        IloExpr sum(env);
        for (int i = 0; i < inp.nC; i++)
            for (int j = 0; j < inp.ri[i]; j++)
                sum += (inp.w[i][j][k]*x_ilo[i][j]);

        model.add(sum + Omega*Q_ilo <= inp.R[k]);
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


double defineRobustBertsimas(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj, double Omega, double sigma, double u)
{
    IloEnv env = model.getEnv();

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

    // robust constraint
    cout << "Omega and sigma are :" << Omega << " " << sigma << endl;
    cout << "... with u = " << u << endl;
    double CC = sqrt(sigma)*Omega*sqrt(u)/2.0;
    double coeff = 0.0;
    for (int k = 0; k < inp.nR; k++)
    {
        IloExpr sum(env);
        for (int i = 0; i < inp.nC; i++)
            for (int j = 0; j < inp.ri[i]; j++)
            {
                coeff = sqrt(sigma)*Omega/(2.0*sqrt(u)) + inp.w[i][j][k];
                //cout << "coeff(" <<i<<","<<j<<") is " << coeff << endl;
                sum += (coeff*x_ilo[i][j]);
            }

        model.add(sum + CC <= inp.R[k]);
        sum.end();	
    }

}


double defineRobustDet(IloModel & model, IloCplex & cplex, TwoD & x_ilo,
        IloObjective & obj, double Omega, double SSsigma, double u)
{
    IloEnv env = model.getEnv();

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

    // robust constraint
    /*
       double CC = sqrt(sigma)*Omega*sqrt(u)/2.0;
       double coeff = 0.0;
       for (int k = 0; k < inp.nR; k++)
       {
       IloExpr sum(env);
       for (int i = 0; i < inp.nC; i++)
       for (int j = 0; j < inp.ri[i]; j++)
       {
       coeff = sqrt(sigma)*Omega/(2.0*sqrt(u)) + inp.w[i][j][k];
       sum += (coeff*x_ilo[i][j]);
       }

       model.add(sum + CC <= inp.R[k]);
       sum.end();	
       }
       */

    // parameter "u" is used to model uncertainty level
    // Note: "u" is no longer used
    for (int k = 0; k < inp.nR; k++)
    {
        IloExpr sum(env);
        for (int i = 0; i < inp.nC; i++)
            for (int j = 0; j < inp.ri[i]; j++)
                sum += (inp.w[i][j][k]*x_ilo[i][j]);

        sum += Omega*sqrt(SSsigma);

        model.add(sum <= inp.R[k]);
        sum.end();	
    }


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
