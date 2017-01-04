/*
 * SampleDecoder.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */
#include <iostream>
#include <iomanip>
#include <sstream>
#include "SampleDecoder.h"
#include <stdio.h>
#include <stdlib.h>
#include <RInside.h>
#include <cassert>

SampleDecoder::SampleDecoder()  { }
SampleDecoder::~SampleDecoder() { }

using namespace std;

double best_score = 999999.9;

extern int nStudents;
extern int nPeriods;
extern int nClusters;
extern int nSeats;
extern int nAttributes;
int natPos    = 9;
int genderPos = 0;
extern int nWG;
extern int  ** attr;   //!< attributes of students (read from disk fil
extern int   * ub;	   //!< number of seats of each type (E, B, F)
extern int  ** clusters;    //!< list of seat of each type (E, B, F)
extern int  ** matrixChart; //!< seating chart using a matrix structure
extern int  ** s;	   //!< seating chart using a matrix structure
extern int  ** y;	   //!< seat assigned in period t to student i
extern int ** s_best_brkga;	//!< best assignment found by brkga
extern int ** y_best_brkga;	//!< best assignment found by brkga

extern int nrowsChart;		//!< seating chart info : nr rows
extern int ncolsChart;		//!< seating chart info : nr cols
extern int *** listCluster; //!< list of students aassigned to each cluster (in each period)
extern int   * totInCluster;
extern int withForbidden; //!< defines whether a forbidden matrix is given
extern int  ** F;	  //!< incompatibilities matrix

extern std::vector < std::vector < int> > N_setReduced;
extern std::vector < std::vector < int> > nationalities;
extern std::vector < std::vector < std::vector < int > > > workgroups;
extern std::map < int, std::vector <int> > genderMap;

double compute_fitness(int ** s, int **y, int ** x, std::vector < std::vector < int> > & N_set, int ** counting, int ** countingReduced);

void update_best(int ** s, int ** y, double score);

double SampleDecoder::decode(const std::vector< double >& chromosome, RInside & R, int ** x, std::vector < std::vector < int> > & N_set, int ** counting, int ** countingReduced) const 
{

#ifdef AAA
  cout << "================================================" << endl;
  cout <<"Chromosome is :: ";
  for (unsigned j = 0; j < chromosome.size(); j++)
    cout << "    " << chromosome[j];
  cout << endl;
  cout << "================================================" << endl;
  cout << "Current Cluster assignment INSIDE ;; " << endl;
  for (int i = 0; i < nStudents; i++)
    {
      for (int t = 0; t < nPeriods; t++)
	cout << "stud(" << i << "," << t << ") = " << x[i][t];
      cout << endl;
    }

#endif

  int ** used = new (int*[nClusters]);   
  for (int c = 0; c < nClusters; c++)
    {
      used[c] = new int[ub[c]];
      for (int j = 0; j < ub[c]; j++)
	used[c][j] = 0;
    }
  
  // encoding (chromosome of length "n"):
  for (int t = 0; t < nPeriods; t++)
    {
      for (int i = 0; i < nSeats; i++)
	s[t][i] = -1;

      for (int j = 0; j < nStudents; j++)
	  y[t][j] = -1;

      // assigning students per cluster
      int progrChrom = 0;
      for (int c = 0; c < nClusters; c++)
	{
#ifdef M_DEBUG
	  cout << "CLUSTER " << c << endl;
	  cout << "UB[" << c << "] = " << ub[c] << endl;
#endif
	  for (int k = 0; k < totInCluster[c]; k++)
	    {
	      int stud_i = listCluster[t][c][k];	 
	      if (x[stud_i][t] != c) cout << " ## ## ## ## ERROR c1 " << endl; 
	  
#ifdef M_DEBUG
	      cout << "student " << stud_i << " in cluster " << c << " in period " << t << endl;
	      cout << "student " << stud_i << " random number r :: " << chromosome[progrChrom] << endl;
#endif
	  
	      //int posR = floor(chromosome[progrChrom]*(ub[c]-k+1));
	      int posR = floor(chromosome[progrChrom]*(ub[c]-k));
	      //	    cout << "  relative position is " << posR << endl;
	      assert(posR >= 0 && posR < ub[c]);
	      int progr = -1;
	      int pos;
	      for (pos = 0; pos < ub[c]; pos++)
		{
		  if (used[c][pos] == 0)
		    progr++;
		  if (progr == posR)
		    break;
		}
	      // 
	      int seat_i = clusters[c][pos];
	      //cout << "student " << stud_i << " is assigned to seat " << seat_i << endl;
	      used[c][pos] = 1;
	      if (s[t][seat_i] != -1) 
		{
		  cout << "error here -----------------------------------" << endl;
		  exit(145);
		}
			
	      s[t][seat_i] = stud_i;
	      y[t][stud_i] = seat_i;
	      progrChrom++;
	    }
#ifdef M_DEBUG
	  int aka;
	  cin >> aka;
#endif

	}

#ifdef M_DEBUG
      cout << "Assignment for period " << t << " ########## " << endl;
      cout << endl << endl;
      cout << "=========================" << endl;
      cout << "** ** Seating chart ** **" << endl;
      cout << "=========================" << endl;
      cout << endl;
      for (int t = 0; t < nPeriods; t++)
	{
	  cout << endl << endl << " ------------------------------PERIOD " << t <<"-------------------------------" << endl;
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
#endif

      // reset data structure for relative position
      for (int c = 0; c < nClusters; c++)
	for (int j = 0; j < ub[c]; j++)
	  used[c][j] = 0;
      
    }

  // compute objective function using the same criterion used in cplex
  // this allows a fair comparison
    double score = compute_fitness(s, y, x, N_set, counting, countingReduced);
  
    if (score < best_score)
	update_best(s,y, score);

  // REM. : It is minimizing...
  //double score = 1.0 - fit;
  //double score = fit;
   
  //cout << "################ score of this chromosome is :: " << score << endl;
  return score;
}

void update_best(int ** s, int ** y, double score)
{
    best_score = score;
    for (int t = 0; t < nPeriods; t++)
    {
	for (int i = 0; i < nSeats; i++)
	    s_best_brkga[t][i] = s[t][i];

	for (int j = 0; j < nStudents; j++)
	    y_best_brkga[t][j] = y[t][j];
    }
}

/// Compute fitness value of a given assigment
double compute_fitness(int ** s, int ** y, int ** x, std::vector < std::vector < int> > & N_set, int ** counting, int ** countingReduced)
{
#ifdef M_DEBUG
  // check that each student is assigned to one seat only
  int * studDist = new int[nStudents];
  for (int t = 0; t < nPeriods; t++)
  {
       for (int i = 0; i < nStudents; i++)
	    studDist[i] = -1;

       for (int j = 0; j < nSeats; j++)
       {
	    int stud = s[t][j];

	    if (stud == -1) continue; // empty seat

	    assert(y[t][stud] == j);

	    if (studDist[stud] != -1)
	    {
		 cout << "Error double assignment with student " << stud << endl;
		 exit(124);
	    }
	    studDist[stud] = j;
       }

       // check that each student has an assignment
       for (int i = 0; i < nStudents; i++)
	    if (studDist[i] == -1)
	    {
		 cout << "Error student " << i << "not assigned " << endl;
		 exit(125);
	    }
  }
       
#endif

  // check incompatibilities among students (forbidden matrix)
  int n_forbidden = 0;
  if (withForbidden)
  {
      for (int i = 0; i < nStudents; i++)
      {
	  for (int j = i+1; j < nStudents; j++)
	  {
	      if (F[i][j] == 1)
	      {
		  //cout << "students " << i << " and " << j << " are forbidden " << endl;
		  for (int t = 0; t < nPeriods; t++)
		  {
		      int seat_i = y[t][i]; // get the seat of student i

		      for (unsigned p = 0; p < N_set[seat_i].size(); p++)
		      {
			  int seat_neigh = N_set[seat_i][p];
			  if (s[t][seat_neigh] == j)
			  {
#ifdef M_DEBUG
			      // violation of incompatibility
			      cout << "Students " << i << " and " << j
				   << " are violating incompatibility in"
				   << " period " << t << " :: " << endl;
			      cout << " *** " << i << " is in " << seat_i << endl;
			      cout << " *** " << j << " is in " << seat_neigh << endl;
#endif
			      
			      n_forbidden++;
			  }	
		      }
		  }
	      }
	  }
      }
  }

   int infeasNeigh  = 0;
   int infeasNat    = 0;
   int infeasGender = 0;
   int infeasWG     = 0;

   // verify that there is no repetition of neighborhs
   //int ** counting = new (int*[nStudents]);
   for (int i = 0; i < nStudents; i++)
   {
       //counting[i] = new int[nStudents];
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
#ifdef M_DEBUG
	    cout << "WARNING :: Students " << i << " and " << j << " repeat neighborhood " << endl;
#endif
	 }
 			  

   // count number of infeasibilities due to nationalities
   for (int i = 0; i < nStudents-1; i++)
       for (unsigned j = 0; j < nationalities[i].size(); j++)
       {
	   int student_j = nationalities[i][j];
	   if (counting[i][student_j] > 0)
	   {
	       infeasNat++;
#ifdef M_DEBUG
	       cout << "WARNING :: Students " << i << " and " << student_j << " are neighbors and have same nationality (" << attr[i][natPos] << ") " << endl; 
	       assert(attr[i][natPos] == attr[student_j][natPos]);
#endif
	   }
       }
   
   //int ** countingReduced = new (int*[nStudents]);
   for (int i = 0; i < nStudents; i++)
   {
       //countingReduced[i] = new int[nStudents];
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
#ifdef M_DEBUG
	    cout << "WARNING :: Students " << student_i << " and " << student_j << " are neighbors and women (" << attr[student_i][genderPos] << ") " << endl; 
#endif
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
#ifdef M_DEBUG
		  cout << "WARNING :: Students " << i << " and " << student_j << " are neighbors and same WG (" << attr[i][wgPos] << ") in period " << t << endl; 
#endif
		  //assert(attr[i][wgPos] == attr[student_j][wgPos]);
		  break;
	       }
	    }
	 }
      }
   }

#ifdef M_DEUBG
   cout << " ====================" << endl;
   cout << "| Summary Statistics |" << endl;
   cout << " ====================" << endl;
   cout << " *  Nr. Neighbors Infeasibilities \t :: " << infeasNeigh << endl;
   cout << " *  Nr. Nationality Infeasibilities \t :: " << infeasNat << endl;
   cout << " *  Nr. Gender Infeasibilities \t\t :: " << infeasGender << endl;
   cout << " *  Nr. WG Infeasibilities \t\t :: " << infeasWG << endl;
#endif
   
   double score = 2.0*(double)n_forbidden + (double)infeasNeigh + 0.1*(double)infeasNat + (double)infeasGender + 0.1*(double)infeasWG;

   return score;
}
