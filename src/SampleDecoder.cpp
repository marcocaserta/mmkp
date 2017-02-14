/*
 * SampleDecoder.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 *      Adapted : marco caserta (for the mmkp problem)
 */



#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <fstream>

#include <sstream>
#include "SampleDecoder.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cassert>

SampleDecoder::SampleDecoder()  { }
SampleDecoder::~SampleDecoder() { }

using namespace std;


/// FUNCTIONS DEFINITION ==========================
/// END FUNCTIONS DEFINITION ==========================

double select_corridor_width_base(double r)
{
   if (r < 0.20)
        return 0.10;
    else if(r < 0.40)
        return 0.20;
    else if(r < 0.50)
        return 0.50;
    else if(r < 0.80)
        return 0.75;
   else
        return 0.9;
/*
     if (r < 0.25)
        return 0.75;
    else if(r < 0.50)
        return 0.85;
    else if(r < 0.75)
        return 0.90;
    else
        return 0.95;
   if (r < 0.25)
        return 0.25;
    else if(r < 0.50)
        return 0.50;
    else if(r < 0.75)
        return 0.75;
    else
        return 0.9;*/
}

int select_nSol_base(double r)
{
 /*    if (r < 0.25)
        return 3;
    else if(r < 0.50)
        return 5;
    else if(r < 0.75)
        return 7;
    else 
        return 10;*/
   if (r < 0.25)
        return 10;
    else if(r < 0.50)
        return 12;
    else if(r < 0.75)
        return 14;
    else 
        return 16;
}

double select_zero_base(double r)
{
    if (r < 0.25)
        return 0.6;
    else if(r < 0.50)
        return 0.7;
    else if (r < 0.75)
        return 0.8;
    else
        return 0.85;

 /*    if (r < 0.25)
        return 0.4;
    else if(r < 0.50)
        return 0.6;
    else if (r < 0.75)
        return 0.7;
    else
    return 0.8;*/
    return 0.8;
}

double select_one_base(double r)
{

    if (r < 0.20)
        return 0.0;
    else if(r < 0.40)
        return 0.05;
    else if (r < 0.60)
        return 0.10;
    else if (r < 0.8)
        return 0.15;
    else
        return 0.25;
}

int select_cut(double r)
{
    return 0;

    if (r < 0.5)
        return 1;
    else
        return 0;
}

double SampleDecoder::decode(const std::vector< double >& chromosome) const 
{
    // decoding (chromosome of length "n"):
    double corridorWidthBase = select_corridor_width_base(chromosome[0]);
    int nSolBase             = select_nSol_base(chromosome[1]);
    double propFixed0        = select_zero_base(chromosome[2]);
    double propFixed1        = select_one_base(chromosome[3]);
    int add_z_cut            = select_cut(chromosome[4]);


    string sBase = "./bin/mmkp -f ../data/I13.txt ";
    stringstream s1;
    s1 << corridorWidthBase;
    stringstream s2;
    s2 << nSolBase;
    stringstream s3;
    s3 << propFixed0;
    stringstream s4;
    s4 << propFixed1;
    stringstream s5;
    s5 << add_z_cut;

    string commLine =  sBase + " -c " + s1.str() + " -n " + s2.str() 
        + " -z " + s3.str() + " -u " + s4.str() + " -a " + s5.str();

    const char * cLine = commLine.c_str();
    double outCL = system(cLine);

    string temp;
    double score;

    ifstream fsol("solution.txt", ios::in);
    fsol >> temp >> score;
    fsol.close();


    // REM. : It is minimizing...
    cout << "################ score of this chromosome is " << score << endl;

    double result = system("cat solution.txt >> summary.txt");

    return -score;
}

