#include <string.h>
#include <ctype.h>
#include <algorithm>
#include <string>

#ifndef NUCSTN_H
#define NUCSTN_H

using namespace std;

const int MAXNUMBEROFNEIGHBORS = 100;
const int NUMBER_OF_PARAMETERS = 200;
const long int FACTORIALS[17] = {1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,6227020800,87178291200,1307674368000,20922789888000};
const double PI = 3.14159265;
const int nilvalue = 33333333;
#define MY_NEIGHMASK 0x3FFFFFFF

namespace NUCSTN{

class CMolecule {
  private:
    //The position of the particle
  public:
    CMolecule();
    virtual ~CMolecule();
    double posx,posy,posz;
    double potential;
    
    int id;
    int neighbors[MAXNUMBEROFNEIGHBORS];
    double neighbordist[MAXNUMBEROFNEIGHBORS];
    //JR: also save the distance in cartesian and spherical coordinates
    double n_diffx[MAXNUMBEROFNEIGHBORS];
    double n_diffy[MAXNUMBEROFNEIGHBORS];
    double n_diffz[MAXNUMBEROFNEIGHBORS];
    double n_r[MAXNUMBEROFNEIGHBORS];
    double n_phi[MAXNUMBEROFNEIGHBORS];
    double n_theta[MAXNUMBEROFNEIGHBORS];

    int n_neighbors;
    int cell;

    double realQ4[9],imgQ4[9];
    double realQ6[13],imgQ6[13];
    double arealQ4[9],aimgQ4[9];
    double arealQ6[13],aimgQ6[13];
    
    double Q6,Q4,W4,W6,AQ6,AQ4,AW4,AW6;
    double frenkelnumber;
    double avq6q6;
    //belongs to which cluster
    int belongsto;
    int issolid;
    int structure;
};

//change this to work without lammps neighbour classes
//now all my functions would go here
void read_ParticleFile(CMolecule *, string , int , double (&)[3]);
double get_absDistance(CMolecule *, int ,int , double(&)[3], double&, double&, double&);
void get_AllNeighborsAndDistances(CMolecule *, int,  double, double(&)[3] );
void convert_SphericalCoordinates(double , double , double , double &, double &, double &);
double PLM(int , int , double );
void YLM(int , int, double , double , double &, double &);
void QLM(int ,int ,double ,double,double &, double &);
void calculate_complexQLM_6(CMolecule *, int, double, double(&)[3] );
double get_NumberFromBond(CMolecule *, int,int );
void calculate_frenkelNumbers(CMolecule *, int, double );
int gather_clusters(CMolecule *, int , double , int );


//end of NUCSIZE namespace        
}

#endif