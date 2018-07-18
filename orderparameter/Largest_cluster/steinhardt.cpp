#include "nucstn.h"
#include <iostream>
#include <cmath>
//--------------------------------------------------------------------------------------------
double NUCSTN::get_absDistance(CMolecule *molecules, int ti ,int tj, double *box,double &diffx, double &diffy, double &diffz)
{
  double abs,boxx,boxy,boxz;


  boxx = box[0];
  boxy = box[1];
  boxz = box[2];

  diffx = molecules[tj].posx - molecules[ti].posx;
  diffy = molecules[tj].posy - molecules[ti].posy;
  diffz = molecules[tj].posz - molecules[ti].posz;
  //nearest image
  if (diffx >  boxx/2.0) {diffx = diffx - boxx;};
  if (diffx < -boxx/2.0) {diffx = diffx + boxx;};
  if (diffy >  boxy/2.0) {diffy = diffy - boxy;};
  if (diffy < -boxy/2.0) {diffy = diffy + boxy;};
  if (diffz >  boxz/2.0) {diffz = diffz - boxz;};
  if (diffz < -boxz/2.0) {diffz = diffz + boxz;};
  
  abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
  return abs;
}


void NUCSTN::get_AllNeighborsAndDistances(CMolecule *molecules, int nmax, double cutoff, double *box )
{

        double nd,d;
        double diffx,diffy,diffz;
        double r,theta,phi;
        nd = cutoff;
        int tii,tjj,ti,tj;
        int nop;
        nop = nmax;
        int ccount[nop];
        int jnum;   
        for (ti = 0;ti<nop;ti++)
        {
                ccount[ti]=0;
                //initially set the cluster to their ids
                molecules[ti].belongsto = molecules[ti].id;

                for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++)
                {
                        molecules[ti].neighbors[tn] = nilvalue;
                        molecules[ti].neighbordist[tn] = -1.0;
                }
        }

        for (ti=0; ti<nop; ti++)
        {
                
                for (tj=ti+1; tj<nop; tj++)
                {
                        
                        d = get_absDistance(molecules,ti,tj,box,diffx,diffy,diffz); 
                        if (d < nd) {
                                molecules[ti].neighbors[ccount[ti]] = tj; 
                                molecules[ti].neighbordist[ccount[ti]] = d; 
                                molecules[ti].n_diffx[ccount[ti]] = diffx;
                                molecules[ti].n_diffy[ccount[ti]] = diffy;
                                molecules[ti].n_diffz[ccount[ti]] = diffz;
                                convert_SphericalCoordinates(diffx, diffy, diffz, r, phi, theta);
                                molecules[ti].n_r[ccount[ti]] = r;
                                molecules[ti].n_phi[ccount[ti]] = phi;
                                molecules[ti].n_theta[ccount[ti]] = theta;
                                molecules[tj].neighbors[ccount[tj]] = ti;
                                molecules[tj].neighbordist[ccount[tj]] = d;
                                molecules[tj].n_diffx[ccount[tj]] = -diffx;
                                molecules[tj].n_diffy[ccount[tj]] = -diffy;
                                molecules[tj].n_diffz[ccount[tj]] = -diffz;
                                convert_SphericalCoordinates(-diffx, -diffy, -diffz, r, phi, theta);
                                molecules[tj].n_r[ccount[tj]] = r;
                                molecules[tj].n_phi[ccount[tj]] = phi;
                                molecules[tj].n_theta[ccount[tj]] = theta;
                                ccount[tj] +=1;
                        }
                }
        }
        for (int ti=0; ti<nop; ti++){
                molecules[ti].n_neighbors = ccount[ti];
        }

}


void NUCSTN::convert_SphericalCoordinates(double x, double y, double z, double &r, double &phi, double &theta)
{
        r = sqrt(x*x+y*y+z*z);
        theta = acos(z/r);
        phi = atan2(y,x);
}

double NUCSTN::PLM(int l, int m, double x)
{
        double fact,pll,pmm,pmmp1,somx2;
        int i,ll;
        pll = 0.0;
//        if (m < 0 || m > l || fabs(x) > 1.0)
//                cerr << "impossible combination of l and m" << "\n";
        pmm=1.0;
        if (m > 0){
                somx2=sqrt((1.0-x)*(1.0+x));
                fact=1.0;
                for (i=1;i<=m;i++){
                        pmm *= -fact*somx2;
                        fact += 2.0;
                }
        }
  
        if (l == m)
                return pmm;
        else{
                pmmp1=x*(2*m+1)*pmm;
                if (l == (m+1))
                        return pmmp1;
                else{
                        for (ll=m+2;ll<=l;ll++){
                                pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
                                pmm=pmmp1;
                                pmmp1=pll;
                        }
                        return pll;
                }
        }
}


void NUCSTN::YLM(int l, int m, double theta, double phi, double &realYLM, double &imgYLM)
{
        double factor;
        double m_PLM;
        m_PLM = PLM(l,m,cos(theta));
        factor = ((2.0*double(l) + 1.0)*FACTORIALS[l-m]) / (4.0*PI*FACTORIALS[l+m]);
        realYLM = sqrt(factor) * m_PLM * cos(double(m)*phi);
        imgYLM  = sqrt(factor) * m_PLM * sin(double(m)*phi);
}


void NUCSTN::QLM(int l,int m,double theta,double phi,double &realYLM, double &imgYLM )
{
        realYLM = 0.0;
        imgYLM = 0.0;
        if (m < 0) {
                YLM(l, abs(m), theta, phi, realYLM, imgYLM);
                realYLM = pow(-1.0,m)*realYLM;
                imgYLM = pow(-1.0,m)*imgYLM;
        }
        else
        {
                YLM(l, m, theta, phi, realYLM, imgYLM);
        }
}


void NUCSTN::calculate_complexQLM_6(CMolecule *molecules, int nmax, double cutoff, double *box)
{
        //nn = number of neighbors
        int nn,nop;
        double realti,imgti;
        double realYLM,imgYLM;
        nop = nmax;
        get_AllNeighborsAndDistances(molecules, nmax, cutoff, box );
        for (int ti= 0;ti<nop;ti++)
        {
                nn = molecules[ti].n_neighbors;
                for (int mi = -6;mi < 7;mi++)
                {
                        realti = 0.0;
                        imgti = 0.0;
                        for (int ci = 0;ci<nn;ci++)
                        {
                                QLM(6,mi,molecules[ti].n_theta[ci],molecules[ti].n_phi[ci],realYLM, imgYLM);
                                realti += realYLM;
                                imgti += imgYLM;
                        }
                        realti = realti/(double(nn));
                        imgti = imgti/(double(nn));
                        molecules[ti].realQ6[mi+6] = realti;
                        molecules[ti].imgQ6[mi+6] = imgti;
                }
        }

}


//before this function, data needs to be exchanged
double NUCSTN::get_NumberFromBond(CMolecule *molecules, int ti,int tj)
{
        double sumSquareti,sumSquaretj;
        double realdotproduct,imgdotproduct;
        double connection;
        sumSquareti = 0.0;
        sumSquaretj = 0.0;
        realdotproduct = 0.0;
        imgdotproduct = 0.0;

        for (int mi = 0;mi < 13;mi++)
        {
                sumSquareti += molecules[ti].realQ6[mi]*molecules[ti].realQ6[mi] + molecules[ti].imgQ6[mi] *molecules[ti].imgQ6[mi];
                sumSquaretj += molecules[tj].realQ6[mi]*molecules[tj].realQ6[mi] + molecules[tj].imgQ6[mi] *molecules[tj].imgQ6[mi];
                realdotproduct += molecules[ti].realQ6[mi]*molecules[tj].realQ6[mi];
                imgdotproduct  += molecules[ti].imgQ6[mi] *molecules[tj].imgQ6[mi];
        }
        connection = (realdotproduct+imgdotproduct)/(sqrt(sumSquaretj)*sqrt(sumSquareti));
        return connection;
}


void NUCSTN::calculate_frenkelNumbers(CMolecule *molecules, int nmax, double threshold)
{
        int frenkelcons;
        double scalar;
        int nop = nmax;
        for (int ti= 0;ti<nop;ti++)
        {
                frenkelcons = 0;
                molecules[ti].avq6q6 = 0.0;
                for (int c = 0;c<molecules[ti].n_neighbors;c++)
                {
                        scalar = get_NumberFromBond(molecules,ti,molecules[ti].neighbors[c]);
                        if (scalar > threshold) frenkelcons += 1;
                        molecules[ti].avq6q6 += scalar;
                }
                molecules[ti].frenkelnumber = frenkelcons;
                molecules[ti].avq6q6 /= molecules[ti].n_neighbors;
                //cout<<this->molecules[ti].frenkelnumber<<""<<"\n";
        } 
}


int NUCSTN::gather_clusters(CMolecule *molecules, int nmax, double avgthreshold, int maxcons)
{
    //int* clusterid = new int[Size]
    int ti,tj;
    int change,done;

    //this loop runs recursilvely until all ids are assigned
    while(1){
        done = 1;
        for(ti=0;ti<nmax;ti++){
            for(tj=0;tj<molecules[ti].n_neighbors;tj++){
                if (molecules[ti].belongsto == molecules[tj].belongsto) continue;
                if ((molecules[ti].avq6q6 > avgthreshold) && ( molecules[ti].frenkelnumber > maxcons) ){
                    if ((molecules[tj].avq6q6 > avgthreshold) && ( molecules[tj].frenkelnumber > maxcons) ){
                        molecules[ti].belongsto = molecules[tj].belongsto = min(molecules[ti].belongsto,molecules[tj].belongsto);
                        done = 0;
                    }       
                }
            }
        }
        if (done) break;
    }
//at this point the clusters should have ids
//now make a new array with the nmax size - and count clusters
    int* clusterfreq = new int[nmax];
    int dummy_id,max;

    for(ti=0;ti<nmax;ti++){
        clusterfreq[ti] = 0;
    }

    for(ti=0;ti<nmax;ti++){
        //ids start from one, so maybe we do id -1
        dummy_id = molecules[ti].id -1;
        clusterfreq[dummy_id] += 1;    
    }

    //the freqs are recorded - find the maximum one
    max = 0;
    for(ti=0;ti<nmax;ti++){
        if (clusterfreq[ti] > max) max=clusterfreq[ti];    
    }

    delete [] clusterfreq;
    return max;

}


NUCSTN::CMolecule::CMolecule()
{

}

NUCSTN::CMolecule::~CMolecule()
{
    
}

//
//----------------------------------------------------------------------------------------------------