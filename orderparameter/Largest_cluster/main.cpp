#include <iostream>
#include "nucstn.h"
#include <sstream>
#include <string>

int main(int argc, char *argv[]){

	string filename="";
	int natoms=0;
	int csize,maxcons;
	double cutoff,threshold,avgthreshold;

	threshold = 0.6;
	avgthreshold = 0.5;
	cutoff = 3.27;
	maxcons = 7;

	if (argc > 2){
    
    	filename = string(argv[1]);
    	natoms = atoi(argv[2]);
	}

	//just testing so enter the number of atoms on the terminal
	NUCSTN::CMolecule *molecules = new NUCSTN::CMolecule[natoms];
  	double box[3];
  	//I should check at this point if this works
  	NUCSTN::read_ParticleFile(molecules,filename,natoms,box);
  	NUCSTN::calculate_complexQLM_6(molecules, natoms, cutoff, box);
  	NUCSTN::calculate_frenkelNumbers(molecules, natoms, threshold);
  	csize = NUCSTN::gather_clusters(molecules,natoms,avgthreshold,maxcons);

  	cout<<csize<<endl;

  	return 0;


}
