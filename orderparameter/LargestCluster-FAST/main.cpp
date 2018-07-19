///////////////////////////////////////////////////////////////////////////////////////////
//
//
// Structure Analysis Code v 1.0
// by Wolfgang Lechner, Amsterdam 2010
//
// wolfgang.lechner@gmail.com
//
//
//
// This source code is licensed under the Academic Free License
// Copyright by Wolfgang Lechner, 2010
// In addition to that I would like you to cite
// Wolfgang Lechner and Christoph Dellago, J. Chem. Phys. 129, 114707 (2008)
// and send me a short email about your project when using this code.
//
// The package shows how to use the averaged version of the local bond order parameters
// in order to determine the local structure of points or particles distributed in
// three dimensional space described in the Paper mentioned above.
// This is useful for a wide range of applications including particle simulations to
// investigate nucleation, determine the crystal structure in metals, and the study of
// glasses.
//
// This package was originally part of a particle simulation to study the nucleation
// of colloidal particles. If you are interested in a collaboration on this field feel
// free to email me about your interests.
// Also, if you have any questions or found a bug do not hesitate writing me!
//
//
// modified version of the code 
// by Jutta Rogal, New York 2017
//
// jutta.rogal@rub.de
//
// reading different formats of input files and changing the output
//
///////////////////////////////////////////////////////////////////////////////////////////

#include "molecular_system.h"
#include "molecule.h"
#include "parameter.h"
#include <ctime>
#include <sstream>
#include <string>
#include <cstdlib>

int main(int argc, char *argv[])
{
	// The program reads a lammps dump file as input
	// parameters are hard coded in parameter.cpp
	string filename="";
	if (argc > 1)
	{
		filename = string(argv[1]);
		//t_parameter = atoi(argv[2]);
	}
	else{
		cout<<"Must provide input file name to program."<<endl;
		cout<<"Exiting program..."<<endl;
		exit(1);
	}
	// The main obejct is created. It hold all the functions and data used in the analysis.
	CMolecularSystem *m_MolSys = new CMolecularSystem;
	// The parameterfile is read, currently parameter are NOT READ BUT SET TO A FIXED VALUE! (parameter.cpp)
	m_MolSys->parameter->readParameter();
	
	// System is initalized, memory allocated, read in particle positions and box
	m_MolSys->readParticleFile(filename);   


	// calls Q6 and frenkelnumber within function and gets largest Cluster
	m_MolSys->calculate_largestClusterparameter_Full();

	//Free the memory.
	m_MolSys->deleteMolecules();

	return 0;
}
