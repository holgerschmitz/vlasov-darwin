/***************************************************************************
                          potential.cpp  -  description
                             -------------------
    begin                : Fri Jun 2 2000
    copyright            : (C) 2000 by Holger Schmitz
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
// potential.cxx


#include "potential.h"
#include "diaghelper.h"

// ----------------------------------------------------------------------
// Potential

void Potential::AddSpecies(ESVlasovSpecies* pS) { 	
    species.push_back(WESVlasovSpecies(pS));
}

/// Get parameters Nx and dx and initialize the potential grid
void Potential::Init () {

    n0 = -1;

    Nx[0] = GlGridSize_x;
    Nx[1] = GlGridSize_y;
    dx[0] = GlGridSpace_x;
    dx[1] = GlGridSpace_y;

    cerr << "Grid Spacing is " << dx << endl;
    cerr << "Grid Size is " << Nx << endl;


    // resize grid
	den.resize(Nx.Data());
    Pot.resize(Nx.Data());
    Ex.resize(Nx.Data());
    Ex.clear();
    Ey.resize(Nx.Data());
    Ey.clear();

    int gssx = GlGridSize_x+1, gssy = GlGridSize_y+1;
    pois = new Poisson;
    pois->resize(
        PositionD(0.0,0.0), 
        PositionD(gssx*dx[0], gssy*dx[1]), 
        PositionI(gssx,gssy),
        Poisson::periodic, 
        Poisson::periodic
    );

    PositionI GrSize(gssx,gssy);

    In.resize(GrSize.Data());
    
    DiagField.Init(this);

    cerr << "Done Potential " << Nx << endl;
}

bool Potential::Execute (double timestep) {
    
	int i;
	double dF;
	ScalarField tmp;
//    PositionI GS(GlGridSize_x,GlGridSize_y);
	tmp.resize(Nx.Data());
	
	// initialise
	den.clear();

    // iterate through the species
	for (int s = species.size() - 1; s >= 0; s--) {
		ESVlasovSpecies* pS = species[s];

		dF = pS->GetCharge();
        
//        cerr << "Creating Density\n";
		pS->MakeRho();      //request to make particle density...
        
        // ... and get it
        for (int j=0; j<GlGridSize_y; ++j) 
            for (int i=0; i<GlGridSize_x; ++i) 
		        tmp(i,j) = dF*pS->Rho()(i,j);    
                
		den += tmp;         // and add
	}

    write_Scalar(den, "den.out");
    
    for (int j=0; j<GlGridSize_y; ++j) 
        for (int i=0; i<GlGridSize_x; ++i) { 
            In(i,j) = -(den(i,j)+n0);
//            std::cerr << "den " << i << " " << j << " " << In(i,j) << std::endl;
        }

//    write_Scalar(In, "In.out");
    
    pois->solve(In,Pot);
    
    write_Scalar(Pot, "Pot.out");

	for (int i=1; i<(GlGridSize_x-1); ++i) 
        for (int j=0; j<GlGridSize_y; ++j) 
    	    Ex(i,j) = (Pot(i-1,j) - Pot(i+1,j)) / (2*dx[0]);

	for (int i=0; i<GlGridSize_x; ++i) 
        for (int j=1; j<(GlGridSize_y-1); ++j) 
    	    Ey(i,j) = (Pot(i,j-1) - Pot(i,j+1)) / (2*dx[1]);


    // Perform wrapping
    for (int j=0; j<GlGridSize_y; ++j) {
        Ex(0,j) = Ex(GlGridSize_x-2,j);
	    Ex(GlGridSize_x-1,j) = Ex(1,j);
        
//          Ey(0,j) = Ey(GlGridSize_x-2,j);
//  	    Ey(GlGridSize_x-1,j) = Ey(1,j);
    }

    // Perform wrapping
    for (int i=0; i<GlGridSize_x; ++i) {
//          Ex(i,0) = Ex(i,GlGridSize_y-2);
//  	    Ex(i,GlGridSize_y-1) = Ex(i,1);
        
        Ey(i,0) = Ey(i,GlGridSize_y-2);
	    Ey(i,GlGridSize_y-1) = Ey(i,1);
    }

    write_Scalar(Ex, "Ex.out");
    write_Scalar(Ey, "Ey.out");

    DiagField.Execute();


//    Task::Execute();
  	return false;
    
}

void write_Scalar(ScalarField &Field, const char* fname) {
    ofstream O(fname);
    for (int i=0; i<=Field.getHigh(0); ++i) {
        for (int j=0; j<=Field.getHigh(1); ++j)
            O   << i << " "  << j << " "  
                << Field(i,j) << endl;
        O << endl;
    }
}

void write_Scalar(ScalarField &Field, ostream& O) {
    for (int i=0; i<=Field.getHigh(0); ++i) {
        for (int j=0; j<=Field.getHigh(1); ++j)
            O   << i << " "  << j << " "  
                << Field(i,j) << endl;
        O << endl;
    }
}



