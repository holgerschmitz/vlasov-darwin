// -*- C++ -*-
// $Id$


#include "potential.h"
#include "diaghelper.h"

// ----------------------------------------------------------------------
// Potential

void Potential::AddSpecies(EFieldForce* pS) { 	
    species.push_back(pEFieldForce(pS));
}

/// Get parameters Nx and dx and initialize the potential grid
void Potential::Init () {

//    n0 = -1;
    n0 = Parameters::instance().bgDensity();

    LBound = Process::instance().getBoundary().scalarLow();
    HBound = Process::instance().getBoundary().scalarHigh();

    dx[0] = Parameters::instance().gridSpace_x();
    dx[1] = Parameters::instance().gridSpace_y();

    std::cout << "Grid Spacing is " << dx << std::endl;
    std::cout << "Grid Size is " << LBound << " to " << HBound << std::endl;


    // resize grid
    den.resize(LBound.Data(),HBound.Data());
    Pot.resize(LBound.Data(),HBound.Data());
    Ex.resize(LBound.Data(),HBound.Data());
    Ex.clear();
    Ey.resize(LBound.Data(),HBound.Data());
    Ey.clear();

    int gssx = HBound[0]-LBound[0], gssy = HBound[1]-LBound[1];
    pois = new Poisson;
    pois->resize(
        PositionD(0.0,0.0), 
        PositionD(gssx*dx[0], gssy*dx[1]), 
        PositionI(gssx,gssy),
        Poisson::periodic, 
        Poisson::periodic
    );

    In.resize(LBound.Data(),HBound.Data());
    
    DiagField.Init(this);

}

void Potential::Execute () {
    
    int lx0 = LBound[0], lx1 = LBound[0]+1;
    int ly0 = LBound[1], ly1 = LBound[1]+1;
    int mx0 = HBound[0], mx1 = HBound[0]-1;
    int my0 = HBound[1], my1 = HBound[1]-1;

	int i;
	double dF;
	ScalarField tmp;

	tmp.resize(LBound.Data(),HBound.Data());
	
	// initialise
	den.clear();

    // iterate through the species
	for (int s = species.size() - 1; s >= 0; s--) {
		EFieldForce *pS = species[s];

		dF = pS->getCharge();
        
        DistMomentRho *distRho = pS->getDerivedRho();
        
        // ... and get it
        ScalarField &rho =distRho->getRho();
        for (int j=ly0; j<=my0; ++j) 
            for (int i=lx0; i<=mx0; ++i) 
		        tmp(i,j) = dF*rho(i,j);    
                
		den += tmp;         // and add
	}
        
    for (int j=ly0; j<=my0; ++j) 
        for (int i=lx0; i<=mx0; ++i) {
            In(i,j) = -(den(i,j)+n0);
//            std::cerr << "den " << i << " " << j << " " << In(i,j) << std::endl;
        }

    pois->solve(In,Pot);
    
	for (int i=lx1; i<=mx1; ++i) 
        for (int j=ly0; j<=my0; ++j) 
    	    Ex(i,j) = (Pot(i-1,j) - Pot(i+1,j)) / (2*dx[0]);

	for (int i=lx0; i<=mx0; ++i) 
        for (int j=ly1; j<=my1; ++j) 
    	    Ey(i,j) = (Pot(i,j-1) - Pot(i,j+1)) / (2*dx[1]);


    // Perform wrapping
    for (int j=ly0; j<=my0; ++j) {
        Ex(lx0,j) = Ex(mx1,j);
        Ex(mx0,j) = Ex(lx1,j);
        
    }

    // Perform wrapping
    for (int i=lx0; i<=mx0; ++i) {
        Ey(i,ly0) = Ey(i,my1);
        Ey(i,my0) = Ey(i,ly1);
    }

    DiagField.Execute();    
}

ScalarField &Potential::GetByName(const std::string& name) {
  if ("Pot"==name)
    return Pot; 
  else  if ("Ex"==name)
    return Ex; 
  else
    return Ey;
}

void write_Scalar(ScalarField &Field, const char* fname, double offset) {
    ofstream O(fname);
    for (int i=Field.getLow(0); i<=Field.getHigh(0); ++i) {
        for (int j=Field.getLow(1); j<=Field.getHigh(1); ++j)
            O   << i << " "  << j << " "  
                << Field(i,j)-offset << endl;
        O << endl;
    }
}

void write_Scalar(ScalarField &Field, ostream& O, double offset) {
    for (int i=Field.getLow(0); i<=Field.getHigh(0); ++i) {
        for (int j=Field.getLow(1); j<=Field.getHigh(1); ++j)
            O   << i << " "  << j << " "  
                << Field(i,j)-offset << endl;
        O << endl;
    }
}



