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
#include "darwin.h"

// ----------------------------------------------------------------------
// Potential

void Darwin::AddSpecies(DarwinVlasovSpecies* pS) { 	
    species.push_back(WDarwinVlasovSpecies(pS));
}

/// Get parameters Nx and dx and initialize the potential grid
void Darwin::Init () {

    n0 = -1;
    csc = 0.1;
    
    Nx[0] = GlGridSize_x;
    Nx[1] = GlGridSize_y;
    dx[0] = GlGridSpace_x;
    dx[1] = GlGridSpace_y;

	dV = GlVolumeQuant;

    cerr << "Grid Spacing is " << dx << endl;
    cerr << "Grid Size is " << Nx << endl;


    // resize grid
	den.resize(Nx.Data());
	om2.resize(Nx.Data());
    
	jx.resize(Nx.Data());
	jy.resize(Nx.Data());
	jz.resize(Nx.Data());

	sx.resize(Nx.Data());
	sy.resize(Nx.Data());
	sz.resize(Nx.Data());

	vxx.resize(Nx.Data());
	vxy.resize(Nx.Data());
	vxz.resize(Nx.Data());
	vyy.resize(Nx.Data());
	vyz.resize(Nx.Data());
	vzz.resize(Nx.Data());

    Pot.resize(Nx.Data());
    Az.resize(Nx.Data());

    Ex.resize(Nx.Data());
    Ey.resize(Nx.Data());
    Ez.resize(Nx.Data());
    
    Bx.resize(Nx.Data());
    By.resize(Nx.Data());
    Bz.resize(Nx.Data());
    
	den.clear();
	om2.clear();
    
	jx.clear();
	jy.clear();
	jz.clear();

	sx.clear();
	sy.clear();
	sz.clear();

    Pot.clear();
    Az.clear();

    Ex.clear();
    Ey.clear();
    Ez.clear();
    
    Bx.clear();
    By.clear();
    Bz.clear();

    int gssx = GlGridSize_x+1, gssy = GlGridSize_y+1;
    pois = new Poisson;
    pois->resize(
        PositionD(0.0,0.0), 
        PositionD(gssx*GlGridSpace_x, gssy*GlGridSpace_y), 
        PositionI(gssx,gssy),
        Poisson::periodic, 
        Poisson::periodic
    );

    helmh = new Helmholtz();
    
    PositionI GrSize(gssx,gssy);

    In.resize(GrSize.Data());
    Lambda.resize(GrSize.Data());
    Out.resize(GrSize.Data());
    
    cerr << "Done Darwin " << Nx << endl;
}

bool Darwin::Execute (double timestep) {
    
	int i;
	double dF, dF2;
    
	ScalarField tmp;
    
    PositionI GS(GlGridSize_x,GlGridSize_y);
	tmp.resize(GS.Data());
	
//    cerr << "Initialising Density\n";
	// initialise
    
    /* *************************************
     *
     */

    /* *************************************
     *  First ask the species to create the density and the current density
     *  and add them all up.
     */
    
	den.clear();
	om2.clear();
    
	jx.clear();
	jy.clear();
	jz.clear();

	sx.clear();
	sy.clear();
	sz.clear();

    vxx.clear();
	vxy.clear();
	vxz.clear();
	vyy.clear();
	vyz.clear();
	vzz.clear();
	
    VelocityD jt;
    FixedArray<double,6> vvt;
    
    // iterate through the species
	for (int s = species.size() - 1; s >= 0; s--) {
		DarwinVlasovSpecies* pS = species[s];

		dF = pS->GetCharge()/dV;
        dF2 = dF * pS->GetCharge();
        
//        cerr << "Creating Density\n";
		pS->MakeRho();          //request to make particle density.
        
        // ... and get it
        for (int j=0; j<GlGridSize_y; ++j) 
            for (int i=0; i<GlGridSize_x; ++i) {
		        den(i,j) += dF*pS->Rho()(i,j);
                om2(i,j) += dF2*pS->Rho()(i,j);
                
                jt = pS->getJ(i,j);
                jx(i,j) += dF*jt[0];
                jy(i,j) += dF*jt[1];
                jz(i,j) += dF*jt[2];
                
                sx(i,j) += dF2*jt[0];
                sy(i,j) += dF2*jt[1];
                sz(i,j) += dF2*jt[2];
                
                vvt =  pS->getVVTens(i,j);
                vxx(i,j) +=  dF*vvt[0];
                vxy(i,j) +=  dF*vvt[1];
                vxz(i,j) +=  dF*vvt[2];
                vyy(i,j) +=  dF*vvt[3];
                vyz(i,j) +=  dF*vvt[4];
                vzz(i,j) +=  dF*vvt[5];
        }
    }

    write_Scalar(den, "den.out");
    write_Scalar(om2, "om2.out");
    write_Scalar(jx, "jx.out");
    write_Scalar(jy, "jy.out");
    write_Scalar(jz, "jz.out");
    write_Scalar(sx, "sx.out");
    write_Scalar(sy, "sy.out");
    write_Scalar(sz, "sz.out");
    write_Scalar(vxx, "vxx.out");
    write_Scalar(vxy, "vxy.out");
    write_Scalar(vxz, "vxz.out");
    write_Scalar(vyy, "vyy.out");
    write_Scalar(vyz, "vyz.out");
    write_Scalar(vzz, "vzz.out");
    
    /* *************************************
     *  With the charge density we can first calculate the 
     *  scalar potential Pot. 
     *  By differentiating we get the longitudinal electric field.
     *  This is stored provisionally in the Ex and Ey fields.
     *  Remember: Ez = 0 in two dimensions
     */

    for (int j=0; j<GlGridSize_y; ++j) 
        for (int i=0; i<GlGridSize_x; ++i) { 
            In(i,j) = -(den(i,j)+n0);
//            std::cerr << "den " << i << " " << j << " " << In(i,j) << std::endl;
        }
    
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
    }

    // Perform wrapping
    for (int i=0; i<GlGridSize_x; ++i) {
        Ey(i,0) = Ey(i,GlGridSize_y-2);
	    Ey(i,GlGridSize_y-1) = Ey(i,1);
    }


    /* *************************************
     *  Now we want to calculate the magnetic field B.
     *  In two dimensions we calculate the z-component of the
     *  vector potential A_z which leads to the magnetic field
     *  components B_x and B_y.
     *  The component B_z can be calculated directly from rot j_z.
     */

    /* *************************************
     * The z--component of the vector potential.
     */
    
    pois->solve(jz,Az);

    
    /* *************************************
     * ... resulting in Bx and By
     */

	for (int i=0; i<GlGridSize_x; ++i) 
        for (int j=1; j<(GlGridSize_y-1); ++j) 
    	    Bx(i,j) = (Az(i,j+1) - Az(i,j-1)) / (2*dx[1]);

    for (int i=0; i<GlGridSize_x; ++i) {
        Bx(i,0) = Bx(i,GlGridSize_y-2);
	    Bx(i,GlGridSize_y-1) = Bx(i,1);
    }
    
	for (int i=1; i<(GlGridSize_x-1); ++i) 
        for (int j=0; j<GlGridSize_y; ++j) 
    	    By(i,j) = (Az(i-1,j) - Az(i+1,j)) / (2*dx[0]);

    for (int j=0; j<GlGridSize_y; ++j) {
        By(0,j) = By(GlGridSize_x-2,j);
	    By(GlGridSize_x-1,j) = By(1,j);
    }

    /* *************************************
     * The z--component of  -rot(j).
     */

	for (int i=1; i<(GlGridSize_x-1); ++i) 
        for (int j=1; j<(GlGridSize_y-1); ++j) 
                In(i,j) = - (jy(i+1,j) - jy(i-1,j)) / (2*dx[0])
                          + (jx(i,j+1) - jx(i,j-1)) / (2*dx[1]);

    for (int i=0; i<GlGridSize_x; ++i) {
        In(i,0) = In(i,GlGridSize_y-2);
	    In(i,GlGridSize_y-1) = In(i,1);
    }

    for (int j=0; j<GlGridSize_y; ++j) {
        In(0,j) =In(GlGridSize_x-2,j);
	    In(GlGridSize_x-1,j) = In(1,j);
    }

    pois->solve(In,Bz);
    
     
    /* *************************************
     *  The transverse electric field has to be calculated individually
     *  for each component. 
     */
    
    // This is the same for every component
	for (int i=0; i<GlGridSize_x; ++i) 
        for (int j=0; j<GlGridSize_y; ++j)
            Lambda(i,j) = csc*om2(i,j);
           
    // x-component of transverse electric field
	for (int i=1; i<(GlGridSize_x-1); ++i) 
        for (int j=1; j<(GlGridSize_y-1); ++j) {
            In(i,j) = csc*(
                (    vxx(i+1,j) - vxx(i-1,j)
                    +vxy(i+1,j) - vxy(i-1,j)
                    +vxz(i+1,j) - vxz(i-1,j)
                ) / (2*dx[0])                       /// grad (rho <vv>)
                +om2(i,j)*Ex(i,j)                   /// om2*E
                +sy(i,j)*Bz(i,j)-sz(i,j)*By(i,j)    /// q rho <v> x B
            );
        }
    
    for (int i=0; i<GlGridSize_x; ++i) {
        In(i,0) = In(i,GlGridSize_y-2);
	    In(i,GlGridSize_y-1) = In(i,1);
    }

    for (int j=0; j<GlGridSize_y; ++j) {
        In(0,j) =In(GlGridSize_x-2,j);
	    In(GlGridSize_x-1,j) = In(1,j);
    }
    
    helmh->solve(Out,In,Lambda);
    
	for (int i=0; i<GlGridSize_x; ++i) 
        for (int j=0; j<GlGridSize_y; ++j)
            Ex(i,j) += Out(i,j);
    
    // y-component of transverse electric field
	for (int i=1; i<(GlGridSize_x-1); ++i) 
        for (int j=1; j<(GlGridSize_y-1); ++j) {
            In(i,j) = csc*(
                (    vxy(i,j+1) - vxy(i,j-1)
                    +vyy(i,j+1) - vyy(i,j-1)
                    +vyz(i,j+1) - vyz(i,j-1)
                ) / (2*dx[1])                       /// grad (rho <vv>)
                +om2(i,j)*Ey(i,j)                   /// om2*E
                +sz(i,j)*Bx(i,j)-sx(i,j)*Bz(i,j)    /// q rho <v> x B
            );
        }
    
    for (int i=0; i<GlGridSize_x; ++i) {
        In(i,0) = In(i,GlGridSize_y-2);
	    In(i,GlGridSize_y-1) = In(i,1);
    }

    for (int j=0; j<GlGridSize_y; ++j) {
        In(0,j) =In(GlGridSize_x-2,j);
	    In(GlGridSize_x-1,j) = In(1,j);
    }
    
    helmh->solve(Out,In,Lambda);
    
	for (int i=0; i<GlGridSize_x; ++i) 
        for (int j=0; j<GlGridSize_y; ++j)
            Ey(i,j) += Out(i,j);
    
    // z-component of transverse electric field
	for (int i=1; i<(GlGridSize_x-1); ++i) 
        for (int j=1; j<(GlGridSize_y-1); ++j) {
            In(i,j) = csc*(sx(i,j)*By(i,j)-sy(i,j)*Bx(i,j));/// q rho <v> x B
        }
    
    for (int i=0; i<GlGridSize_x; ++i) {
        In(i,0) = In(i,GlGridSize_y-2);
	    In(i,GlGridSize_y-1) = In(i,1);
    }

    for (int j=0; j<GlGridSize_y; ++j) {
        In(0,j) =In(GlGridSize_x-2,j);
	    In(GlGridSize_x-1,j) = In(1,j);
    }
    
    helmh->solve(Out,In,Lambda);
    
	for (int i=0; i<GlGridSize_x; ++i) 
        for (int j=0; j<GlGridSize_y; ++j)
            Ez(i,j) = Out(i,j);
    
    /* *************************************
     *  This should be it!
     *  Now we can write out the diagnostics
     */
    
    write_Scalar(Ex, "Ex.out");
    write_Scalar(Ey, "Ey.out");
    write_Scalar(Ez, "Ez.out");
    write_Scalar(Bx, "Bx.out");
    write_Scalar(By, "By.out");
    write_Scalar(Bz, "Bz.out");
    write_Scalar(Az, "Az.out");
    write_Scalar(Pot,"Phi.out");
 
//    DiagField.Execute();


//    Task::Execute();
  	return false;
    
}

