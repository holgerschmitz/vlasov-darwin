// -*- C++ -*-
// $Id$

#include "wrapvlasov.h"
#include <complex>

#ifndef BERNSTEIN_H
#define BERNSTEIN_H

/** @brief Initialises the Vlasov distribution function with a
 *  Maxwellian plus a perturbation which excites Bernstein waves
 *
 *  The parameters k and omega have to be set through the input.
 *  They need to correspond to a single Bernstein mode. Use code in
 *  helpers/berntein to calculate values
 *
 *  The wave vector lies in xy-direction and the magnetic field is assumed
 *  to lie in the z-direction
 */
class VlasovBernsteinInit : public VlasovInitialiser {
  protected:
    /// density in particles\f$m^{-{\rm dim}}\f$
    double N;
    /// thermal velocity
    VelocityD v_th;
    /// The perturbation wave vector k
    PositionI wavevec;
    /// The perturbation wave vector used in the formula for the distribution function
    double distk;
    /// The perturbation frequency used in the formula for the distribution function
    double omega;
    /// The perturbation amplitude
    double amp;
    /// The integration accuracy
    int accuracy;

    /// pointer to the owning VlasovSpecies class
    ForceFieldBase *pVlasov; 
  public:
	  /// Default constructor
    VlasovBernsteinInit();
	  /// Destructor
    virtual ~VlasovBernsteinInit();

    /** @brief Perform distribution initialisation.
     *  The Maxwellian distribution will be written into the
     *  dist parameter.
     */
    virtual void initialise(ForceFieldBase *pVlasov);
  protected:
    /// Create the parameter map for the readable object
    virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
  private:
    typedef std::complex<double> Complex;
    Complex ComplexCore(double vperp, double phi, int n);
    double DistFactor(double vperp, int n);
    Complex Distribution1(double vperp, double phi);
    Complex Distribution1Integral(double vxmin, double vymin, double vxmax, double vymax, int steps);
};

#endif
