// -*- C++ -*-

#include <string>
#include "config.h"
#include "numeric.h"
#include "potential.h"
#include "boundary.h"

#ifndef WRAPVLASOV_H
#define WRAPVLASOV_H

class Potential;

template<class ForceField>

/** @brief This class wraps the Vlasov advancers in a SpeciesBase class
 *  so that Opar can use it.
 *
 *  The abstract methods MakeRho and Rho are redefined to integrate over
 *  the velocity space and return the resulting density.
 *
 *  The template class ForceField must declare a method "Force" that
 *  takes a PositionI as argument an returns a VelocityD
 *
 *  The ForceField template class additionally has to supply the Init, 
 *  and MakeParamMap method. The Init mehod has to take a scaling factor
 *  of type double whith which to multiply the forcefield.
 *
 */
class VlasovSpecies : public ForceField {
protected:
  PositionD GridSize; ///< the size of the numerical grid in the space direction
  PositionI GridSizeI; ///< the size of the numerical grid in the space direction
  VelocityD VelSize,  ///< the size of the numerical grid in the velocity directions
    VelSizeH; ///< half of VelSize
  PositionD BoxRange; ///< the physical extent of the spatial grid
  PositionD dx;       ///< the physical space between grid points

  double f_infty;
protected:
  VlasovDist Distribution;
  /// The Runge Kutta Temporaries
  VlasovDist T1, T2;
  int cycle, RKState;
  VlasovDist TempDist;
protected:
  ScalarField gRho;   ///< contains the particle density
  ScalarField EKin;   ///< contains the kinetic energy when calculated
  ScalarField Jx, Jy, Jz;   
  ScalarField Vxx, Vxy, Vxz;   
  ScalarField Vyy, Vyz, Vzz;   
  double Charge;      ///< The Charge of the particles
  double Mass;        ///< The Mass of the particles
  double dt;          ///< The timestep (in s)
  VelocityD VRange;   ///< range in velocty space
        
  typename ForceField::FieldType *pPot;
        
  Boundary *boundary;
        
  int t;
public:
  /// Default constructor
  VlasovSpecies (Boundary *boundary_);
  /// Destructor
  virtual ~VlasovSpecies ();

  virtual void Init ();
  virtual void Execute ();

  void MakeRho();
  void MakeJs();
  
  ScalarField &getJx() { return Jx; }
  ScalarField &getJy() { return Jy; }
  ScalarField &getJz() { return Jz; }
        
  /// Return the Charge
  double GetCharge () const { return Charge; };

  /// Return the time step
  double Getdt () const { return dt; };

  PositionD position(const PositionI &Xi) {
    return BoxRange*Xi/GridSize;
  }
        
  int gridSize(int i) {
    return GridSizeI[i];
  }
        
  double deltaX() {
    return BoxRange[0]/GridSize[0];
  }
        
  double deltaY() {
    return BoxRange[1]/GridSize[1];
  }
        
  double deltaVx() {
    return VRange[0]/VelSizeH[0];
  }
        
  double deltaVy() {
    return VRange[1]/VelSizeH[1];
  }
        
  double deltaVz() {
    return VRange[2]/VelSizeH[2];
  }
        
  /**@brief Returns the velocity in the physically
   * normalized units.
   * (NOT in grid points per timestep!!)
   */
  VelocityD velocity(const VelocityI &Vi) {
    return (VRange*(VelocityD(-0.5,-0.5,-0.5)+Vi))/VelSizeH - VRange;
  }
  /// Return a reference to the density grid
  ScalarField& Rho () { return gRho; };
                  
  void resize(PhasePositionI &low, PhasePositionI &high);

  void write(ostream &);
  void writeYVySlice(int x,int vx, int vz, string fname);
  void writeVxVySlice(int t, int x, int y, int vz, ostream &Slice);

  VelocityD getJ(int i, int j);
  FixedArray<double,6> getVVTens(int i, int j);
  ScalarField &KineticEnergy();
  double TotalEnergy();
        
        
  bool master() {
    return boundary->master();
  }
  int procnum() {
    return boundary->procnum();
  }
        
private:
  void advance(double timestep);
  void advanceSpace_x(double timestep);
  void advanceSpace_y(double timestep);
  void advanceVel_x(double timestep);
  void advanceVel_y(double timestep);
  void advanceVel_z(double timestep);
        
  void advanceStepA(double timestep);
  void advanceStepB(double timestep);
  void advanceStepFull(double timestep);


  double epsilonLeft(double fj, double fjp);
  double epsilonRight(double fj, double fjm);
  void InterpolationInitStep(const VlasovDist &Dist);
  double densityError();
  void correctDensityError(double err);
  double interpolateX    (const PositionI &Xi,
                          const VelocityI &Vi,
                          int j,
                          double alpha,
                          const VlasovDist &Dist);
                            
  double interpolateY    (const PositionI &Xi,
                          const VelocityI &Vi,
                          int j,
                          double alpha,
                          const VlasovDist &Dist);
                            
  double interpolateVx   (const PositionI &Xi,
                          const VelocityI &Vi,
                          int j,
                          double alpha,
                          const VlasovDist &Dist);
        
  double interpolateVy   (const PositionI &Xi,
                          const VelocityI &Vi,
                          int j,
                          double alpha,
                          const VlasovDist &Dist);

  double interpolateVz   (const PositionI &Xi,
                          const VelocityI &Vi,
                          int j,
                          double alpha,
                          const VlasovDist &Dist);
        
  int bound_plus()   { return 1;}
  int bound_minus()  { return 1;}
}; // VlasovSpecies


class Potential;

/** @brief Implements the electrostatic force for plugging into the 
 *  VlasovSpecies
 */
class EFieldForce {
private:
  /// Pointer to the potential
  Potential* pE;
  /// Scaling constant
  double dttx;
  ScalarField FEngy;
public:
  EFieldForce() {};

  ScalarField &GetEx();
  ScalarField &GetEy();
        
  ScalarField &FieldEnergy();

  VelocityD Force(const PositionI &Pos, 
                  const VelocityD &Vel,
                  double dt);
  void Init(double dttx_, double vnorm_, Potential* pPot_);
        
  typedef Potential FieldType;
};

/** @brief Implements the electrostatic force with a constant magnetic
 *  field for plugging into the VlasovSpecies
 */
class EBFieldForce {
private:
  Potential* pPot;
  /** @brief Magnetic Field normalized by 
   *  \f$\frac{v_A}{c}\frac{1}{B_0}\f$.
   */
  VelocityD B;
  double dttx;
  double vnorm;
public:
  EBFieldForce() {};
        
  ScalarField &GetEx();
  ScalarField &GetEy();
        
  VelocityD Force(const PositionI &Pos, 
                  const VelocityD &Vel,
                  double dt);
  void Init(double dttx_, double vnorm_, Potential* pPot_);
  void setBField(VelocityD B_);

  typedef Potential FieldType;
};

class Darwin;

/** @brief Implements a force field that plugs into the concrete advancers
 *  taking the force from the Opar electric field.
 *
 *  Implements the Force method for the Advancer and the Init and MakeParamMap
 *  method to work together with the VlasovSpecies.
 */
class EMDarwinForce {
public:
  Darwin* pFields;
  double dttx;
  double vnorm;
  ScalarField FEngy;
public:
  EMDarwinForce() {};

  ScalarField &GetEx();
  ScalarField &GetEy();
  ScalarField &GetEz();
  ScalarField &GetBx();
  ScalarField &GetBy();
  ScalarField &GetBz();
        
  ScalarField &FieldEnergy();
        
  VelocityD Force(const PositionI &Pos, 
                  const VelocityD &Vel,
                  double dt);
  void Init(double dttx_, double vnorm_, Darwin* pFields_);

  typedef Darwin FieldType;
};

void write_Distribution(VlasovDist &, const char *);

#include "wrapvlasov_temp.cpp"

#endif
