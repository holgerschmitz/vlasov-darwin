#include <string>
#include <map>

#ifndef DERIVEDFIELDS_H
#define DERIVEDFIELDS_H

#include "stlpwrapper.h"
#include "boundary.h"

class ForceFieldBase;

class DistributionDerivedField {
  protected:
      Boundary *boundary;
  public:
      DistributionDerivedField(Boundary *boundary_);
      virtual ~DistributionDerivedField() {}
      virtual void calc(ForceFieldBase &) = 0;
      virtual const char* name() = 0;
      virtual ScalarField& getField(std::string name) = 0;
};

typedef PtrWrapper<DistributionDerivedField> pDistributionDerivedField;

class DerivedFieldsContainer {
  private:
      typedef std::map<std::string,pDistributionDerivedField> MapType;
      MapType derivedFields;
  public:
      void add(pDistributionDerivedField field);
      void update(ForceFieldBase &dist);
      /// Returns a derived field by its name or NULL if it is not found
      pDistributionDerivedField getField(std::string);
};


class DistMomentRho : public DistributionDerivedField {
  private:
      /// contains the particle density
      ScalarField Rho;
  public:
      DistMomentRho(Boundary *boundary_);
      void calc(ForceFieldBase &vlasov);
      const char* name() {return "rho";}
      ScalarField& getRho() { return Rho; }
      ScalarField& getField(std::string name);
}; 
   
class DistMomentVelocities : public DistributionDerivedField {
  private:
      /// The current densities in x,y and z-direction
      ScalarField Jx, Jy, Jz;
      /// The mixed second order velocity moments of the distribution
      ScalarField Vxx, Vxy, Vxz;   
      ScalarField Vyy, Vyz, Vzz;
  public:
      DistMomentVelocities(Boundary *boundary_);
      void calc(ForceFieldBase &vlasov);
      const char* name() {return "vels";}

      /** Return a reference to the grid containing the 
       *  x--component of the current density
       */
      ScalarField &getJx() { return Jx; }

      /** Return a reference to the grid containing the 
       *  y--component of the current density
       */
      ScalarField &getJy() { return Jy; }

      /** Return a reference to the grid containing the 
       *  z--component of the current density
       */
      ScalarField &getJz() { return Jz; }

      /// Returns the vector of the current density at a grid position
      VelocityD getJ(int i, int j);

      /** @brief Returns the values of the second order velocity moments
       *  in a six-component vector
       */
      FixedArray<double,6> getVVTens(int i, int j);
      ScalarField& getField(std::string name);
};    

#endif
