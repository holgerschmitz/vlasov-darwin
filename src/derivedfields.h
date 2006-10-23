#include <string>
#include <map>

#ifndef DERIVEDFIELDS_H
#define DERIVEDFIELDS_H

#include "scalarfield.h"
#include "stlpwrapper.h"
#include "boundary.h"

  /** @file derivedfields.h
	* @brief derived field classes 
	*
	* Implements density distribution, its integral and other derived fields
	*/


class ForceFieldBase;

///abstract base class for derived fields
class DistributionDerivedField {
  protected:
	///boundary conditions
      Boundary *boundary;
  public:
	//---------------------------------------------------------
	//constructor & destructor
	///constructor
      DistributionDerivedField(Boundary *boundary_);
	///virtual destructor, since this is a base class
      virtual ~DistributionDerivedField() {}
	//---------------------------------------------------------
	/** @brief prototype for the calc function
	  *
	  * purely virtual function, needs to be overwritten by each derived class
	  */
      virtual void calc(ForceFieldBase &) = 0;
	/** @brief prototype for the name function
	 *
	 * purely virtual function, needs to be overwritten by each derived class
	 */
      virtual const char* name() = 0;
	/** @brief prototype for the return-by-name function
	 *
	 * purely virtual function, needs to be overwritten by each derived class
	 */
      virtual ScalarField& getField(std::string name) = 0;
};

///wrapped pointer to DistributionDerivedField
typedef PtrWrapper<DistributionDerivedField> pDistributionDerivedField;

///container class for derived fields
class DerivedFieldsContainer {
  private:
	  ///typedef for a map of strings against derived fields
      typedef std::map<std::string,pDistributionDerivedField> MapType;
      ///map for  storing pointers to DerivedField objects
      MapType derivedFields;
      ///boundary of the field
      Boundary *boundary;
  public:
	/// constructor, initializes boundary
      DerivedFieldsContainer(Boundary *boundary_);
	/// add a field to container
      pDistributionDerivedField add(pDistributionDerivedField field);
      /// calls calc(dist) for each derived field 
      void update(ForceFieldBase &dist);
      /// Returns a derived field by its name or NULL if it is not found
      pDistributionDerivedField getField(std::string);
};

///the particle density distribution
class DistMomentRhoBase : public DistributionDerivedField {
  protected:
      /// contains the particle density
      ScalarField Rho;
  public:
	///constructor
      DistMomentRhoBase(Boundary *boundary_);
      ///returns name, which is "rho"
      const char* name() {return "rho";}
      ///returns the scalarfield Rho
      ScalarField& getRho() { return Rho; }
	/// returns the scalarfield Rho by name
      ScalarField& getField(std::string name);
};

/** @brief The density integral in first order. 
 * 
 *  Note that the density is
 *  obtained simply by summing over the distribution. Since not the 
 *  distribution itself but the INTEGRAL over one cell is stored, this
 *  first order summation is actually EXACT!
 */
class DistMomentRhoOne : public DistMomentRhoBase {
  public:
	  ///constructor, initializes boundary
      DistMomentRhoOne(Boundary *boundary_) 
        : DistMomentRhoBase(boundary_) {}
    /** Calculate the density integral from the distribution function.
      * The region of the electromagnetic fields is on
      * grid cell smaller in every direction
      * It's important that the density is only calculated on the
      * inner cells, since ScalarFieldReduce simply adds all the densities of
      * the processes.
	*/
      void calc(ForceFieldBase &vlasov);
};

///current density and second order velocity moments
class DistMomentVelocitiesBase : public DistributionDerivedField {
  protected:
      ScalarField 		Jx, 	///< The x-component of the current densitiy
					Jy, 	///< The y-component of the current densitiy
					Jz; 	///< The z-component of the current densitiy
      
      ScalarField 		Vxx, 	///< mixed second order velocity moment of the distribution
					Vxy,  ///< mixed second order velocity moment of the distribution
					Vxz;  ///< mixed second order velocity moment of the distribution
      ScalarField 		Vyy, 	///< mixed second order velocity moment of the distribution
					Vyz,  ///< mixed second order velocity moment of the distribution
					Vzz;	///< mixed second order velocity moment of the distribution
  public:
	///constructor, initializes boundary
      DistMomentVelocitiesBase(Boundary *boundary_);
      ///returns name, which is "vels"
      const char* name() {return "vels";}

      /** Return a reference to the grid containing the 
       *  x-component of the current density
       */
      ScalarField &getJx() { return Jx; }

      /** Return a reference to the grid containing the 
       *  y-component of the current density
       */
      ScalarField &getJy() { return Jy; }

      /** Return a reference to the grid containing the 
       *  z-component of the current density
       */
      ScalarField &getJz() { return Jz; }

      /// Returns the vector of the current density at a grid position
      VelocityD getJ(int i, int j);

      /** @brief Returns the values of the second order velocity moments
	 *
       *  in a six-component vector
       */
      FixedArray<double,6> getVVTens(int i, int j);
	///get field by name
      ScalarField& getField(std::string name);
};    

/** @brief velocity distribution in first order 
  *
  * since this information is stored in the velocity space numerical grid it is exact.
  */
class DistMomentVelocitiesOne : public DistMomentVelocitiesBase {
  public:
	  ///constructor, initializes boundary
      DistMomentVelocitiesOne(Boundary *boundary_) 
        : DistMomentVelocitiesBase(boundary_) {}
	
     /** Calculate the velocity integral from the distribution function
       */
      void calc(ForceFieldBase &vlasov);
};

/** @brief velocity distribution in second order
  *
  * used for the diagonal elements of the pressure tensor, where the velocity appears in second order
  *
  * \f$ P_{ii} = \int v_i  v_i f(\vec v) d^3 v \f$
  *
  * Appart from the diagonal elements identical to the first order velocity class,
  * since only the non-mixed second order moments have to be evaluated using the higher order
  * scheme
  */
class DistMomentVelocitiesTwo : public DistMomentVelocitiesBase {
  public:
	///constructor, initializes boundary
      DistMomentVelocitiesTwo(Boundary *boundary_) 
        : DistMomentVelocitiesBase(boundary_) {}
	/// Calculate velocites from the distribution function
      void calc(ForceFieldBase &vlasov);
};

/// scalar field, containing the heat flux
class DistMomentHeatFluxBase : public DistributionDerivedField {
  protected:
      /// scalarfield, containing the heat flux
      ScalarField HFluxX;
  public:
	///constructor
      DistMomentHeatFluxBase(Boundary *boundary_);
	///name of this field
      const char* name() {return "hflux";}
	///accessor method return heat flux 
      ScalarField& getFluxX() { return HFluxX; }
	///return field by name
      ScalarField& getField(std::string name);
};

/// heat flux in first order, used for diagnostics
class DistMomentHeatFluxOne : public DistMomentHeatFluxBase {
  public:
	///constructor
      DistMomentHeatFluxOne(Boundary *boundary_) 
        : DistMomentHeatFluxBase(boundary_) {}
	/// Calculate the heat flux from the distribution function
      void calc(ForceFieldBase &vlasov);
};

///first order density distribution
typedef DistMomentRhoOne DistMomentRho;
///second order velocity distribution
typedef DistMomentVelocitiesTwo DistMomentVelocities;
/// first order heat flux
typedef DistMomentHeatFluxOne DistMomentHeatFlux;

#endif
