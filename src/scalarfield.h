// -*- C++ -*-
// $Id$

#include "index.h"
#include "matrix.h"
#include "vlasov.h"

#ifndef SCALARFIELD_H
#define SCALARFIELD_H

/** A class defining scalar fields in the simulation.
 * The scalar field needs to carry some additional information about the
 * parity. A typedef is therefore not sufficient.
 * 
 * The parity defines how the field behaves under the parity (mirror) transformation.
 * Eg the electric field components are odd while the magnetic fields components are
 * even.
 * 
 * For boundary conditions we also need to know which field component the scalar 
 * field resebles. The perpendicular component behaves differently than the parallel
 * components.
 */
class ScalarField : public NumMatrix<double, DIMENSION>
{
  private:
      /// The parity value (+1 or -1)
      int parity;
      /// The component index 0,1,2 stands for x, y and z.
      int component;
  public:
      //static const int ScalarComponent;
      //static const int XComponent;
      //static const int YComponent;
      //static const int ZComponent;
      enum {ScalarComponent=0, XComponent=1, YComponent=2, ZComponent=3};
      static const int EvenParity;
      static const int OddParity;
  public:
      /// Accessor method setting the parity
      void setParity(int parity_) { parity=parity_; }
      /// Accessor method setting the parity and the component
      void setParity(int parity_, int component_) { parity=parity_; component=component_; }

      /// Accessor method getting the parity
      int getParity() { return parity; }

      /// Accessor method setting the component
      void setComponent(int component_) { component=component_; }
      /// Accessor method getting the component
      int getComponent() { return component; }

      /// Set the parity to +1
      void setEven() { parity=1; }
      /// Set the parity to -1
      void setOdd() { parity=-1; }
      /// Check if the parity is +1
      bool isEven() { return 1==parity; }
      /// Check if the parity is -1
      bool isOdd() { return -1==parity; }
};

#endif
