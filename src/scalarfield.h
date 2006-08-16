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
  public:
      typedef enum 
        {ScalarComponent=0, XComponent=1, YComponent=2, ZComponent=3, 
         XYComponent=4, XZComponent=5, YZComponent=6} ComponentType;
            
      typedef enum 
        {BX, BY ,BZ, EX, EY, EZ, JX, JY, JZ, RHO, PHI, AX, AY, AZ, THETA, unspecified} FieldComponent;
      static const int EvenParity;
      static const int OddParity;

  private:
      /// The parity value (+1 or -1)
      int parity;
      /// The component index 0,1,2 stands for x, y and z.
      ComponentType component;

      FieldComponent ftype;
  public:
      /// Accessor method setting the parity
      void setParity(int parity_) 
        { parity=parity_; ftype = unspecified;}
        
      /// Accessor method setting the parity and the component
      void setParity(int parity_, ComponentType component_) 
        { parity=parity_; component=component_; ftype = unspecified;}

      /// Accessor method getting the parity
      int getParity() { return parity; }

      /// Accessor method setting the component
      void setComponent(ComponentType component_)
        { component=component_; ftype = unspecified;}
        
      /// Accessor method getting the component
      ComponentType getComponent() { return component; }
      FieldComponent getFieldType() { return ftype; }

      void setFieldType(FieldComponent type)
      {
        ftype = type;
        switch (ftype)
        {
          case BX:
            component = XComponent;
            parity = EvenParity;
            break;
          case BY:
            component = YComponent;
            parity = EvenParity;
            break;
          case BZ:
            component = ZComponent;
            parity = EvenParity;
            break;
            
          case EX:
          case JX:
          case AX:
            component = XComponent;
            parity = OddParity;
            break;
          case EY:
          case JY:
          case AY:
            component = YComponent;
            parity = OddParity;
            break;
          case EZ:
          case JZ:
          case AZ:
            component = ZComponent;
            parity = OddParity;
            break;

          case RHO:
          case PHI:
          case THETA:
            component = ScalarComponent;
            parity = EvenParity;
            break;
            
          default: 
            break;
        }
      }
      

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
