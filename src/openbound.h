#ifndef OPENBOUND_H
#define OPENBOUND_H

/** Replacement for the VlasovDist when the scheme is used for the open
 *  boundary conditions.
 *  Wraps a distribution and presents only the boundary values
 */
class BoundMatrix  
{
  public:
    typedef enum {top, bottom, left, right} DirectionType;
  private:
    VlasovDist &dist;
    DirectionType direction;

    int lo[5];
    int hi[5];
    int dims[5];
    
  public:
    typedef enum {top, bottom, left, right} DirectionType;
    
    /** constructur */
    BoundMatrix(VlasovDist &dist_, DirectionType direction_=top);

    void setDirection(BoundMatrix::DirectionType direction_);

    //======================================================    
    //===== Methods implementing matrix operations
    //======================================================    
    
    /** */
    const int* getLow() const
    {
      return lo;
    }
    
    /** */
    const int* getHigh() const
    {
      return hi;
    }
    
    /** */
    const int* getDims() const
    {
      return dims;
    }
    
    
    /** */
    int getLow(int i) const
    {
      return lo[i];
    }
    
    /** */
    int getHigh(int i) const
    {
      return hi[i];
    }
    
    /** */
    int getDims(int i) const
    {
      return dims[i];
    }

    /** index operator, writing */
    T& operator()(int i, int j, int k, int l, int m)
    {
      return dist(i,j,k,l,m);
    }
    
    /** index operator, reading */
    T  operator()(int i, int j, int k, int l, int m) const
    {
      return dist(i,j,k,l,m);
    }

    //======================================================    
    //===== Methods implementing boundary operations
    //====================================================== 
    
    void copyBoundary();
    
  private:
    void copyVelocity(int fromi, int fromj, int toi, int toj);

};


class OpenBoundForce {
  protected:
    /// the physical space between grid points
    PositionD dx;
    VelocityD dv;
    ForceFieldBase &base;
    BoundMatrix Distribution;
  private:
    double density;
    VelocityD current;
    NumMatrix<VelocityD, DIMENSION> veloffset;
  public:
    OpenBoundForce
      (
        ForceFieldBase &base_
      );

    void setDirection(BoundMatrix::DirectionType direction);
    
        
    //======================================================    
    //===== Methods needed by the Scheme  
    //======================================================    
    
    
    /** @brief The grid spacing in v_x--direction in physical
     *  dimensions
     */
    double deltaVx() {
      return dv[0];
    }

    /** @brief The grid spacing in v_y--direction in physical
     *  dimensions
     */
    double deltaVy() {
      return dv[1];
    }

    /** @brief The grid spacing in v_z--direction in physical
     *  dimensions
     */
    double deltaVz() {
      return dv[2];
    }

    /**@brief Returns the velocity in the physically
     * normalized units.
     * (NOT in grid points per timestep!!)
     */
    VelocityD velocity(const VelocityI &Vi) {
      return base.velocity(Vi);
    }

    VelocityD velocity(const VelocityD &Vi) {
      return base.velocity(Vi);
    }

    /** The force at one position given the velocity.
     *  Actually the displacement in the velocity space is returned
     */
    VelocityD Force(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

    /** The forces individually by component. Calculation is the same as in 
     *  the Force method but all operations not needed by the x--component
     *  are omitted. This method is therefore slightly faster.
     */
    double ForceX(const PositionI &Pos, 
                  const VelocityD &Vel,
                  double dt);

    /** The forces individually by component. Calculation is the same as in 
     *  the Force method but all operations not needed by the y--component
     *  are omitted. This method is therefore slightly faster.
     */
    double ForceY(const PositionI &Pos, 
                  const VelocityD &Vel,
                  double dt);

    /** The forces individually by component. Calculation is the same as in 
     *  the Force method but all operations not needed by the z--component
     *  are omitted. This method is therefore slightly faster.
     */
    double ForceZ(const PositionI &Pos, 
                  const VelocityD &Vel,
                  double dt);
                  
                  
    //======================================================    
    //===== Methods for calculating boundary conditions
    //======================================================    
    
    void setMoments(double density_, VelocityD current_);
    void initialize();

  private:
    
    double integrateDensity(int i, int j);
    VelocityD integrateCurrent(int i, int j)
};

// may convert this to templated class
class OpenBound : public PosFluxCons3rdOrder<OpenBoundForce>
{
  public:
    typedef BoundMatrix::DirectionType DirectionType;
    
    OpenBound(ForceFieldBase& base, DirectionType direction);
    void apply();
};

#endif
