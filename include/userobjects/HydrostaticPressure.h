#ifndef HYDROSTATICPRESSURE_H
#define HYDROSTATICPRESSURE_H

#include "GeneralUserObject.h"

// Forward Declarations
class HydrostaticPressure;

template<>
InputParameters validParams<HydrostaticPressure>();


/**
 * The HydrostaticPressure user object returns:
 * p = p(h, \vec{q}) where \vec{q} = (q_x, q_y) is the momentum vector
 * and various partial derivatives.
 *
 * We derive from UserObject so we can have a validParams()
 * function and be registered in the factory.
 */
class HydrostaticPressure : public GeneralUserObject
{
public:
  // Constructor
  HydrostaticPressure(const std::string & name, InputParameters parameters);

  // Destructor
  virtual ~HydrostaticPressure();

  /**
   * Called when this object needs to compute something.
   */
  virtual void execute() {};

  /**
   * Called before execute() is ever called so that data can be cleared.
   */
  virtual void initialize(){};

  /**
   * Finalize.  This is called _after_ execute() and _after_ threadJoin()!  This is probably where you want to do MPI communication!
   */
  virtual void finalize() {};

  // The interface for derived HydrostaticPressure objects to implement...
  virtual Real pressure(Real h=0., RealVectorValue vector_q=0.) const;

  // The derivative of pressure wrt h
  virtual Real dp_dh(Real h=0., RealVectorValue vector_q=0.) const;

  // The derivative of pressure wrt x-momentum (q_x)
  virtual Real dp_dqx(Real h=0., RealVectorValue vector_q=0.) const;

  // The derivative of pressure wrt y-momentum (q_y)
  virtual Real dp_dqy(Real h=0., RealVectorValue vector_q=0.) const;

  // Sound speed squared
  virtual Real c2(Real h=0., RealVectorValue vector_q=0.) const;

protected:
  Real _gravity;


};


#endif // HYDROSTATICPRESSURE_H

