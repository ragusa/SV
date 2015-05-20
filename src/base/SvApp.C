#include "SvApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

// kernels
#include "SV_Continuity.h"
#include "SV_Continuity_Explicit.h"
#include "SV_Momentum.h"
#include "SV_Momentum_Explicit.h"
#include "SV_ArtificialViscFlux.h"
#include "SV_ArtificialViscFlux_Explicit.h"

// auxkernels
#include "BathymetryAux.h"
#include "PressureAux.h"
#include "EntropyAux.h"
#include "EntropyFluxAux.h"
#include "VelocityAux.h"
#include "NormVectorAux.h"

// ics
#include "StepIC.h"

// bcs
#include "SVSetWaterHeightInletBC.h"
#include "SVSetWaterHeightOutletBC.h"
#include "SVSetWaterVelocityInletBC.h"
#include "SVSetWaterVelocity.h"
#include "SolidWallBC.h"

// userobject: eos
#include "HydrostaticPressure.h"

// materials
#include "ComputeViscCoeff.h"

// bathymetry function
#include "Bathymetry1D_step.h"


template<>
InputParameters validParams<SvApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

SvApp::SvApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  SvApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  SvApp::associateSyntax(_syntax, _action_factory);
}

SvApp::~SvApp()
{
}

void
SvApp::registerApps()
{
  registerApp(SvApp);
}

void
SvApp::registerObjects(Factory & factory)
{
  // kernels
  registerKernel(SV_Continuity);
  registerKernel(SV_Continuity_Expl);
  registerKernel(SV_Momentum);
  registerKernel(SV_Momentum_Expl);
  registerKernel(SV_ArtificialViscFlux);
  registerKernel(SV_ArtificialViscFlux_Expl);

  // auxkernels
  registerAux(BathymetryAux);
  registerAux(PressureAux);
  registerAux(EntropyAux);
  registerAux(EntropyFluxAux);
  registerAux(VelocityAux);
  registerAux(NormVectorAux);

  // ics
  registerInitialCondition(StepIC);

  // bcs
  registerBoundaryCondition(SVSetWaterHeightInletBC);
  registerBoundaryCondition(SVSetWaterHeightOutletBC);
  registerBoundaryCondition(SVSetWaterVelocityInletBC);
  registerBoundaryCondition(SVSetWaterVelocity);
  registerBoundaryCondition(SolidWallBC);

  // userobject: eos
  registerUserObject(HydrostaticPressure);

  // materials
  registerMaterial(ComputeViscCoeff);

  // functions
  registerFunction(Bathymetry1D_step);
}

void
SvApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
