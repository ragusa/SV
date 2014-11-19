#include "SvApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

template<>
InputParameters validParams<SvApp>()
{
  InputParameters params = validParams<MooseApp>();
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
}

void
SvApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
