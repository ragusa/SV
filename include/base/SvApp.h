#ifndef SVAPP_H
#define SVAPP_H

#include "MooseApp.h"

class SvApp;

template<>
InputParameters validParams<SvApp>();

class SvApp : public MooseApp
{
public:
  SvApp(InputParameters parameters);
  virtual ~SvApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* SVAPP_H */
