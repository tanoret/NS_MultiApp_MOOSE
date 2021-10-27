#include "AirfoilAppApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
AirfoilAppApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

AirfoilAppApp::AirfoilAppApp(InputParameters parameters) : MooseApp(parameters)
{
  AirfoilAppApp::registerAll(_factory, _action_factory, _syntax);
}

AirfoilAppApp::~AirfoilAppApp() {}

void
AirfoilAppApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"AirfoilAppApp"});
  Registry::registerActionsTo(af, {"AirfoilAppApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
AirfoilAppApp::registerApps()
{
  registerApp(AirfoilAppApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
AirfoilAppApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  AirfoilAppApp::registerAll(f, af, s);
}
extern "C" void
AirfoilAppApp__registerApps()
{
  AirfoilAppApp::registerApps();
}
