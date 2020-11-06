#ifndef PROBE_ENERGIECALCULATION_H
#define PROBE_ENERGIECALCULATION_H

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include <map>
#include <unordered_map>
#include "toyChemUtil.hh"

double energie_calculation (std::string molecule,map<string,double>& look_up_map);//(std::unordered_map<std::string, ggl::chem::Molecule*>& probe_smilesstrings);

#endif
