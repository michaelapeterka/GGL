#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>

#include <iostream>
#ifndef PROBE_ENERGIECALCULATION_H
#define PROBE_ENERGIECALCULATION_H

#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include <map>

using namespace std;


double energie_calculation(std::map<std::string, double>&);

#endif
