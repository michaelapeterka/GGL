#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include "probe_energiecalculation.h"
#include "toyChemUtil.hh"
#include <map>
#include<boost/algorithm/string.hpp>

using namespace std;
using namespace OpenBabel;
using namespace boost::algorithm;

double energie_calculation(std::string molecule,std::map<string,double>& look_up_map)
{
	//canonicalization of the input smiles
	//1. stringstream for input_smile
	   string input = molecule;
           stringstream input_smiles;
 
           //fill stringstream
	   input_smiles << input << endl;

	//2. stringstream for output
	   stringstream output_smiles;
	
	//3. Set up OpenBabel format converter
	  OBConversion conv(&input_smiles,&output_smiles);
	  if(!conv.SetInAndOutFormats("smi","can")){
		  cerr << "Can not open output file" << endl; 
          }
	   
	  //convert molecule
	  int n = conv.Convert();
	  cout << n << "molecules converted" << endl;
         
        //4. Energy value from look-up table
	  //use boost library to trim (delete the whitespace before and after the canonical string to get the same size)
 	  string canonical = output_smiles.str();
	  string new_canonical = trim_copy(canonical);

	  //iterate through map
	   double value = 0.0;
	   map<string,double>::iterator iten = look_up_map.find(new_canonical);
	   if(iten != look_up_map.end())
	   {
		   value = iten ->second;
		   cout << "GOT IT" << iten -> first << " : " << iten->second << endl;
           }
	   else
	   { 
		   cerr << "sorry, no such molecule in look-up table" << endl; 
	   }	    
return value; 
}
