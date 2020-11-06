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

double energie_calculation(std::string molecule,std::map<string,double>& look_up_map)//(std::unordered_map<std::string, ggl::chem::Molecule*>& probe_smiles)
{
	//canonicalization of the input smiles
	//1. stringstream for input_smile
	   string input = molecule;
           stringstream input_smiles;
 
           //fill stringstream
	   input_smiles << input << endl;

	//2. stringstream for output
	   stringstream output_smiles;

	/*declare map
	map<string,double> SMILES_energy;
	string SMILES_string;*/
	
	//3. Set up OpenBabel format converter
	  OBConversion conv(&input_smiles,&output_smiles);
	  if(!conv.SetInAndOutFormats("smi","can")){
		  cerr << "Can not open output file" << endl; 
          }
	   
	  //convert molecule
	  int n = conv.Convert();
	  cout << n << "molecules converted" << endl; /*OBConversion conv;conv.SetInFormat("smi");*/
         
        //4. Energy value from look-up table
	  //use boost library to trim (delete the whitespace before and after the canonical string to get the same size)
 	  string canonical = output_smiles.str();
	  string new_canonical = trim_copy(canonical);

	  //iterate through map
	   double value = 0.0;
	   map<string,double>::iterator iten = look_up_map.find(new_canonical);
	   if(iten != look_up_map.end()){
		   value = iten ->second;
		   cout << "GOT IT" << iten -> first << " : " << iten->second << endl;
           }
	   else{ 
		   cerr << "sorry, no such molecule in look-up table" << endl; 
	   }	   
	/*set up OpenBabel Mol
	OBMol mol;

	//set up forcefield
	OBForceField *FF = OBForceField::FindType("MMFF94");
	if(!FF){
	cout << "can not set up forcefield"<< endl;
	}	

	//read the stringstream
	while(input_smiles >> SMILES_string){
		if(!conv.ReadString(&mol,SMILES_string)){
			cout << "can not read SMILES" << SMILES_string << endl;
			return 1;
		}
		if(!FF ->Setup(mol)){
			cout << "can not set up forcefield" << endl;
			return 1;
		}
		//iterator
		map<string,double>::iterator iter;
		iter = SMILES_energy.find(SMILES_string);
		if(iter == SMILES_energy.end()){
			SMILES_energy[SMILES_string] = FF->Energy();
		}
		mol.Clear();
	}
	/*cout <<"--------------Energyvalues in Funktion Energyberechnung------------" << endl;
	for(map<string,double>::iterator i = SMILES_energy.begin();i !=SMILES_energy.end(); ++i)
	{
		cout << i->first <<":"<< i->second << endl;
	
	}
	cout <<"------------------ENDE EXTERNE FUNKTION--------------" << endl;*/ 
	return value; //SMILES_energy[SMILES_string]; 
}
  /*summe initialisieren fuer den return-wert

  double sum = 0.0;
  double sum1 = 0.0;
  //declare the iterator for the map smiles_energy
  std::map<std::string, double>::iterator iter;
//als allererstes inhalt der map in einen stringstream schreiben, da sonst obenergy nichts berechnen kann
//stringstream input_smile zum speichern und lesen des uebergebenen parameters fuer die Energieberechnung(smilestrings); "fuellen" mit den SMILES-Strings (also keys) der map.
  
  std::string input;
  input = educts;
  std::stringstream input_smiles;

	  input_smiles << input <<std::endl; 
	 
	

//2. stringstream fuer output_smiles
  
  std::stringstream output_smiles;
  //std::ofstream file_save_smiles_energy("smilesWithEnergy.smi", std::ios::app );//| std::ios::out);

  
 //3. declare map

  std::map<std::string, double> smiles_energy;
  std::string SmilesString;
  double value_energy;

 //4. set up openbabel format converter

 OpenBabel::OBConversion conv (&input_smiles, &output_smiles);//!!!!!!!!!!!!konvertierung von input smile in ein canonical 

       if(!conv.SetInAndOutFormats("smi","can"))
         {
           std::cerr << "Cannot open the output file" << std::endl;
         }

       // molecules are converted  !!!!!!!!!!
       int n = conv.Convert();

       std::cout << n << "molecules converted"<< std::endl;
       //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


//6. starting calculation obenergy
//set up the openbabel OBMol

 OpenBabel::OBMol mol;

 //declaration forcefield

 OpenBabel::OBForceField *FF = OpenBabel::OBForceField::FindType("MMFF94");
//not get the forcefield
if (!FF)
{
std::cerr << "Cannot find the forcefield" << std::endl << std::endl;
}


//read the file line by line and setup forcefield

 while (output_smiles >> SmilesString)
  {
   if (!conv.ReadString(&mol, SmilesString))
    {
      std::cout << "Cannot read SMILES" << SmilesString<< std::endl;
      return 1;
    }

   if (!FF->Setup(mol)) //setup forcefield for mol
    {
     std::cout << "Cannot set up forcefield" << std::endl;
                                                                return 1;
    }

   //print the energy and the unit for the molecules
   //looks if the SmileString is already in the map or not
   //save the smilesstring and the value in a file
   //!!!!achtung neu!!!!//if(smiles_energy.find(SmilesString) == smiles_energy.end())//if the string is not in the map

   iter = smiles_energy.find(SmilesString);
   if(iter == smiles_energy.end())
    {
      //std::cout << "key value pair not present in the map so calculate the energy for this smilestring" << std::endl;
      value_energy = FF->Energy();
       smiles_energy[SmilesString] = value_energy;//FF ->Energy();
    }
   mol.Clear();
  }//!!!!!!!

 return smiles_energy[SmilesString];//(sum);
}*/
 
                      

