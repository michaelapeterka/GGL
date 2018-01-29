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
#include <unordered_map>

double energie_calculation(std::unordered_map<std::string, ggl::chem::Molecule*>& probe_smiles)
{

  //summe initialisieren fuer den return-wert

  double sum = 0.0;
  double sum1 = 0.0;
  //declare the iterator for the map smiles_energy
  std::map<std::string, double>::iterator iter;
//als allererstes inhalt der map in einen stringstream schreiben, da sonst obenergy nichts berechnen kann
//stringstream input_smile zum speichern und lesen des uebergebenen parameters fuer die Energieberechnung(smilestrings); "fuellen" mit den SMILES-Strings (also keys) der map.
  

  std::stringstream input_smiles;
  
  for (SMILES_container::const_iterator it = probe_smiles.begin(); it != probe_smiles.end(); ++it)
	{
	  input_smiles << it->first <<std::endl; 
	}
	

//2. stringstream fuer output_smiles
  
  std::stringstream output_smiles;
  //std::ofstream file_save_smiles_energy("smilesWithEnergy.smi", std::ios::app );//| std::ios::out);

  
 //3. declare map

  std::map<std::string, double> smiles_energy;
  std::string SmilesString;
  double value_energy;

 //4. set up openbabel format converter

 OpenBabel::OBConversion conv (&input_smiles, &output_smiles);

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
      std::cout << "key value pair not present in the map so calculate the energy for this smilestring" << std::endl;
      value_energy = FF->Energy();
       smiles_energy[SmilesString] = value_energy;//FF ->Energy();
    }
    else
    {
      std::cout << "Key-value pair is present in the map" << endl;
      value_energy = iter->second;
      // oder : value_energy = smiles_energy[SmilesString];
    }
    sum += value_energy;    
   mol.Clear();
  }//!!!!!!!

 return (sum);
}
 
                      

