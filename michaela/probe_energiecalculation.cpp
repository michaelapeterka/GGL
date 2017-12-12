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
#include <map>

double energie_calculation(std::map<std::string, double>& probe_smiles)
{


//als allererstes inhalt der map in einen stringstream schreiben, da sonst obenergy nichts berechnen kann
//stringstream input_smile zum speichern und lesen des uebergebenen parameters fuer die Energieberechnung(smilestrings); "fuellen" mit den SMILES-Strings (also keys) der map.
  

  std::stringstream input_smiles;
  
  for (std::map<std::string, double>::iterator it = probe_smiles.begin(); it != probe_smiles.end(); ++it)
	{
	  input_smiles << it->first <<std::endl; 
	}
	

//2. stringstream fuer output_smiles
  
  std::stringstream output_smiles;
  std::ofstream file_save_smiles_energy("smilesWithEnergy.smi", std::ios::app );//| std::ios::out);

  
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

   if(smiles_energy.find(SmilesString) == smiles_energy.end())//if the string is not in the map
     {
      smiles_energy[SmilesString] = FF ->Energy();
      file_save_smiles_energy <<SmilesString << " " <<FF->Energy()<<std::endl;
     }
   //hier fehlt noch das der string genommen wird und was er dann macht, wenn er schon in der map ist. ansonsten nimm den wert aus der map
   mol.Clear();

  }
 // looking for the highest value in the map and put it out on the screen just to check

double max_energy = 0.0;
std::string max_string;

 for (std::map<std::string, double>::iterator it = smiles_energy.begin(); it!=smiles_energy.end();++it)
   {
     //std::cout << "key found: " << " " << it -> first << " " << "associated value" << it -> second << std::endl;
     if (it -> second > max_energy)
     {
       max_energy = it -> second;
       max_string = it -> first;
     }
   }

 //std::cout << "Highest Energy: Smiles: " << max_string << "\n Energy: " << max_energy << std::endl;
 
 
 file_save_smiles_energy.close();
 
                      
return max_energy;
}
