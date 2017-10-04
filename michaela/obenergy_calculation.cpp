//obenergy calculation including file-manipulation of the input-file

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <utility>

int main(int argc, char **argv)
{
//1. open the inputfile
  std::ifstream input_file("mols.smi");
  

//2. open output_file for writing the canonical output
  std::ofstream Output_SMILES("Output_SMILES.can");
  std::ofstream file_save_smiles_energy("smilesWithEnergy.smi", std::ios::app );//| std::ios::out);

  //open the manipulated_inputfile for writing
  std::ifstream file_manipulated("Output_SMILES.can");
  
 //3. declare map

  std::map<std::string, double> smiles_energy;
  std::string SmilesString;
  double value_energy;

  

  //4. if files can not be open
 if (!input_file.is_open())
     {
       std::cerr << "Cannot open the file" << std::endl;
     }

 
 if (! Output_SMILES.is_open())
 {
   std::cerr << "Cannot open the output file" << std::endl;
 }


 //5. set up openbabel format converter

       OpenBabel::OBConversion conv(&input_file,&Output_SMILES);
       
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

 while (file_manipulated >> SmilesString)
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
   //looks if the SmileString is in the map already or not
   //save the smilesstring and the value in a file
   if(smiles_energy.find(SmilesString) == smiles_energy.end())
     {
       value_energy = FF->Energy();
       smiles_energy[SmilesString] = value_energy;
      file_save_smiles_energy <<SmilesString << " " <<FF->Energy()<<std::endl;
     }
   mol.Clear();

  }
 //ausgabe zur Kontrolle

 for (std::map<std::string, double>::iterator it = smiles_energy.begin(); it!=smiles_energy.end();++it)
   {
     std::cout << "key found: " << " " << it -> first << " " << "associated value" << it -> second << std::endl;
   }

 file_manipulated.close();
 file_save_smiles_energy.close();
 Output_SMILES.close();
 input_file.close();
  return 0;
}
