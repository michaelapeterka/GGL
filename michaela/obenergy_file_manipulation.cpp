//obenergy calculation including file-manipulation of the input-file

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

int main(int argc, char **argv)
{

  //open the inputfile
  std::ifstream file("HCN5.smi");
  std::string smiles;

//if file can not be open
if (!file.is_open())
     {

       std::cerr << "Cannot open the file" << std::endl;
     }
       
     else{
       //!!!!!
//open output_file for writing 
       std::ofstream Output_SMILES("Output_SMILES.smi");

      
       //!!!!! reading input file line per line 
       while (file>>smiles)
	 { //searching in std::string smiles for not needed characters
 	  for (int i = 0; i<smiles.size();)
 	    {
 	      if (smiles[i] == ':' ||  smiles[i] == '0' || smiles[i] == '1' || smiles[i] == '2' || smiles[i] == '3' || smiles[i] == '4' || smiles[i] == '5' || smiles[i] == '6' || smiles[i] == '7' || smiles[i] == '8' || smiles[i] == '9')
 		{
                 smiles.erase(i,1); // erase the not wanted character at position i otherwise move on to the next position
 		}
	      else
		{
 		  i++;
		}
            }//write the manipulated/changed strings into the new file called output. 
//std::cout << smiles << std::endl<< std::endl;
	      Output_SMILES << smiles << std::endl;
	}
       Output_SMILES.close();
         }


//!!!!zusammenfuehren von der manipulation der

//open the manipulated file

 std::ifstream file_manipulated("Output_SMILES.smi");

 if (!file_manipulated.is_open())
   {
     std::cerr << "Cannot open the file" << std::endl;
   }

 std::string line;

 //setup the openbabel OPMOL

 OpenBabel::OBMol mol;

//set up the openbabel format converter

OpenBabel::OBConversion obconv;
obconv.SetInFormat("smi");

//now get the forcefield

OpenBabel::OBForceField *FF = OpenBabel::OBForceField::FindType("MMFF94");

if (!FF)
  {
std::cerr << "Cannot find the forcefield" << std::endl << std::endl;
}

//read the file line by line and setup forcefield

//std::cout << "Filename" << file << std::endl;

while (std::getline(file_manipulated,line))
  {
   if (!obconv.ReadString(&mol, line))
    {
      std::cout << "Cannot read SMILES" << line << std::endl;
      return 1;
    }

    if (!FF->Setup(mol)) //setup forcefield for mol
    {
     std::cout << "Cannot set up forcefield" << std::endl;
     return 1;
    }

//Print the Energy and the unit for the molecule

std::cout << line << std::endl;
std::cout << "Molecule total energy: " << FF->Energy() << " " << FF->GetUnit() << std::endl;

mol.Clear();
}
      file_manipulated.close();
       file.close();
  
   return 0;
 }

//zusammenfuehrun ENDE !!!!!!
