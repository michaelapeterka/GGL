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

  //open the inputfile
  std::ifstream input_file("molecuel_try.smi");
  std::string smiles;

//if file can not be open
if (!input_file.is_open())
     {

       std::cerr << "Cannot open the file" << std::endl;
     }
       
//open output_file for writing 
 std::ofstream Output_SMILES("Output_SMILES.can",std::ios::app | std::ios::out);

       if (! Output_SMILES.is_open())
	 {
	   std::cerr << "Cannot open the output file" << std::endl;
	 }

//set up openbabel format converter

       OpenBabel::OBConversion conv(&input_file,&Output_SMILES);
       
       if(!conv.SetInAndOutFormats("smi","can"))
	 {
	   std::cerr << "Cannot open the output file" << std::endl;
	 }

       // molecules are converted  !!!!!!!!!!     
       int n = conv.Convert();

       std::cout << n << "molecules converted"<< std::endl;
       //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       //!!!!calculation energy

       //open the manipulated file

std::ifstream file_manipulated("Output_SMILES.can");
std::string key;
double value;
std::map<std::string, double> smiles_energy;
std::string line;

 if (!file_manipulated.is_open())
   {
     std::cerr << "Cannot open the file" << std::endl;
   }

 //set up the openbabel OPMol

 OpenBabel::OBMol mol;

 //declaration forcefield

 OpenBabel::OBForceField *FF = OpenBabel::OBForceField::FindType("MMFF94");
//not get the forcefield
if (!FF)
  {
std::cerr << "Cannot find the forcefield" << std::endl << std::endl;
}


//read the file line by line and setup forcefield
//open file to write in the SMILES-String and the Energy

 std::ofstream SMILES_Energy("SMILES_Energy");

 
 while (file_manipulated >> line)
  {
   if (!conv.ReadString(&mol, line))
    {
      std::cout << "Cannot read SMILES" << line << std::endl;
      return 1;
    }

   if (!FF->Setup(mol)) //setup forcefield for mol
    {
     std::cout << "Cannot set up forcefield" << std::endl;
     return 1;
    }

   //print the energy and the unit for the molecules

   //std::cout << line << " : " << FF->Energy() << FF->GetUnit() << std::endl;

   SMILES_Energy <<line << " " <<FF->Energy()<<std::endl;
   mol.Clear();

  }

 //open for reading the file and to insert the values to the map

 std::ifstream file_SMILES_Energy("SMILES_Energy");
 
 if (!file_SMILES_Energy.is_open())
   {
     std::cout << "Cannot open the File" << line << std::endl;
   }


 //!!schaue ob das element schon in der map ist

 std::map<std::string, double>::iterator it;
 typedef std::map<std::string, double>::iterator it_type;
 std::pair<it_type,bool> controll;
 
 
 //!!!!!!!!!!!!
 //insert the values of the file to the map
 while (file_SMILES_Energy >> key >> value)
   {
        controll=smiles_energy.insert(make_pair(key,value));

    if (controll.second)
     {
       it = controll.first;  //insertion succeeded;save iterator
     }
     // else
    // {
    //    std::cout<<"already exist in map"<< std::endl;
    //   }
   
     //smiles_energy[key]=value; 
    // smiles_energy.insert(std::pair<std::string, double>(key,value));//make_pair(key,value));
     
   }

 
 //Ausgabe zur Kontrolle

 for (std::map<std::string, double>::iterator it = smiles_energy.begin();it !=smiles_energy.end(); ++it){
   std::cout << it->first << " " << it ->second << std::endl;
 }
 //Ausgabe Kontrolle Ende
 
 SMILES_Energy.close();
 Output_SMILES.close();
 input_file.close();
       
	 return 0;

	 
}
       
