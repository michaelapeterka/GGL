//modiefied code from openbabel.org/api/2.3.0/obforcefield_energy_8cpp-example.shtml

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/shared_ptr.h>
#include <openbabel/forcefield.h>

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace OpenBabel;

//shared_ptr = shared pointer; vorinstallierte  in openbabel

//********Function to read a molecule from the file

std::shared_ptr<OBMol> GetMol(const std::string &filename) 
{
  //1. create the OBMol object

  std::shared_ptr<OBMol> mol(new OBMol);

    //2. Create the OBConerversion object

  OBConversion conv;
    
    //base class for fileformats
    //FormatFromExt -> searches registered foramts for an ID the same as the file extension
  OBFormat *format = conv.FormatFromExt(filename.c_str());//

    if(!format || !conv.SetInFormat(format))//sets the input format from an id 
    {
      std::cerr << "Can not find file input format" << filename << std::endl;
      return mol;
    }

   //3.Open the file

	std::ifstream file(filename.c_str());

	//look if it can be open or not

	if (!file)
	  {
	    std::cout << "Can not open the file for reading" << std::endl;
	    return mol;
	  }

   //4. Read the molecule (or not)

	if (!conv.Read(mol.get(), &file))//reads ab object of a class derived from OBBase -- the input stream must be specified and 
	  {
	    std::cout << "Molecule from file can not be read" << std::endl;
	    return mol;
	  }
	return mol;
	
}
//********End Function read molecule from the file

//********Beginn Main function

  int main(int argc, char **argv)
  {//under construction!!!! not ready yet!!!!!!check the number of parameters and tell the user how to run the program
    //not happy with this solution!!!!
    if (argc < 2)
    {
      //std::ifstream file;
      //std::string line;
	
	
      //std::cout << "Enter filename: " << std::endl;
	std::cout << "Enter filename: " << argv[0] << "filename" << std::endl;
	return 1;
    }
	/*file.open(line.c_str());

	if (file.is_open())

	  {
	    while (std::getline(file,line))
	      {std::cout << "molecule: " << line << std::endl;
	      }
	      file.close();
	  }
	    else
	      {
		std::cerr << "Cannot open the file" << std::endl;
	      }
	return 0;
		//}

		//Create the OBConversion Object

		OBConversion conv;

		//BaseClass for fileformats

		OBFormat *format = conv.FormatFromExt(line.c_str());

		if (!format || !conv.SetInFormat(format))
		  {

		  }*/

    //Read file

    std::shared_ptr<OBMol> mol = GetMol(argv[1]);

    //now get the forcefield!!!
    OBForceField *FF = OBForceField::FindType("MMFF94");//oder MMFF94etc.
    if (!FF)
      {
	std::cout << "Can not find the forcefield" << std::endl;
	return 1;
      }

    //setup forcefield

    if(!FF->Setup(*mol))// set up the force field for mol
      {
	std::cout << "Can not set up forcefield"<< std::endl;
	return 1;
      }

    //Print the Energy and the unit for the molecule
    //Change: zuerst das molekuel, fuer das die energie berechnet wird ausgeben, dann erst die dazugehoerige energie
    
    std::cout << "Molecule total energy: " << FF->Energy() << " " << FF->GetUnit()<< std::endl;

    return 0;

  }
  
//********End Main function
