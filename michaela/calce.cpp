/* Last changed Time-stamp: <2017-09-07 01:29:04 xtof> */

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>

#include <fstream>
#include <string>

using namespace OpenBabel;

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Usage " << argv[0] << " filename" << std::endl;
    return 1;
  }

  // open file with one SMILES per line
  std::ifstream file (argv[1]);
  if (!file) {
    std::cerr << "Cannot open the file " << argv[1] << std::endl;
    return 1;
  }

  std::string line;
  OBMol mol;
  // setup an openbabel format converter
  OBConversion obconv;
  obconv.SetInFormat("smi");

  // setup a force field
  OBForceField *FF = OBForceField::FindType("MMFF94");
  if (!FF) {
    std::cout << "Can not find the forcefield" << std::endl;
    return 1;
  }

  // read file line by line
  while (std::getline(file, line)) {

    if (!obconv.ReadString(&mol, line)) {
      std::cerr << "Could not read SMILES " << line << std::endl;
      return 1;
    }

    if (!FF->Setup(mol)) {
      std::cout << "Can not set up forcefield" << std::endl;
      return 1;
    }

    std::cout
      << line << " " << FF->Energy() << " " << FF->GetUnit() << std::endl;

    mol.Clear();
  }
    
  file.close();

  return 0;
}
