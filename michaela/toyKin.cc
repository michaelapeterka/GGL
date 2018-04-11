

#include "toyChemUtil.hh"

#include <fstream>
#include <sstream>
#include <utility>
#include <algorithm>
#include <iterator>
#include <cstring>

#include <boost/lexical_cast.hpp>


#include "biu/OptionParser.hh"


#include "sgm/SGM_vf2.hh"
#include "sgm/Graph_boost.hh"
#include "sgm/PA_OrderCheck.hh"


#include "ggl/MR_ApplyRule.hh"
#include "ggl/Graph_GML_writer.hh"
#include "ggl/GS_stream.hh"

#include <ggl/chem/SMILESparser.hh>
#include "ggl/chem/GS_SMILES.hh"
#include "ggl/chem/SMILESwriter.hh"

#include "ggl/chem/Reaction.hh"
#include "ggl/chem/MR_Reactions.hh"
#include "ggl/chem/RRC_ArrheniusLaw.hh"
#include "ggl/chem/EC_MoleculeDecomposition.hh"
#include "RRC_QuantumMechanics.hh"

#include "ggl/chem/AP_NSPDK.hh"
#include "ggl/chem/AP_disabled.hh"

#include "version.hh"

//!!inkludieren der header datei fuer die Berechnung der Energie
#include "probe_energiecalculation.h"

//////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <ctime>
#include <cmath>

size_t
getRandomNumber( const size_t maxExluding )
{
	// get random number in float interval [0,1] exluding the upper bound 1
	double range01 = (double)rand() / (1.0 + (double)RAND_MAX);

	// get value within integer range < maxExcluding
	return (size_t) floor( ( range01 * (double)maxExluding )) ;
}

//Aenderung Gillespie Beginn
//!!!Funktionsdefinition fuer die erstellung der allgemeinen smiles
//gibt den allgemein smile zurueck ohne labelung um die anzahl der molecule population zu definieren zurueck

std::string extract_generic_smile(string metabolite)
{
  bool in_brackets = false;
  bool after_colon = false;
  string generic_smile;

  for (int i = 0; i < metabolite.size(); ++i)
    {
      char c = metabolite[i];

      if(c =='[')
	{
	  in_brackets = true;
	  continue;
	}
      if(c ==']')
	{
	  in_brackets = false;
	  after_colon = false;
	  continue;
	}
      if(in_brackets && c == ':')
	{
	  after_colon = true;
	  continue;
	}
      if (after_colon)
	{
	  continue;
	}

      /*if(c =='#')
	{
	  in_brackets = false;
	  after_colon = false;
	  continue;
	  }*/

      generic_smile.push_back(c);
    }
 return generic_smile;
}

double partial_sum(vector<double>& a_my, double index)
{
	double sum = 0.0;
	for (int i = 0; i <= index; ++i)
	{
	 sum += a_my.at(i);
	}
	return sum;
}
//Änderung fuer Gillespie ENDE

////////////////////////////////////////////////////////////////////////////////
#if HAVE_UNORDERED_MAP > 0
#include <unordered_map>
#define USED_HASH_MAP_DESCRIPTION "std::unordered_map"

#elif HAVE_TR1_UNORDERED_MAP > 0
#include <tr1/unordered_map>
#define USED_HASH_MAP_DESCRIPTION "std::tr1::unordered_map"

#elif HAVE_GNU_HASH_MAP > 0
#include <ext/hash_map>
#define USED_HASH_MAP_DESCRIPTION "__gnu_cxx::hash_map"

#else
#include <map>
#define USED_HASH_MAP_DESCRIPTION "std::map"

#endif


//////////////////////////////////////////////////////////////////////////

bool
isWithoutAtomClass( sgm::Graph_Interface & graph );

bool
allWithAtomClass( sgm::Graph_Interface & graph );

size_t
extractAtomClass( ggl::chem::Molecule_Graph & mol, std::set<size_t> usedAtomClasses );

size_t
setAtomClass( ggl::chem::Molecule & mol, const std::set<size_t> usedAtomClasses, const size_t nextClassLabelOption );

//////////////////////////////////////////////////////////////////////////

void
initAllowedArguments( biu::OptionMap & allowedArgs, std::string &infoText );  


//////////////////////////////////////////////////////////////////////////

void
printMoleculeGraph( std::ostream& out, ggl::chem::Molecule & m );

//////////////////////////////////////////////////////////////////////////

void
printMoleculeGraphs( std::ostream& out, SMILES_container & smiles );

//////////////////////////////////////////////////////////////////////////

void
printMoleculeGraphsGML( std::ostream& out, SMILES_container & smiles, const bool withSpaces );

//////////////////////////////////////////////////////////////////////////

void
giveRuleGMLExample();

//////////////////////////////////////////////////////////////////////////

void
giveGraphGMLExample();

//////////////////////////////////////////////////////////////////////////

CompareStringPointer compareStringPointer;

//////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv ) {
	using namespace std;
	using namespace ggl;
	using namespace chem;
	
	//////////////////////////////////////////////////////////////
	// data to fill
	//////////////////////////////////////////////////////////////
	int exitValue = 0;
	
	enum InfoMode {OUT_SILENT, OUT_NORMAL, OUT_VERBOSE};
	InfoMode infoMode = OUT_NORMAL;
	
	std::ostream* out = &std::cout;
	std::ofstream* outFile = NULL;
	
	ReactionRateCalculation*  rateCalc = NULL;
	
	AromaticityPerception * aromaticityPrediction = NULL;

	size_t iterations = 0;
	 // rule pattern for each number of connected components
	RulePatternMap rulePattern;
	  // the supported molecule groups that can be abbreviated in rule and
	  // molecule input
	ggl::chem::GroupMap moleculeGroups;
	
	// central idea of the implementation !! dont change !!!
	const bool ignoreAtomClass = true;

	
	//////////////////////////////////////////////////////////////
	// parameter parsing and checking
	//////////////////////////////////////////////////////////////
	
	biu::OptionMap allowedArgs;	//< map of allowed arguments
	std::string infoText;		//< info string of the program
	initAllowedArguments(allowedArgs,infoText);	// init
	
		// parse programm arguments	
	biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc, 
													argv, infoText);
		// arguments parseable and all mandatory arguments given
	if (opts.argExist("help")) {
		opts.coutUsage();
		return 0;
	}
	if (opts.argExist("ruleExample")) {
		giveRuleGMLExample();
		return 0;
	}
	if (opts.argExist("graphExample")) {
		giveGraphGMLExample();
		return 0;
	}
	if (opts.argExist("version")) {
		giveVersion();
		std::cout	<<"\n used SMILES_container : " <<USED_HASH_MAP_DESCRIPTION <<"\n"
					<<std::endl;
		return 0;
	}
	if (opts.noErrors()) {
		  // rules are officially not mandatory.. we have to request manually
		if (!opts.argExist("rules") || opts.getStrVal("rules").size() == 0 ) {
			std::cerr <<"\n\n\tERROR : needed argument not given : 'rules' check usage\n\n";
			return -1;
		}
		  // molecules are officially not mandatory.. we have to request manually
		if (	(!opts.argExist("smiles") || opts.getStrVal("smiles").size()==0) 
			&&	(!opts.argExist("mols") || opts.getStrVal("mols").size()==0) ) 
		{
			std::cerr <<"\n\n\tERROR : needed argument not given :"
				<<" neither 'smiles' nor 'mols' is present, check usage\n\n";
			return -1;
		}
		if (opts.getBoolVal("v")) {
			infoMode = OUT_VERBOSE;
		}
	} else {
		return -1;
	}

	SMILES_container c1;
	SMILES_container c2;
	SMILES_container& targetSmiles = c1;
	SMILES_container& producedSmiles = c2;
	
	try { 


		//////////////////////////////////////////////////////////////
		// set up iteration parameters
		//////////////////////////////////////////////////////////////
		
		if (opts.getIntVal("iter") < 1) {
			std::ostringstream oss;
			oss	<<"number of rule application iterations (" 
			<<opts.getIntVal("iter") <<") has to be at least 1";
			throw ArgException(oss.str());
		}
		iterations = (size_t)opts.getIntVal("iter");
		
		
		//////////////////////////////////////////////////////////////
		// set up random number generator
		//////////////////////////////////////////////////////////////

		if (opts.getIntVal("srand") < 0) {
			// initialize with system time
			const unsigned int randSeed = time(NULL);
			srand( randSeed );
			std::cerr <<"# seed init = "<<randSeed<<std::endl;
		} else {
			// initialize with provided seed value
			srand( (unsigned int) opts.getIntVal("srand") );
		}

		//////////////////////////////////////////////////////////////
		// check for multiple use of STDIN as input parameter value
		//////////////////////////////////////////////////////////////
		
		{
			int conflicts = 1;
			conflicts *= (opts.getStrVal("rules").compare("STDIN")==0?2:1);
			conflicts *= (opts.getStrVal("smiles").compare("STDIN")==0?3:1);
			conflicts *= (opts.getStrVal("rules").compare("STDIN")==0?5:1);
			
			if (conflicts > 5) {
				std::ostringstream oss;
				oss <<"cannot read ";
				switch (conflicts) {
					case 6: oss <<"RULES and SMILES"; break;
					case 10: oss <<"RULES and MOLECULES"; break;
					case 15: oss <<"SMILES and MOLECULES"; break;
					case 30: oss <<"RULES, SMILES and MOLECULES"; break;
				}
				oss <<" from STDIN!";
				throw ArgException(oss.str());
			}
		}
		
		//////////////////////////////////////////////////////////////
		// set up output streams
		//////////////////////////////////////////////////////////////
		
		  // set output stream
		if (opts.getStrVal("out").size() == 0) {
			throw ArgException("no output file given");
		} else if (opts.getStrVal("out").compare("STDOUT") != 0) {
			outFile = new std::ofstream(	opts.getStrVal("out").c_str()
											, std::ofstream::out );
			if (!outFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open output file '" <<opts.getStrVal("out") <<"'";
				throw ArgException(oss.str());
			}
			out = outFile;
		}


		//////////////////////////////////////////////////////////////
		// setup aromaticity perception model
		//////////////////////////////////////////////////////////////

		if (infoMode == OUT_VERBOSE) {
			(*out) <<"\n # LOAD AROMATICITY MODELS ..."; out->flush();
		}
		switch(opts.getCharVal("aromaticity")) {
		case 'o' :
		case 'O' :
			aromaticityPrediction = new AP_NSPDK("OpenBabel:2013");
			assert( ((AP_NSPDK*)aromaticityPrediction)->getModel() != NULL );
			break;
		case 'm' :
		case 'M' :
			aromaticityPrediction = new AP_NSPDK("Marvin:general:2013");
			assert( ((AP_NSPDK*)aromaticityPrediction)->getModel() != NULL );
			break;
		case 'n' :
		case 'N' :
			aromaticityPrediction = new AP_disabled();
			break;
		default:
				std::ostringstream oss;
				oss	<<"aromaticity perception type '" <<opts.getCharVal("aromaticity") <<"'"
					<<" is not supported";
				throw ArgException(oss.str());
		}
		if (infoMode == OUT_VERBOSE) {
			(*out) <<" DONE\n"; out->flush();
		}


		
		//////////////////////////////////////////////////////////////
		// parse groups from input
		//////////////////////////////////////////////////////////////


		if (opts.argExist("groups")) {
			if (infoMode == OUT_VERBOSE) {
				(*out) <<"\n # PARSE GROUPS ..."; out->flush();
			}
			  // check if input specifications are compatible
			if ( opts.argExist("smiles") && opts.getStrVal("groups") == opts.getStrVal("smiles") ) {
				std::ostringstream oss;
				oss	<<"Cannot read both groups and smiles from '" <<opts.getStrVal("groups") <<"'";
				throw ArgException(oss.str());
			}
			if ( opts.argExist("mols") && opts.getStrVal("groups") == opts.getStrVal("mols") ) {
				std::ostringstream oss;
				oss	<<"Cannot read both groups and mols from '" <<opts.getStrVal("groups") <<"'";
				throw ArgException(oss.str());
			}
			if ( opts.argExist("rules") && opts.getStrVal("groups") == opts.getStrVal("rules") ) {
				std::ostringstream oss;
				oss	<<"Cannot read both groups and rules from '" <<opts.getStrVal("groups") <<"'";
				throw ArgException(oss.str());
			}

			  // parse molecule group definitions that can be abbreviated in
			  // rules and molecules
			parseGroups( opts.getStrVal("groups"), moleculeGroups );

			if (infoMode == OUT_VERBOSE) {
				(*out) <<" DONE\n"; out->flush();
			}
		}

		//////////////////////////////////////////////////////////////
		// parse Rules from input
		//////////////////////////////////////////////////////////////
		
		if (infoMode == OUT_VERBOSE) {
			(*out) <<"\n # PARSE RULES ..."; out->flush();
		}
		std::vector<ggl::chem::ChemRule> rules;
		  // parse all rules for given input
		parseRules( opts.getStrVal("rules"), rules, moleculeGroups );

		if (infoMode == OUT_VERBOSE) {
			(*out) <<" DONE\n"; out->flush();
		}
		
		if (rules.empty()) {
			(*out)	<<"\n PROBLEM : no rules found in given input!\n"
						<<std::endl;
			return -1;
		}
		
		if (infoMode == OUT_VERBOSE) {
			(*out) <<"\n # CHECK RULES ..."; out->flush();
		}
		{
			  // check if all rules are valid
			bool allRulesOK = true;
			// check if rules do not contain atom class information (not allowed)
			for (size_t r=0; r<rules.size(); ++r) {
				// check if nodes of leftside pattern are without atom class
				ggl::chem::LeftSidePattern leftRule( rules.at(r) );
				allRulesOK = isWithoutAtomClass( leftRule );
				// check if node label contains class separator (not allowed)
				if ( ! allRulesOK ) {
					(*out) <<"\n PROBLEM : rule " <<(r+1) <<" '"
							<<rules.at(r).getID()
							<<"' : left side contains an atom class label separator ':', which is not allowed\n";
					break;
				}
				// check if nodes of rightside pattern are without atom class
				ggl::chem::RightSidePattern rightRule( rules.at(r) );
				allRulesOK = isWithoutAtomClass( rightRule );
				// check if node label contains class separator (not allowed)
				if ( ! allRulesOK ) {
					(*out) <<"\n PROBLEM : rule " <<(r+1) <<" '"
							<<rules.at(r).getID()
							<<"' : right side contains an atom class label separator ':', which is not allowed\n";
					break;
				}
			}
			// additional chemical sanity checks of rules
			if (allRulesOK && !opts.argExist("noRuleCheck")) {
				for (size_t r=0; r<rules.size(); ++r) {
					size_t conStatus = rules.at(r).isConsistent();
					if( conStatus != ggl::chem::ChemRule::C_Consistent ) {
						allRulesOK = false;
						(*out) <<"\n PROBLEM : rule " <<(r+1) <<" '"
								<<rules.at(r).getID() <<"' is not chemically correct"
									" or contains unsupported properties:\n";
						rules.at(r).decodeConsistencyStatus( conStatus, (*out) );
					}
				}
			}
			if (!allRulesOK) {
				return -1;
			}
		}
		if (infoMode == OUT_VERBOSE) {
			(*out) <<" DONE\n"; out->flush();
		}
		
		
		  // generate left side pattern of each rule needed for its application
		  // and store separate rule lists based on the component number of
		  // their left side pattern
		for (size_t r=0; r<rules.size(); ++r) {
			ggl::chem::LeftSidePattern* pattern = new ggl::chem::LeftSidePattern(rules[r]);
			  // store in the pattern list according to the component number
			size_t compNumber = pattern->getFirstOfEachComponent().size();
			rulePattern[compNumber].push_back( pattern );
		}
		
		//////////////////////////////////////////////////////////////
		// parse molecules from input
		//////////////////////////////////////////////////////////////
		
		if (infoMode == OUT_VERBOSE) {
			(*out) <<"\n # PARSE INPUT ..."; out->flush();
		}
		size_t nextAtomClass = 1;
		if (opts.argExist("smiles")) {
			// parse molecules from SMILES and add atom numbering
			nextAtomClass = parseSMILES( opts.getStrVal("smiles"), targetSmiles, moleculeGroups, nextAtomClass );
		}
		if (opts.argExist("mols")) {
			// parse molecules from GML and add atom numbering
			nextAtomClass = parseMolGML( opts.getStrVal("mols"), targetSmiles, moleculeGroups, nextAtomClass );
		}
		if (infoMode == OUT_VERBOSE) {
			(*out) <<" DONE\n"; out->flush();
		}
		
		if (!opts.argExist("noInputCorrection")) {
			if (infoMode == OUT_VERBOSE) {
				(*out) <<"\n # CORRECT INPUT ..."; out->flush();
			}
			  // correct aromaticity and adjacent protons of input molecules
			nextAtomClass = correctInputMolecules( targetSmiles, producedSmiles, aromaticityPrediction, true, nextAtomClass );
			  // clear input molecules
			for (SMILES_container::iterator it=targetSmiles.begin(); it!=targetSmiles.end() ; ++it ) {
				delete it->second;
			}
			targetSmiles.clear();
			  // setup target molecules
			std::swap( producedSmiles, targetSmiles );
			if (infoMode == OUT_VERBOSE) {
				(*out) <<" DONE\n"; out->flush();
			}
		}

		if (infoMode == OUT_VERBOSE) {
			(*out) <<"\n # CHECK INPUT ..."; out->flush();
		}
		{
			  // check if all rules are valid
			bool allMolOK = true;

			// set of all atom classes already present within the input
//			std::set<size_t> usedAtomClasses;

			//////////////////////////////////////////////////////////////
			// get used atom class labels and check for redundancy
			//////////////////////////////////////////////////////////////

//			// extract used atom class information within all molecules
//			for (SMILES_container::iterator it=targetSmiles.begin(); allMolOK && it!=targetSmiles.end() ; ++it ) {
//				ggl::chem::Molecule_Graph molGraph( *(it->second) );
//				size_t redundantLabel = extractAtomClass( molGraph, usedAtomClasses );
//				if (redundantLabel > 0) {
//					allMolOK = false;
//					(*out) <<"\n PROBLEM : molecule '"
//							<<it->first <<"' contains the redundant class label '"
//							<<redundantLabel<<"', which is not allowed\n";
//				}
//			}

			if (allMolOK) {

				//////////////////////////////////////////////////////////////
				// number atoms without class label (excluding present labels)
				//////////////////////////////////////////////////////////////

				// extract used atom class information within all molecules
//				size_t nextClassLabel = 1;
//				for (SMILES_container::iterator it=targetSmiles.begin(); allMolOK && it!=targetSmiles.end() ; ++it ) {
//					nextClassLabel = setAtomClass( *(it->second), usedAtomClasses, nextClassLabel );
//				}

				//////////////////////////////////////////////////////////////
				// ensure all atoms are numbered within class label
				//////////////////////////////////////////////////////////////

				for (SMILES_container::iterator it=targetSmiles.begin(); allMolOK && it!=targetSmiles.end() ; ++it ) {
					// check if molecule is without atom class information
					ggl::chem::Molecule_Graph molGraph( *(it->second) );
					allMolOK = allWithAtomClass( molGraph );
					if (!allMolOK) {
						(*out) <<"\n PROBLEM : molecule '"
								<<it->first <<"' misses some atom classes, which is not allowed\n";
					}
				}
			}

			if (allMolOK && !opts.argExist("noInputCheck")) {
				for (SMILES_container::iterator it=targetSmiles.begin(); it!=targetSmiles.end() ; ++it ) {
					size_t conStatus = ggl::chem::MoleculeUtil::isConsistent( *(it->second) );
					if( conStatus != ggl::chem::MoleculeUtil::C_Consistent ) {
						allMolOK = false;
						(*out) <<"\n PROBLEM : molecule '"
								<<it->first <<"' is not chemically correct"
									" or contains unsupported properties:\n";
						ggl::chem::MoleculeUtil::decodeConsistencyStatus( conStatus, (*out) );
					}
				}
			}
			if (!allMolOK) {
				return -1;
			}
		}
		if (infoMode == OUT_VERBOSE) {
			(*out) <<" DONE\n"; out->flush();
		}


		//////////////////////////////////////////////////////////////
		// print to stream if in verbose mode
		//////////////////////////////////////////////////////////////

		if (infoMode == OUT_VERBOSE) {

			  // print groups
			(*out) <<"\n ######## PARSED AND AVAILABLE GROUPS #######\n";
			printGroups( *out, moleculeGroups );

			(*out) <<"\n ######## PARSED RULES #######\n";
			printRules( *out, rules );
			
			(*out) <<"\n ###### PARSED MOLECULES #####\n\n";
			printSMILES( *out, targetSmiles );
			(*out) <<std::endl; out->flush();
		}
		
		//////////////////////////////////////////////////////////////
		// perform rule application iterations
		//////////////////////////////////////////////////////////////

		//!!!Definition fuer Energy_calculation produceSmiles/targetSmiles
		double energy_producedSmiles = 0.0;
		double energy_targetSmiles = 0.0;

		//metabolites_in_producedReactions
		std::map<string,int> metabolites_in_Reactions;

		//how many metabolites in produced REcations
		std::map<string,int> value;

		//!!!Ende
		
		  // set up graph matcher
		sgm::SGM_vf2 sgm;
		
		  // the list of all possible reaction instances
		MR_Reactions::Reaction_Container producedReactions;

			
		  // make target molecules virtuelly the ones produced in the last iteration
		std::swap(producedSmiles, targetSmiles);

		
		  // perform iterations
		for (size_t it=0; it<iterations; ++it ) {

			// generate reactions possible for new molecules from the last iteration
			if (producedSmiles.size() > 0) {

//				std::cout <<"## new mols : \n";
//				for (SMILES_container::const_iterator m=producedSmiles.begin(); m!=producedSmiles.end(); m++) {
//					std::cout <<"\t\t\t"<<m->first<<" "<<m->second<<std::endl;
//				}


				// all molecules that can be produced by a reaction involving a molecule from producedSmiles
				SMILES_container toFill;
				// apply all rules but ensure that always at least one molecule from producedSmiles is involved
				applyRules(rulePattern
						, targetSmiles, producedSmiles
						, toFill, producedReactions
						, rateCalc
						, true
						, *aromaticityPrediction
						, ignoreAtomClass
						, true /* enforceUniqueAtomMatch */
						);



				//!!!versuch calculation energy produced smiles 
				energy_producedSmiles = energie_calculation(producedSmiles);
					cout << "deltaEnergy = " << energy_producedSmiles << endl;
				//!!Versuch Ende
					
				// merge targetSmiles and producedSmiles into targetSmiles
				targetSmiles.insert(producedSmiles.begin(), producedSmiles.end());
				producedSmiles.clear();

				//!!!!Versuch Beginn calculation energy target Smiles
				energy_targetSmiles =  energie_calculation(targetSmiles);
		        	cout << "Delta Energy : " << energy_targetSmiles << endl;
			       //!!Versuch Ende

									    
				// clear reaction product container, since most will not be needed
				for (SMILES_container::iterator it=toFill.begin(); it!=toFill.end() ; ++it ) {
					delete it->second;
				}
				toFill.clear();
			}

			// pick a reaction TODO change accordingly!!

			//fuer die nummerierung der Regeln
			int rule_index = 0;
			std::map<string,int> rule_indexs;
			vector<string> rule_names_ids;
			string rule_id;

			//fuer die Berechnung von h_my
			vector<double>h_my;
			double h_my_multiplicator = 1.0;

			//molecular population level -> Anzahl der molekuelpopulationen eruieren (in targetSmile)
                          std::map<string,int> molecular_population_level;

                          for (SMILES_container::const_iterator ts=targetSmiles.begin(); ts!= targetSmiles.end(); ++ts)
                            {
                              std::string generic_smile = extract_generic_smile(ts->first);

                              if(molecular_population_level.find(generic_smile) == molecular_population_level.end())
                                {
                                  molecular_population_level[generic_smile] = 1;
                                }
                                else
                                {
                                    molecular_population_level[generic_smile]++;
                                }

                            }

                          //just to check -> temporary
                           cout << "!!!!!!!!!!!!!!!Anzahl der Smiles:" << endl;

                           for(map<string,int>::iterator im = molecular_population_level.begin();im != molecular_population_level.end();++im)
                             {
                               cout << "Smile: " << im -> first << " how often: " << im -> second << endl;
                             }

			// get number of reaction instances for each rule
			typedef std::map< std::string, size_t > RuleHist;
			RuleHist ruleId2reactions;
			for (MR_Reactions::Reaction_Container::const_iterator r=producedReactions.begin(); r != producedReactions.end(); r++) {
				// check if not seen so far
				if (ruleId2reactions.find(r->rule_id) == ruleId2reactions.end()) {
					ruleId2reactions[r->rule_id] = 1;
					//nummerierung der Regeln!!!!!!!!!!!!!!!!!!1
					rule_indexs[r->rule_id] = rule_index;
					++rule_index;
					rule_names_ids.push_back(r->rule_id);
					//Ende
				} else {
					// increase counter
					ruleId2reactions[r->rule_id]++;
				}
				double h_my_multiplicator = 1.0;
				//um die molecular_population_levels der jeweiligen reaction zuzuordnen
				for (Reaction::Metabolite_Container::const_iterator me = r->metabolites.begin();me != r -> metabolites.end();++me)
				  {
					  string generic_smile = extract_generic_smile(*m);
					  h_my_multiplicator = h_my_multiplicator*molecular_population_level[generic_smile];
					  
				   /* if(metabolites_in_Reactions.find(*me) == metabolites_in_Reactions.end())
				      {
					metabolites_in_Reactions[*me] = 1;
				      }
				    else
				      {
					metabolites_in_Reactions[*me] ++;
				      }
*/
				  }
				//versuch 1: in der Schleife:
				h_my.at(rule_indexs[r->rule_id]) = h_my_multiplicator; 
			}
			 							
			//2. versuch: ohne connection zu rule_id --> dazu muss aber hy_multiplicator ausserhalb der schleife sein. --> tendiere eher zu version 1 anstatt version 2 
			h_my.push_back(h_my_multiplicator);
						
   			//!!!!nur temporaer zur ueberpruefung
			  /*for (std::map<string,int>::iterator iti = metabolites_in_Reactions.begin();iti!=metabolites_in_Reactions.end();++iti)
			    {
			      std::cout << "Metabolite: " << iti -> first << " how often: " << iti -> second << endl;
			    }							     
		       */
			  

			  //!!!!!Beginn Aenderung fuer Gillespie_Algorithmus!!!!!
			  //1. Teil: Initialisierung

			//Anzahl der Reaktionen M
			int M = ruleId2reactions.size();  //std::cout << "ruleID2reactions_groesse_ausgabe : " << M << endl;
			
			//number of reactants N = targetSmiles
			int N = targetSmiles.size();
			    
			//reaction counter n --> probably not necessary because of the for-loop
			int n = 0;

			//maximale Anzahl der Reaktionen --> also probably not necessary because of the for-loop
			int n_max = 100;

			//Zeit t
			double t = 0.0;
			
			//!!!Kalkulation der Energiedifferenz fuer die reaction rates c_mu und die dazugehoerigen Instanzen -> deltaE = sum_energy(producedSmiles) - sum_energy(targetSmiles)
			/*double energy_producedSmiles = energie_calculation(producedSmiles);
			double energy_targetSmiles = energie_calculation(targetSmiles);
			double deltaE = energy_producedSmiles - energy_targetSmiles;*/
			//double deltaE = energie_calculation(producedSmiles) - energie_calculation(targetSmiles);
			double deltaE = energy_producedSmiles - energy_targetSmiles;
			    std::cout<<"Delta E :::::::" << deltaE <<std::endl;
			
		       
			//vector c_mu -> reaction rates (of size M)
			vector<double> c_mu(M);

			//c_mu mit den Energiewerten fuellen
			  for(int i = 0 ; i<M; ++i)
			    {
			      c_mu.at(i) = deltaE;
			    }
			  
		       /* !!!!!nach oben verschoben vor ruleHist!!!!molecular population level -> die Anzahl der unterschiedlichen molekuelpopulationen anschauen
			  std::map<string,int> molecular_population_level;

			  for (SMILES_container::const_iterator ts=targetSmiles.begin(); ts!= targetSmiles.end(); ++ts)
			    {
			      std::string generic_smile = extract_generic_smile(ts->first);

			      if(molecular_population_level.find(generic_smile) == molecular_population_level.end())
				{
				  molecular_population_level[generic_smile] = 1;
				}
				else
				{
				    molecular_population_level[generic_smile] ++;
				}

			    }

			  //just to check -> temporary
			   cout << "!!!!!!!!!!!!!!!Anzahl der Smiles:" << endl;
  
			   for(map<string,int>::iterator im = molecular_population_level.begin();im != molecular_population_level.end();++im)
			     {
			       cout << "Smile: " << im -> first << " how often: " << im -> second << endl;
			     }
       */

		       //2. Teil: Beginn des Algorithmus

			   //Generieren der Zufallszahlen
			   double r0 = double(rand())/(RAND_MAX);
			   double r1 = double(rand())/(RAND_MAX);
			   double r2 = double(rand())/(RAND_MAX);
			   
			   //berechnung von h_my
			   //wurde oben schon gemacht. --> kommt drauf an, wo man die while schleife setzt --> wenn man sie setzt
			   
			 //Berechnung der propensity a_my -> i = 1,....,Anzahl der Reaktionen für jede Reaktion
			   std::vector<double> a_my;

			     for(size_t a = 0; a < h_my.size();++a)
			       {
				 a_my.push_back(h_my(a)*c_my(a));
			       }

		       //Berechnung der Summe ueber a_my

			   double a0 = 0.0;

			   for (size_t i = 0; i < a.size();++j)
			     {
			       a0 = a0 + a.at(i);
			     }

		      //Berechnen der Zeit, wann die naechste Reaktion stattfinden wird.

			   doube tau = 0.0;
			   tau = (1/a0)*log(1/r0);
			   
		     //Berechnen welche Rule genommen/passieren wird (welche Reaktion wird genommen)

			   for (int j=0; j<M;j++)
			     {
			       if(partial_sum(a_my,j-1) < (r2*a0) && (r2*a0) <= partial_sum(a_my,j)
				 {
				   rule = j;
				   break;
				 }
			     }
		    //jetzt wird ermittelt welche Regel zu der ermittelten Zahl gehoert

		    //das zuweisen welche Regel genommen wird, um damit dann weiterarbeite zu koennen und zu schauen, wieviele instances es gibt. 
			  rule_id = rule_names_ids.at(rule);
			  
		    //--> jetzt weiss ich welche Regel genommen wird. 
			       //Jetzt moechte ich die Instanz dazu ermitteln
			       //iterieren ueber producedReactions um die ermittelte rule_id mit der rule_id von producedReactions zu vergleichen.
			       //mittels string comparison
			
			vector<MR_Reactions::Reaction_Container::const_iterator> instance_index;       
			for (MR_Reactions::Reaction_Container::const_iterator pr=producedReactions.begin(); pr != producedReactions.end(); pr++)
			  {
			    if (rule_id.compare(pr->first) == 0) //(Regel die aus dem Gillespie ermittelt wurde))
			      {
				      cout << "ist gleich" << endl;
				      //rausspeichern der einträge auf die pr gerade zeigt!!!!
				      instance_index.push_back(pr);
			      }
			  }
			 
			/*
			  !!!!!!!!!!!
			  Berechnung des Boltzmannfaktors: exp(-deltaE/(kB * T))
			  !!!!!!!!!!!
                        */
			    std::cout << "===========Start calculate Boltzmannfaktor" << std::endl;
			    
			    //Berechnung der Energie der einzelnen Instanzen
			    //da die Energie logischerweise die gleiche ist, wie fuer die genommen Regel kann hier mit der schon berechneten Energie deltaE weitergearbeitet werden --> deltaE
			    
			    //3. kBT in kcal/mol at room temperature (298°K)
			    double kBT;
			    kBT = 0.5924;

			    //4.Berechnung eines einzelnen Boltzmanfaktors BF
			    double BF;
			    BF =exp(deltaE/kBT);
			    
			    //Aufsummieren aller Boltzmannfaktoren (soviele wie es Instanzen gibt) zu Z; Z = Summe zur Normierung;
			    double Z = 0.0;

			    for (int z = 0; z < instance_index.size(); ++z)
			      {
				Z += BF;
			      }
			    
			    //Berechnen der Wahrscheinlichkeit
			    vector<double> p_BF;

			    for(int k=0; k < instance_index.size();++k)
			      {
				p_BF.push_back(k*BF/Z);
			      }

			   //Berechnung welche Instanz genommen wird

			    double sum1 = 0.0;
			    int instance = 0;
			  	
			  	for(int i = 0; i < p.size(); i++)
				{
				  sum1 += p.at(i);
				   if(sum1 > r1)
				    {
				      instance= i;
			              break;
				    }	
				}	

			    //die gewahelte Instance wird pickedReaction zugewiesen.
			    // die pickedReaction bekommt die Instanz zugewiesen, die in producedReactions vorhanden ist.
			    // --> die Reaction Instanz, die genommen wird!!!!
			    //
				MR_Reactions::Reaction_Container pickedReaction;
				pickedReaction = instance_index.at(instance);

				
				/*
				  ich weiss jetzt welche instance ich nehme und welche regel
				  
				*/

				
			//**********************************************************************************************************************************************************************
			//!!!!!begin aenderung just for me to check --> not important to the projekt
			//!!! ausgabe producedREactions
			std::cout <<" ===============Produced Reactions" << std::endl;
			for (MR_Reactions::Reaction_Container::const_iterator itt = producedReactions.begin(); itt!=producedReactions.end();itt++)
			{

			  std::cout << "Produced Reactions" << *itt << std::endl;
			  //	std::cout << "Produced Reactions metabolites" << itt->metabolites << endl;
			}
			//!!!!Ausgabe von producedSmiles/targetSmiles
			std::cout << "================ProducedSmiles" << std::endl;
			for (SMILES_container::const_iterator ittry =producedSmiles.begin(); ittry != producedSmiles.end(); ittry++)
			{
			  //std::cout << std::endl << std::endl;
			  std::cout << "Ausgabe ProducedSmiles" << ittry->first << " " << ittry -> second << std::endl;
			}					 
			std::cout << "================ TargetSmiles" << std::endl; 
			for (SMILES_container::const_iterator iteratortry =targetSmiles.begin(); iteratortry != targetSmiles.end(); iteratortry++)
			{
			  //std::cout << std::endl << std::endl;
			  std::cout << "Ausgabe TargetSmiles" << iteratortry->first << std::endl;
			}
			std::cout << "===========ruleId2Reactions" << std::endl; 
			    for (std::map<std::string,size_t>::const_iterator iterator1 = ruleId2reactions.begin(); iterator1 !=ruleId2reactions.end(); iterator1++)
			    {
				
			   std::cout << " Ausgabe ID: key: " << iterator1->first <<" value: "<< iterator1->second << std::endl;
			   }
		
				//!!!!! aenderung ende */

			    
			    //========Ende aenderung gillespie algorithm

			    
			// print stats on molecules and reactions
			(*out)	<<"\n# " <<(it) <<". iteration :"
					<<"\t molecules = " <<targetSmiles.size()
					<<"\t reactions = " <<producedReactions.size() <<"\n";

			if (infoMode == OUT_VERBOSE) {
				// print number of reactions per rule
				for( RuleHist::const_iterator r=ruleId2reactions.begin(); r != ruleId2reactions.end(); r++) {
					(*out) <<"\t\t" <<r->second<<"\t" <<r->first <<"\n";
				}
				(*out)	<<std::endl;
				// print all reactions
				for (MR_Reactions::Reaction_Container::const_iterator r=producedReactions.begin(); r != producedReactions.end(); r++) {
					(*out) <<"\t\t"<<*r<<"\n";
				}
				(*out) <<std::endl;
			}

			// stop processing if no reactions possible
			if (producedReactions.empty()) {
				break;
			}

			// pick a random rule
			size_t pickedRule = getRandomNumber( ruleId2reactions.size() );
			// get hist data of picked rule
			RuleHist::const_iterator pickedRuleHist = ruleId2reactions.begin();
			while( pickedRule > 0 ) {
				pickedRuleHist++;
				pickedRule--;
			}

			//!!!beginn aenderung
			std::cout << "=============pickedRule (int number)" << std::endl;
			std::cout << pickedRule << std::endl;


			//!!ende aenderung 
			// pick a reaction for this rule
			size_t pickedReactionNumber = getRandomNumber( pickedRuleHist->second );
			// find according reaction (TODO could be done via binary search since reactions are ordered by ruleID)
			MR_Reactions::Reaction_Container::const_iterator pickedReaction = producedReactions.end();
			for (MR_Reactions::Reaction_Container::const_iterator r=producedReactions.begin(); r != producedReactions.end(); r++) {
				// check if reaction for picked rule
				if (r->rule_id == pickedRuleHist->first) {
					// check if this is ithe reaction to report
					if (pickedReactionNumber == 0) {
						// store picked reaction
						pickedReaction = r;
						// stop search
						break;
					}
					pickedReactionNumber--;
				}
			}

				 //!!beginn aenderung!!!
				 std::cout << "========= pickedReactionNumber" << std::endl;
				 std::cout << pickedReactionNumber << std::endl;

				 std::cout <<"=======producedReaction nach der Wahl fuer die rule" << std::endl;

				 for (MR_Reactions::Reaction_Container::const_iterator iti = producedReactions.begin(); iti!=producedReactions.end();iti++)
				   {

				     std::cout << "Produced Reactions: " << *iti << std::endl; 
				   }
										    std::cout << "======Ende====Ausgabe=====" << std::endl;

				
				 //!!ende aenderung!!!

			// print information for picked reaction
			(*out)	<<"!!!!!  "<<(*pickedReaction) <<"\n"
					<<std::endl;

			
			std::vector<SMILES_container::value_type::second_type> used_educts;
			// remove educts from targetSmiles
			for (Reaction::Metabolite_Container::const_iterator m = pickedReaction->metabolites.begin(); m != pickedReaction->metabolites.end(); m++) {
				// get molecule
				SMILES_container::iterator metabolite = targetSmiles.find( *m );
				assert( metabolite != targetSmiles.end() ); // should never happen...
				// delete molecule graph

				
				delete( metabolite->second );
				// remove metabolite from targetSmiles
				targetSmiles.erase( metabolite );
			}

			// store products for next iteration in producedSmiles
			for (Reaction::Product_Container::const_iterator p = pickedReaction->products.begin(); p != pickedReaction->products.end(); p++) {
				if ( producedSmiles.find( *p ) == producedSmiles.end() ) {
					  // parse product SMILES to graph representation
					std::pair<ggl::chem::Molecule,int> result
						= ggl::chem::SMILESparser::parseSMILES( *p );
					  // check parsing result
					assert( result.second == -1 ); // parsing should always succeed!
					  // store molecule for next iteration
					producedSmiles[ *p ] = new ggl::chem::Molecule(result.first);
				}
			}

			// remove all reactions that share metabolites with picked reaction
			Reaction::Metabolite_Container consumedMetabolites(pickedReaction->metabolites);
			MR_Reactions::Reaction_Container::const_iterator r = producedReactions.begin();
			while (r != producedReactions.end()) {
				bool keepCurReaction = true;
				for (Reaction::Metabolite_Container::const_iterator m = consumedMetabolites.begin(); keepCurReaction && m != consumedMetabolites.end(); m++) {
					// check if no consumed metabolite is part of the metabolites of this reaction
					keepCurReaction = (r->metabolites.find(*m) == r->metabolites.end());
				}
				if ( ! keepCurReaction ) {
					// memorize element to be deleted
					MR_Reactions::Reaction_Container::const_iterator toRemove = r;
					// go to next reaction
					r++;
					// remove checked reaction (should keep r valid)
					producedReactions.erase( toRemove );
				} else {
					// go to next reaction
					r++;
				}
			}
			//!!!!Ausgabe ProducedSmiles!!!!!!
			std::cout << "================ProducedSmiles" << std::endl;
			for (SMILES_container::const_iterator itryi =producedSmiles.begin(); itryi != producedSmiles.end(); itryi++)
			{
			  //std::cout << std::endl << std::endl;
			  std::cout << "Ausgabe ProducedSmiles" << itryi->first << std::endl;
			}

		        //!!!!Ausgabe TargetSmiles
		        std::cout << "================ TargetSmiles" << std::endl; 
			for (SMILES_container::const_iterator iteratotry1 =targetSmiles.begin(); iteratotry1 != targetSmiles.end(); iteratotry1++)
			{
			  //std::cout << std::endl << std::endl;
			  std::cout << "Ausgabe TargetSmiles" << iteratotry1->first << std::endl;
			  }

		       //!!!Ausgabe deltaE-versuch
			    double sum_energy_producedSmiles = energie_calculation(producedSmiles);
			    std::cout <<"Summe producedSmiles" << sum_energy_producedSmiles << std::endl;
			    
			    double sum_energy_targetSmiles = energie_calculation(targetSmiles);
										 std::cout <<"Summe targtSmiles" << sum_energy_targetSmiles << std::endl;	       
			    double deltaE1 = sum_energy_producedSmiles - sum_energy_targetSmiles;
			    std::cout << "=======Ausgabe deltaE========"<<std::endl;
									  std::cout << deltaE1 << std::endl;
									  double kb = 8.61733*(pow(10,-5));
									  double T = 273.15;
									  double Bf = exp(-deltaE1/kb*T);
			    std::cout << "=======Ausgabe Boltzmannfaktor=======" << std::endl; 									  
			    std::cout << Bf << std::endl;


		 //erhoehung der Zeit durch tau
			    t += tau;
			    
		//erhoehung von n wenn noetig
		// n += 1;
			    
			    std::cout << "==============Beginn naechste Iteration===========" << std::endl;
										      
		} // end rule application iteration loop
		
		
		//////////////////////////////////////////////////////////////
		// print final set of molecules
		//////////////////////////////////////////////////////////////

		// print stats on molecules and reactions
		(*out)	<<"\n# final molecules :\n\n";
		for (SMILES_container::const_iterator m=targetSmiles.begin();m!=targetSmiles.end();m++) {
			(*out) <<m->first <<"\n";
		}
		for (SMILES_container::const_iterator m=producedSmiles.begin();m!=producedSmiles.end();m++) {
			(*out) <<m->first <<"\n";
		}


	
	} catch (std::exception& ex) {
		(*out) <<"\n\n ERROR : " <<ex.what() <<"\n"<<std::endl;
		exitValue = -1;
	}
	
	
	//////////////////////////////////////////////////////////////
	// garbage collection
	//////////////////////////////////////////////////////////////

	if (aromaticityPrediction != NULL)	delete aromaticityPrediction;
	if (rateCalc != NULL)			delete rateCalc;
	out = &std::cout;
	if (outFile != NULL)			delete outFile;
	for (RulePatternMap::iterator pat = rulePattern.begin();
			pat != rulePattern.end(); ++pat)
	{
		for (size_t p=0; p<pat->second.size(); ++p)
			delete pat->second[p];
	}
	for (SMILES_container::iterator it=targetSmiles.begin(); it!=targetSmiles.end() ; ++it ) {
		delete it->second;
	}
	return exitValue;
}

//////////////////////////////////////////////////////////////////////////

/*!
 * Initialises allowed parameters and their default values and prepares an
 * information string for the tool.
 */
void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText )  
{
	infoText = "\n"
		"Reads a list of molecules and chemical rules and does a chemical "
		"reaction simulation.\n"
		"\n"
		"Note, neither molecules nor rules are allowed to specify atom class label information!"
		"\n"
		"Rules have to be in GML format (use '-ruleExample' for an example)."
		"\n"
		"It is possible to specify the molecules in GML format as well."
		;
	
	allowedArgs.push_back(biu::COption(	
							"rules", true, biu::COption::STRING, 
							"Chemical rules input : 'STDIN' when to read from standard input, or a ':'-separated list of file names (use -ruleExample for a rule sample)"));
	allowedArgs.push_back(biu::COption(	
							"smiles", true, biu::COption::STRING, 
							"SMILES input : 'STDIN' when to read from standard input, or a ':'-separated list of file names."));
	allowedArgs.push_back(biu::COption(	
							"mols", true, biu::COption::STRING, 
							"Molecules input in GML format : 'STDIN' when to read from standard input, or a ':'-separated list of file names."));
	allowedArgs.push_back(biu::COption(	
							"groups", true, biu::COption::STRING,
							"Predefined molecule groups abbreviated in rules and molecules in GML format : 'STDIN' when to read from standard input, or a ':'-separated list of file names."));
	allowedArgs.push_back(biu::COption(
							"iter", true, biu::COption::INT, 
							"Number of rule application iterations",
							"1"));
	allowedArgs.push_back(biu::COption(
							"srand", true, biu::COption::INT,
							"Seed value for random number generator initialization. If -1, the current time is used.",
							"-1"));
	allowedArgs.push_back(biu::COption(	
							"out", true, biu::COption::STRING, 
							"Output file name or 'STDOUT' when to write to standard output",
							"STDOUT"));
	allowedArgs.push_back(biu::COption(
							"aromaticity", true, biu::COption::CHAR,
							"The aromaticity perception model to be used : (M)arvin general model, (O)penBabel model, or (N)o aromaticity perception.",
							"N"));
	allowedArgs.push_back(biu::COption(
							"noInputCorrection", true, biu::COption::BOOL,
							"Dont correct the input molecules (aromaticity perception, proton filling, ...)"));
	allowedArgs.push_back(biu::COption(
							"noInputCheck", true, biu::COption::BOOL,
							"Dont check the input molecules for consistency (atom/bond label, ...)"));
	allowedArgs.push_back(biu::COption(
							"noRuleCheck", true, biu::COption::BOOL, 
							"Dont check the rules for consistency"));
	allowedArgs.push_back(biu::COption(	
							"help", true, biu::COption::BOOL, 
							"Displays help on all parameters"));
	allowedArgs.push_back(biu::COption(	
							"ruleExample", true, biu::COption::BOOL, 
							"Displays an example for the chemical reaction GML encoding"));
	allowedArgs.push_back(biu::COption(	
							"graphExample", true, biu::COption::BOOL, 
							"Displays an example for the molecule graph GML encoding"));
	allowedArgs.push_back(biu::COption(	
							"v", true, biu::COption::BOOL, 
							"Verbose output"));
	allowedArgs.push_back(biu::COption(	
							"version", true, biu::COption::BOOL, 
							"Version information"));
}

//////////////////////////////////////////////////////////////////////////

void
printMoleculeGraph( std::ostream& out, ggl::chem::Molecule & m )
{
	sgm::Graph_boost<ggl::chem::Molecule> g(m);
	for (size_t i=0; i<g.getNodeNumber(); ++i ) {
		out <<std::setw(6) <<i <<" ("<<g.getNodeLabel(i) <<")  --> ";
		for (sgm::Graph_Interface::OutEdge_iterator e = g.getOutEdgesBegin(i),
				eEnd = g.getOutEdgesEnd(i); e != eEnd; ++e)
		{
			out <<" | " <<e->getToIndex() <<" (" <<e->getEdgeLabel() <<")";
		}
		out <<" |\n";
	}

}

//////////////////////////////////////////////////////////////////////////

void
printMoleculeGraphs( std::ostream& out, SMILES_container & smiles )
{
	for (	SMILES_container::const_iterator s = smiles.begin(); 
			s!= smiles.end(); ++s )
	{
		  // parse SMILES to graph
		std::pair<ggl::chem::Molecule,int> result 
			= ggl::chem::SMILESparser::parseSMILES(s->first);
		  // check parsing result
		if (result.second != -1) {
			std::ostringstream oss;
			oss	<<"printMoleculeGraphs : parsing error in SMILES string '"
				<<s->first
				<<"' at position " <<result.second;
			throw ArgException(oss.str());
		}
		printMoleculeGraph(out, result.first);
		out <<"\n";
	}

}


//////////////////////////////////////////////////////////////////////////

void
printMoleculeGraphsGML( std::ostream& out, SMILES_container & smiles, const bool withSpaces  )
{
	size_t i=0;
	for (	SMILES_container::const_iterator s = smiles.begin(); 
			s!= smiles.end(); ++s )
	{
		  // parse SMILES to graph
		std::pair<ggl::chem::Molecule,int> result 
			= ggl::chem::SMILESparser::parseSMILES(s->first);
		  // check parsing result
		if (result.second != -1) {
			std::ostringstream oss;
			oss	<<"printMoleculeGraphs : parsing error in SMILES string "
				<<s->first
				<<" at position " <<result.second;
			throw ArgException(oss.str());
		}
		out <<((i==0)?"":"\n")<<"# result graph " <<i <<"\n";
		ggl::Graph_GML_writer::write(out, result.first, withSpaces);
		i++;
	}
	out <<std::endl;

}


//////////////////////////////////////////////////////////////////////////

void
giveRuleGMLExample()
{
	std::cout <<"\n\
====== RULE TO ENCODE =======================\n\
\n\
  Diels-Alder reaction:\n\
\n\
  0(C)   5(C) = 4(C)               0(C) - 5(C) - 4(C) \n\
                                                      \n\
   ||                     ==>       |             |   \n\
                                                      \n\
  1(C) - 2(C) = 3(C)               1(C) = 2(C) - 3(C) \n\
\n\
======= RULE IN GML =========================\n\
\n\
rule [\n\
 ruleID \"Diels-Alder reaction\"\n\
 context [\n\
   node [ id 0 label \"C\" ]\n\
   node [ id 1 label \"C\" ]\n\
   node [ id 2 label \"C\" ]\n\
   node [ id 3 label \"C\" ]\n\
   node [ id 4 label \"C\" ]\n\
   node [ id 5 label \"C\" ]\n\
 ]\n\
 left [\n\
   edge [ source 0 target 1 label \"=\" ]\n\
   edge [ source 1 target 2 label \"-\" ]\n\
   edge [ source 2 target 3 label \"=\" ]\n\
   edge [ source 4 target 5 label \"=\" ]\n\
   constrainNoEdge [ source 0 target 4 ]\n\
   constrainNoEdge [ source 3 target 5 ]\n\
 ]\n\
 right [\n\
   edge [ source 0 target 1 label \"-\" ]\n\
   edge [ source 0 target 5 label \"-\" ]\n\
   edge [ source 1 target 2 label \"=\" ]\n\
   edge [ source 2 target 3 label \"-\" ]\n\
   edge [ source 3 target 4 label \"-\" ]\n\
   edge [ source 4 target 5 label \"-\" ]\n\
 ]\n\
]\n\
\n\
=============================================\n"
                <<std::endl;
}


//////////////////////////////////////////////////////////////////////////

void
giveGraphGMLExample()
{
	std::cout <<"\n\
====== MOLECULE TO ENCODE =======================\n\
\n\
  SMILES :  'C=CC(C)=C' \n\
\n\
======= MOLECULE IN GML =========================\n\
\n\
# molecule C=CC(C)=C \n\
graph [\n\
  node [ id 0 label \"C\" ]\n\
  node [ id 1 label \"C\" ]\n\
  node [ id 2 label \"C\" ]\n\
  node [ id 3 label \"C\" ]\n\
  node [ id 4 label \"C\" ]\n\
  edge [ source 1 target 0 label \"=\" ]\n\
  edge [ source 2 target 1 label \"-\" ]\n\
  edge [ source 3 target 2 label \"-\" ]\n\
  edge [ source 4 target 2 label \"=\" ]\n\
]\n\
\n\
==============================================\n"
                <<std::endl;
}

//////////////////////////////////////////////////////////////////////////

bool
isWithoutAtomClass( sgm::Graph_Interface & graph )
{
	for (size_t i=0; i<graph.getNodeNumber(); i++) {
		// check if node label contains class separator (not allowed)
		if ( graph.getNodeLabel(i).find(":") != std::string::npos ) {
			// not allowed!
			return false;
		}
	}
	// no invalid atom label found
	return true;
}

//////////////////////////////////////////////////////////////////////////

bool
allWithAtomClass( sgm::Graph_Interface & graph )
{
	for (size_t i=0; i<graph.getNodeNumber(); i++) {
		// check if node label contains class separator (required)
		if ( graph.getNodeLabel(i).find(":") == std::string::npos
				// ensure there is something behind the separator
			|| graph.getNodeLabel(i).find(":")+1 >= graph.getNodeLabel(i).size() )
		{
			// not allowed!
			return false;
		}
	}
	// no invalid atom label found
	return true;
}

//////////////////////////////////////////////////////////////////////////

size_t
extractAtomClass( ggl::chem::Molecule_Graph & mol, std::set<size_t> usedAtomClasses )
{
	for (size_t i=0; i<mol.getNodeNumber(); i++) {
		// ensure the class label was not known
		const size_t classLabel = ggl::chem::MoleculeUtil::getClass( mol.getNodeLabel(i) );
		if ( classLabel > 0 && !usedAtomClasses.insert( classLabel ).second ) {
			// return redundant class label
			return classLabel;
		}
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////
