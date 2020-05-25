

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

//RandomNumberGenerator
double random_number01()
{
	return (double(rand())/(RAND_MAX));
}


//create generic smiles
//gives back the generic smile for calculating the number of moleculs 

std::string extract_generic_smile(string molecule){
  bool in_brackets = false;
  bool after_colon = false;
  string generic_smile;

  for (int i = 0; i < molecule.size(); ++i){
      char c = molecule[i];

      if(c =='['){
	  in_brackets = true;
	  continue;
	}
      if(c ==']'){
	  in_brackets = false;
	  after_colon = false;
	  continue;
	}
      if(in_brackets && c == ':'){
	  after_colon = true;
	  continue;
	}
      if (after_colon){
	  continue;
	}
      /*if(c =='#'){
	  in_brackets = false;
	  after_colon = false;
	  continue;
	  }*/
      generic_smile.push_back(c);
    }
 return generic_smile;
}

//calculate the energy difference of products and metabolites including boltzmann: exp(-deltaE/kBT) 
//stand 03.05. --> only calculate the difference because boltzmann does not work correct at the moment!
double energy_per_instance(double sum_products, double sum_metabolites) {
 //boltzmannconstant
 double kB;
 //absolut temperature
 double T;
 //boltzmannfactor
 double kBT;
//calculting boltzmann 
double boltzmann=0.0; 
//storing the energy_difference
 double dE = 0.0;

 dE = sum_products - sum_metabolites; 
 //for later, when boltzmann is working
 //boltzmann = exp((-dE)/kBT
 //
return dE; 
}

//version 03.05.2020 Beginn
//calculate the energy of metabolites and products and return a multimap 
multimap<double,vector<string> > rule_dE(vector<string> metabolites, vector<string> products){
 multimap<double,vector<string > > rule_;
 double sum_metabolites = 0.0;
 double sum_products = 0.0;
 double dE = 0.0;

 //calculate the energy of the metabolites
 for(int i= 0; i<metabolites.size();++i){
	 //tno
	 cout <<"--------temporaere ausgabe ther metabolites ------------" << metabolites.at(i) << endl; 
  sum_metabolites += energie_calculation(metabolites.at(i));
 }
 //calculate the energy of the products
 for(int j= 0; j<products.size();++j){
  sum_products += energie_calculation(products.at(j));
 }
 //calculate the energy difference and do it statistically correct 
 dE = energy_per_instance(sum_products,sum_metabolites);

 //put all together in one multimap
 rule_.insert(make_pair(dE,metabolites));

 return rule_;
}

//add all energy values per instance for one rule and norm it 

map<string,double> ruleID_calc(multimap<string,multimap<double,vector<string> > > ruleID){
 map<string,double> ruleID_total;
 //add all values of the multimap ruleID
 for(multimap<string,multimap<double,vector<string> > >::iterator iter=ruleID.begin();iter!=ruleID.end();++iter){
  if(ruleID_total.find(iter->first)==ruleID_total.end()){
   for(multimap<double,vector<string> >::iterator ii=(*iter).second.begin();ii!=(*iter).second.end();++ii){
    ruleID_total[iter->first] = (*ii).first;
   }
  }
  else{
   for(multimap<double,vector<string> >::iterator ii=(*iter).second.begin();ii!=(*iter).second.end();++ii){
    ruleID_total[iter->first]+=(*ii).first;
   }
  }
 }
cout << "in unterfunktion---------------------rule_ID_total - add all values of the multimap ruleID" << endl;
for(map<string,double>::iterator i = ruleID_total.begin();i!=ruleID_total.end();++i)
{
  cout << i->first <<":" << i->second << endl; 
}
 //add all values together to Z_total
 double Z_total = 0.0;
 for(map<string,double>::iterator it=ruleID_total.begin();it!=ruleID_total.end();++it){
  Z_total += it->second;
 }
//temporary
cout <<"Z_total" << Z_total<< endl;

 //divide all entries of ruleID_total by Z_total to get percentage per rule between 0 and 1 (all together it has to be 1)
 //store percentage in a map<string,double>
 map<string,double> ruleID_norm;

 for(map<string,double>::iterator i = ruleID_total.begin();i!=ruleID_total.end();++i){
  if(ruleID_norm.find(i->first)==ruleID_norm.end()){
   ruleID_norm[i->first] = (i->second)/Z_total;
  }
  /*else{
   ruleID_norm[i->first] +=(i->second)/Z_total;
  }*/
 }
 return ruleID_norm;
}

//calculate the reactionrate per ruleID
double reaction_rates_calc(vector<string>targetSmile,map<string,int> mpl_ts){
 //for store the generic_smiles
  map<string,double> _generic_smile;
  //calculate how often a metabolite exist per reactioninstance per rule
  for(int i= 0; i < targetSmile.size();++i){
   _generic_smile[extract_generic_smile(targetSmile.at(i))]+=1;
  }

  //calculate the sole reactionrates per rule 
  double v1 = 1.0;
  double v2 = 1.0;
  double v3 = 1.0;
  double value = 0.0; 

 for(map<string,double>::iterator i = _generic_smile.begin();i!=_generic_smile.end();++i){
  if(i->second == 1){
   v1= v1*mpl_ts[i->first];
  }
  if(i->second ==2){
   v2=(mpl_ts[i->first]*(mpl_ts[i->first]-1))/2;
  }
  if(i->second==3){
   v3=(mpl_ts[i->first]*(mpl_ts[i->first]-1)*(mpl_ts[i->first]-2))/6; 
  } 
 }
 value = v1*v2*v3;
 return value;
}

//calculate the rule rates per rule
map<string,double> rule_rates_calc(map<string,double> _ruleID_percentage,map<string,double> _ruleID_reactionRates){
  //store the rule rates with rule_id
  map<string,double> a_my;

  for(map<string,double>::iterator i = _ruleID_percentage.begin();i!=_ruleID_percentage.end();++i){
    map<string,double>::iterator j = _ruleID_reactionRates.find(i->first);
    if(j!=_ruleID_reactionRates.end()){
      a_my[i->first] = i->second * j->second; 
    }
  }
return a_my;
}

//calculate the total rule rate
double rule_rate_total_calc(map<string,double> _rule_rates){
  double a0 = 0.0;

  for(map<string,double>::iterator iter = _rule_rates.begin();iter!=_rule_rates.end();++iter){
   a0 += iter->second;
  }
  return a0;
}

//calculate the ruleID according to reactionrate and energy -- gillespie part
string calculate_ruleID(map<string,double> _rule_rates, double _rule_rate_total, double _r2){
	double sum_a0 = 0.0;
	string taken_ruleID;
  
 for(map<string,double>::iterator iter= _rule_rates.begin();iter!=_rule_rates.end();++iter){
	 sum_a0+=iter->second;
	 if(sum_a0 > _rule_rate_total*_r2){
	  taken_ruleID = iter->first;
	  cout << "--------------ruleid-taken-----------"<< iter->first << endl;
	 break; 
	 }
 }
 return taken_ruleID;
}

//calculate the reaction which will be taken for the rule_id
int reaction_taken(vector<double> _instance_energy)
{
 //sum all values up
 srand(time(NULL));
 double r2 = random_number01();
 double r3 = random_number01();
  
 double sum_instance = 0.0;

 for(int i = 0; i<_instance_energy.size();++i){
	 sum_instance += _instance_energy.at(i);
 }
cout <<"summe instance : " << sum_instance << endl; 
 //norm the entries
 vector<double> norm_instance_energy;
 for(int j= 0; j<_instance_energy.size();++j){
	 norm_instance_energy.push_back((_instance_energy.at(j))/sum_instance);
 }
 //get to know the reaction instance
 double sum_ = 0.0;
 int r = 0;

 for(int i = 0; i <norm_instance_energy.size();++i){
	 sum_ += norm_instance_energy.at(i);
	 if(sum_ > r3){
	  r = i;
	  break;
	 }
 }
 cout << "Diese instanz wurde genommen" << r  << endl; 
 return r; 
}

//calculate the ruleID and the metabolites according to the ruleID
map<string,vector<string> > calc_ruleID_met(multimap<string,multimap<double,vector<string> > > _ruleID, string _ruleID_taken){
	vector<double> instance_energy;
	vector<string> vec_tmp;
	int rp = 0; 

 //iterators for equal range
 typedef multimap<string,multimap<double,vector<string> > >::iterator ip;
 pair<ip,ip> result = _ruleID.equal_range(_ruleID_taken);

 //iterator outermap
 multimap<string,multimap<double,vector<string> > >::iterator om;
 //iterator innermap
 multimap<double,vector<string> >::iterator im;
//iteraotr vector
vector<string>::iterator vecin;

//store energyvaluesi
for(om = result.first;om!=result.second;++om){
 for(im = (*om).second.begin();im!=(*om).second.end();++im){
  instance_energy.push_back((*im).first);
 }
}

 //give the vector of energies to anoter function for calculating the reaction instance
 rp = reaction_taken(instance_energy);

 //pick the metabolites to the ruleID_taken and store all in a map
 
 //laufvariable
 int l = 0; 

 for(om = result.first;om!=result.second;++om){
	 cout << om->first << endl;
	 if(rp == l ){
	   for(im = (*om).second.begin();im!=(*om).second.end();++im){
	     for(vecin = (*im).second.begin();vecin!=(*im).second.end();++vecin){
		      vec_tmp.push_back(*vecin);
		      cout << *vecin << endl;  
	      }
	     break;
	    }
	 }
	   ++l;
 }
map<string,vector<string> > _ruleID_i;
_ruleID_i[_ruleID_taken]= vec_tmp;

return _ruleID_i; 
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

		//Changes for Gillespie 
		//Calculate the molecular_population_level in targetSmile
		
		std::map<std::string,int> mpl_ts;
		std::map<std::string,int> mpl_ps;
		
               	 //initialize the time t
                   double t = 0.0;
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
					
				// merge targetSmiles and producedSmiles into targetSmiles
				targetSmiles.insert(producedSmiles.begin(), producedSmiles.end());
				producedSmiles.clear();
				    
				// clear reaction product container, since most will not be needed
				for (SMILES_container::iterator it=toFill.begin(); it!=toFill.end() ; ++it ) {
					delete it->second;
				}
				toFill.clear();
			}

			// pick a reaction TODO change accordingly!!
			//Step 1

		cout << "-------------Ausgabe TargetSmiles--------------"<<endl ; 

			for(SMILES_container::const_iterator tts = targetSmiles.begin();tts!=targetSmiles.end();tts++){
				cout << tts->first <<   endl; 
			}
			//calculate the molcule population via TargetSmiles : how many molecules are present at them moment 
			for (SMILES_container::const_iterator its = targetSmiles.begin();its != targetSmiles.end();its++){
				string generic_smile = extract_generic_smile(its->first);

				if (mpl_ts.find(generic_smile) == mpl_ts.end()){
					mpl_ts[generic_smile] = 1;
				}
				else{
					mpl_ts[generic_smile]++;
				}
			}

			cout <<"----------mpl_ts------------" << endl;

			for (std::map<std::string,int>::const_iterator itera = mpl_ts.begin();itera != mpl_ts.end();itera++){
				cout << itera -> first << " : " << itera -> second << endl;	
			}

		       //generate random numbers
                           double r0 = random_number01();
			   cout << " r0 : " << r0 << endl; 
                           double r1 = random_number01();
			   cout << " r1 : " << r1 << endl;  
			
			//NEUE VERSION VOM 03.05.2020 
			vector<string> metabolites_all;	
			vector<string> products_all;
			//multimap for store the deltaE and metabolites(-->metabolites are unique)
			typedef multimap<double,vector<string> > dE_met;
			dE_met Emet;
		        multimap<string,dE_met> ruleID;	//store all ruleID-energyvalue-metabolites
			map<string,double> ruleID_percentage;//stores the percentage per ruleID (c_my)
			map<string,double> ruleID_reactionRates;//molecule combination (h_my)
			map<string,double> rule_rates; //a_my
			double rule_rate_total=0.0; //a0;
			string rule_taken;
			double dt = 0.0; //random waiting time
			//NEUE VERSION VOM 03.05.2020 ENDE
			// get number of reaction instances for each rule
			typedef std::map< std::string, size_t > RuleHist;
			RuleHist ruleId2reactions;
			for (MR_Reactions::Reaction_Container::const_iterator r=producedReactions.begin(); r != producedReactions.end(); r++) {
				//ausgabe von producedreaction -- selbst geschrieben
				cout<< "---------------ausgabe r------------"<< endl; 
                                 cout << (*r) << endl;

				 cout <<"-------------ausgabe r ende ------------" << endl;  
				// check if not seen so far  
				  if (ruleId2reactions.find(r->rule_id) == ruleId2reactions.end()) {
					ruleId2reactions[r->rule_id] = 1;
					//nummerierung der Regeln!!!!!!!!!!!!!!!!!!
				        cout <<"-----------ruleId2reactions-----------" << endl;
				       for (map<string,size_t>::iterator rulerator= ruleId2reactions.begin();rulerator != ruleId2reactions.end();++rulerator)
				            {
						    cout << rulerator -> first << " : " << rulerator -> second << endl; 
				            }		    
				       cout << "----------ruleId2reactions-Ende----------" << endl;
					//metabolites
					for(Reaction::Metabolite_Container::const_iterator me1 = r-> metabolites.begin(); me1 != r->metabolites.end();++me1){
						//all metabolites which occur in one reaction instance to store for adding later on
						metabolites_all.push_back(*me1);
					}		
                                        //products
					for(Reaction::Metabolite_Container::const_iterator pr1 = r->products.begin(); pr1 != r->products.end(); ++pr1){
						//all products which occur in one reaction instance to store for adding later on 
						products_all.push_back(*pr1); 
					} 
				
				       //fill the innermap of ruleID with energydifference and metabolites(metabolites because of uniqueness)
                                         Emet = rule_dE(metabolites_all,products_all);
				       
				      //fill the outermap ruleID
				         ruleID.insert(make_pair(r->rule_id,Emet));

				      //calculate the reactionrates according to the molecules occur targetsmiles
				        ruleID_reactionRates[r->rule_id]=reaction_rates_calc(metabolites_all,mpl_ts);
					//just for representation --> temporary
					cout <<"-----------------------reation_rate_per rule ---------------------" << endl; 
					for(map<string,double>::iterator i = ruleID_reactionRates.begin();i!=ruleID_reactionRates.end();++i){
						cout << i->first <<": " << i->second << endl; 
				        }		
			                metabolites_all.clear();	
				        products_all.clear(); 	
					//Ende
				  }
				  
				  else {
					// increase counter
					ruleId2reactions[r->rule_id]++;

					//Version4.5.2020
					//metabolites 
					for(Reaction::Metabolite_Container::const_iterator me2 = r-> metabolites.begin();me2!=r->metabolites.end();++me2){
						//all metabolites which occur in one reaction instance to store for adding later on
                                                metabolites_all.push_back(*me2); 
					}
                                       
                                  	//products
				        for(Reaction::Metabolite_Container::const_iterator pr2 = r->products.begin();pr2 != r->products.end();++pr2){
						 //all products which occur in one reaction instance to store for adding later on
                                                products_all.push_back(*pr2); 
					}

					//fill the innermap of ruleID with energydifference and metabolites(metabolites because of uniqueness)
                                         Emet = rule_dE(metabolites_all,products_all);

                                      //fill the outermap ruleID
                                         ruleID.insert(make_pair(r->rule_id,Emet));

					metabolites_all.clear();
					products_all.clear();
					}				 
			}//for-loop end
			//latest version
			//
			cout <<"----------------------------ruleID-ENERGYVALUE-METABOLITES-------------"<<endl; 
		        for(multimap<string,dE_met>::iterator irule=ruleID.begin();irule!=ruleID.end();++irule){
                                          cout << irule->first << ":";
                                           for(multimap<double,vector<string> >::iterator idE=(*irule).second.begin();idE!=(*irule).second.end();++idE){
                                            cout <<(*idE).first << ":";
                                            for(vector<string>::iterator iv=(*idE).second.begin();iv!=(*idE).second.end();++iv){
                                             cout << (*iv) << endl;
                                            }
                                           }
                                        }

	               //contains the ruleID and the percentage
		        ruleID_percentage = ruleID_calc(ruleID);

			//just for output - temporary
			cout << "-------------ruleID_percentage--------------";
			for (map<string,double>::iterator i= ruleID_percentage.begin();i!=ruleID_percentage.end();++i){
			 cout << i->first <<":" << i->second<< endl; 
			}
//calculate the rule_rates for each rule
rule_rates = rule_rates_calc(ruleID_percentage,ruleID_reactionRates);
//just for output - temporary
cout << "----------------rule_rate----------------"<< endl; 
for(map<string,double>::iterator i=rule_rates.begin();i!=rule_rates.end();++i){
	cout <<i->first << ": " << i->second << endl;
}

//calculate the total rule_rate
rule_rate_total = rule_rate_total_calc(rule_rates);
cout<<"----------rule_rate_total--------------" << rule_rate_total << endl;

//calculate dt (=waiting time=a0) 
dt = (-log(r0))/rule_rate_total;
//calculate the ruleID
rule_taken = calculate_ruleID(rule_rates,rule_rate_total,r1); 
//ausgabe der genommenen ruleID - just temporary
cout <<"-------------------ruleID which is taken---------------" << rule_taken << endl;

//calculate the reaction_instance of the rule ID in an extern function- store ruleID and belonging metabolites!
map<string,vector<string> > ruleID_metabolites;
ruleID_metabolites = calc_ruleID_met(ruleID,rule_taken);

//temporary output of map
cout << "ruleID and belonging metabolites of the instance---------------------------" << endl; 
for(map<string,vector<string> >::iterator i = ruleID_metabolites.begin();i!=ruleID_metabolites.end();++i)
{
	cout << i->first <<": " <<  endl;
  for(vector<string>::iterator j = (*i).second.begin();j!=(*i).second.end();++j)
  {
    cout << (*j) << endl; 
  }
}
cout <<"RUleID and Metabolites END ---------------------" << endl; 

//latest version end    
			    //========Ende aenderung gillespie algorithm======================

			    
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
										    
			std::cout << "======Ende====Ausgabe=====" << std::endl;

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

		     
		 //update of the time
			    t = t+dt;
	    
	        std::cout << "==============Beginn naechste Iteration===========" << std::endl; 						
		mpl_ts.clear();		
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
