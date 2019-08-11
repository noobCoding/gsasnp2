#include <iostream>
#include "boost/program_options.hpp" 
//#include "boost/filesystem.hpp" 
#include "gsa.h"

using namespace std;
using boost::math::normal; // typedef provides default type is double.
using std::left; using std::showpoint; using std::noshowpoint;
using std::setw; using std::setprecision;
using std::numeric_limits;

typedef struct SetInfo {
	string setname;
	int genecount;
	int setsize;
	double zscore;
	double pvalue;
	double qvalue;
	double adjZscore;
	double adjPval;
	double adjQval;
} setinfo;

namespace
{
	const size_t ERROR_IN_COMMAND_LINE = 1;
	const size_t SUCCESS = 0;
	const size_t ERROR_UNHANDLED_EXCEPTION = 2;

} // namespace 

//#define handle_error(msg) \  do { perror(msg); exit(EXIT_FAILURE); } while (0)

#define ADJACENT_CORRELATION	"data/EUR_Adjacent_correlation"
#define SNP_LOC_MAP					"data/rsloc_hg19"
#define GENE_LOC_FILE				"data/hg19GeneList"
#define DEFAULT_SNP_FILE			"data/DIAGRAM"
#define DEFAULT_GENE_FILE			"data/DIAGRAMgene"
#define DEFAULT_MAP_FILE			"data/db19_20k"
#define DEFAULT_SET_FILE			"data/c2.cp.v5.2.symbols.gmt"
#define DEFAULT_OUTPUT				"DIAGRAM_zscore_result.txt"
#define MINSETSTR					"10"
#define MAXSETSTR					"200"
#define MINSET					10
#define MAXSET					200
#define TRUNCATION				0.25
#define CORRELATION				0.5



#define OUTPUT_HEADER		"Set\tSize\tCount\tz-score\tAdj. z-score\tAdj. p-value\tAdj. q-value\tList of genes\n"
#define USER_USAGE			"Basic usage:\n\n\tgsasnp2 [--input/-i] <snp/gene_list_file> [--pathway/-p] <set/pathway_list_file> [--snpgene/-s] <0/1> <additional parameters.....>\n"

//		MAIN	///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
	// definition
	string	inputFile = DEFAULT_SNP_FILE;
	string	setFile = DEFAULT_SET_FILE;
	string	geneMapFile = DEFAULT_MAP_FILE;
	string	adjGeneFile = ADJACENT_CORRELATION;
	string	outputFile = DEFAULT_OUTPUT;
	//double	truncation_level = TRUNCATION;
	double	correlation = CORRELATION;
	int minSetSize = MINSET, maxSetSize = MAXSET;
	bool	convert2symbol = false;	
	bool	inputType = SNP_INPUT;

	// initialization input

	try
	 {
		 std::string appName = argv[0]; //boost::filesystem::basename(argv[0]);		
		 std::vector<std::string> sentence;

		 /// Define and parse the program options
		 ///
		 namespace po = boost::program_options;
		 po::options_description desc("Valid options");
		 desc.add_options()
			 ("help,?", "Print a brief user usage")			
			("pathway,p", po::value<std::string>()->required(), "Path to input set/pathway list file. Necessary!")
			("input,i", po::value<std::string>()->required(), "Path to input SNP/GENE file. Necessary!")
			("snpgene,s", po::value<bool>()->required(), "This parameter indicates that input file is  SNP (0) or GENE (1) p-value.")
			("output,o", po::value<std::string>()->required()->default_value(DEFAULT_OUTPUT), "Path to output result. Default is 'adjusted_zscore_result.txt'.")
			("genemap,g", po::value<std::string>()->required()->default_value(DEFAULT_MAP_FILE), "Gene map file.")
			("adj,a", po::value<std::string>()->required()->default_value(ADJACENT_CORRELATION), "Inter-gene genotype correlation by race file.")			
			("minset", po::value<int>()->required()->default_value(MINSET), "Minimum set size (number of genes). Default minset = 10")
			("maxset", po::value<int>()->required()->default_value(MAXSET), "Minimum set size (number of genes). Default maxset = 200")					
			("symbol,b", po::value<bool>()->required()->default_value(false), "Convert gene ID in gene set/pathway map to symbols.");
		
		
		 po::variables_map vm;
		 try
		 {
			 po::store(po::parse_command_line(argc, argv, desc), vm);
		
			 /// --help option ///
			 if (vm.count("help"))
			 {
				 std::cout << USER_USAGE << std::endl << std::endl;
				 cout << desc << endl;
				 return SUCCESS;
			 }
			 po::notify(vm); // throws on error, so do after help in case 
							  //there are any problems 
		 }
		 catch (boost::program_options::required_option& e)
		 {
			 std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
			 std::cout << USER_USAGE << std::endl << std::endl;
			 cout << desc << endl;
			 return ERROR_IN_COMMAND_LINE;
		 }
		 catch (boost::program_options::error& e)
		 {
			 std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
			 std::cout << USER_USAGE << std::endl << std::endl;
			 cout << desc << endl;
			 return ERROR_IN_COMMAND_LINE;
		 }
		
		if (vm.count("input"))		inputFile = vm["input"].as<string>();
		if (vm.count("snpgene"))	inputType = vm["snpgene"].as<bool>();
		if (vm.count("pathway"))	setFile = vm["pathway"].as<string>();
		if (vm.count("output"))		outputFile = vm["output"].as<string>();// , cout << "Update output file: " << outputFile << endl;
		if (vm.count("map"))		geneMapFile = vm["map"].as<string>();// , cout << "Update gene map file: " << geneMapFile << endl;
		if (vm.count("adj"))		adjGeneFile = vm["adj"].as<string>();// , cout << "Update adjacent gene file: " << adjGeneFile << endl;
		if (vm.count("minset"))		minSetSize = vm["minset"].as<int>();// , cout << "Update min size of gene set/pathway: " << minSetSize << endl;
		if (vm.count("maxset"))		maxSetSize = vm["maxset"].as<int>();// , cout << "Update max size of gene set/pathway: " << maxSetSize << endl;		
		if (vm.count("symbol"))		convert2symbol = vm["symbol"].as<bool>();// , cout << "Update max size of gene set/pathway: " << maxSetSize << endl;
	 }
	 catch (std::exception& e)
	 {
		 std::cerr << "Unhandled Exception reached the top of main: "
			 << e.what() << ", application will now exit. Please note: " << std::endl;		
		 return ERROR_UNHANDLED_EXCEPTION;
	 } 

	// analyzing	
	gsa* test = new gsa();
	test->minSetSize = minSetSize;
	test->maxSetSize = maxSetSize;

	string	hgFile = "";
	string	geneListFile = "";
	int	hgVersion = 1;		// 0 : hg18, 1 : hg19, 2 : hg38

	if (geneMapFile.find("hg18") != string::npos) hgVersion = 0;
	if (geneMapFile.find("hg19") != string::npos) hgVersion = 1;
	if (geneMapFile.find("hg38") != string::npos) hgVersion = 2;
	
	switch (hgVersion)
	{
	case 0:
		hgFile = "data/rsloc_hg18";
		geneListFile = "data/hg18GeneList";
		break;
	case 1:
		hgFile = "data/rsloc_hg19";
		geneListFile = "data/hg19GeneList";
		break;
	case 2:
		hgFile = "data/rsloc_hg38";
		geneListFile = "data/hg38GeneList";
		break;
	default:
		hgFile = "data/rsloc_hg19";
		geneListFile = "data/hg19GeneList";
		break;
	}
	
	cout << "*** Welcome to GSA-SNP2 ***\n\n";
	time_t start = time(0);
	cout << "* Loading input data from file " << inputFile << "...";
	if (inputType == SNP_INPUT)		
		test->snpPvalueLoad(inputFile); // change file path	
	else // GENE_INPUT
		test->genePvalueLoad(inputFile);	

	cout << ".............. finished.\n";

	cout << "* Loading gene map from " << geneMapFile << "...";
	test->snpGeneMapGenerate(geneMapFile, inputType);
	cout << ".............. finished.\n";

	cout << "* Loading gene sets/pathways " << setFile << "...";
	test->setGeneLoad(setFile, convert2symbol);
	cout << ".............. finished.\n";

	if (inputType == SNP_INPUT) {
		cout << "* Loading gene information from " << adjGeneFile << "...";
		test->geneLocMapGenerate(geneListFile);
		test->adjacentGeneMapLoad(adjGeneFile);
		cout << ".............. finished.\n";
		cout << "* Loading SNP information...";
		test->snpLocMapGenerate(hgFile);
		cout << ".............. finished.\n";
	}
	
	time_t stop = time(0);
	cout << " - Loading time: " << difftime(stop, start) << " seconds\n\n";

	//////////////////////////////////////////////////////////////////////////
	cout << "* Processing data...";
	start = time(0);
	test->geneScoring(inputType); 
	if (inputType == SNP_INPUT) {
		test->geneFilter(correlation); 
	}

	test->adjustedGeneScoring(inputType);
	cout << ".............. finished.\n";
	cout << "* Adjusting scores...";
	test->adjustedPvalueBH();
	cout << ".............. finished.\n";

	stop = time(0);
	cout << " - Processing time: " << difftime(stop, start) << " seconds\n\n";
	//////////////////////////////////////////////////////////////////////////	
	//const int szSetInput = test->setGlobal.size();
	map<double, setinfo> mySet;
	double myepsilon = 0.000000001;
	int k = 0; // k * epsilon

	cout << "* Building output structures...";
	for (unsigned int i = 0; i < test->setGlobal.size(); i++)
	{
		setinfo tmp;
		tmp.setname = test->setGlobal[i];
		tmp.genecount = test->setGeneGlobal[tmp.setname].size();
		tmp.setsize = test->setGlobalActualSize[tmp.setname];
		tmp.zscore = test->zsOrigin[tmp.setname];	
		tmp.adjZscore = test->zsAdjusted[tmp.setname];
		tmp.adjPval = test->pvAdjusted[tmp.setname];
		tmp.adjQval = test->fdrAdjusted[tmp.setname];

		if (mySet.find(tmp.adjQval) != mySet.end())
			mySet[tmp.adjQval] = tmp;  // sorting based on  adjusted zscore		
		else
			mySet[tmp.adjQval + ++k * myepsilon] = tmp;  // sorting based on  original zscore		
	}
	cout << ".............. finished.\n";

	// printing output
	ofstream fou(outputFile);
	fou << OUTPUT_HEADER;
	map<double, setinfo>::iterator rit;
	cout << "* Writing output to file...";
	for (rit = mySet.begin(); rit != mySet.end(); ++rit)
	{
		// Insert the first item
		fou << rit->second.setname << "\t";
		fou << rit->second.setsize << "\t";
		fou << rit->second.genecount << "\t";
		fou << rit->second.zscore << "\t";		
		fou << rit->second.adjZscore << "\t";
		fou << rit->second.adjPval << "\t";
		fou << rit->second.adjQval << "\t";

		k = 1; // k * epsilon
		map<double, string> geneScore; // sort gene score & present them
		for (auto& vi : test->setGeneGlobal[rit->second.setname]) {
			double tmpscore = test->adjGeneScore[test->geneRefinedIndex[vi]];
			if (geneScore.find(tmpscore) != geneScore.end()) { // exist score
				geneScore[tmpscore + ++k * myepsilon] = vi;
			}
			else {
				geneScore[tmpscore] = vi;
			}
		}

		map<double, string>::reverse_iterator ri;
		for (ri = geneScore.rbegin(); ri != geneScore.rend(); ri++)
		{
			fou << ri->second.c_str() << "(" << ri->first << "); ";
		}	
		fou <<endl;	
	}
	fou.close();
	cout << ".............. finished.\n";

	////////////////////////////////////////////////////////////////////////// Finishing	
	delete test; test = 0;
	cout << "\nResult is saved at: " << outputFile << "\n";
	//cout << "\nAll done! Press any key to exit." << endl;
	//cin.get();
	return EXIT_SUCCESS;
}
