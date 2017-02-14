#pragma once
#include <cmath>
#include <vector>
#include <thread>
#include <chrono>
#include <random>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <memory>
#include <map>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <windows.h>
#include <process.h>
#include <boost/math/distributions/normal.hpp> // for normal_distribution
#include <boost/math/distributions.hpp>
#include <afxwin.h>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/policies/policy.hpp>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <functional>
#include "regression.h"
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/layered/SugiyamaLayout.h>
#include <ogdf/energybased/SpringEmbedderFR.h>


using namespace std;
using boost::math::normal; // typedef provides default type is double.
using std::setw;
using std::setprecision;
using std::numeric_limits;

//typedef vector<uint32_t> data_t;

#define kGene	2
#define SIGNIFICANT_THRESHOLD_OF_PVALUE 0.0001
#define M_PI           3.14159265358979323846  /* pi */
#define SNP_INPUT		0
#define GENE_INPUT		1

typedef struct GeneLocation {
	string geneName;
	int	chromosome;
	int start;
	int end;
	virtual ~GeneLocation() {};
}	genelocation;

typedef	struct ItsAGene {
	string	geneName;		// gene name	~ gene1Name		
	string	snpL, snpM, snpR;
	virtual ~ItsAGene() {};
}	itsagene;

typedef struct AdjacentGene {
	itsagene	geneA;
	itsagene	geneB;
	double		corr[9] = { 0. }; //cor_LeLe, cor_LeMi, cor_LeRi, cor_MiLe, cor_MiMi, cor_MiRi, cor_RiLe, cor_RiMi, cor_RiRi;
	virtual ~AdjacentGene() {};
}	adjacentpair;

struct netparam {
	vector<double> *adjGeneScore;
	map<string, vector<string>> *setGeneRefined;
	map<string, double> *fdrAdjusted;
	map<string, int> *geneRefinedIndex;
	map<string, map<string, double>> *genenet;
	map<string, vector<string>> *geneSnpGlobal;
	map<string, int> *lgRsInputIndex;
	vector<double> *lgRsInput;
	map<string, int> *geneInputIndex;

	virtual ~netparam() {};
};


struct initparam {
	bool inputType;			// 0: SNP_INPUT, 1: GENE_INPUT
	bool convert2Symbol;	// 0: default NOT, 1: convert		// for pathway file only
	int	hgVersion;		// 0 : hg18, 1 : hg19, 2 : hg38

	int  minSetSize;
	int  maxSetSize;
	string inputFile;
	string snpGeneMapFile;
	string setFile;			// pathway file
	string adjacentGeneFile;
	string networkFile;		// default: STRING_NETWORK
	double geneCutoff;

	virtual ~initparam() {};
};


struct setinfo {
	string setname;
	int genecount;
	int setsize;
	double zscore;
	double adjZscore;
	double adjPval;
	double adjQval;

	virtual ~setinfo() {};
};


struct resparam {
	vector<double> adjGeneScore;
	map<string, vector<string>> setGeneRefined;
	map<string, double> fdrAdjusted;
	map<string, int> geneRefinedIndex;
	map<string, map<string, double>> genenet;
	map<string, vector<string>> geneSnpGlobal;
	map<string, int> lgRsInputIndex;
	vector<double> lgRsInput;
	map<string, int> geneInputIndex;

	virtual ~resparam() {};
};


class gsa
{
private:

public:
	map<string, vector<string>> setGeneGlobal;
	map<string, vector<string>> geneSnpGlobal;
	map<string, int> geneInputIndex;
	map<string, int> setGlobalActualSize;
	map<string, int> lgRsInputIndex;
	vector<double> lgRsInput;

	map<string, double> geneInputScore;
	vector<string> geneInput;
	vector<string> setGlobal;

	map<string, double> zsOrigin;
	map<string, double>	zsAdjusted;
	map<string, double>	pvAdjusted;
	map<string, double> fdrAdjusted;

	vector<double> orgGeneScore;
	vector<double>		adjGeneScore;
	map<string, vector<string>> setGeneRefined;
	vector<string> geneRefined;
	map<string, int> geneRefinedIndex;
	map<string, string> bestSnpGene; // <gene, best-snp>	

	map<string, int> snpLocMap; // <snp, position>
	map<string, string> famiGene; // <gene, metaGene>
	map<string, adjacentpair> geneInOrder; // <geneA, geneA-geneB-relation>		
	map<string, genelocation> geneLoc; // <gene, start_loc>
	map<double, int>	sorted_pval_index;

	map<string, map<string, double>> genenet;
	int minSetSize, maxSetSize;

public:
	gsa()
	{
		minSetSize = 10;
		maxSetSize = 200;
	};

	~gsa() {

	};

	void refreshAll() {
		if (!geneInputScore.empty()) geneInputScore.clear(), map<string, double>().swap(geneInputScore);
		if (!genenet.empty()) genenet.clear(), map<string, map<string, double>>().swap(genenet);
		if (!sorted_pval_index.empty()) sorted_pval_index.clear(), map<double, int>().swap(sorted_pval_index);
		if (!setGeneGlobal.empty()) setGeneGlobal.clear(), map<string, vector<string>>().swap(setGeneGlobal);
		if (!setGlobalActualSize.empty()) setGlobalActualSize.clear(), map<string, int >().swap(setGlobalActualSize);
		if (!geneSnpGlobal.empty()) geneSnpGlobal.clear(), map<string, vector<string>>().swap(geneSnpGlobal);
		if (!geneInputIndex.empty()) geneInputIndex.clear(), map<string, int>().swap(geneInputIndex);
		if (!lgRsInputIndex.empty()) lgRsInputIndex.clear(), map<string, int>().swap(lgRsInputIndex);
		if (!lgRsInput.empty()) lgRsInput.clear(), vector<double>().swap(lgRsInput);
		if (!geneInput.empty()) geneInput.clear(), vector<string>().swap(geneInput);
		if (!setGlobal.empty()) setGlobal.clear(), vector<string>().swap(setGlobal);
		if (!zsOrigin.empty()) zsOrigin.clear(), map<string, double>().swap(zsOrigin);
		if (!zsAdjusted.empty()) zsAdjusted.clear(), map<string, double>().swap(zsAdjusted);
		if (!pvAdjusted.empty()) pvAdjusted.clear(), map<string, double>().swap(pvAdjusted);
		if (!fdrAdjusted.empty()) fdrAdjusted.clear(), map<string, double>().swap(fdrAdjusted);
		if (!adjGeneScore.empty()) adjGeneScore.clear(), vector<double>().swap(adjGeneScore);
		if (!geneInOrder.empty()) geneInOrder.clear(), map<string, adjacentpair>().swap(geneInOrder); // <geneA, geneA-geneB-relation>	
		if (!snpLocMap.empty()) snpLocMap.clear(), map<string, int>().swap(snpLocMap); // <snp, position>
		if (!geneLoc.empty()) geneLoc.clear(), map<string, genelocation>().swap(geneLoc); // <gene, start_loc>
		if (!setGeneRefined.empty()) setGeneRefined.clear(), map<string, vector<string>>().swap(setGeneRefined);
		if (!geneRefined.empty()) geneRefined.clear(), vector<string>().swap(geneRefined);
		if (!geneRefinedIndex.empty()) geneRefinedIndex.clear(), map<string, int>().swap(geneRefinedIndex);
		if (!bestSnpGene.empty()) bestSnpGene.clear(), map<string, string>().swap(bestSnpGene); // <gene, best-snp>
		if (!famiGene.empty()) famiGene.clear(), map<string, string>().swap(famiGene); // <gene, best-snp>				
	}

public:

	/////////////////////////////////////////////////////////////////////////////////////////////
	void Initializing(initparam initParam) {
		snpPvalueLoad(initParam.inputFile);
		snpGeneMapGenerate(initParam.snpGeneMapFile, initParam.inputType);
		setGeneLoad(initParam.setFile, initParam.convert2Symbol);

		string	hgFile = "";
		string	geneListFile = "";
		switch (initParam.hgVersion)
		{
		case 0:
			hgFile = "data\\rsloc_hg18";
			geneListFile = "data\\hg18GeneList";
			break;
		case 1:
			hgFile = "data\\rsloc_hg19";
			geneListFile = "data\\hg19GeneList";
			break;
		case 2:
			hgFile = "data\\rsloc_hg38";
			geneListFile = "data\\hg38GeneList";
			break;
		default:
			hgFile = "data\\rsloc_hg19";
			geneListFile = "data\\hg19GeneList";
			break;
		}
		snpLocMapGenerate(hgFile);
		geneLocMapGenerate(geneListFile);
		familyGeneLoad();
		adjacentGeneMapLoad(initParam.adjacentGeneFile);
	}

	void Analyzing(initparam initParam)// resparam &resParam) 
	{
		geneFilter();
		geneScoring(initParam.inputType); // original
		adjustedGeneScoring(initParam.inputType); // original
		adjustedPvalueBH();

		/////////////////////////////////////////////////////////////////////////
	/*	resParam.adjGeneScore = adjGeneScore;
		resParam.fdrAdjusted = fdrAdjusted;
		resParam.geneInputIndex = geneInputIndex;
		resParam.geneRefinedIndex = geneRefinedIndex;
		resParam.geneSnpGlobal = geneSnpGlobal;
		resParam.lgRsInput = lgRsInput;
		resParam.lgRsInputIndex = lgRsInputIndex;
		resParam.setGeneRefined = setGeneRefined;
		resParam.genenet = genenet;*/
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	// initialization & pre-processing

	// for all available genes, load the corresponding partial network
	void networkLoading(string network_file = "data\\STRING_NETWORK.txt") // Sora's update
	{
		ifstream fin(network_file);
		if (!fin) {
			cerr << "Cannot access specified file: " << network_file << endl;
			EXIT_FAILURE;
		}

		string inpLine;
		while (getline(fin, inpLine)) {
			istringstream ss(inpLine);
			map<string, double> tmp;
			string genea, geneb;
			double weight;
			ss >> genea >> geneb >> weight;
			weight = weight / 1000; // re-scale to 0 ~ 1
			if (geneInputIndex.find(genea) != geneInputIndex.end() && geneInputIndex.find(geneb) != geneInputIndex.end()) {	// both genes are under consideration
				if (genenet.find(genea) == genenet.end()) {	// not yet recognized before
					tmp[geneb] = weight; // add new edge of ogdf::Graph
					genenet[genea] = tmp; // add new vertex of ogdf::Graph
					tmp.clear();	// reset tmp for next new edge
				}
				else {
					if (genenet[genea].find(geneb) == genenet[genea].end()) {// this is a new edge ~ linkage genea-geneb has not been recognized before
						genenet[genea][geneb] = weight;
					}
					else {	// edge genea-geneb has been recognized before, update if there is a better weight (bigger)
						if (weight > genenet[genea][geneb]) {
							genenet[genea][geneb] = weight;
						}
					}
				}
			} /// end if both genes
		} /// end while		
		fin.close();
	}

	//
	double findMedianOfTopMax(map<double, int> sorted_pval, double topPercentile = 0.10) {
		//	ofstream fou("sorted_logp");
		int asize = (int)round(0.5 * topPercentile * sorted_pval.size());
		map<double, int>::reverse_iterator it = sorted_pval.rbegin();
		for (int i = 0; i < asize; i++, it++) {}
		return it->first;
	}

	//////////////	Truncation 
	double find4threshold(double truncation_threshold) {
		//	ofstream fou("sorted_logp");
		int asize = (int)round(truncation_threshold * sorted_pval_index.size());
		map<double, int>::iterator it = sorted_pval_index.begin();
		for (int i = 0; i < asize; i++, it++) {
			//	fou << it->first << "\t" << it->second << endl;
		}
		//fou.close();
		return it->first;
	}

	void adjustedGeneScoring(bool inputType = SNP_INPUT) // adjusted z-score , default = less FP
	{
		double threshold = 0.25;
		vector<double> newGeneScore;
		vector<double> newGeneSize;
		map<int, vector<double>> valHistogram;
		double upperGS = 0.;
		double lowerGS = 100.;
		vector<double> tmpTop;

		vector<double>().swap(newGeneScore);
		vector<double>().swap(newGeneSize);

		ofstream fou;
		for (auto& vi : geneInput)		// scan for all genes
		//for (auto& vi : geneRefined)		// scan for all genes
		{
			double tmpScore;
			if (inputType == SNP_INPUT) {
				tmpScore = lgRsInput[lgRsInputIndex[bestSnpGene[vi]]];
			}
			else if (inputType == GENE_INPUT) {
				tmpScore = geneInputScore[vi];
			}

			newGeneScore.push_back(tmpScore);	//ordered by geneInput	
			double gPCASize = 0;
			if (gPCASize <= 0)	gPCASize = geneSnpGlobal[vi].size();

			newGeneSize.push_back(gPCASize);
			valHistogram[(int)gPCASize].push_back(tmpScore);
		}
		//fou.close();

		//////////////////////////////////////////////////////////////////////////
		// regression
		//fou.open("lowerGS");
		double outlierBound = 0;
		map<int, double> sizeToSigma;
		map<int, double> sizeToBound;
		bool firstpin = true;
		vector<double> duopin;
		int countpin = 0;
		int firstpinkey = 0;

		for (auto mi : valHistogram) {
			countpin++;
			if (!firstpin || countpin == valHistogram.size()) {
				for (auto vi : mi.second) duopin.push_back(vi);

				double mu = meanv(duopin);	// m
				double si = sdv(duopin);	// s
				double nSigma = mu + 3 * si;

				for (auto vi = valHistogram[firstpinkey].begin(); vi != valHistogram[firstpinkey].end(); vi++) {
					sizeToSigma[firstpinkey] = nSigma;
					if (*vi > nSigma) {	// x >= m + 3s  # outlier - '=' is a must coz of in case single size value, si == 0, mu = single value
											//fou << *vi << "\t" << mu << "\t" << si << "\t" << nSigma << endl;					
											//	eid.push_back(idx);
					}
					else {
						if (upperGS < *vi)
							//fou << "upper: " << upperGS << endl,
							upperGS = *vi;
					}
					//idx++;
				}

				for (auto vi = mi.second.begin(); vi != mi.second.end(); vi++) {
					sizeToSigma[mi.first] = nSigma;
					if (*vi >= nSigma) {	// x >= m + 3s  # outlier - '=' is a must coz of in case single size value, si == 0, mu = single value
											//fou << *vi << "\t" << mu << "\t" << si << "\t" << nSigma << endl;					
					//	eid.push_back(idx);
					}
					else {
						if (upperGS < *vi)
							//fou << "upper: " << upperGS << endl,
							upperGS = *vi;
					}
					//idx++;
				}
				firstpin = true;
			}
			else {
				duopin = mi.second;
				firstpin = false;
				firstpinkey = mi.first;
			}
		}

		if (upperGS > 8.0)	upperGS = 8.0; // sorted_pval_index.end()->first;
		outlierBound = upperGS; // fix - maximum value in lower GS partition
		lowerGS = sorted_pval_index.begin()->first;

		vector<double> topGS;
		vector<double> topSize;
		vector<double> cutoffGS;
		vector<double> cutoffSize;

		bool solved = false;
		vector<double> lastFive;
		double avgGS; // desired cutoff gene score
		int iteration = 0;
		vector<double> coef;
		do {
			avgGS = (upperGS + lowerGS) / 2;

			coef.clear(); vector<double>().swap(coef);
			topGS.clear(); vector<double>().swap(topGS);
			topSize.clear(); vector<double>().swap(topSize);
			cutoffGS.clear(); vector<double>().swap(cutoffGS);
			cutoffSize.clear(); vector<double>().swap(cutoffSize);

			for (int i = 0; i < (int)newGeneScore.size(); ++i) {
				if (newGeneScore[i] > avgGS) {
					if (newGeneScore[i] < outlierBound) {
						topGS.push_back(newGeneScore[i]);
						topSize.push_back(newGeneSize[i]);
					}
				}
				else {
					cutoffGS.push_back(newGeneScore[i]);
					cutoffSize.push_back(newGeneSize[i]);
				}
			}

			vcRegression myreg;
			myreg.Polynomial(topSize, topGS, 1, coef);
			iteration++;

			if (-1e-10 < coef[0] && coef[0] < 1e-10) // coef ~ 0
				solved = true;
			else {
				if (coef[0] < 0) {
					upperGS = avgGS;
				}
				else {
					lowerGS = avgGS;
				}
			}

			if (lastFive.size() < 5) {	// not enough five members
				lastFive.push_back(coef[0]);
			}
			else {
				double tmp = lastFive[lastFive.size() - 1];
				bool allsame = true;
				for (int i = (int)lastFive.size() - 2; i > (int)lastFive.size() - 5; i--) {
					if (tmp != lastFive[i])
						allsame = false;
				}

				if (!allsame && iteration <= 5000) {
					lastFive.push_back(coef[0]);
				}
				else {
					solved = true;
				}
			}
		} while (!solved);

		int maxSize = 0, minSize = 1;
		double maxGS = 0., minGS = 100;

		for (int i = 0; i < (int)cutoffGS.size(); i++) {
			//	fou << cutoffSize[i] << "\t" << cutoffGS[i] << endl;
			if (maxSize < cutoffSize[i]) maxSize = (int)cutoffSize[i];
			if (maxGS < cutoffGS[i]) maxGS = cutoffGS[i];
		}
		//fou.close();

		//fou.open("upperGS");
		for (int i = 0; i < (int)topGS.size(); i++) {
			//	fou << topSize[i] << "\t" << topGS[i] << endl;
			if (maxSize < topSize[i]) maxSize = (int)topSize[i];
		}
		//fou.close();

		//// regression for lower part with optimal cutoff
		//// 1 - logarithmic estimation
		//vcRegression gsreg;
		//coef.clear(); vector<double>().swap(coef);
		//gsreg.Logarithmic(cutoffSize, cutoffGS, coef);
		//fou.open("coef2");
		//fou << coef[0] << "\t" << coef[1] << endl;
		//fou << "y = " << coef[0] << " * ln(x) + " << coef[1] << endl;	
		//fou.close();
		//	
		//fou.open("logall.txt");
		//double p_threshold = find4threshold(threshold);
		//for (int i = 0; i < (int)newGeneScore.size(); i++) {
		//	if (newGeneScore[i] < p_threshold)	newGeneScore[i] = 0.;
		//	if (newGeneScore[i] > outlierBound)	newGeneScore[i] = outlierBound;
		//	double trendline = coef[0] * log(newGeneSize[i]) + coef[1];
		//	fou << newGeneSize[i] << "\t" << trendline << endl;
		//	if (adjzMethod == 0) {		// default = LESS FALSE POSITIVE
		//		newGeneScore[i] -= trendline;
		//	}
		//	else {						// MORE POWER
		//		if (newGeneScore[i] < trendline)
		//			newGeneScore[i] = trendline;
		//	}
		//}
		//fou.close();
		//adjGeneScore = newGeneScore;

		//// 2 - monotonic cubic spline		
		int nCutoff = (int)cutoffGS.size();
		vector<int> segment;
		map<int, int> sizeHistogram;
		for (auto vi : cutoffSize)
			sizeHistogram[(int)vi]++;

		int numSize = sizeHistogram.size();
		int count = 0;
		bool once = false;
		double percentile = 0.05;
		double ks = percentile;
		segment.push_back(0); // for knots 0
		countpin = 1;
		int knot = 0;
		int countk = 1;
		for (auto mi : sizeHistogram) {
			// 2^k ~ 
			knot = (int)pow(2, countk) - 1;
			if (countpin == knot) {
				segment.push_back(mi.first);
				countk++;
			}

			// Fibonacci:
			/*knot = fiboA + fiboB;
			if (countpin == knot) {
			segment.push_back(mi.first);
			knot = fiboB;
			fiboB = fiboA + knot;
			fiboA = knot;
			}		*/

			// count pins			
			//if (countpin % 15 == 1 || countpin == numSize) segment.push_back(mi.first);
			countpin++;

			//if (count < nCutoff * ks) { // count percentile
			//	count += mi.second;
			//	if (count > nCutoff * ks)
			//	{
			//		segment.push_back(mi.first);
			//		ks += percentile;
			//	}
			//}
			//else {
			//	ks += percentile;
			//}
		}

		//////////////////////////////////////////////////////////////////////////
		int nSeg = (int)segment.size();
		vector<myPoint> knots(nSeg + 1); // gene sizes for last knot
		vector<int> bincount(nSeg + 1);
		int range = 3;
		for (int i = 1; i < nSeg; i++)		knots[i].x = segment[i], knots[i].y = 0., bincount[i] = 0; // 1, 2, 3, ...(nSeg - 1)
		knots[0].x = 0.00000000001;	// the first knot
		knots[0].y = 0.00000000001;
		bincount[0] = 1;
		knots[nSeg].x = maxSize;	// the last knot
		knots[nSeg].y = maxGS;
		bincount[nSeg] = 1;
		int tmpsize = 0, yOneCount = 0;

		for (int i = 0; i < (int)cutoffSize.size(); i++) {
			tmpsize = (int)cutoffSize[i];
			for (int k = 1; k < nSeg; k++) {
				if (segment[k] - range / 2 < tmpsize && tmpsize <= segment[k] + range / 2) {
					knots[k].y += cutoffGS[i];
					bincount[k]++;
				}
			}
		}
		for (int i = 1; i < nSeg; i++) {
			knots[i].y /= (double)bincount[i];
		}

		vector<myPoint> upperKnots; // knots for upper line
		vector<myPoint> lowerKnots; // knots for lower line

									// upper knots
		upperKnots.push_back(knots[0]); // 1st knot
		double previousScore = knots[0].y;
		for (int i = 1; i < nSeg; i++) {
			if (knots[i].y > previousScore) {
				upperKnots.push_back(knots[i]);
				previousScore = knots[i].y;
			}
		}
		upperKnots.push_back(knots[nSeg]); // last knot

										   // lower knots
		lowerKnots.push_back(knots[nSeg]); // last knot
		previousScore = knots[nSeg].y;
		for (int i = nSeg - 1; i > 0; i--) {
			if (knots[i].y < previousScore) {
				lowerKnots.insert(lowerKnots.begin(), knots[i]);
				previousScore = knots[i].y;
			}
		}
		lowerKnots.insert(lowerKnots.begin(), knots[0]); // 1st knot

		vcRegression gsreg, upperReg, lowerReg;
		upperReg.monotoneCubicSplineInit(upperKnots);
		lowerReg.monotoneCubicSplineInit(lowerKnots);

		//predicted median knots
		vector<myPoint> mediumKnots; // knots for median line
		mediumKnots.push_back(knots[0]); // 1st knot
		myPoint tmp;
		for (int i = 1; i < nSeg; i++) {
			tmp.x = knots[i].x;
			tmp.y = 0.5 * (upperReg.monotoneCubicSplineInterp(tmp.x) + lowerReg.monotoneCubicSplineInterp(tmp.x));
			mediumKnots.push_back(tmp);
		}
		mediumKnots.push_back(knots[nSeg]); // last knot
		gsreg.monotoneCubicSplineInit(mediumKnots);

		double p_threshold = 0.;
		if (inputType == SNP_INPUT)
			p_threshold = find4threshold(threshold);

		for (int i = 0; i < (int)newGeneScore.size(); i++) {
			if (newGeneScore[i] < p_threshold)	newGeneScore[i] = 0.;
			if (newGeneScore[i] > outlierBound)	newGeneScore[i] = outlierBound;
			double trendline = gsreg.monotoneCubicSplineInterp(newGeneSize[i]);

			//if (adjzMethod == 0) 
			{		// LESS FALSE POSITIVE = default
				newGeneScore[i] -= trendline;
			}
			//else {						// MORE POWER
			//	if (newGeneScore[i] < trendline)
			//		newGeneScore[i] = trendline;
			//}
		}
		adjGeneScore = newGeneScore;

		//////////////////////////////////////////////////////////////////////////////////
		// for set		
		if (inputType == GENE_INPUT) {
			setGeneRefined = setGeneGlobal;
			geneRefined = geneInput;
			geneRefinedIndex = geneInputIndex;
		}

		double muPop = meanv(adjGeneScore);
		double sigPop = sdv(adjGeneScore);

		for (auto& mi : setGeneRefined) {
			vector<double> genePvalue;

			for (auto& vi : mi.second) {
				double curGeneScore = adjGeneScore[geneInputIndex[vi]];
				genePvalue.push_back(curGeneScore);
			}
			zsAdjusted[mi.first] = zsScoring(genePvalue, muPop, sigPop, (int)adjGeneScore.size());
			pvAdjusted[mi.first] = z2p(zsAdjusted[mi.first]);
		}

	}

	//////////////////////////////////////////////////////////////////////	
	void geneFilter(double inputGeneCorrelation = 0.1)
	{
		//ofstream fou("genefilter");		
		// find representative gene if exist
		map<string, bool> geneRefinedtmp;
		for (auto gi : geneInput) {
			geneRefinedtmp[gi] = 1;
		}

		int newIdx = 0;
		for (auto gi : geneRefinedtmp) {
			geneRefined.push_back(gi.first);
			geneRefinedIndex[gi.first] = newIdx++;
		}

		// adjacent gene filter
		for (auto sg : setGeneGlobal)
		{
			// set-gene map: 1st phase: representative assigning
			// 2nd phase: sorting gene in chromosome order				
			map<int, string> refinedgeneinorder; // <gene-loc, genelocation> sorted ascending
			int aliasLoc = 1;
			//for (auto gi : sg.second) {// gene list
			for (int i = 0; i < (int)sg.second.size(); i++) // gene list	
			{
				string gi = sg.second[i];

				//////////////////////////////////////////////////////////////////////
				if (famiGene.find(gi) != famiGene.end()) { // there is a representative for this gene
					if (bestSnpGene[famiGene[gi]] == bestSnpGene[gi]) // represented by the same SNP
					{
						bool famiFound = false; // = true if they are in the same set/pathway
						for (int j = 0; j < (int)sg.second.size(); j++) {
							if (j != i) {
								string fakeFami = sg.second[j];
								if (famiGene[gi].compare(fakeFami) == 0) { // identical, they are in the same set
									famiFound = true;
									break;
								}
							}
						}
						if (famiFound) {
							//fou << "familiy gene removed: " << gi << endl; 
							continue;
						}	// not count gi						
					}
				}

				//////////////////////////////////////////////////////////////////////////					
				if (refinedgeneinorder.find(geneLoc[gi].start) == refinedgeneinorder.end()) // same start
					refinedgeneinorder[geneLoc[gi].start] = gi;
				else
				{
					//	fou << "start duplicated! " << refinedgeneinorder[geneLoc[gi].start] << "\t" << gi << endl;
					string previousGene = refinedgeneinorder[geneLoc[gi].start]; // gene which has the same start
					int previousLength = geneLoc[previousGene].end - geneLoc[previousGene].start;
					int giLength = geneLoc[gi].end - geneLoc[gi].start;

					if (giLength < previousLength) { // the shorter, the sooner (<) - later (>)
						geneLoc[previousGene].start = geneLoc[gi].start + aliasLoc; // refined gene location
						refinedgeneinorder[geneLoc[gi].start + aliasLoc] = previousGene;
						refinedgeneinorder[geneLoc[gi].start] = gi;
					}
					else {
						geneLoc[gi].start = geneLoc[gi].start + aliasLoc; // refined gene location
						refinedgeneinorder[geneLoc[gi].start + aliasLoc] = gi;
					}

					aliasLoc++;
				}
				//fou << geneLoc[gi].start << "\t" << gi << endl;				
			}

			vector<string> sortedRefinedGene;
			for (auto li : refinedgeneinorder) {
				sortedRefinedGene.push_back(li.second);
				//fou << li->second << endl;
			}

			// 3rd: checking correlation of 2 adjacent genes			
			vector<string> finalRefinedGene(0);
			for (int i = 0; i < (int)sortedRefinedGene.size() - 1; i++)
			{
				finalRefinedGene.push_back(sortedRefinedGene[i]);			// always keep current gene
				if (geneInOrder.find(sortedRefinedGene[i]) != geneInOrder.end()) // there is some neighbor
				{
					if (geneInOrder[sortedRefinedGene[i]].geneB.geneName.compare(sortedRefinedGene[i + 1]) == 0) // equal~ adjacent gene detected
					{
						// find best SNP of geneA, geneB	
						string	snpA = bestSnpGene[sortedRefinedGene[i]];
						string	snpB = bestSnpGene[sortedRefinedGene[i + 1]];
						snpA.replace(0, 2, "");
						snpB.replace(0, 2, "");

						// located snpA, snpB against (LL, LM, LR, ML, MM, MR, RL, RM, RR)
						adjacentpair abPair = geneInOrder[sortedRefinedGene[i]];
						int locA = snpLocMap[snpA]; // they should be valid, God please again!!
						int locB = snpLocMap[snpB];
						if (locA != 0 && locB != 0) {
							const int	corrLoc[3][3] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
							//	fou << "locA: " << locA << "\t" << "locB: " << locB << endl;

							int minidxA = 0;
							int minDiffA = abs(locA - snpLocMap[abPair.geneA.snpL]);
							if (minDiffA > abs(locA - snpLocMap[abPair.geneA.snpM]))		minidxA = 1, minDiffA = abs(locA - snpLocMap[abPair.geneA.snpM]);
							if (minDiffA > abs(locA - snpLocMap[abPair.geneA.snpR]))		minidxA = 2;

							int minidxB = 0;
							int minDiffB = abs(locB - snpLocMap[abPair.geneB.snpL]);
							if (minDiffB > abs(locB - snpLocMap[abPair.geneB.snpM]))		minidxB = 1, minDiffB = abs(locB - snpLocMap[abPair.geneB.snpM]); ;
							if (minDiffB > abs(locB - snpLocMap[abPair.geneB.snpR]))		minidxB = 2;

							double corrAB = abPair.corr[corrLoc[minidxA][minidxB]];
							if (corrAB >= inputGeneCorrelation) { // high correlation detected, remove B							
								//fou << "removed: " << sortedRefinedGene[i + 1] << "\t" << lgRsInput[lgRsInputIndex[bestSnpGene[sortedRefinedGene[i + 1]]]] << endl;								
								i++; // step over geneB, move to next gene
							}
							else
							{
								if (i == (int)sortedRefinedGene.size() - 2)		   // if current gene is the last gene, continue/ bypass;
									finalRefinedGene.push_back(sortedRefinedGene[i + 1]);			// always keep current gene
							}
						}
						else
							if (i == (int)sortedRefinedGene.size() - 2)		   // if current gene is the last gene, continue/ bypass;
								finalRefinedGene.push_back(sortedRefinedGene[i + 1]);			// always keep current gene
					}
					else
						if (i == (int)sortedRefinedGene.size() - 2)		   // if current gene is the last gene, continue/ bypass;
							finalRefinedGene.push_back(sortedRefinedGene[i + 1]);			// always keep current gene
				}
				else
					if (i == (int)sortedRefinedGene.size() - 2)		   // if current gene is the last gene, continue/ bypass;
						finalRefinedGene.push_back(sortedRefinedGene[i + 1]);			// always keep current gene				
			}	/// end of for (int i = 0; i < (int)sortedRefinedGene.size() - 1; i++)

				// 4th: update setGeneRefined
			setGeneRefined[sg.first] = finalRefinedGene;
			/*if (sg.second.size() != setGeneRefined[sg.first].size()) {
				countrefinedSet++;*/
				//fou << "Gene removed: " << sg.first << "\t" << sg.second.size() - finalRefinedGene.size() << endl;
		//}
		}	/// end of for (auto sg : setGeneGlobal)		

	}	/// end of function geneFilter

	void familyGeneLoad(string filename = "data\\FamilyGeneList")
	{
		ifstream fin(filename);
		if (!fin.is_open())
		{
			cerr << "Unavailable data: " << filename << " is missed!" << endl;
			EXIT_FAILURE;
		}

		string inpLine;								// remove 1 meta lines
		for (int j = 0; j < 1; ++j)
		{
			getline(fin, inpLine);
		}

		while (getline(fin, inpLine))
		{
			string agene, metagene, dumbdata;
			stringstream iss(inpLine);
			iss >> dumbdata >> agene >> dumbdata >> dumbdata >> dumbdata >> metagene;
			famiGene[agene] = metagene;
		}
		fin.close();
	}

	void geneLocMapGenerate(string filename = "data\\hg19GeneList") {
		ifstream fin(filename);
		if (!fin.is_open())
		{
			cerr << "Not all data are available: " << filename << " is missed!" << endl;
			EXIT_FAILURE;
		}

		string inpLine;
		// remove 1 meta lines
		for (int j = 0; j < 1; ++j)
		{
			getline(fin, inpLine);
		}

		while (getline(fin, inpLine))
		{
			genelocation agene;
			stringstream iss(inpLine);
			iss >> agene.geneName >> agene.chromosome >> agene.start >> agene.end;
			geneLoc[agene.geneName] = agene;
		}
		fin.close();
	}

	void adjacentGeneMapLoad(string filename = "data\\EUR_Adjacent_correlation") {
		ifstream fin(filename);
		if (!fin.is_open())
		{
			cerr << "Not all data are available: " << filename << " is missed!" << endl;
			EXIT_FAILURE;
		}

		string inpLine;
		// remove 1 meta lines
		for (int j = 0; j < 1; ++j)
		{
			getline(fin, inpLine);
		}

		while (getline(fin, inpLine))
		{
			adjacentpair apair;
			stringstream iss(inpLine);
			iss >> apair.geneA.geneName >> apair.geneB.geneName
				>> apair.geneA.snpL >> apair.geneA.snpM >> apair.geneA.snpR
				>> apair.geneB.snpL >> apair.geneB.snpM >> apair.geneB.snpR;
			for (int i = 0; i < 9; i++) {
				iss >> apair.corr[i];
			}
			geneInOrder[apair.geneA.geneName] = apair;
		}
		fin.close();
	}

	/*const string allChromosome[22] = { "chr22.txt", "chr21.txt", "chr20.txt" ,
		"chr19.txt", "chr18.txt", "chr17.txt", "chr16.txt", "chr15.txt", "chr14.txt", "chr13.txt", "chr12.txt", "chr11.txt",
		"chr10.txt", "chr9.txt", "chr8.txt", "chr7.txt", "chr6.txt", "chr5.txt", "chr4.txt", "chr3.txt", "chr2.txt", "chr1.txt" };*/

	void snpLocMapGenerate(string filename = "data\\rsloc_hg19") {
		ifstream fin;
		string inpLine, curSnp;
		fin.open(filename);

		if (!fin.is_open())
		{
			cerr << "Unavailable data: " << filename << " is missed!" << endl;
			EXIT_FAILURE;
		}

		snpLocMap.clear();
		std::map<string, int>().swap(snpLocMap);
		////////////////////////////////////////////////////////			

		while (getline(fin, inpLine))
		{
			string chromosome, curSnp;
			int position;
			stringstream iss(inpLine);
			//iss >> chromosome >> position >> curSnp;
			iss >> curSnp >> position;
			//curSnp.replace(0, 2, "");	// remove "rs"
			snpLocMap[curSnp] = position;
			//fou << curSnp << "\t" << position << "\n";
		}
		////////////////////////////////////////////////////////
		fin.close();
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	void geneScoring(bool inputType = SNP_INPUT)
	{
		vector<double> newGeneScore;

		if (inputType == SNP_INPUT) {
			for (auto& vi : geneInput)		// scan for all genes
			{
				double tmpScore = -MAXINT32;
				for (auto& svi : geneSnpGlobal[vi]) // for corresponding list of SNPs
				{
					if (lgRsInputIndex.find(svi) != lgRsInputIndex.end()) /// some not welcome snp may appear!!!
					{
						// find the best (greatest) log score				
						if (lgRsInput[lgRsInputIndex[svi]] > tmpScore)	// if corresponding log score is greater 
						{
							tmpScore = lgRsInput[lgRsInputIndex[svi]];
							bestSnpGene[vi] = svi;
						}
					}
				}
				newGeneScore.push_back(tmpScore);	//ordered by geneInput		
			}
		}
		else if (inputType == GENE_INPUT) {
			for (auto& vi : geneInput)		// scan for all genes
			{
				double tmpScore = geneInputScore[vi];
				newGeneScore.push_back(tmpScore);	//ordered by geneInput			
			}
		}
		//////////////////////////////////////////////////////////////////////////////////
		// for set
		orgGeneScore = newGeneScore;///geneProfile[0];

		double muPop = meanv(orgGeneScore);
		double sigPop = sdv(orgGeneScore);
		for (auto& mi : setGeneGlobal)
		{
			vector<double> genePvalue;
			for (auto& vi : mi.second)
			{
				double curGeneScore = orgGeneScore[geneInputIndex[vi]];
				if (geneInputIndex.find(vi) != geneInputIndex.end()) // only genes included in sets ~ USEFUL GENES
				{
					genePvalue.push_back(curGeneScore);
				}
			}
			zsOrigin[mi.first] = zsScoring(genePvalue, muPop, sigPop);
		}
	}

	void setGeneLoad(string setFile = "data\\c2.cp.v5.2.symbols.gmt", bool convert2symbol = 0)
	{
		string set_gene_GO = "data\\Gene Ontology";
		string set_gene_c5all = "data\\c5.all.v5.0.symbols.gmt";
		string set_gene_KEGG = "data\\KEGG";


		ifstream fin(setFile);
		if (!fin.is_open())
		{
			cerr << "Cannot access specified file: " << setFile << endl;
			EXIT_FAILURE;
		}

		auto_ptr<map<string, string>> ensemblID(new map<string, string>);	//  <ensemblID, symbol>
		auto_ptr<map<string, string>> entrezID(new map<string, string>);		//  <entrezID, symbol>
		auto_ptr<map<string, string>> uniswissprotID(new map<string, string>);	// <uniswissprotID, symbol>
		if (convert2symbol) { // load conversion table

			string id2gene = "data\\id2gene";  // convert other ID type to symbol type
			ifstream ftable(id2gene);
			if (!ftable.is_open())
			{
				cerr << "Cannot access conversion table: " << id2gene << endl;
				EXIT_FAILURE;
			}

			string aline;
			getline(ftable, aline);		// header: ensembl_gene_id >>	entrezgene	>>	external_gene_name	>>	uniprot_swissprot
			while (getline(ftable, aline)) {
				string ensembl, entrez, genesymbol, uniswiss;
				stringstream tmp(aline);
				tmp >> ensembl >> entrez >> genesymbol >> uniswiss;

				(*ensemblID)[ensembl] = genesymbol;

				if (entrez.compare("NA") != 0) { // available 
					(*entrezID)[entrez] = genesymbol;
				}

				if (uniswiss.length() > 0) { // available 
					(*uniswissprotID)[uniswiss] = genesymbol;
				}
			}
			ftable.close();
		}

		string inpLine;
		while (getline(fin, inpLine))
		{
			string tmpSetName;
			int realSetSize = 0;

			if (inpLine.find("hsa") != string::npos) // KEGG: remove 
			{
				inpLine.replace(inpLine.find_first_of("\t"), 1, " : ");
			}

			// pick set name
			tmpSetName = inpLine.substr(0, inpLine.find_first_of("\t"));
			inpLine.replace(0, inpLine.find_first_of("\t") + 1, "");	// remove setGlobal name from current stream

			if (inpLine.find("http") != string::npos)	// c5all: remove link to group "http://www....."
			{
				inpLine.replace(0, inpLine.find_first_of("\t") + 1, "");
			}

			vector<string> term;
			string s, tmps;
			while (inpLine.find("\t") != string::npos)
			{
				// pick gene terms
				tmps = inpLine.substr(0, inpLine.find_first_of("\t"));
				realSetSize++;

				/// check conversion
				s = tmps;
				if (convert2symbol) {
					if ((*ensemblID).find(tmps) != (*ensemblID).end()) { // ensemblID
						s = (*ensemblID)[tmps];
					}

					if ((*entrezID).find(tmps) != (*entrezID).end()) { // entrezID
						s = (*entrezID)[tmps];
					}

					if ((*uniswissprotID).find(tmps) != (*uniswissprotID).end()) { // uniswissprotID
						s = (*uniswissprotID)[tmps];
					}
				}

				//////////////////////////////////////////////////////////////////////////				
				if (geneInputIndex.find(s) != geneInputIndex.end()) // available genes
				{
					term.push_back(s);	// geneGlobal terms
				}

				inpLine.replace(0, inpLine.find_first_of("\t") + 1, "");
			}

			if (inpLine.length() > 0)
			{
				realSetSize++;
				if (geneInputIndex.find(inpLine) != geneInputIndex.end()) // available genes
				{
					term.push_back(inpLine);	// the last term having no "\t"-trait
				}
			}

			setGlobalActualSize[tmpSetName] = realSetSize;

			// set trunking minSetSize, maxSetSize
			if ((minSetSize <= (int)term.size()) && ((int)term.size() <= maxSetSize))
				//if ((minSetSize <= realSetSize) && (realSetSize <= maxSetSize))
			{
				setGlobal.push_back(tmpSetName);
				setGeneGlobal[tmpSetName] = term; // geneGlobal list of setGlobal[i]
			}

			tmpSetName.clear();
			realSetSize = 0;
			term.clear();				// clear current terms			
		}

		fin.close();
	}

	void snpGeneMapGenerate(string snp_gene_file = "data\\db19_20k", bool inputType = SNP_INPUT)
	{
		auto_ptr<map<string, vector<string>>> snpGeneGlobal(new map<string, vector<string>>);
		ifstream fin(snp_gene_file);
		if (!fin.is_open())
		{
			cerr << "Cannot access specified file: " << snp_gene_file << endl;
			EXIT_FAILURE;
		}

		if (inputType == SNP_INPUT) {
			string inpLine;
			while (getline(fin, inpLine))
			{
				while (inpLine.find(",") != string::npos)
				{
					inpLine.replace(inpLine.find_first_of(","), 1, "\t");
				}
				vector<string> term;
				string rs;
				string  s;

				rs = inpLine.substr(0, inpLine.find_first_of("\t"));
				rs.insert(0, "rs");
				inpLine.replace(0, inpLine.find_first_of("\t") + 1, "");

				if (lgRsInputIndex.find(rs) != lgRsInputIndex.end()) // available rsInput				
				{
					while (inpLine.find("\t") != string::npos)
					{
						s = inpLine.substr(0, inpLine.find_first_of("\t"));
						term.push_back(s);	// geneGlobal terms				
						inpLine.replace(0, inpLine.find_first_of("\t") + 1, "");
					}
					if (inpLine.length() > 0)
						term.push_back(inpLine);	// the last term having no "\t"-trait;
													//snpGeneGlobal.push_back(term); // geneGlobal list of rsGlobal[i]
					(*snpGeneGlobal)[rs] = term;
					term.clear();				// clear current terms			
				}
			}
			fin.close();

			/// with each gene, scan for available SNP, which means available p-val/ -log(p-val)
			/// if there's no SNP available, remove gene	
			// all genes include at least ONE snp because we only read in available snp gene lists
			for (auto& mi : *snpGeneGlobal)
			{
				for (auto& vi : mi.second)	// gene list
				{
					geneSnpGlobal[vi].push_back(mi.first);

					if (geneInputIndex.find(vi) == geneInputIndex.end()) // not yet included
					{
						geneInput.push_back(vi);
						geneInputIndex[vi] = geneInput.size() - 1;
					}
				}
			}
		}
		else if (inputType == GENE_INPUT)
		{
			string inpLine;
			while (getline(fin, inpLine))
			{
				while (inpLine.find(",") != string::npos)
				{
					inpLine.replace(inpLine.find_first_of(","), 1, "\t");
				}
				vector<string> term;
				string rs;
				string agene;

				rs = inpLine.substr(0, inpLine.find_first_of("\t"));
				rs.insert(0, "rs");
				inpLine.replace(0, inpLine.find_first_of("\t") + 1, "");

				while (inpLine.find("\t") != string::npos)
				{
					agene = inpLine.substr(0, inpLine.find_first_of("\t"));
					if (geneInputIndex.find(agene) != geneInputIndex.end()) { // gene input
						geneSnpGlobal[agene].push_back(rs);
					}
					inpLine.replace(0, inpLine.find_first_of("\t") + 1, "");
				}
				if (inpLine.length() > 0)
					if (geneInputIndex.find(inpLine) != geneInputIndex.end()) // gene input
						geneSnpGlobal[inpLine].push_back(rs);	// the last term having no "\t"-trait;												
			}
			fin.close();
		}
	}

	void genePvalueLoad(string gene_pval_file = "data\\DIAGRAMgene")
	{
		ifstream fin(gene_pval_file);
		if (!fin.is_open())
		{
			cerr << "Cannot access specified file: " << gene_pval_file << endl;
			EXIT_FAILURE;
		}

		string curGene, tmpG;
		string pval;
		double tmpP;

		auto_ptr<map<string, string>> ensemblID(new map<string, string>);	//  <ensemblID, symbol>
		auto_ptr<map<string, string>> entrezID(new map<string, string>);		//  <entrezID, symbol>
		auto_ptr<map<string, string>> uniswissprotID(new map<string, string>);	// <uniswissprotID, symbol>
		string id2gene = "data\\id2gene";  // convert other ID type to symbol type
		ifstream ftable(id2gene);
		if (!ftable.is_open())
		{
			cerr << "Cannot access conversion table: " << id2gene << endl;
			EXIT_FAILURE;
		}

		string aline;
		getline(ftable, aline);		// header: ensembl_gene_id >>	entrezgene	>>	external_gene_name	>>	uniprot_swissprot
		while (getline(ftable, aline)) {
			string ensembl, entrez, genesymbol, uniswiss;
			stringstream tmp(aline);
			tmp >> ensembl >> entrez >> genesymbol >> uniswiss;

			(*ensemblID)[ensembl] = genesymbol;

			if (entrez.compare("NA") != 0) { // available 
				(*entrezID)[entrez] = genesymbol;
			}

			if (uniswiss.length() > 0) { // available 
				(*uniswissprotID)[uniswiss] = genesymbol;
			}
		}
		ftable.close();

		//ofstream tfo("DIAGRAMgene");
		while (fin >> tmpG) {
			fin >> tmpP;
			if (tmpP > 0) {
				if ((*ensemblID).find(tmpG) != (*ensemblID).end()) { // ensemblID
					curGene = (*ensemblID)[tmpG];
				}
				else
					if ((*entrezID).find(tmpG) != (*entrezID).end()) { // entrezID
						curGene = (*entrezID)[tmpG];
					}
					else
						if ((*uniswissprotID).find(tmpG) != (*uniswissprotID).end()) { // uniswissprotID
							curGene = (*uniswissprotID)[tmpG];
						}
						else
							curGene = tmpG; // symbol

				geneInput.push_back(curGene);
				geneInputIndex[curGene] = geneInput.size() - 1;
				geneInputScore[curGene] = -log10(tmpP);
				//tfo << curGene << "\t" << pow(10, -tmpP) << endl;
			}
		}
		fin.close();
		//tfo.close();
	}

	void snpPvalueLoad(string snp_pval_file = "data\\DIAGRAM")//example_height_snp")
	{
		ifstream fin(snp_pval_file);
		if (!fin.is_open())
		{
			cerr << "Cannot access specified file: " << snp_pval_file << endl;
			EXIT_FAILURE;
		}

		string curSnp;
		string pval, dumbdata;
		double tmpP;

		while (fin >> curSnp) {
			//	identify p-value
			fin >> pval;
			if (tmpP = atof(pval.c_str())) {
				if (tmpP <= 0) {
					cerr << "Invalid p-value!!!" << endl;
					EXIT_FAILURE;
				}

				tmpP = -log10(tmpP);
				lgRsInput.push_back(tmpP);		// list of -log(p-value)
				lgRsInputIndex[curSnp] = lgRsInput.size() - 1; // tmp indexing
				/*if (sorted_pval_index.find(tmpP) != sorted_pval_index.end()) {
					tmpP += countesp * myesp;
					countesp++;
				}*/
				sorted_pval_index[tmpP] = lgRsInput.size() - 1;
			}
		}
		fin.close();
	}

	string setSubNetDraw(netparam *param, string setName, double cutoff) {
		string filePath = "network";
		cutoff = -log10(cutoff);
		vector<string> minList;
		map<string, bool> tmpnode;
		map<string, map<string, double>> subnet;

		if ((*(param->setGeneRefined)).find(setName) != (*(param->setGeneRefined)).end()) { // exist
			vector<string> geneList = (*(param->setGeneRefined))[setName];
			map<string, int>	nodeDegree;
			map<string, double> curNode;

			for (int i = 0; i < (int)geneList.size(); i++) {
				if ((*(param->geneRefinedIndex)).find(geneList[i]) != (*(param->geneRefinedIndex)).end()) {
					if ((*(param->adjGeneScore))[(*(param->geneInputIndex))[geneList[i]]] >= cutoff) {
						tmpnode[geneList[i]] = 1;
						for (int j = 0; j < (int)geneList.size(); j++) {
							if (j != i)
								if ((*(param->adjGeneScore))[(*(param->geneInputIndex))[geneList[j]]] >= cutoff) // too strict
									if ((*(param->geneRefinedIndex)).find(geneList[j]) != (*(param->geneRefinedIndex)).end()) {
										if ((*(param->genenet))[geneList[i]].find(geneList[j]) != (*(param->genenet))[geneList[i]].end()) { // having a linkage between genes								
											curNode[geneList[j]] = (*(param->genenet))[geneList[i]][geneList[j]];
											nodeDegree[geneList[i]] ++;
											nodeDegree[geneList[j]] ++;
											tmpnode[geneList[j]] = 1;
										} // else nothing happen											
									}
						}
						subnet[geneList[i]] = curNode;
						curNode.clear();
					}
				}
			}

			for (auto tn : tmpnode) {
				minList.push_back(tn.first);
			}
			//////////////////////////////////////////////////////////////////////////						
			const int minDegree = 2;
			vector<int> eid;
			int idx = 0;
			for (auto vi : geneList) {
				if (nodeDegree[vi] < minDegree) {
					eid.push_back(idx);
				}
				idx++;
			}
			while (eid.size() > 0) {
				geneList.erase(geneList.begin() + eid.back());
				eid.pop_back();
			}

			////	SNP-GENE INFORMATION	//////////////////////////////////////////////////////////////////////
			const string geneTag1 = "<div id='myModal";
			const string geneTag2 = "' class='modal'>  <div class='modal-content'><div class='modal-header'><span id='close";
			const string geneTag3 = "' class='close'>&times;</span><h2>Gene-SNPs Reference</h2></div><div class='modal-body'><table><col width='200'><col width='auto'><tbody><tr><td><p><b>Gene<b></p></td><td><a href = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene=";
			const string geneTag4 = "'><b>";
			const string geneTag5 = "</b></a><i>("; // +gene score
			const string geneTag6 = ")</i></td></tr><tr><td colspan='2'>&nbsp; </td></tr><tr><td><p><b>SNPs</b> <i>(p_value &le; 0.05) </i> </p></td> <td><p>";

			const string snpTag1 = "<a href = 'https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=";
			const string snpTag2 = "'>";
			const string snpTag3 = "</a><i>(";
			const string snpTag4 = ")</i>&nbsp; ";

			const string subGeneTag1 = "</p></td></tr><tr><td colspan='2'>&nbsp; </td></tr><tr><td><p><b>Neighboring genes </b><i>(interaction weight)</i></p></td> <td><p>";
			const string subGeneTag2 = "<a href = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene="; // + gene +
			const string subGeneTag3 = "'><b>"; // + gene +
			const string subGeneTag4 = "</b></a><i>("; // + gene score +
			const string subGeneTag5 = ")</i>&nbsp; ";

			const string geneTagEnd = "</p></td></tr></tbody></table></div><div class='modal-footer'><h4>" + setName + "</h4></div></div></div></br>";

			string snpGeneInfo = "";
			string agene = "";

			for (auto& vi : minList)		// scan for all genes in minimal List
			{
				double genescore = pow(10, -((*(param->adjGeneScore))[(*(param->geneInputIndex))[vi]]));
				string strGeneScore = to_string(genescore);
				agene = geneTag1 + vi + geneTag2 + vi + geneTag3 + vi + geneTag4 + vi + geneTag5 + strGeneScore + geneTag6;

				map<double, string> snpval;
				double esp = 0.000000001;
				int k = 0;
				for (auto& svi : (*(param->geneSnpGlobal))[vi]) // for corresponding list of SNPs
				{
					if ((*(param->lgRsInputIndex)).find(svi) != ((*(param->lgRsInputIndex))).end()) /// some not welcome snp may appear!!!
					{
						if ((*(param->lgRsInput))[(*(param->lgRsInputIndex))[svi]] >= cutoff) {// significant SNPs
							if (snpval.find((*(param->lgRsInput))[(*(param->lgRsInputIndex))[svi]]) == snpval.end()) {// not exist
								snpval[(*(param->lgRsInput))[(*(param->lgRsInputIndex))[svi]]] = svi;
							}
							else {
								snpval[(*(param->lgRsInput))[(*(param->lgRsInputIndex))[svi]] + ++k * esp] = svi;
							}
						}
					}
				}

				for (map<double, string>::reverse_iterator sv = snpval.rbegin(); sv != snpval.rend(); sv++) {
					double tp = pow(10, -sv->first);
					string strTp = to_string(tp);
					agene += snpTag1 + sv->second + snpTag2 + sv->second + snpTag3 + strTp + snpTag4;
				}
				agene += subGeneTag1;

				//////////////////////////////////////////////////////////////////////////
				map<double, string> geneweight;
				k = 0;

				for (auto& svi : minList) { // find neighbor in current minimal list
					if (svi.compare(vi) != 0) // not identical
					{
						if (subnet[vi].find(svi) != subnet[vi].end()) // having link from vi to svi
						{
							if (geneweight.find(subnet[vi][svi]) != geneweight.end()) // not exist
								geneweight[subnet[vi][svi]] = svi;
							else {
								geneweight[subnet[vi][svi] + ++k * esp] = svi;
							}
						}
						else
							if (subnet[svi].find(vi) != subnet[svi].end()) // having link from vi to svi
							{
								if (geneweight.find(subnet[svi][vi]) != geneweight.end()) // not exist
									geneweight[subnet[svi][vi]] = svi;
								else {
									geneweight[subnet[svi][vi] + ++k * esp] = svi;
								}
							}
					}
				}

				for (map<double, string>::reverse_iterator gw = geneweight.rbegin(); gw != geneweight.rend(); gw++) {
					string strEdgeWeight = to_string(gw->first);
					agene += subGeneTag2 + gw->second + subGeneTag3 + gw->second + subGeneTag4 + strEdgeWeight + subGeneTag5;
				}
				//////////////////////////////////////////////////////////////////////////
				agene += geneTagEnd;
				snpGeneInfo += agene;
			}

			/////////////////////////////////////////////////////////
			//double curGeneScore = orgGeneScore[geneInputIndex[<geneName>]];	 // NODE WEIGHT
			//filePath = setSubNetBuilding("geneset_test", geneList, subnet, orgGeneScore, geneInputIndex);
			ogdf::Graph G;
			ogdf::GraphAttributes GA(G, ogdf::GraphAttributes::nodeGraphics |
				ogdf::GraphAttributes::nodeStyle |
				ogdf::GraphAttributes::nodeLabel |
				ogdf::GraphAttributes::edgeGraphics |
				ogdf::GraphAttributes::edgeArrow |
				ogdf::GraphAttributes::edgeStyle);
			map <string, ogdf::node> PROTEIN_MAP;

			graph_construct(minList, subnet, PROTEIN_MAP, G, GA);
			node_colorize(G, GA, *(param->adjGeneScore), *(param->geneInputIndex), PROTEIN_MAP);
			node_setsize(G, GA);
			edge_adjust(G, GA);
			ogdf::SugiyamaLayout SL;
			SL.call(GA);
			ogdf::SpringEmbedderFR SE;
			SE.call(GA);

			ogdf::GraphIO::SVGSettings svgs = ogdf::GraphIO::SVGSettings();
			svgs.fontSize(15);
			svgs.fontColor("#000000");
			filePath = setName + ".html";
			ofstream fou(filePath);
			ogdf::GraphIO::drawSVGe(GA, fou, svgs, snpGeneInfo, setName);
			if (fou.is_open())	fou.close();

			PROTEIN_MAP.clear();
			map <string, ogdf::node>().swap(PROTEIN_MAP);
		}
		// finishing			
		subnet.clear();
		map<string, map<string, double>>().swap(subnet);
		//////////////////////////////////////////////////////////////////////////
		return filePath;
	}

	string coreNetDraw(netparam *param, string setName, double qvalcutoff, double cutoff)
	{
		string filePath = "network";
		cutoff = -log10(cutoff);
		vector<string> minList(0);
		map<string, bool> tmpnode;

		for (auto& si : *(param->setGeneRefined)) {
			if ((*(param->setGeneRefined)).find(si.first) != (*(param->setGeneRefined)).end()) // exist
				if ((*(param->fdrAdjusted))[si.first] <= qvalcutoff)
				{
					vector<string> geneList = (*(param->setGeneRefined))[si.first];
					for (auto& vi : geneList)
						if ((*(param->geneRefinedIndex)).find(vi) != (*(param->geneRefinedIndex)).end())
							if ((*(param->adjGeneScore))[(*(param->geneInputIndex))[vi]] >= cutoff)
								tmpnode[vi] = 1;
				}
		}

		for (auto ti : tmpnode)		minList.push_back(ti.first);

		map<string, map<string, double>> subnet;
		map<string, int>	nodeDegree;
		map<string, double> curNode;
		for (int i = 0; i < (int)minList.size(); i++) {
			for (int j = 0; j < (int)minList.size(); j++) {
				if ((*(param->genenet))[minList[i]].find(minList[j]) != (*(param->genenet))[minList[i]].end()) { // having a linkage between genes								
					curNode[minList[j]] = (*(param->genenet))[minList[i]][minList[j]];
					nodeDegree[minList[i]] ++;
				}
			}
			subnet[minList[i]] = curNode;
			curNode.clear();
		}

		//////////////////////////////////////////////////////////////////////////						
		const int minDegree = 1;
		vector<int> eid;
		int idx = 0;
		for (auto vi : minList) {
			if (nodeDegree[vi] < minDegree) {
				eid.push_back(idx);
			}
			idx++;
		}
		while (eid.size() > 0) {
			minList.erase(minList.begin() + eid.back());
			eid.pop_back();
		}

		////	SNP-GENE INFORMATION	//////////////////////////////////////////////////////////////////////		
		string strcutoff = to_string(pow(10, -cutoff));
		const string geneTag1 = "<div id='myModal";
		const string geneTag2 = "' class='modal'>  <div class='modal-content'><div class='modal-header'><span id='close";
		const string geneTag3 = "' class='close'>&times;</span><h2>Gene-SNPs Reference</h2></div><div class='modal-body'><table><col width='200'><col width='auto'><tbody><tr><td><p><b>Gene<b></p></td><td><a href = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene=";
		const string geneTag4 = "'><b>";
		const string geneTag5 = "</b></a><i>("; // +gene score
		const string geneTag6 = ")</i></td></tr><tr><td colspan='2'>&nbsp; </td></tr><tr><td><p><b>SNPs</b> <i>(p_value &le; " + strcutoff + ") </i> </p></td> <td><p>";

		const string snpTag1 = "<a href = 'https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=";
		const string snpTag2 = "'>";
		const string snpTag3 = "</a><i>(";
		const string snpTag4 = ")</i>&nbsp; ";

		const string subGeneTag1 = "</p></td></tr><tr><td colspan='2'>&nbsp; </td></tr><tr><td><p><b>Neighboring genes </b><i>(interaction weight)</i></p></td> <td><p>";
		const string subGeneTag2 = "<a href = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene="; // + gene +
		const string subGeneTag3 = "'><b>"; // + gene +
		const string subGeneTag4 = "</b></a><i>("; // + gene score +
		const string subGeneTag5 = ")</i>&nbsp; ";

		const string geneTagEnd = "</p></td></tr></tbody></table></div><div class='modal-footer'><h4>" + setName + "</h4></div></div></div></br>";

		string snpGeneInfo = "";
		string agene = "";

		for (auto& vi : minList)		// scan for all genes in minimal List
		{
			double genescore = pow(10, -((*(param->adjGeneScore))[(*(param->geneInputIndex))[vi]]));
			string strGeneScore = to_string(genescore);
			agene = geneTag1 + vi + geneTag2 + vi + geneTag3 + vi + geneTag4 + vi + geneTag5 + strGeneScore + geneTag6;
			map<double, string> snpval;
			double esp = 0.000000001;
			int k = 0;
			for (auto& svi : (*(param->geneSnpGlobal))[vi]) // for corresponding list of SNPs
			{
				if ((*(param->lgRsInputIndex)).find(svi) != ((*(param->lgRsInputIndex))).end()) /// some not welcome snp may appear!!!
				{
					if ((*(param->lgRsInput))[(*(param->lgRsInputIndex))[svi]] >= cutoff) {// significant SNPs
						if (snpval.find((*(param->lgRsInput))[(*(param->lgRsInputIndex))[svi]]) == snpval.end()) {// not exist
							snpval[(*(param->lgRsInput))[(*(param->lgRsInputIndex))[svi]]] = svi;
						}
						else {
							snpval[(*(param->lgRsInput))[(*(param->lgRsInputIndex))[svi]] + ++k * esp] = svi;
						}
					}
				}
			}

			for (map<double, string>::reverse_iterator sv = snpval.rbegin(); sv != snpval.rend(); sv++) {
				double tp = pow(10, -sv->first);
				string strTp = to_string(tp);
				agene += snpTag1 + sv->second + snpTag2 + sv->second + snpTag3 + strTp + snpTag4;
			}
			agene += subGeneTag1;
			//////////////////////////////////////////////////////////////////////////

			map<double, string> geneweight;
			k = 0;
			for (auto& svi : minList) { // find neighbor in current minimal list
				if (svi.compare(vi) != 0) // not identical
				{
					if (subnet[vi].find(svi) != subnet[vi].end()) // having link from vi to svi
					{
						if (geneweight.find(subnet[vi][svi]) != geneweight.end()) // not exist
							geneweight[subnet[vi][svi]] = svi;
						else {
							geneweight[subnet[vi][svi] + ++k * esp] = svi;
						}
					}
					else
						if (subnet[svi].find(vi) != subnet[svi].end()) // having link from svi to vi
						{
							if (geneweight.find(subnet[svi][vi]) != geneweight.end()) // not exist
								geneweight[subnet[svi][vi]] = svi;
							else {
								geneweight[subnet[svi][vi] + ++k * esp] = svi;
							}
						}
				}
			}

			for (map<double, string>::reverse_iterator gw = geneweight.rbegin(); gw != geneweight.rend(); gw++) {
				string strEdgeWeight = to_string(gw->first);
				agene += subGeneTag2 + gw->second + subGeneTag3 + gw->second + subGeneTag4 + strEdgeWeight + subGeneTag5;
			}
			//////////////////////////////////////////////////////////////////////////
			agene += geneTagEnd;
			snpGeneInfo += agene;
		}

		/////////////////////////////////////////////////////////
		//double curGeneScore = orgGeneScore[geneInputIndex[<geneName>]];	 // NODE WEIGHT
		//filePath = setSubNetBuilding("geneset_test", geneList, subnet, orgGeneScore, geneInputIndex);
		ogdf::Graph G;
		ogdf::GraphAttributes GA(G, ogdf::GraphAttributes::nodeGraphics |
			ogdf::GraphAttributes::nodeStyle |
			ogdf::GraphAttributes::nodeLabel |
			ogdf::GraphAttributes::edgeGraphics |
			ogdf::GraphAttributes::edgeArrow |
			ogdf::GraphAttributes::edgeStyle);
		map <string, ogdf::node> PROTEIN_MAP;

		graph_construct(minList, subnet, PROTEIN_MAP, G, GA);
		node_colorize(G, GA, *(param->adjGeneScore), *(param->geneInputIndex), PROTEIN_MAP);
		node_setsize(G, GA);
		edge_adjust(G, GA);

		ogdf::SugiyamaLayout SL;
		SL.call(GA);
		ogdf::SpringEmbedderFR SE;
		SE.call(GA);

		ogdf::GraphIO::SVGSettings svgs = ogdf::GraphIO::SVGSettings();
		svgs.fontSize(15);
		svgs.fontColor("#000000");
		filePath = setName + ".html";
		ofstream fou(filePath);
		ogdf::GraphIO::drawSVGe(GA, fou, svgs, snpGeneInfo, setName);
		if (fou.is_open())			fou.close();

		PROTEIN_MAP.clear();
		map <string, ogdf::node>().swap(PROTEIN_MAP);
		// finishing
		subnet.clear();
		map<string, map<string, double>>().swap(subnet);
		//////////////////////////////////////////////////////////////////////////
		return filePath;
	}

	//////////////////////////////////////////////////////////////////////////
	string commonNetworkLoading(string network_file = "data\\HIPPIE_NETWORK.txt") // Sora's update
	{
		map<string, map<string, double>> genenet2;
		map<string, map<string, double>> cmnet;

		ifstream fin(network_file);
		if (!fin) {
			cerr << "Cannot access specified file: " << network_file << endl;
		}

		string inpLine;
		while (getline(fin, inpLine))
		{
			istringstream ss(inpLine);
			map<string, double> tmp;
			string genea, geneb;
			double weight;
			ss >> genea >> geneb >> weight;
			weight = weight; // re-scale to 0 ~ 1
			if (geneInputIndex.find(genea) != geneInputIndex.end() && geneInputIndex.find(geneb) != geneInputIndex.end()) {	// both genes are under consideration
				if (genenet2.find(genea) == genenet2.end()) {	// not yet recognized before
					tmp[geneb] = weight; // add new edge of ogdf::Graph
					genenet2[genea] = tmp; // add new vertex of ogdf::Graph
					tmp.clear();	// reset tmp for next new edge
				}
				else {
					if (genenet2[genea].find(geneb) == genenet2[genea].end()) {// this is a new edge ~ linkage genea-geneb has not been recognized before
						genenet2[genea][geneb] = weight;
					}
					else {	// edge genea-geneb has been recognized before, update if there is a better weight (bigger)
						if (weight > genenet2[genea][geneb]) {
							genenet2[genea][geneb] = weight;
						}
					}
				}
			} /// end if both genes
		} /// end while		
		fin.close();

		//// matching
		double cutoff = 0.01;
		cutoff = -log10(cutoff);
		vector<string> minList(0);
		map<string, bool> tmpnode;

		for (auto si : genenet2) {
			if (geneRefinedIndex.find(si.first) != geneRefinedIndex.end() && // under consideration
				adjGeneScore[geneInputIndex[si.first]] >= cutoff)		// cutoff
				//orgGeneScore[geneInputIndex[si.first]] >= cutoff)		// cutoff
				for (auto vi : si.second)
					if (geneRefinedIndex.find(vi.first) != geneRefinedIndex.end())
						//if (orgGeneScore[geneInputIndex[vi.first]] >= cutoff)
						if (adjGeneScore[geneInputIndex[vi.first]] >= cutoff)
							tmpnode[vi.first] = 1;

			if (tmpnode.size() > 0)
				tmpnode[si.first] = 1; // this node not empty so add it;
		}

		for (auto ti : tmpnode)		minList.push_back(ti.first);

		map<string, int>	nodeDegree;
		map<string, double> curNode;
		for (int i = 0; i < (int)minList.size(); i++) {
			for (int j = 0; j < (int)minList.size(); j++) {
				if (genenet2[minList[i]].find(minList[j]) != genenet2[minList[i]].end())
					if (genenet.find(minList[i]) != genenet.end())
						if (genenet[minList[i]].find(minList[j]) != genenet[minList[i]].end()) { // having a linkage between genes in both STRING & HIPPIE
							curNode[minList[j]] = genenet[minList[i]][minList[j]];
							nodeDegree[minList[i]] ++;
						}
			}

			if (curNode.size() > 0)	cmnet[minList[i]] = curNode;
			curNode.clear();
		}

		ofstream fou("commonnet0.01.txt");
		int nEdge = 0;
		for (auto agene : cmnet) {
			for (auto bgene : agene.second) {
				fou << agene.first << "\t" << bgene.first << "\t" << bgene.second << endl;
				nEdge++;
			}
		}
		fou << cmnet.size() << endl;
		fou << nEdge << endl;
		fou.close();
		////////////////////////////////////////////////////////////////////////////

		////string coreNetDraw(netparam *param, string setName, double qvalcutoff, double cutoff)		
		string filePath = "common_network";

		//////	SNP-GENE INFORMATION	//////////////////////////////////////////////////////////////////////
		//const string geneTag1 = "<div id='myModal";
		//const string geneTag2 = "' class='modal'>  <div class='modal-content'><div class='modal-header'><span id='close";
		//const string geneTag3 = "' class='close'>&times;</span><h2>Gene-SNPs Reference</h2></div><div class='modal-body'><table><col width='200'><col width='auto'><tbody><tr><td><p><b>Gene<b></p></td><td><a href = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene=";
		//const string geneTag4 = "'><b>";
		//const string geneTag5 = "</b></a><i>("; // +gene score
		//const string geneTag6 = ")</i></td></tr><tr><td colspan='2'>&nbsp; </td></tr><tr><td><p><b>SNPs</b> <i>(p_value &le; 0.05) </i> </p></td> <td><p>";

		//const string snpTag1 = "<a href = 'https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=";
		//const string snpTag2 = "'>";
		//const string snpTag3 = "</a><i>(";
		//const string snpTag4 = ")</i>&nbsp; ";

		//const string subGeneTag1 = "</p></td></tr><tr><td colspan='2'>&nbsp; </td></tr><tr><td><p><b>Neighboring genes </b><i>(interaction weight)</i></p></td> <td><p>";
		//const string subGeneTag2 = "<a href = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene="; // + gene +
		//const string subGeneTag3 = "'><b>"; // + gene +
		//const string subGeneTag4 = "</b></a><i>("; // + gene score +
		//const string subGeneTag5 = ")</i>&nbsp; ";

		//const string geneTagEnd = "</p></td></tr></tbody></table></div><div class='modal-footer'><h4> COMMON NETWORK </h4></div></div></div></br>";

		//string snpGeneInfo = "";
		//string agene = "";

		//for (auto& vi : minList)		// scan for all genes in minimal List
		//{
		//	double genescore = pow(10, -(((orgGeneScore))[((geneInputIndex))[vi]]));
		//	string strGeneScore = to_string(genescore);
		//	agene = geneTag1 + vi + geneTag2 + vi + geneTag3 + vi + geneTag4 + vi + geneTag5 + strGeneScore + geneTag6;

		//	map<double, string> snpval;
		//	double esp = 0.000000001;
		//	int k = 0;
		//	for (auto& svi : ((geneSnpGlobal))[vi]) // for corresponding list of SNPs
		//	{
		//		if (((lgRsInputIndex)).find(svi) != (((lgRsInputIndex))).end()) /// some not welcome snp may appear!!!
		//		{
		//			if (((lgRsInput))[((lgRsInputIndex))[svi]] >= cutoff) {// significant SNPs
		//				if (snpval.find(((lgRsInput))[((lgRsInputIndex))[svi]]) == snpval.end()) {// not exist
		//					snpval[((lgRsInput))[((lgRsInputIndex))[svi]]] = svi;
		//				}
		//				else {
		//					snpval[((lgRsInput))[((lgRsInputIndex))[svi]] + ++k * esp] = svi;
		//				}
		//			}
		//		}
		//	}

		//	for (map<double, string>::reverse_iterator sv = snpval.rbegin(); sv != snpval.rend(); sv++) {
		//		double tp = pow(10, -sv->first);
		//		string strTp = to_string(tp);
		//		agene += snpTag1 + sv->second + snpTag2 + sv->second + snpTag3 + strTp + snpTag4;
		//	}
		//	agene += subGeneTag1;

		//	//////////////////////////////////////////////////////////////////////////
		//	map<double, string> geneweight;
		//	k = 0;

		//	for (auto& svi : minList) { // find neighbor in current minimal list
		//		if (svi.compare(vi) != 0) // not identical
		//		{
		//			if (cmnet[vi].find(svi) != cmnet[vi].end()) // having link from vi to svi
		//			{
		//				if (geneweight.find(cmnet[vi][svi]) != geneweight.end()) // not exist
		//					geneweight[cmnet[vi][svi]] = svi;
		//				else {
		//					geneweight[cmnet[vi][svi] + ++k * esp] = svi;
		//				}
		//			}
		//			else
		//				if (cmnet[svi].find(vi) != cmnet[svi].end()) // having link from vi to svi
		//				{
		//					if (geneweight.find(cmnet[svi][vi]) != geneweight.end()) // not exist
		//						geneweight[cmnet[svi][vi]] = svi;
		//					else {
		//						geneweight[cmnet[svi][vi] + ++k * esp] = svi;
		//					}
		//				}
		//		}
		//	}

		//	for (map<double, string>::reverse_iterator gw = geneweight.rbegin(); gw != geneweight.rend(); gw++) {
		//		string strEdgeWeight = to_string(gw->first);
		//		agene += subGeneTag2 + gw->second + subGeneTag3 + gw->second + subGeneTag4 + strEdgeWeight + subGeneTag5;
		//	}
		//	//////////////////////////////////////////////////////////////////////////
		//	agene += geneTagEnd;
		//	snpGeneInfo += agene;
		//}

		///////////////////////////////////////////////////////////
		////double curGeneScore = orgGeneScore[geneInputIndex[<geneName>]];	 // NODE WEIGHT
		////filePath = setSubNetBuilding("geneset_test", geneList, subnet, orgGeneScore, geneInputIndex);
		//ogdf::Graph G;
		//ogdf::GraphAttributes GA(G, ogdf::GraphAttributes::nodeGraphics |
		//	ogdf::GraphAttributes::nodeStyle |
		//	ogdf::GraphAttributes::nodeLabel |
		//	ogdf::GraphAttributes::edgeGraphics |
		//	ogdf::GraphAttributes::edgeArrow |
		//	ogdf::GraphAttributes::edgeStyle);
		//map <string, ogdf::node> PROTEIN_MAP;

		//graph_construct(minList, cmnet, PROTEIN_MAP, G, GA);
		//node_colorize(G, GA, (orgGeneScore), (geneInputIndex), PROTEIN_MAP);
		//node_setsize(G, GA);
		//edge_adjust(G, GA);
		//ogdf::SugiyamaLayout SL;
		//SL.call(GA);
		//ogdf::SpringEmbedderFR SE;
		//SE.call(GA);

		//ogdf::GraphIO::SVGSettings svgs = ogdf::GraphIO::SVGSettings();
		//svgs.fontSize(15);
		//svgs.fontColor("#000000");
		//filePath += ".html";
		//ofstream fou(filePath);
		//ogdf::GraphIO::drawSVGe(GA, fou, svgs, snpGeneInfo, filePath);
		//if (fou.is_open())	fou.close();

		//PROTEIN_MAP.clear();
		//map <string, ogdf::node>().swap(PROTEIN_MAP);
		// finishing					
		//////////////////////////////////////////////////////////////////////////
		return filePath;
	}

	//static void addedge(ogdf::node node1, ogdf::node node2, double value, ogdf::Graph &G, ogdf::GraphAttributes &GA);
	//static void addNode(string label, ogdf::Graph &G, ogdf::GraphAttributes &GA, std::map<string, ogdf::node> &PROTEIN_MAP);
	//static string ith(int col);
	//static string colorize(double weight);
	//static void node_colorize(ogdf::Graph G, ogdf::GraphAttributes &GA, std::vector<double> orgGeneScore, std::map <string, int> geneInputIndex, std::map <string, ogdf::node> PROTEIN_MAP);
	//static double percent(int degree, std::vector<int> degrees);
	//static double percent(double weight, std::vector<double> weights);
	//static std::vector <float> sort_weight(ogdf::Graph G, ogdf::GraphAttributes GA);
	//static void graph_construct(std::vector <string> SubNetGenes, std::map<string, std::map<string, double>> subnet, std::map <string, ogdf::node> &PROTEIN_MAP, ogdf::Graph &G, ogdf::GraphAttributes &GA);
	//static std::vector <int> sort_degree(ogdf::Graph G);
	//static void node_setsize(ogdf::Graph G, ogdf::GraphAttributes &GA);
	//static void edge_adjust(ogdf::Graph G, ogdf::GraphAttributes &GA);
	//static string setSubNetDraw(netparam *param, string setName, double cutoff = 0.05);
	//static string coreNetDraw(netparam *param, string setName = "DIABETES_CORE_GENE_NET", double qvalcutoff = 0.10, double cutoff = 0.001);

	void addedge(ogdf::node node1, ogdf::node node2, double value, ogdf::Graph &G, ogdf::GraphAttributes &GA) {
		ogdf::edge e = G.newEdge(node1, node2);
		GA.arrowType(e) = ogdf::EdgeArrow::eaNone;
		GA.strokeWidth(e) = (float)value;
		GA.strokeColor(e) = ogdf::Color("#21577E");
		// Middleblue "#2D69D2"
		// Deepblue "#21577E"
	}

	string ith(int col) { // integer to hex
		stringstream ss;
		ss << hex << col;
		string result(ss.str());
		if (col < 10) { result = "0" + result; }
		return(result);
	}

	string colorize(double weight) {
		string Red, Green, Blue;
		weight = weight > 0 ? weight : 0;
		if (weight > 5)
		{
			Red = ith(0); Green = ith(138); Blue = ith(69);
		}
		else
			if (weight > 1.3)
			{
				Red = ith((int)round(0 + (77 - 0)*(5 - weight) / (5)));
				Green = ith((int)round(138 + (184 - 138)*(5 - weight) / (5)));
				Blue = ith((int)round(69 + (133 - 69)*(5 - weight) / (5)));
			}
			else
			{
				Red = ith((int)round(0 + (194 - 0)*(5 - weight) / 5));
				Green = ith((int)round(138 + (255 - 138)*(5 - weight) / 5));
				Blue = ith((int)round(69 + (229 - 69)*(5 - weight) / 5));
			}

		string res = "#" + Red + Green + Blue;
		return(res);
	}

	void addNode(string label, ogdf::Graph &G, ogdf::GraphAttributes &GA, map<string, ogdf::node> &PROTEIN_MAP) {
		ogdf::node v = G.newNode();
		PROTEIN_MAP.insert({ label, v });
		GA.label(v) = label;
		GA.shape(v) = ogdf::Shape::shEllipse;
		GA.width(v) = 54;
		GA.height(v) = 54;
		GA.strokeWidth(v) = 2;
	}

	void node_colorize(ogdf::Graph G, ogdf::GraphAttributes &GA, vector<double> orgGeneScore, map <string, int> geneInputIndex, map <string, ogdf::node> PROTEIN_MAP) {
		ogdf::node v;
		cout << "Node_Coloring" << endl;
		forall_nodes(v, G) { // gene weight coloring, -log10(0.05) = 1.3
							 //GA.fillColor(v) = Color("#4C89AE"); // set default_color as 				
			GA.fillColor(v) = ogdf::Color("#C2FFE5"); // set default_color as 						
			GA.fillColor(PROTEIN_MAP[GA.label(v)]) = ogdf::Color(colorize(orgGeneScore[geneInputIndex[GA.label(v)]]));
		}
	}

	double percent(int degree, vector<int> degrees) {
		int pos = find(degrees.begin(), degrees.end(), degree) - degrees.begin();
		return(1 - double(pos + 1) / degrees.size());
	}

	float percent(float weight, vector<float> weights) {
		int pos = find(weights.begin(), weights.end(), weight) - weights.begin();
		return(1 - (float)(pos + 1) / weights.size());
	}

	void graph_construct(vector <string> geneList, map<string, map<string, double>> subnet, map <string, ogdf::node> &PROTEIN_MAP, ogdf::Graph &G, ogdf::GraphAttributes &GA) {
		for (auto iter = geneList.begin(); iter != geneList.end(); ++iter) {
			if (PROTEIN_MAP.count(*iter) == 0) { addNode(*iter, G, GA, PROTEIN_MAP); }
			for (auto it = subnet[*iter].begin(); it != subnet[*iter].end(); ++it) {
				if (PROTEIN_MAP.count(it->first) == 0) { addNode(it->first, G, GA, PROTEIN_MAP); }
				if (it->second != 0) { addedge(PROTEIN_MAP[*iter], PROTEIN_MAP[it->first], it->second, G, GA); }
			}
		}
	}

	vector <int> sort_degree(ogdf::Graph G) {
		vector<int> degrees;
		ogdf::node v;
		forall_nodes(v, G) { degrees.push_back(v->degree()); }
		sort(degrees.begin(), degrees.end());
		degrees.erase(unique(degrees.begin(), degrees.end()), degrees.end());
		return degrees;
	}

	void node_setsize(ogdf::Graph G, ogdf::GraphAttributes &GA) {
		cout << "Node_sizing" << endl;
		vector<int> degrees = sort_degree(G);
		ogdf::node v;
		forall_nodes(v, G) { // default size : 54
			GA.width(v) += (1.0 - percent(v->degree(), degrees)) * 54;
			GA.height(v) += (1.0 - percent(v->degree(), degrees)) * 54;
		}
	}

	vector<float> sort_weight(ogdf::Graph G, ogdf::GraphAttributes GA) {
		vector<float> weights;
		ogdf::edge e;
		forall_edges(e, G) { weights.push_back(GA.strokeWidth(e)); }
		sort(weights.begin(), weights.end());
		weights.erase(unique(weights.begin(), weights.end()), weights.end());
		return weights;
	}

	void edge_adjust(ogdf::Graph G, ogdf::GraphAttributes &GA) {
		cout << "Edge_adjusting" << endl;
		vector <float> weights = sort_weight(G, GA);
		ogdf::edge e;
		forall_edges(e, G) {
			GA.strokeWidth(e) = ((float)1.0 - percent(GA.strokeWidth(e), weights)) * 5 + (float)0.4;
			if (GA.strokeWidth(e) > 2.1) { GA.strokeColor(e) = ogdf::Color("#000000"); } // emphasize edge with Black		
		}
	}

	void adjustedPvalueBH()
		/*- Step-up Hochberg:
		p#(i) = n/i *p(i)
		where (i <= n - 1),
		n = number of p - values
		p(i) = ith smallest p - value
		p#(i) = adjusted value of p(i)
		*/
	{
		auto_ptr<map<double, string>> sort(new map<double, string>);
		//sorting min to max
		for (auto& mi : pvAdjusted)
		{
			(*sort)[mi.second] += mi.first + "\t";
		}

		int countRank = 1;
		int maxRank = pvAdjusted.size();
		double lastFdr = (double)-MAXINT32;

		for (auto& mi2 : (*sort))
		{
			string iss = mi2.second;
			string  set;
			int count;
			double tmp = 0.0;

			count = 0;
			while (iss.find("\t") != iss.npos)
			{
				set = iss.substr(0, iss.find_first_of("\t")); // position of "\t" == count(0~ previous "\t")				

				fdrAdjusted[set] = pvAdjusted[set] * ((double)maxRank / countRank);

				if (fdrAdjusted[set] > 1.0) { fdrAdjusted[set] = 1.0; }
				if (fdrAdjusted[set] < lastFdr) { fdrAdjusted[set] = lastFdr; }
				lastFdr = fdrAdjusted[set];

				iss.replace(0, iss.find_first_of("\t") + 1, "");
				count++;
			}

			if (count > 1)		countRank += count;
			else				countRank++;
		}
	}

	double z2p(double z)
	{
		/**
		* Computes cumulative density function of standard Gaussian using Taylor approximation
		* @param z z-score
		* @return (1 - standard Gaussian cdf(z) )
		*/

		if (isnan(z))
		{
			return 1.0;
		}

		if (z < -8.0)
		{
			return 1.0;
		}

		if (z > 8.0)
		{
			return 0.0;
		}

		double sum = 0.0;
		double term = z;

		for (int i = 3; sum + term != sum; i += 2)
		{
			sum += term;
			term = (term * z * z) / i;
		}

		return (0.5 - sum * exp(-z * z / 2) / sqrt(2 * M_PI));
	}

	double p2z(double p) {
		normal s;
		return -quantile(complement(s, p));
	}

	double meanv(vector<double> xSample)
	{
		double sum = 0.0;
		for (auto& vi : xSample)	sum += vi;
		return sum / xSample.size();
	}

	double sdv(vector<double> xSample)
	{
		double szSample = (double)xSample.size();
		if (szSample == 0)		return 0.0;
		if (szSample == 1)		return 0.0;

		double mean = meanv(xSample);
		double sq_sum = 0.0;

		for (auto& vi : xSample) {
			sq_sum += (vi - mean) * (vi - mean);
		}

		return  sqrt(sq_sum / (szSample - 1));
	}

	double zsScoring(vector<double> xSub, double meanX, double stdevX)
	{
		double szSub = xSub.size();
		if (szSub == 0 || stdevX == 0)
			return -7.9;//szSub += 0.00000001; // no division by 0 ~ smoothing
		return (meanv(xSub) - meanX) / (stdevX / sqrt(szSub));
	}

	double zsScoring(vector<double> xSub, double meanX, double stdevX, int G)
	{
		double szSub = xSub.size();
		if (szSub == 0 || stdevX == 0 || G == 1 || szSub == 1)
			return -8.0;//szSub += 0.00000001; // no division by 0 ~ smoothing		
		return (meanv(xSub) - meanX) / (sqrt((G - szSub) / (G - 1)) * stdevX / sqrt(szSub));
	}
};
