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
#include <map>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <boost/math/distributions/normal.hpp> // for normal_distribution
#include <boost/math/distributions.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <set>
#include <functional>
#include "regression.h"

using namespace std;
using boost::math::normal; // typedef provides default type is double.
using std::setw;
using std::setprecision;
using std::numeric_limits;

#define kGene	2
#define SIGNIFICANT_THRESHOLD_OF_PVALUE 0.0001
#define M_PI           3.14159265358979323846  /* pi */
#define MAXINT32	2147483647
#define SNP_INPUT				0
#define GENE_INPUT				1

typedef struct GeneLocation {
	string geneName;
	int	chromosome;
	int start;
	int end;
}	genelocation;

typedef	struct ItsAGene {
	string	geneName;		// gene name	~ gene1Name		
	string	snpL, snpM, snpR;
}	itsagene;

typedef struct AdjacentGene {
	itsagene	geneA;
	itsagene	geneB;
	double		corr[9] = { 0. }; //cor_LeLe, cor_LeMi, cor_LeRi, cor_MiLe, cor_MiMi, cor_MiRi, cor_RiLe, cor_RiMi, cor_RiRi;
}	adjacentpair;

struct MYPARAM
{
	vector<double> geneNullScore;
	//////////////////////////////////
	vector<double> subNewScore;
	vector<double> orgGeneProfile;
	map<string, double> nsSetScore;
	//////////////////////////
	map<string, int> subHistogram;
	vector<vector<double>> subGeneProfile;
	vector<vector<double>> subSetProfile;

	int nTrial;
	int iThread;
	double pCutoff = SIGNIFICANT_THRESHOLD_OF_PVALUE;
	int snpSize;
	int geneSize;
	int* pivot;
	vector<int>* snpList;
	int* snpN;
	double*	snpValue;
	string* geneName;
	vector<string> setGlobal;
	map <string, int> obSetList;
	vector<int>* geneList;
	int* geneN;
	map<int, int> newIndex;
};

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
	if (szSub == 0 || stdevX == 0 || G == 1)
		return -8.0;//szSub += 0.00000001; // no division by 0 ~ smoothing
	return (meanv(xSub) - meanX) / (sqrt((G - szSub) / (G - 1)) * stdevX / sqrt(szSub));
}

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

	~gsa()
	{
	};

public:

	/////////////////////////////////////////////////////////////////////////////////////////////
	void adjustedPvalueBH()
		/*- Step-up Hochberg:
			p#(i) = n/i *p(i)
			where (i <= n - 1),
			n = number of p - values
			p(i) = ith smallest p - value
			p#(i) = adjusted value of p(i)
		*/
	{
		//////////////////////////////////////////////////////////////////////////
		map<string, double>	ogNorm = pvAdjusted;

		map<double, string> sort;

		//sorting min to max
		for (auto& mi : ogNorm)
		{
			sort[mi.second] += mi.first + "\t";
		}

		int countRank = 1;
		int maxRank = ogNorm.size();
		double lastFdr = (double)-MAXINT32;

		for (auto& mi2 : sort)
		{
			string iss = mi2.second;
			string  set;
			int count;
			
			count = 0;
			while (iss.find("\t") != iss.npos)
			{
				set = iss.substr(0, iss.find_first_of("\t")); // position of "\t" == count(0~ previous "\t")				

				fdrAdjusted[set] = ogNorm[set] * ((double)maxRank / countRank);

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
			if (!firstpin || countpin == (int)valHistogram.size()) {
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

				/*while (eid.size() > 0) {
					mi.second.erase(mi.second.begin() + eid.back());
					eid.pop_back();
				}*/

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
		//		upperGS = findMedianOfTopMax(.5);		
		//fou << lowerGS << endl << upperGS << endl;

		vector<double> topGS;
		vector<double> topSize;
		vector<double> cutoffGS;
		vector<double> cutoffSize;

		bool solved = false;
		vector<double> lastFive;
		double avgGS; // desired cutoff gene score
		int iteration = 0;
		vector<double> coef(2,0);
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

		int maxSize = 0;
		double maxGS = 0.;

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

		//// 2 - monotonic cubic spline		
		vector<int> segment(0);
		map<int, int> sizeHistogram;
		for (auto vi : cutoffSize)
			sizeHistogram[(int)vi]++;

		//double percentile = 0.05;
		//double ks = percentile;
		segment.push_back(0); // for knots 0
		countpin = 1;
		int knot = 0;
		int countk = 1;
		//int fiboA = 0, fiboB = 1;
		for (auto mi : sizeHistogram) {			
			// 2^k - 1 ~ 2(k-1)+1
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
		//	if (countpin % 15 == 1 || countpin == numSize) segment.push_back(mi.first);
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
		int tmpsize = 0;//, yOneCount = 0;

		/*map<int, map<double, int>> size2pin;
		for (int i = 0; i < (int)cutoffSize.size(); i++) {
			tmpsize = (int)cutoffSize[i];
			for (int k = 1; k < nSeg; k++) {
				if (segment[k] - range < tmpsize && tmpsize <= segment[k]) {
					size2pin[k][cutoffGS[i]] = 1;
				}
			}
		}
		for (int i = 1; i < nSeg; i++) {
			knots[i].y = findMedianOfTopMax(size2pin[i], 1.0);
		}*/

		for (int i = 0; i < (int)cutoffSize.size(); i++) {
			tmpsize = (int)cutoffSize[i];
			for (int k = 1; k < nSeg; k++) {
				if (segment[k] - range/2 < tmpsize && tmpsize <= segment[k] + range/2) {
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

			newGeneScore[i] -= trendline;			
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
	void geneFilter(double inputCorrelation = 0.5)
	{
		//int countrefinedSet = 0;
		//int countrefinedGene = 0;
		// find representative gene if exist
		map<string, bool> geneRefinedtmp;
		for (auto gi : geneInput) {
			if (famiGene.find(gi) != famiGene.end()) { // there is a representative for this gene
													   //geneRefinedtmp[famiGene[gi]] = 1;
				geneRefinedtmp[gi] = 1;
			}
			else
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
			//fou << "refinedgeneinorder" << endl;
			map<int, string> refinedgeneinorder; // <gene-loc, genelocation> sorted ascending
			int aliasLoc = 1;
			for (auto gi : sg.second) {// gene list
									   //if (geneRefinedIndex.find(gi) == geneRefinedIndex.end()) // has been replaced by representative gene
									   //if (geneInputIndex.find(gi) == geneInputIndex.end()) // has been replaced by representative gene
									   //{
									   //	/*bool found = false;
									   //	for (auto ti : sg.second) {
									   //		if (ti.compare(famiGene[gi]) == 0) found = true;
									   //		break;
									   //	}
									   //	if (found)	refinedgeneinorder[geneLoc[famiGene[gi]].start] = famiGene[gi];
									   //	else */
									   //	refinedgeneinorder[geneLoc[gi].start] = gi;
									   //	//fou << geneLoc[famiGene[gi]].start << "\t" << famiGene[gi] << endl;
									   //}
									   //else 
				if (refinedgeneinorder.find(geneLoc[gi].start) == refinedgeneinorder.end()) // same start
					refinedgeneinorder[geneLoc[gi].start] = gi;
				else
				{
					//fou << "start duplicated! " << refinedgeneinorder[geneLoc[gi].start] << "\t" << gi << endl;
					string previousGene = refinedgeneinorder[geneLoc[gi].start]; // gene which has the same start
					int previousLength = geneLoc[previousGene].end - geneLoc[previousGene].start;
					int giLength = geneLoc[gi].end - geneLoc[gi].start;

					if (giLength > previousLength) { // the shorter the later
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

			//fou << "sortedRefinedGene" << endl;
			vector<string> sortedRefinedGene;
			for (auto li : refinedgeneinorder) {
				sortedRefinedGene.push_back(li.second);
				//fou << li->second << endl;
			}

			// 3rd: checking correlation of 2 adjacent genes
			//fou << "finalRefinedGene\n";
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
						//fou << "snpA: " << snpA << "\t" << "snpB: " << snpB << endl;

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
							//	fou << minidxA << "\t" << minidxB << "\t" << corrAB << endl;
							if (corrAB >= inputCorrelation) { // high correlation detected, remove B							
																
								i++; // step over geneB, move to next gene
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
				}
				else
					if (i == (int)sortedRefinedGene.size() - 2)		   // if current gene is the last gene, continue/ bypass;
						finalRefinedGene.push_back(sortedRefinedGene[i + 1]);			// always keep current gene				
			}	/// end of for (int i = 0; i < (int)sortedRefinedGene.size() - 1; i++)

				// 4th: update setGeneRefined
			setGeneRefined[sg.first] = finalRefinedGene;			
		}	/// end of for (auto sg : setGeneGlobal)			
			

	}	/// end of function geneFilter

	/*void familyGeneLoad(string filename = "data\\FamilyGeneList")
	{
		ifstream fin(filename);
		if (!fin.is_open())
		{
			cerr << "Unavailable data: " << filename << " is missed!" << endl;
			
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
	}*/

	void geneLocMapGenerate(string filename = "data\\hg19GeneList") {
		ifstream fin(filename);
		if (!fin.is_open())
		{
			cerr << "Not all data are available: " << filename << " is missed!" << endl;
			
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

	void snpLocMapGenerate(string filename = "data\\rsloc_hg19") {
		ifstream fin;
		string inpLine, curSnp;
		fin.open(filename);

		if (!fin.is_open())
		{
			cerr << "Unavailable data: " << filename << " is missed!" << endl;
			
		}
		////////////////////////////////////////////////////////			
		// remove 1 meta lines
		for (int j = 0; j < 1; ++j)
		{
			getline(fin, inpLine);
		}

		while (getline(fin, inpLine))
		{
			string chromosome, curSnp;
			int position;
			stringstream iss(inpLine);
			iss >> curSnp >> position;

			snpLocMap[curSnp] = position;
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
		vector<double> orgGeneScore = newGeneScore;///geneProfile[0];

		double muPop = meanv(orgGeneScore);
		double sigPop = sdv(orgGeneScore);
		map<double, bool> upl;
		for (auto& mi : setGeneGlobal)
		{
			vector<double> genePvalue;
			for (auto& vi : mi.second)
			{
				double curGeneScore = orgGeneScore[geneInputIndex[vi]];
				if (geneInputIndex.find(vi) != geneInputIndex.end()) // only genes included in sets ~ USEFUL GENES
				//if (upl.find(curGeneScore) == upl.end()) 
					//if (curGeneScore != 0)
				{
					genePvalue.push_back(curGeneScore);
					//	upl[curGeneScore] = 1;
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
		string id2gene = "data\\id2gene";  // convert other ID type to symbol type

		map<string, string> ensemblID;	//  <ensemblID, symbol>
		map<string, string> entrezID;		//  <entrezID, symbol>
		map<string, string> uniswissprotID;// <uniswissprotID, symbol>

		ifstream fin(setFile);
		if (!fin.is_open())
		{
			cerr << "Cannot access specified file: " << setFile << endl;
			
		}

		if (convert2symbol) { // load conversion table
			ifstream ftable(id2gene);
			if (!fin.is_open())
			{
				cerr << "Cannot access conversion table: " << id2gene << endl;				
			}

			string aline;
			getline(ftable, aline);		// header: ensembl_gene_id >>	entrezgene	>>	external_gene_name	>>	uniprot_swissprot
			while (getline(ftable, aline)) {
				string ensembl, entrez, genesymbol, uniswiss;
				stringstream tmp(aline);
				tmp >> ensembl >> entrez >> genesymbol >> uniswiss;

				ensemblID[ensembl] = genesymbol;

				if (entrez.compare("NA") != 0) { // available 
					entrezID[entrez] = genesymbol;
				}

				if (uniswiss.length() > 0) { // available 
					uniswissprotID[uniswiss] = genesymbol;
				}
			}
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
					if (ensemblID.find(tmps) != ensemblID.end()) { // ensemblID
						s = ensemblID[tmps];
					}

					if (entrezID.find(tmps) != entrezID.end()) { // entrezID
						s = entrezID[tmps];
					}

					if (uniswissprotID.find(tmps) != uniswissprotID.end()) { // uniswissprotID
						s = uniswissprotID[tmps];
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

			// set trunking minGene, maxGene
			if ((minSetSize <= (int)term.size()) && ((int)term.size() <= maxSetSize))
			{
				setGlobal.push_back(tmpSetName);
				setGeneGlobal[tmpSetName] = term; // geneGlobal list of setGlobal[i]
			}

			tmpSetName.clear();
			realSetSize = 0;
			term.clear();			// clear current terms			
		}

		fin.close();
	}

	void snpGeneMapGenerate(string snp_gene_file = "data\\db19_20k", bool inputType = SNP_INPUT)
	{
		ifstream fin(snp_gene_file);
		if (!fin.is_open())
		{
			cerr << "Cannot access specified file: " << snp_gene_file << endl;
			
		}
		map<string, vector<string>> snpGeneGlobal;
		
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
					(snpGeneGlobal)[rs] = term;
					term.clear();				// clear current terms			
				}
			}
			fin.close();

			/// with each gene, scan for available SNP, which means available p-val/ -log(p-val)
			/// if there's no SNP available, remove gene	
			// all genes include at least ONE snp because we only read in available snp gene lists
			for (auto& mi : snpGeneGlobal)
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
	
		}

		string curGene, tmpG;
		string pval;
		double tmpP;

		map<string, string> ensemblID;	//  <ensemblID, symbol>
		map<string, string> entrezID;		//  <entrezID, symbol>
		map<string, string> uniswissprotID;	// <uniswissprotID, symbol>
		string id2gene = "data\\id2gene";  // convert other ID type to symbol type
		ifstream ftable(id2gene);
		if (!ftable.is_open())
		{
			cerr << "Cannot access conversion table: " << id2gene << endl;
			
		}

		string aline;
		getline(ftable, aline);		// header: ensembl_gene_id >>	entrezgene	>>	external_gene_name	>>	uniprot_swissprot
		while (getline(ftable, aline)) {
			string ensembl, entrez, genesymbol, uniswiss;
			stringstream tmp(aline);
			tmp >> ensembl >> entrez >> genesymbol >> uniswiss;

			(ensemblID)[ensembl] = genesymbol;

			if (entrez.compare("NA") != 0) { // available 
				(entrezID)[entrez] = genesymbol;
			}

			if (uniswiss.length() > 0) { // available 
				(uniswissprotID)[uniswiss] = genesymbol;
			}
		}
		ftable.close();

		//ofstream tfo("DIAGRAMgene");
		while (fin >> tmpG) {
			fin >> tmpP;
			if (tmpP > 0) {
				if ((ensemblID).find(tmpG) != (ensemblID).end()) { // ensemblID
					curGene = (ensemblID)[tmpG];
				}
				else
					if ((entrezID).find(tmpG) != (entrezID).end()) { // entrezID
						curGene = (entrezID)[tmpG];
					}
					else
						if ((uniswissprotID).find(tmpG) != (uniswissprotID).end()) { // uniswissprotID
							curGene = (uniswissprotID)[tmpG];
						}
						else
							curGene = tmpG; // symbol

				geneInput.push_back(curGene);
				geneInputIndex[curGene] = geneInput.size() - 1;
				geneInputScore[curGene] = -log10(tmpP);
				
			}
		}
		fin.close();
	}


	void snpPvalueLoad(string snp_pval_file = "data\\DIAGRAM")//example_height_snp")
	{
		ifstream fin(snp_pval_file);
		if (!fin.is_open())
		{
			cerr << "Cannot access specified file: " << snp_pval_file << endl;
			
		}

		string curSnp;
		string pval, dumbdata;
		double tmpP;
		//double myesp = 0.00000000000000001;
		//int countesp = 1;

		while (fin >> curSnp) {
			//	identify p-value
			fin >> pval;
			tmpP = atof(pval.c_str()); 
			if (tmpP <= 0) {
				cerr << "Invalid p-value!!!" << endl;
				continue;
			}

			
			tmpP = -log10(tmpP);
			lgRsInput.push_back(tmpP);		// list of -log(p-value)
			lgRsInputIndex[curSnp] = lgRsInput.size() - 1; // tmp indexing
			sorted_pval_index[tmpP] = lgRsInput.size() - 1;			
		}
		fin.close();
	}
};
