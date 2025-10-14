#ifndef SIM
#define SIM
#include <iostream>
#include <vector>
#include <utility>   
#include <string>
#include <unordered_map>
#include <random>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include "Tile.h"

using std::string;
using std::vector;
using std::pair;
using std::unordered_map;
using std::ofstream;

class Sim
{
public:
	// R1, R2, R3 is number of reaction pairs for each reaction
	Sim() : gen(rd()), dis(0, 1.0), _timeTrack{ 0.0 } {}
	~Sim() {}

	void nextReacTime(const double& totalProp);
	int getRandIndex(int size);  // Get index for reaction
	int determineReaction(const double& totalProp, const vector<double>& reacProp);
	void performReaction(const int& rxn, Tile& tile);

	void simStep(Tile& tile);
	//void printConc(ofstream& dataFile);
	void printConcToFile(const string& filename);
	void runSim(double maxTime);
	void createPlot();

private:
	vector<vector<pair<string,int>>> _concOverTime;  // Keeps track of concentrations over time
	vector<int> _rxnTrack;  // Keep track of all reactions
	vector<double> _timeTrack;  // Keeps track of time for each reaction
	std::random_device rd;  // Obtain seed 
	std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis;  // (0, 1.0)
	double getRand() { return dis(gen); }  // Computes random number
};

#endif