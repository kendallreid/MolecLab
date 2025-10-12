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

class Sim
{
public:
	// R1, R2, R3 is number of reaction pairs for each reaction
	Sim(int A, int B, int U) : _conc{ A, B, U }, gen(rd()), dis(0, 1.0), _timeTrack{ 0.0 },
		_rxnTrack(), _concOverTime{ {A}, {B}, {U} } {}
	~Sim() {}

	void nextReacTime(const double& totalProp);
	int determineReaction(const double& totalProp, const vector<double>& reacProp);
	void performReaction(const int& rxn, Tile& tile);
	void initConc(Tile& tile);
	void updateConcentrations(int reaction);
	int getRandIndex(int size);  // Get index for reaction

	void simStep(Tile& tile);
	void runSim(double maxTime);

private:
	unordered_map<string, int> _conc;  // Reactant key holds # of pixels as concentration
	vector<vector<int>> _concOverTime;  // Keeps track of concentrations over time
	vector<int> _rxnTrack;  // Keep track of all reactions
	vector<double> _timeTrack;  // Time vector
	std::random_device rd;  // Obtain seed 
	std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis;  // (0, 1.0)
	double getRand() { return dis(gen); }  // Computes random number
};

#endif