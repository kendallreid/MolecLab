#ifndef SIM
#define SIM
#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <random>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include "Tile.h"

using std::string;
using std::vector;
using std::pair;

// REACTIONS //
// R1: A+B=2U
// R2: A+U=2A
// R3: B+U=2B

class Sim
{
public:
	// R1, R2, R3 is number of reaction pairs for each reaction
	Sim(int A, int B, int U) : _conc{ A, B, U }, gen(rd()), dis(0, 1.0), _time{ 0.0 },
		_rxn(), _concOverTime{ {A}, {B}, {U} } {}
	~Sim() {}

	void nextReacTime(const double& totalProp);
	int determineReaction(const double& totalProp, const vector<double>& reacProp);
	void performReaction(const int& rxn, Tile& tile);
	void updateConcentrations(int reaction);
	double getRand() { return dis(gen); }  // Computes random number
	int getRandIndex(int size);  // Get index for reaction

	void simStep(Tile& tile);
	void runSim(double maxTime);

private:
	vector<int> _conc;  // Concentrations vector for A, B, C
	vector<vector<int>> _concOverTime;  // Keeps track of concentrations over time
	vector<int> _rxn;  // Keep track of all reactions
	vector<double> _time;  // Time vector
	std::random_device rd;  // Obtain seed 
	std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis;  // (0, 1.0)
};

#endif