#ifndef TILE
#define TILE
#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <utility>
#include <random>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

using std::string;
using std::pair;
using std::unordered_set;
using std::vector;
using std::fstream;
using std::stringstream;
using std::cout;
using std::endl;

class Tile
{
public:
	// update resize and rules and products manually???
	Tile(const int row = 0, const int col = 0, const double k1 = 1, const double k2 = 1, const double k3 = 1, const double totalProp = 0) :
		_rowSize(row), _colSize(col), _reacRate{ k1, k2, k3 }, _totalProp(totalProp), rule1{ "A", "B" }, rule2{ "A", "U" },
		rule3{ "B", "U" }, prod1{ "U" , "U" }, prod2{ "A", "A" }, prod3{ "B", "B" }
	{
		// Number of reactions for program (update if different)
		_reactantPixelPairPos.resize(3);
		_reacProp.resize(3);
	};
	~Tile() {};

	void populateMatrix();  // Using file, input values to matrix
	void updateSizeParams();  // Calculate size of matrix for later functions
	void findPixelPairs();  // Obtain all pixel pairs
	void populateReacPosVec(pair<int, int> pos1, pair<int, int> pos2);  // Fills vector with all pixel pairs - sorted by reaction numbers 
	bool bothInSet(const string& reac1, const string& reac2, const unordered_set<string>& rule);  // Check if their is a reaction
	void tileSimStep();
	void printMatrix();

	void calcReacProp();
	void calcTotalProp();

	const double& getTotalProp() const { return _totalProp; }
	const vector<double>& getReacProp() const { return _reacProp; }
	const vector<vector<pair<pair<int, int>, pair<int, int>>>>& getPixelPairPos() const { return _reactantPixelPairPos; }

	void updateMatrix(int rxn, int rxnIndex);

	// Possible reactants and their products
	const unordered_set<string> rule1;  // A+B=2U
	const unordered_set<string> rule2;  // A+U=2A
	const unordered_set<string> rule3;  // B+U=2B
	const pair<string, string> prod1;
	const pair<string, string> prod2;
	const pair<string, string> prod3;


private:
	vector<vector<string>> _pixelMatrix;  // Creates a grid of pixels to represent tile - string representation of reactant
	vector<vector<pair<pair<int, int>, pair<int, int>>>> _reactantPixelPairPos;  // Vector holding position of reactants for each pixel pair (pair of coordinate pairs)
	int _rowSize, _colSize, _numReactions;
	double _totalProp;
	vector<double> _reacRate;  // Reaction rates for each possible reaction
	vector<double> _reacProp;  // Propensities for each reaction
	//std::random_device rd;  // Obtain seed 
	//std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
	//std::uniform_real_distribution<> dis;  // (0, 1.0)
};

#endif