#ifndef TILE
#define TILE
#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <random>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

using std::string;
using std::pair;
using std::unordered_set;
using std::unordered_map;
using std::vector;
using std::fstream;
using std::stringstream;
using std::cout;
using std::endl;

struct Reaction {
	unordered_set<string> reactants;
	pair<string, string> products;
	double rate;

	Reaction(const unordered_set<string>& r, const pair<string, string>& p, const double k)
		: reactants(r), products(p), rate(k) {}
};

class Tile
{
public:
	Tile(const string& filename) : _rowSize(0), _colSize(0), _totalProp(0)
	{
		readFromFile(filename);  // Read matrix and reactions
		updateSizeParams();
		initConc();
	}
	~Tile() {};

	void readFromFile(const string& filename);
	void populateMatrix(string& line);  // Using file, input values to matrix
	void populateReactions(string& line);  // Using file, input reactions to list of possible 
	void updateSizeParams();  // Calculate size of matrix for later functions
	void initConc();  // Initilaize concentration values

	void findPixelPairs();  // Obtain all pixel pairs
	void populateReacPosVec(pair<int, int> pos1, pair<int, int> pos2);  // Fills vector with all pixel pairs - sorted by reaction numbers 
	bool bothInSet(const string& reac1, const string& reac2, const Reaction& rxn);  // Check if their is a reaction
	void tileSimStep();
	void printMatrix();

	void calcReacProp();
	void calcTotalProp();

	const double& getTotalProp() const { return _totalProp; }
	const vector<double>& getReacProp() const { return _reacProp; }
	const vector<vector<pair<pair<int, int>, pair<int, int>>>>& getPixelPairPos() const { return _reactantPixelPairPos; }
	const vector<Reaction>& getReactions() const { return _reactions; }
	const vector<vector<string>>& getPixelMatrix() const { return _pixelMatrix; }
	const unordered_map<string, int>& getConc() const { return _conc; }

	void updateMatrix(int rxn, int rxnIndex);
	void updateConc(int rxn);
	const vector<pair<string, int>> concToVector();

private:
	vector<vector<string>> _pixelMatrix;  // Creates a grid of pixels to represent tile - string representation of reactant
	vector<vector<pair<pair<int, int>, pair<int, int>>>> _reactantPixelPairPos;  // Vector holding position of reactants for each pixel pair (pair of coordinate pairs)
	int _rowSize, _colSize, _numReactions;
	double _totalProp;
	vector<double> _reacProp;  // Propensities for each reaction
	vector<Reaction> _reactions;  // List of reactions
	unordered_map<string, int> _conc;  // Key = reactant --> holds # of pixels as concentration
};
#endif