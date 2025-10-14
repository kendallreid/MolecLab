#include "Tile.h"

void Tile::readFromFile(const string& filename)
{
	fstream input(filename);
	if (input.good())
	{
		string line;
		bool readMatrix = false, readReactions = false;
		while (getline(input, line))
		{
			if (line == "---MATRIX---")
			{
				readMatrix = true;
				continue;
			}
			else if (line == "---REACTIONS---")
			{
				readMatrix = false;
				readReactions = true;
				continue;
			}
			if (readMatrix)
			{
				populateMatrix(line);
				continue;
			}
			else if (readReactions)
			{
				populateReactions(line);
				continue;
			}
		}
	}
	else { cout << "File not properly opened" << endl; }
	input.close();
}

void Tile::populateMatrix(string& line)
{
	if (line.empty())  // Check if reached end of matrix input within file
		return;
	vector<string> row;
	stringstream ss(line);
	string reactant;

	while (getline(ss, reactant, ','))
		row.push_back(reactant);  // Add reactant to proper row in matrix
	_pixelMatrix.push_back(row);  // Add rows to matrix 
}

void Tile::populateReactions(string& line)
{
	if (line.empty())  // Check if reached end of reaction input within file
		return;
	// Split into reaction and rate
	stringstream ss(line);
	string reaction, rate;
	getline(ss, reaction, ',');
	getline(ss, rate);

	// Split into reactants and products
	string reactant, product;
	stringstream reactSs(reaction);
	getline(reactSs, reactant, '=');
	getline(reactSs, product);

	// Parse reactants
	unordered_set<string> reactants;  // Set of all reactants
	stringstream reacSs(reactant);
	string val;
	while (getline(reacSs, val, '+'))
		reactants.insert(val);

	// Parse products
	stringstream prodSs(product);
	string first, second;
	getline(prodSs, first, '+');
	getline(prodSs, second);
	pair<string, string> products(first, second);

	// Create and store reaction
	Reaction rxn(reactants, products, stod(rate));
	_reactions.push_back(rxn);
}

void Tile::updateSizeParams()
{
	if (!_pixelMatrix.empty())  // check if empty
	{
		_rowSize = _pixelMatrix.size();
		_colSize = _pixelMatrix[0].size();  // How many elements in row
	}
	if (!_reactions.empty())
	{
		_numReactions = _reactions.size();
		_reactantPixelPairPos.resize(_numReactions);  // Rows = number of reactions
		_reacProp.resize(_numReactions);
	}
}

void Tile::initConc()
{
	for (const auto& row : _pixelMatrix)  // Loop through to get starting concentrations
		for (const auto& reactant : row)
			++_conc[reactant];  // Add 1 to concentration for each reactant found
}

void Tile::findPixelPairs()
{
	for (int row = 0; row < _rowSize; ++row)
	{
		for (int col = 0; col < _colSize; ++col)
		{
			if (col < _colSize - 1)  // Check if in bounds (next column) right
			{
				if (_pixelMatrix[row][col] != _pixelMatrix[row][col + 1])  // If reactants aren't the same
					populateReacPosVec({ row, col }, { row, col + 1 });  // Enter positions of each reactant into pixel pair vector

			}
			if (col > 0)  // Check in bounds (previous column) left
			{
				if (_pixelMatrix[row][col] != _pixelMatrix[row][col - 1])
					populateReacPosVec({ row, col }, { row, col - 1 });

			}
			if (row < _rowSize - 1)  // Check if in bounds (lower row) down
			{
				if (_pixelMatrix[row][col] != _pixelMatrix[row + 1][col])
					populateReacPosVec({ row, col }, { row + 1, col });

			}
			if (row > 0)  // Check in bounds (above row) up
			{
				if (_pixelMatrix[row][col] != _pixelMatrix[row - 1][col])
					populateReacPosVec({ row, col }, { row - 1, col });

			}
		}
	}
}

void Tile::populateReacPosVec(pair<int, int> pos1, pair<int, int> pos2)
{
	string reac1 = _pixelMatrix[pos1.first][pos1.second];
	string reac2 = _pixelMatrix[pos2.first][pos2.second];

	for (int i = 0; i < _reactions.size(); ++i)
	{
		if (bothInSet(reac1, reac2, _reactions[i]))  // If reactants match reaction add to pixel pair list
			_reactantPixelPairPos[i].emplace_back(pos1, pos2);
	}
}

bool Tile::bothInSet(const string& reac1, const string& reac2, const Reaction& rxn)
{
	return rxn.reactants.count(reac1) && rxn.reactants.count(reac2);
}

void Tile::tileSimStep()
{
	for (int i = 0; i < _reactantPixelPairPos.size(); ++i)  // Clear for next possible reactions
	{
		_reactantPixelPairPos[i].clear();
	}
	findPixelPairs();
	calcReacProp();
	calcTotalProp();
}

void Tile::printMatrix()
{
	for (const auto& row : _pixelMatrix)
	{
		for (const auto& pixel : row)
		{
			cout << pixel << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void Tile::calcReacProp()
{
	for (int i = 0; i < _reactions.size(); ++i)  // For each possible reaction, calculate propensity
	{
		_reacProp[i] = _reactantPixelPairPos[i].size() * _reactions[i].rate;
	}
}

void Tile::calcTotalProp()
{
	_totalProp = _reacProp[0] + _reacProp[1] + _reacProp[2];  // Sum _reacProp 
}

void Tile::updateMatrix(int rxn, int rxnIndex)
{
	pair<pair<int, int>, pair<int, int>> pair = _reactantPixelPairPos[rxn][rxnIndex];  // Get position of pair to update. Gives pair<pair, pair>

	for (int i = 0; i < _reactions.size(); ++i)
	{
		if (rxn == i)
		{
			_pixelMatrix[pair.first.first][pair.first.second] = _reactions[rxn].products.first;  // Update first pixel in matrix with new product
			_pixelMatrix[pair.second.first][pair.second.second] = _reactions[rxn].products.second;  // Update second pixel in matrix with new product

			updateConc(rxn);  // Update concentrations of matrix since matrix was changed
		}
	}
}

void Tile::updateConc(int rxn)
{
	for (int i = 0; i < _reactions.size(); ++i)
	{
		if (rxn == i)
		{
			++_conc[_reactions[rxn].products.first];  // Increase product conc
			++_conc[_reactions[rxn].products.second];  // Increase product conc

			for (const auto& reactant : _reactions[rxn].reactants)  // Decrease each reactant conc
			{
				if (_reactions[rxn].reactants.size() == 1)  // If reactants are the same
				{
					_conc[reactant] -= 2;
				}
				else  // Reactants are different
				{
					--_conc[reactant];
				}
			}
		} 
	}
}

const vector<pair<string, int>> Tile::concToVector()
{
	vector<pair<string, int>> concs;
	for (const auto& conc : _conc)  // Add current conc to tracker
	{
		concs.push_back({ conc.first, conc.second });
	}
	return concs;
}

