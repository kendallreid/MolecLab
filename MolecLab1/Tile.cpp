#include "Tile.h"

void Tile::populateMatrix()
{
	fstream input;
	input.open("input.csv");
	if (input.good())
	{
		string line;
		while (getline(input, line))  // Read line of file
		{
			vector<string> row;
			stringstream ss(line);
			string reactant;
			while (getline(ss, reactant, ','))
			{
				row.push_back(reactant);  // Add reactant to proper row in matrix
			}
			_pixelMatrix.push_back(row);  // Add rows to matrix 
		}
		updateSizeParams();
	}
	else { cout << "file not properly opened" << endl; }
	input.close();
}

void Tile::updateSizeParams()
{
	if (!_pixelMatrix.empty())  // check if empty
	{
		_rowSize = _pixelMatrix.size();
		_colSize = _pixelMatrix[0].size();  // How many elements in row
	}
}

void Tile::findPixelPairs()
{
	for (int row = 0; row < _rowSize; ++row)
	{
		for (int col = 0; col < _colSize; ++col)
		{
			if (col < _colSize - 1)  // Check if in bounds (next column)
			{
				if (_pixelMatrix[row][col] != _pixelMatrix[row][col + 1])  // If reactants aren't the same
					populateReacPosVec({ row, col }, { row, col + 1 });  // Enter positions of each reactant into pixel pair vector

			}
			if (col > 0)  // Check in bounds (previous column)
			{
				if (_pixelMatrix[row][col] != _pixelMatrix[row][col - 1])
					populateReacPosVec({ row, col }, { row, col - 1 });

			}
			if (row < _rowSize - 1)  // Check if in bounds (lower row)
			{
				if (_pixelMatrix[row][col] != _pixelMatrix[row + 1][col])
					populateReacPosVec({ row, col }, { row + 1, col });

			}
			if (row > 0)  // Check in bounds (above row)
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

	// Add pixel pair to possible reactions
	if (bothInSet(reac1, reac2, rule1))
		_reactantPixelPairPos[0].emplace_back(pos1, pos2);  // Reaction 1
	if (bothInSet(reac1, reac2, rule2))
		_reactantPixelPairPos[1].emplace_back(pos1, pos2);  // Reaction 2
	if (bothInSet(reac1, reac2, rule3))
		_reactantPixelPairPos[2].emplace_back(pos1, pos2);  // reaction 3
}

bool Tile::bothInSet(const string& reac1, const string& reac2, const unordered_set<string>& rule)
{
	return rule.count(reac1) && rule.count(reac2);
}

void Tile::tileSimStep()
{
	for (int i = 0; i < _reactantPixelPairPos.size(); ++i)
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
	_reacProp[0] = _reactantPixelPairPos[0].size() * _reacRate[0];  // Number of reaction pairs * k1
	_reacProp[1] = _reactantPixelPairPos[1].size() * _reacRate[1];  // Number of reaction pairs * k2
	_reacProp[2] = _reactantPixelPairPos[2].size() * _reacRate[2];  // Number of reaction pairs * k3
}

void Tile::calcTotalProp()
{
	_totalProp = _reacProp[0] + _reacProp[1] + _reacProp[2];  // Sum _reacProp 
}

void Tile::updateMatrix(int rxn, int rxnIndex)
{
	pair<pair<int, int>, pair<int, int>> pair = _reactantPixelPairPos[rxn][rxnIndex];  // Get position of pair to update

	switch (rxn)
	{
	case 0:  // A+B=2U
		// Update matrix after reaction occurs
		_pixelMatrix[pair.first.first][pair.first.second] = prod1.first;
		_pixelMatrix[pair.second.first][pair.second.second] = prod1.second;
		break;
	case 1:  // A+U=2A
		_pixelMatrix[pair.first.first][pair.first.second] = prod2.first;
		_pixelMatrix[pair.second.first][pair.second.second] = prod2.second;
		break;
	case 2:  // B+U=2B
		_pixelMatrix[pair.first.first][pair.first.second] = prod3.first;
		_pixelMatrix[pair.second.first][pair.second.second] = prod3.second;
		break;
	}
}
