#include "Sim.h"

void Sim::nextReacTime(const double& totalProp)
{
	double dt = (-log(getRand())) / totalProp;  // Get next reaction time interval
	_timeTrack.push_back(dt + _timeTrack[_timeTrack.size() - 1]);  // Add current time plus dt to find next reaction time
}

int Sim::getRandIndex(int size)
{
	std::uniform_int_distribution<> dist(0, size);
	return dist(gen);
}

int Sim::determineReaction(const double& totalProp, const vector<double>& reacProp)
{
	double stepReac = getRand() * totalProp;  // Random value (0 - 1.0) to determine reaction
	double cumulativeSum = 0.0;

	for (int j = 0; j < reacProp.size(); ++j)  // Size is number of reactions possible
	{
		cumulativeSum += reacProp[j];
		if (cumulativeSum >= stepReac)  // Reaction j occurred
		{
			_rxnTrack.push_back(j);  // Add reaction to rxn vector
			return j;
		}
	}
}

// Updates matrix by performing reaction
void Sim::performReaction(const int& rxn, Tile& tile)
{
	int rxnIndex = getRandIndex(tile.getPixelPairPos()[rxn].size() - 1);  // Find a random index for specific reaction
	tile.updateMatrix(rxn, rxnIndex);  // Update matrix pixels
}

void Sim::simStep(Tile& tile)
{
	tile.tileSimStep();  // Pixel pairs and propensities
	if (tile.getTotalProp() != 0)  // Check if any possible reactions left
	{
		nextReacTime(tile.getTotalProp());
		int reaction = determineReaction(tile.getTotalProp(), tile.getReacProp());
		performReaction(reaction, tile);
		
		vector<pair<string, int>> concs;
		for (const auto& conc : tile.getConc())  // Add current conc to tracker
		{
			concs.push_back({ conc.first, conc.second });
		}
		_concOverTime.push_back(concs);
	}
}

void Sim::runSim(double maxTime)
{
	Tile tile("input.csv");

	tile.tileSimStep();  //NEED???????????????

	// Open data file for plotting values
	std::ofstream dataFile("concentrations.txt");
	dataFile << "Time - Reactants" << std::endl;

	while (_timeTrack[_timeTrack.size() - 1] < maxTime && tile.getTotalProp() > 0)  // Time not run out & reactions still possible
	{
		tile.printMatrix();  ///////////////////VISUALS////////////////////////////
		simStep(tile);  // Run reaction 

		// Log current time and concentrations to file
		double currentTime = _timeTrack[_timeTrack.size() - 1];
		dataFile << currentTime << " - ";
		if (!_concOverTime.empty())
		{
			for (const auto& conc : _concOverTime[_concOverTime.size() - 1])  // Print each reactants concentration
			{
				dataFile << "(" << conc.first << ", " << conc.second << ") ";
			}
		}
		dataFile << endl;
	}
	dataFile.close();
	//std::system("gnuplot -p -e \"plot 'concentrations.txt' using 1:2 with lines title 'A', '' using 1:3 with lines title 'B', '' using 1:4 with lines title 'U'\"");
}
