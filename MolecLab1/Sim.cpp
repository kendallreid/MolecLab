#include "Sim.h"

void Sim::nextReacTime(const double& totalProp)
{
	double dt = (-log(getRand())) / totalProp;  // Get next reaction time interval
	_timeTrack.push_back(dt + _timeTrack[_timeTrack.size() - 1]);  // Add current time plus dt to find next reaction time
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

// Updates matrix
void Sim::performReaction(const int& rxn, Tile& tile)
{
	int rxnIndex = getRandIndex(tile.getPixelPairPos()[rxn].size() - 1);  // Find a random index for specific reaction
	tile.updateMatrix(rxn, rxnIndex);  // Update matrix pixels
}

void Sim::initConc(Tile& tile)
{
	const auto& matrix = tile.getPixelMatrix();
	for (const auto& row : matrix)  // Loop through to get starting concentrations
		for (const auto& reactant : row)
			++_conc[reactant];  // Add 1 to concentration for each reactant found
}

void Sim::updateConcentrations(int reaction)
{
	switch (reaction)  // Based on which reaction
	{
	case 0:  // A+B=2U
		--_conc[0];  // Decrement
		--_conc[1];  // Decrement
		_conc[2] = _conc[2] + 2;  // Add 2
		break;

	case 1:  // A+U=2A
		++_conc[0];
		--_conc[2];
		break;

	case 2:  // B+U=2B
		++_conc[1];
		--_conc[2];
		break;
	}
}

int Sim::getRandIndex(int size)
{
	std::uniform_int_distribution<> dist(0, size);
	return dist(gen);
}

void Sim::simStep(Tile& tile)
{
	tile.tileSimStep();
	if (tile.getTotalProp() != 0)
	{
		// Determine time of next reaction
		nextReacTime(tile.getTotalProp());
		// Determine which reaction
		int reaction = determineReaction(tile.getTotalProp(), tile.getReacProp());
		// Execute Reaction
		performReaction(reaction, tile);
		updateConcentrations(reaction);

		_concOverTime[0].push_back(_conc[0]);  // Store concentrations
		_concOverTime[1].push_back(_conc[1]);
		_concOverTime[2].push_back(_conc[2]);
	}
}

void Sim::runSim(double maxTime)
{
	Tile tile("input.csv");
	// Calculate propensities to have starting values
	/*tile.populateMatrix();
	tile.updateSizeParams();*/
	tile.tileSimStep();

	// Open data file for plotting values
	std::ofstream dataFile("concentrations.txt");
	dataFile << "Time A B U" << std::endl;

	while (_time[_time.size() - 1] < maxTime && tile.getTotalProp() > 0)  // Time not run out & reactions still possible
	{
		tile.printMatrix();
		simStep(tile);  // Run reaction 

		// Log current time and concentrations to file
		double currentTime = _time[_time.size() - 1];
		dataFile << currentTime << " " << _conc[0] << " " << _conc[1] << " " << _conc[2] << std::endl;
	}
	dataFile.close();
	std::system("gnuplot -p -e \"plot 'concentrations.txt' using 1:2 with lines title 'A', '' using 1:3 with lines title 'B', '' using 1:4 with lines title 'U'\"");
}
