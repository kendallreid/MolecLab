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
		_concOverTime.push_back(tile.concToVector());  // Add new concentrations to tracker
	}
}

void Sim::printConcToFile(ofstream& dataFile)
{
	double currentTime = _timeTrack[_timeTrack.size() - 1];
	dataFile << currentTime << " ";
	if (!_concOverTime.empty())
	{
		for (const auto& conc : _concOverTime[_concOverTime.size() - 1])  // Print each reactants concentration
		{
			dataFile << conc.second << " ";
		}
	}
	dataFile << endl;
}

void Sim::runSim(double maxTime)
{
	Tile tile("input.csv");

	_concOverTime.push_back(tile.concToVector());  // Starting concentrations
	tile.tileSimStep();  // Calculates starting propensities so sim can run

	// Open data file for plotting values
	ofstream dataFile("concentrations.txt");
	dataFile << "Time ";
	for (const auto& reactant : _concOverTime[0])  // Print all reactants in header
	{
		dataFile << reactant.first << " ";
	}
	dataFile << endl;

	while (_timeTrack[_timeTrack.size() - 1] < maxTime && tile.getTotalProp() > 0)  // Time not run out & reactions still possible
	{
		tile.printMatrix();  ///////////////////VISUALS////////////////////////////
		printConcToFile(dataFile);
		simStep(tile);  // Run reaction 
	}
	dataFile.close();
	createPlot();
}

void Sim::createPlot()
{
	std::ofstream gp("plot.gp");
	gp << "set key autotitle columnhead\n";
	gp << "set xlabel 'Time'\n";
	gp << "set ylabel 'Concentration'\n";
	gp << "set datafile separator ' '\n";  // <-- use space as separator
	gp << "plot for [col=2:*] 'concentrations.txt' using 1:col with lines lw 2\n";
	gp.close();

	std::system("gnuplot -p plot.gp");
}
