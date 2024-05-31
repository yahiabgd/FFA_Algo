#pragma once


#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>

#define PI 3.14159265




class originalFFA
{
public:

	originalFFA(int epoch=1000, int popSize=100, float gamma=0.001, float beta_base=2,
		float alpha=0.2, float alpha_damp=0.99, float  delta=0.05, int exponent=2 );
// we took ub and lb as vectors because the dimension might be dim > 1
// but in general the values are the same in lb[i] ,ex : lb[0]=lb[1]=...=lb[n]=-100
//so lb[i] means the bound of axe[i] and the same with the ub[i]
//	@param[in] vector - lb is the lower bound where the fireflies can move
//	@param[in] vector - ub is the uper bound where the fireflies can move
	//initial position: we generate a random positionnig of firelies to start the algorithme
	//PositionINIT  = lb+(ub-lb)*randvector(dimension) => the position is restricted by lb and ub
//  @param[in] func - function that calculate the fitness value(brightness) of each firefly of a given dimension 
	//brightness is propotional to attractivness.
	//fitness[i] = SUM(x[i] -1)^2 ; i[1,dimension] ; xi = func(position)
	//with xi is the value calculated(returned) by param func(objective function)
//  Update firlies positions:
//  for each firefly we update position 
//  we amend position based on the fitness value(brithness) of firfly:
	//  move towards brighter firfly : newPosition(t+1) = CurrentPosition(t) + Attractivness*(xj(t)-xi(t)) +stepSize
		//newPosition - vector : The position of the firefly after the update
		//CurrentPosition - vector : The position of the firefly before the update
		//Attractivness : Beta * exp(-lambda*(r^2))
			//r : the distance between two fireflies i at position x[i] and j at position x[j] : sqrt(SUM((xi-xj)^2))
			//beta : beta0 * exp(-d_gamma * r^2)); beat0 = 1
		//stepSize : d_alpha*randomvector(dimension)
		//if the brightness is the same move randomlly
	std::vector<double> solve(std::vector<double> &lb, std::vector<double> &ub, std::function<double(std::vector<double>)> func,const std::string& minmax);
	void afficherPopulation() const;
private:
	int d_epoch;
	int d_popSize;
	int d_exponent;
	float d_gamma, d_beta_base, d_alpha, d_alpha_damp, d_delta;
	std::vector<std::vector<double>> d_population;
	float d_dyn_alpha;
	std::vector<std::vector<double>> d_gbestHistory;
};

