#include "FFA.h"

#include <random>

bool compareFitness(std::function<double(std::vector<double>)> func ,std::vector<double> pop1, std::vector<double> pop2, const std::string& minmax)
{
	bool resMin{ func(pop1) < func(pop2) };
	if (minmax == "min")
		//return func(pop1) < func(pop2);
		return resMin;
	else
		//return func(pop1) > func(pop2);
		return !resMin;
}

std::vector<double> get_gBest(std::vector<std::vector<double>> population, std::function<double(std::vector<double>)> func,const std::string& minmax)
{
	//std::vector<double> fitness(population.size());
	//fitness[0] = func(population[0]);
	double gbestIdx = 0;
	for (int i = 1; i < population.size(); ++i)
	{

		//fitness[i] = func(population[i]);
		if ( compareFitness( func, population[i], population[gbestIdx], minmax) )
			gbestIdx = i;
	}
	return population[gbestIdx];
}


void originalFFA::afficherPopulation() const
{
	for (int i = 0; i < d_population.size(); ++i)
	{
		for (int j = 0; j < d_population[i].size(); ++j)
		{
			std::cout<<d_population[i][j]<<" | ";
		}
	std::cout << "\n";
	}
	std::cout << "\n";

}



//random double entre 0 et 1
double randomDouble()
{
	return rand() * 1.0 / RAND_MAX;
}
originalFFA::originalFFA(int epoch, int popSize,float gamma, float beta_base,float alpha , float alpha_damp , float  delta , int exponent )
	: d_epoch{epoch} , d_popSize{popSize} , d_gamma{gamma}, d_beta_base{beta_base}, d_alpha{alpha}, d_alpha_damp{alpha_damp}, d_delta{delta}
	, d_exponent{ exponent }, d_population(popSize), d_dyn_alpha{d_alpha}, d_gbestHistory()
{}

std::vector<double> originalFFA::solve(std::vector<double>& lb, std::vector<double>& ub, std::function<double(std::vector<double>)> func,const std::string& minmax)
{

	//std::default_random_engine generator;
	//std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
	//comparer la dimension du upper/lower bound
	if (lb.size() != ub.size() )
		throw std::runtime_error("la dimension de lower and upper bound est differente");
	//comparer la dimension du population  avec lower/ upper bound


	//generation de la population initial
	for (auto& x : d_population) {
		x.resize(lb.size());
		for (size_t i = 0; i < x.size(); ++i) {
			x[i] = lb[i] + ( randomDouble() * (ub[i] - lb[i]) );
		}
	}

	//affichage pour le test
	//std::cout << "population 0 "<< '\n';
	//afficherPopulation();
	//end affichage

	for (int i = 1; i <= d_epoch; ++i)
	{
		//affichage pour le test
		/*std::cout << "population " << i << '\n';
		afficherPopulation();*/
		//end affichage

		for (int idx = 0; idx < d_popSize; ++idx)
		{
			//On peut remplacer agent avec d_population de idx
			std::vector<double> agent=d_population[idx];

			//calcule la distance maximale
			double dmax =sqrt(lb.size());

			//popChild stock tous les mouvements possible de idx
			//puis � la fin de la boucle on prend la meilleur => local_gbest
			std::vector<std::vector<double>> popChild;

			for (int j = 0  ; j < d_popSize; ++j)
			{

				//if fitness value de la population j  > fitness value de la population idx
				if(compareFitness(func,d_population[j], agent,minmax))
				{
					//double rij = distance(agent,d_population[j])/ dmax;
					//calcule du rayon entre i et idx
					double rij{ 0.0 };
					for (size_t k = 0; k < agent.size(); ++k)
					{
							rij += std::pow(agent[k] - d_population[j][k], 2);
					}
					rij = std::sqrt(rij / (dmax * dmax));

					//calcule du beta
					double beta = d_beta_base * exp(-d_gamma * pow(rij,d_exponent));

					//calcule du mutation vector
					//std::vector<double> mutationVector = vectorMultCoef(d_delta ,generateVector(0,1,d_dimension));
					std::vector<double> mutationVector(agent.size());
					for (size_t k = 0; k < agent.size(); ++k)
					{
						mutationVector[k] = d_delta * (randomDouble()-0.5 );
					}

					//calcule du temp
					//std::vector<double> temp= vectorMultVector(vectorMoinsVector(d_population[j], agent), generateVector(0, 1, d_dimension));
					std::vector<double> temp(agent.size());
					for (size_t k = 0; k < agent.size(); ++k) {
						temp[k] = (d_population[j][k] - agent[k]) * randomDouble();
					}

					//calcule de la nouvelle position
					//  move towards brighter firfly : newPosition(t+1) = CurrentPosition(t) + Attractivness*(xj(t)-xi(t)) + stepSize
					//std::vector<double> newPos = vectorPlusVector(vectorPlusVector(vectorMultCoef(d_dyn_alpha, mutationVector), vectorMultCoef(beta, temp)) , agent);
					std::vector<double> newPos(agent.size());
					for (size_t k = 0; k < agent.size(); ++k) {
						newPos[k] = agent[k] + d_alpha * mutationVector[k] + beta * temp[k];
						// verification des bounds d la nouvelle position
						if (newPos[k] < lb[k]) newPos[k] = lb[k];
						if (newPos[k] > ub[k]) newPos[k] = ub[k];
					}

					//agent = newPos;

					//ajouter la nouvelle solution calculer au tableau des solutions possible pour idx : popChild
					popChild.push_back(newPos);
				}
			}

			//remplacer les agents qui n'ont pas subi de changement par des agents aleatoire
			//if (popChild.size() < d_popSize)
			//{
			//	for (int k = popChild.size(); k < d_popSize ; ++k)
			//		popChild.push_back(randVector(lb,ub,d_dimension));
			//}
			for (size_t j = popChild.size(); j < d_popSize; ++j) {
				std::vector<double> randVect(agent.size());
				for (size_t k = 0; k < agent.size(); ++k) {
					randVect[k] = lb[k] + randomDouble() * (ub[k] - lb[k]);
				}
				popChild.push_back(randVect);
			}

			//local_best la meilleur solution que la firfly[idx] pour faire
			std::vector<double> localBest = get_gBest(popChild,func, minmax);

			//Comparer la fitness de la localBest avec l'agent
			if (compareFitness(func, localBest, agent, minmax))
				//changer l'agent avec la meilleur solution dans le tableau devant idx
				d_population[idx] = localBest;
		}

		//il existe dans le code python
		//d_population.push_back(get_gBest(d_population,func));

		d_dyn_alpha = d_alpha_damp * d_alpha;

		//we don't have history of g_best and g_worst so we don't need to put it inside for epoch loop
		// because for every iteration in this loop g_best get better
		auto g_best = get_gBest(d_population, func, minmax);
		d_gbestHistory.push_back(g_best);

		////Affichage de history de g_best
		//std::cout << "History g_best :\n";
		//for (int i = 0; i < d_gbestHistory.size(); ++i)
		//{
		//	for (int j = 0; j < d_gbestHistory[i].size(); ++j)
		//	{
		//		std::cout << d_gbestHistory[i][j] << " | ";
		//	}
		//	std::cout << "\n";
		//}
		//std::cout << "\n";
		////end affcihage history g_best


		std::cout << "Epoch " << i  << ", best objective = " << func(g_best) << std::endl;//", best solution = [";
		//for (size_t i = 0; i < g_best.size(); ++i) {
		//	std::cout << g_best[i] << (i == g_best.size() - 1 ? "]" : ", ");
		//}
		//std::cout << std::endl;
	}

	//return la meilleur solution g_best de la dernier epoch
	return get_gBest(d_population,func,minmax);

}
