#include "FFA.h"
#include "FFA_v2.h"

#include "Timer.h"
void stat(std::vector<double>results,int it){
    double sum = 0.0;
    double sq_sum = 0.0;
    for (int i = 0; i < it; ++i) {
        std::cout<<results[i]<<std::endl;
        sum += results[i];
        sq_sum += results[i] * results[i];
    }
    double moyenne = sum / it;
    double ecartType = std::sqrt(sq_sum / 10 - moyenne * moyenne);

    // Print the mean and standard deviation of the results
    std::cout << "Moyenne = " << moyenne << std::endl;
    std::cout << "�cart type = " << ecartType << std::endl;
}

double fonctionObjective(std::vector<double> x)
{
	double res=0;
	res = pow(x[0], 3) + pow(x[1], 2) + pow(x[2], 4);
	return res;
}
void testFFA()
{
	originalFFA ffa{ 20,50 };
	std::vector<double> lb = { -1,5,-7 }, ub{ 4,10,-4 };
	std::function<double(std::vector<double>)> funObj = fonctionObjective;
	std::vector<double> g_best=ffa.solve(lb,ub , funObj,"min");


	std::cout << "best solution = [ ";
	for (size_t i = 0; i < g_best.size(); ++i) {
		std::cout << g_best[i] << (i == g_best.size() - 1 ? "]" : ", ");
	}
	std::cout << std::endl;
	std::cout << "best solution : " << funObj(g_best) << "\n";
}

// Define the Rosenbrock function as the objective function
double fonctionRosnebork(const std::vector<double>& x) {
	double sum = 0.0;
	for (size_t i = 0; i < x.size() - 1; ++i) {
		double diff = x[i + 1] - x[i] * x[i];
		sum += 100.0 * diff * diff + (1.0 - x[i]) * (1.0 - x[i]);
	}
	return sum;
}
void testRosenbork()
{
    int it = 1;

	int dim = 30;
	int popSize = 30;
    std::vector<double> results(it);
//     std::vector<vector<double>> resultsfinal(it,vector<double>(dim) );
	for (int i = 0;i<it;i++){
       originalFFA ffa{2000,popSize };
        std::vector<double> lb(dim, -10), ub(dim, 10);
        std::function<double(std::vector<double>)> funObj = fonctionRosnebork;
        std::vector<double> g_best = ffa.solve(lb, ub, funObj,"min");


        std::cout << "best solution = [ ";
        for (size_t i = 0; i < g_best.size(); ++i) {
            std::cout << g_best[i] << (i == g_best.size() - 1 ? "]" : ", ");
        }
        std::cout << std::endl;
        std::cout << "best solution : " << funObj(g_best) << "\n";
        results[i] = funObj(g_best);
        //resultsfinal[i]=g_best;
	}
	stat(results,it);
//	for(int i =0;i<resultsfinal.size();i++){
//        for(int j =0;j<resultsfinal[i].size();j++){
//            std::cout<<resultsfinal[i][j]<<"  |  "
//        }
//	}


}


double fonctionRastrigin(std::vector<double> x)
{
	const double A = 10.0;
	int n = x.size();
	double sum = A * n;

	for (int i = 0; i < n; ++i) {
		sum += (x[i] * x[i] - A * cos(2.0 * M_PI * x[i]));
	}

	return sum;
}
void testRastrigin()
{
	//Timer t("testRastrigin");
	int dim = 30;
	int popSize = 30;


    // Run the optimization multiple times
    int it =20;
    std::vector<double> results(it);
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < it; ++i) {
        //originalFFA(int epoch=1000, int popSize=100, float gamma=0.001, float beta_base=2,
        //	float alpha=0.2, float alpha_damp=0.99, float  delta=0.05, int exponent=2 );
        //10000, 30, 0.001, 0.8351, 0.15009, 0.99, 0.05, 2
        originalFFA ffa{10000 ,popSize, 0.001, 0.8351, 0.2, 0.999, 0.05, 2 };
        std::vector<double> lb(dim,-5.12), ub(dim,5.12);
        std::function<double(std::vector<double>)> funObj = fonctionRastrigin;
        std::vector<double> g_best = ffa.solve(lb, ub, funObj,"min");
        std::cout << "best solution = [ ";
        for (size_t i = 0; i < g_best.size(); ++i) {
            std::cout << g_best[i] << (i == g_best.size() - 1 ? "]" : ", ");
        }
        std::cout << std::endl;
        std::cout << "best solution : " << funObj(g_best) << "\n";

        results[i] = funObj(g_best);
    }

	stat(results,it);

}

double fonctionMichalewicz(const std::vector<double>& x) {
	int n = x.size();
	double m = 10;
	double sum = 0.0;

	for (int i = 0; i < n; ++i) {
		sum += -sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / M_PI), 2 * m);
	}

	return sum;
}
void testMichalewicz()
{
	int dim = 30;
	int popSize = 30;
	originalFFA ffa{ 1000,popSize, 0.001, 0.8351, 0.15009, 0.99, 0.05, 2};
	std::vector<double> lb(dim, 0), ub(dim, PI);
	std::function<double(std::vector<double>)> funObj = fonctionMichalewicz;
	std::vector<double> g_best = ffa.solve(lb, ub, funObj, "min");

	std::cout << "best solution = [ ";
	for (size_t i = 0; i < g_best.size(); ++i) {
		std::cout << g_best[i] << (i == g_best.size() - 1 ? "]" : ", ");
	}
	std::cout << std::endl;
	std::cout << "best solution : " << funObj(g_best) << "\n";
}

/// <summary>
/// test of version 2 FFA
/// </summary>
void testFFA_V2()
{
	// Define the problem
	std::vector<double> lb = { -1,5,-7 }, ub{ 4,10,-4 };

	// Run the optimization
	FFA_V2 ffa(20, 50, 0.001, 2.0, 0.2, 0.99, 0.05, 2);
	std::vector<double> best_x = ffa.solve(fonctionObjective, lb, ub);

	// Print the best solution
	std::cout << "Best solution = [";
	for (size_t i = 0; i < best_x.size(); ++i) {
		std::cout << best_x[i] << (i == best_x.size() - 1 ? "]" : ", ");
	}
	std::cout << std::endl;
	std::cout << "best solution : " << fonctionObjective(best_x) << "\n";
}

void testRosenbrockFFA_V2()
{
	int dim = 30;
	// Define the problem
	std::vector<double> lb(dim,-10), ub(dim,10);

	// Run the optimization
	FFA_V2 ffa(2000, 30, 0.001, 2.0, 0.2, 0.99, 0.05, 2);
	std::vector<double> best_x = ffa.solve(fonctionRosnebork, lb, ub);

	// Print the best solution
	std::cout << "Best solution = [";
	for (size_t i = 0; i < best_x.size(); ++i) {
		std::cout << best_x[i] << (i == best_x.size() - 1 ? "]" : ", ");
	}
	std::cout << std::endl;
	std::cout << "best solution : " << fonctionRosnebork(best_x) << "\n";

}

void testRastriginFFA_V2()
{
	// Define the problem
	int dim = 30;
	std::vector<double> lb(dim, -5.12), ub(dim, 5.12);
	// Run the optimization
	FFA_V2 ffa(10000, 30, 0.001, 0.8351, 0.15009, 0.99, 0.05, 2);
	std::vector<double> best_x = ffa.solve(fonctionRastrigin, lb, ub);

	// Print the best solution
	std::cout << "Best solution = [";
	for (size_t i = 0; i < best_x.size(); ++i) {
		std::cout << best_x[i] << (i == best_x.size() - 1 ? "]" : ", ");
	}
	std::cout << std::endl;
	std::cout << "best solution : " << fonctionRastrigin(best_x) << "\n";


}
void testMichalewiczFFA_V2()
{
	// Define the problem
	int dim = 15;
	std::vector<double> lb(dim, 0), ub(dim, PI);
	// Run the optimization
	FFA_V2 ffa(100, 30, 0.001, 2.0, 0.2, 0.99, 0.05, 2);
	std::vector<double> best_x = ffa.solve(fonctionMichalewicz, lb, ub);

	// Print the best solution
	std::cout << "Best solution = [";
	for (size_t i = 0; i < best_x.size(); ++i) {
		std::cout << best_x[i] << (i == best_x.size() - 1 ? "]" : ", ");
	}
	std::cout << std::endl;
	std::cout << "best solution : " << fonctionMichalewicz(best_x) << "\n";


}

double schwefelFunction(const std::vector<double>& x) {
    double sum = 0.0;
    int d = x.size();

    for (int i = 0; i < d; ++i) {
        sum += x[i] * sin(sqrt(fabs(x[i])));
    }

    return 418.9829 * d - sum;
}
void testschwefelFunction(){
int dim = 30;
	int popSize = 30;
	//originalFFA(int epoch=1000, int popSize=100, float gamma=0.001, float beta_base=2,
        //	float alpha=0.2, float alpha_damp=0.99, float  delta=0.05, int exponent=2 );
	originalFFA ffa{ 10000, 30, 0.1, 0.5, 0.8, 0.99, 0.05, 2};
	std::vector<double> lb(dim, -500), ub(dim, 500);
	std::function<double(std::vector<double>)> funObj = schwefelFunction;
	std::vector<double> g_best = ffa.solve(lb, ub, funObj, "min");

	std::cout << "best solution = [ ";
	for (size_t i = 0; i < g_best.size(); ++i) {
		std::cout << g_best[i] << (i == g_best.size() - 1 ? "]" : ", ");
	}
	std::cout << std::endl;
	std::cout << "best solution : " << funObj(g_best) << "\n";
}



int main()
{
	srand((unsigned int)std::time(0));
	try
	{


		//{

		//	std::cout << "Execution :" << i << "\n";
			//testFFA();
			// std::cout<<std::endl;
			//testFFA_V2();
			//testRosenbork();
			// std::cout<<std::endl;
			testRosenbrockFFA_V2();
			//testRastrigin();

			 //std::cout<<std::endl;
			//testRastriginFFA_V2();
			//testMichalewicz();
//			 std::cout<<std::endl;
//			testMichalewiczFFA_V2();
            //testschwefelFunction();
			//std::cout << "**************************************************\n";
		//}

	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}

	return 0;
}
