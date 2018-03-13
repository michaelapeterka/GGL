#include<iostream>
#include<stdexcept>
#include<vector>
#include<string>
#include<ctime>
#include<stdlib.h>
#include<cmath>

using namespace std;

double partial_sum(vector<double>& a_my, double index)
{
	double sum = 0.0;
	for(int i = 0; i <= index;++i)
	{
		sum += a_my.at(i); 
		
	} 
	return sum; 
}
int main()
{
	//Zeit
	double t = 0.0; 
	
	//Anfangspopulation
	
	double population_A = 0;
	double population_B = 0;
	double population_AB = 5;  
	
	//Anzahl der Reaktionen
	double M = 2; 
	
	//Population zum bestimmten Zeitpunkt t+tau
	
	double X; 
	//rate
	vector<double> c_mu;
	c_mu.push_back(0.01);
	c_mu.push_back(0.01);
	
	
	//reaction_counter
	double n=0;
	
	//maximale Anzahl an Reaktionen
	double n_max = 10; 
	
	//zufallszahlengenerator initialisieren
	srand(time(NULL));
	
	while (n<n_max)
	{
	
	
	cout << "n: " << n << endl << endl; 
	
	cout << "t: " << t << endl; 
	 cout << "Population:" << endl;  
	 cout << "A : "<< population_A << endl;
	 cout << "B : "<< population_B << endl;
	 cout << "AB: "<< population_AB << endl; 
	 
	 
	 
	//1)Generieren der Zufalsszahlen
	
	  double r0 = double(rand())/(RAND_MAX);
	  double r1 = double(rand())/(RAND_MAX);
	  double r2 = double(rand())/(RAND_MAX);
	  
	  
	  //Berechnen von h_my --> Anzahl der Reaktanten-Kombis die im aktuellen Zustand möglich sind
	  
	  vector<double> h_mu;
	  h_mu.push_back(population_A*population_B);
	  h_mu.push_back(population_AB);
	  //Berechnen der propensity function a_i - i = 1,....,Anzahl der Reaktionen für jede Reaktion 
	  
	  /*double a1 = 0.0;
	  double a2 = 0.0; 
	  a1 = population_A*population_B*c1;
	  a2 = population_AB*c2;*/
	  
	  vector<double> a;
	  for (int ii = 0; ii < h_mu.size(); ++ii)
	  {
	  	a.push_back(h_mu.at(ii)*c_mu.at(ii)); 
	  }

	 //Berechnen der Summe der propensity functions a_i;
	 
	 double a0 = 0.0;
	 for (int j = 0; j < a.size(); ++j)
	 {
	 	a0 = a0 + a.at(j);
	}
	
	 
	  
	 
	 //Berechnen der zeit wann die nächste reaktion stattfinden wird
	 
	 double tau = 0.0;
	 tau = (1/a0)*log(1/r1);
	 
	 //Berechne welche Reaktion passieren wird
	 double reaction = 0.0;
	
	
	double sum = 0.0;
	
	/*for (int i = 0; i < M ; i++)
	{
		sum = sum + a.at(i);
		if ((r2*a0) < sum)
		{
			reaction = i;
			cout << "reaction taken: " << reaction << endl; 
			break; 
		}
	}*/
	
	for (int j = 0; j < M ; ++j)
	  {
	  	
		if (partial_sum(a,j-1) < (r2*a0) && (r2*a0)<=partial_sum(a,j))
		{
				reaction = j; 
				cout << "Das ist jetzt die regel die genommen wird: " << j << endl; 
				break; 	
		}
		
	  }	
	
	
	if (reaction == 0)
	{
		population_A -= 1;
		population_B -= 1;
		population_AB += 1;
	}
	
	if (reaction == 1)
	{
		population_A += 1;
		population_B += 1;
		population_AB -= 1; 
	}
	
	t = t+ tau; 
	n = n+1; 
	}
	return 0; 
}
