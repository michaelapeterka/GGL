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
	//1 Initzialisierung
	int M = 2;
	int N = 3; 
	
	//cy --> Reaktionsrate die die Reaktions repräsentiert
	vector<double> cy;
	cy.push_back(0.2);
	cy.push_back(0.1);
	
	//Startwert der Population der 3 Moleküle
	vector<double> pop_start;
	pop_start.push_back(0);
	pop_start.push_back(0);
	pop_start.push_back(5);
	
	/*for (int i = 0; i < pop_start.size();i++)
	{
		cout << "pop_start: " << i << ": " << pop_start.at(i) << endl; 
	}*/
	//die Zeit t
	double t = 0;
	
	//Maximale Anzahl an Interationen ( = wie oft sollen die Reaktionen maximal passieren?)
	double iter_max = 15; 
	
	//iteration counter
	double n = 0;
	//2 Beginn des Algorithmus
	
	//berechne die einzelnen Wahrscheinlichkeitsraten ob ein Dimer gebildet wird oder zerstört wird
	srand(time(NULL));
	  
	  
	  
	//gibt es einen Dimer-Formationsrate bzw. eine Dimerdestruktorrate
	while(n < iter_max)
	{
	
	vector <double> a; 
	a.push_back(cy.at(0)*pop_start.at(0)*pop_start.at(1));
	a.push_back(cy.at(1)*pop_start.at(2));
	
	//die Totale_Reaktionsrate bilden
	double a_total = 0.0;
	for (int i = 0; i<a.size();i++)
	{
		a_total += a.at(i); 
	}
	cout << "a_total : " << a_total <<endl;
	//
	
	
	//3 berechnen von tau und my 
	//bzw. welche Reaktion wird genommen
	double r0 = double(rand())/(RAND_MAX);
	  double r1 = double(rand())/(RAND_MAX);
	  double r2 = double(rand())/(RAND_MAX);//+(M-1);
	  cout << r2 << endl; 
	//Die Wahrscheinlichkeit, dass A+B --> AB wird
	//bzw. AB --> A + B 
	//3 tau und my berechnen --> welche Reaktion wird genommen ??? 
	  
	  double tau = (1/a_total)*(log(1/r1)); 
	double my; 
	
	//a_tot multipliziert mit der randomzahl
	
	double r_multipliziert_mit_a_tot = r2*a_total;
	cout << "r2*a_total: " << r_multipliziert_mit_a_tot << endl; 
	double summe; 
	
	
	  for (int j = 0; j < M ; ++j)
	  {
	  	
		if (partial_sum(a,j-1) < (r2*a_total) && (r2*a_total)<=partial_sum(a,j))
		{
				my = j; 
				cout << "Das ist jetzt die regel die genommen wird: " << j;
				break; 	
		}
		
	  }	
	if (my == 0)
	 {
	 	pop_start.at(2) = pop_start.at(2)+1;
		pop_start.at(1) = pop_start.at(1)-1;
		pop_start.at(0) = pop_start.at(0)-1;
	  } 
	  
	  if (my == 1)
	  {
	  	pop_start.at(2) = pop_start.at(2)-1;
		pop_start.at(1) = pop_start.at(1)+1;
		pop_start.at(0) = pop_start.at(0)+1;
	  }
	
	n = n+1; 
	cout << " ! " << endl; 
	for (int k = 0; k< pop_start.size(); k++)
	{
		cout << pop_start.at(k) << endl; 
	}
	t += tau;
	
	cout << " t:  " << t <<  endl; 
   cout << "!" << endl; 
   }
	return 0; 
}

