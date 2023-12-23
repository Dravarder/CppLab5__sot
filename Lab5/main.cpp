#include <iostream>
#include <complex>
#define _SILENCE_NONFLOATING_COMPLEX_DEPRECATION_WARNING

import Math;
using namespace std;

int main()
{
	cout << "Vvedite deistvitelnuy chast complexnogo chisla:";
	double re;
	cin >> re;
	cout << "Vvedite mnimuy chast complexnogo chisla:";
	double im;
	cin >> im;
	Math::Complex z(re, im);
	cout << "Resultat funkcii na complexnih chislah:" << Math::f(z) << std::endl;
	cout << "Vvedite chislitel racionalnogo chisla:";
	double nom;
	cin >> nom;
	cout << "Vvedite znamenatel racionalnogo chisla:";
	double denom;
	cin >> denom;
	Math::Rational r(nom, denom);
	cout << "Resultat funkcii na racionalnih chislah:" << Math::f(r) << endl;
	cout << "Enter the number: ";
	double u;
	cin >> u;
	cout << "The result of function:" << Math::f(u) << endl;
	return 0;
}