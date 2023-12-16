#include <iostream>
#include <complex>
#define _SILENCE_NONFLOATING_COMPLEX_DEPRECATION_WARNING
import Math;
using namespace std;
//using namespace Math;
int main()
{
	setlocale(LC_ALL, "ru");
	complex<double> aC(5, -5);
	complex<double> bC(5.6, -2.3);
	Complex a(5, -5);
	Complex b(5.6, -2.3);
	Rational num1(1, 2);
	Rational num2(3, 4);
	double c = -12;
	cout << "Мой класс complex : a = " << a << " Real = " << a.Re() << " Imag = " << a.Im() << " Mod = " << a.Mod() << " Arg = " << a.Arg() << endl;
	cout << "Встроенный класс complex: : a = " << aC << " Real = " << aC.real() << " Imag = " << aC.imag() << " Mod = " << abs(aC) << " Arg = " << arg(aC) << endl;
	cout << "Мой класс complex: a+b = " << a << " + " << b << " = " << a + b << endl;
	cout << "Встроенный класс complex: a+b = " << aC << " + " << bC << " = " << aC + bC << endl;
	cout << "Мой класс complex: a-b = " << a << " + " << b << " = " << a - b << endl;
	cout << "Встроенный класс complex: a-b = " << aC << " + " << bC << " = " << aC - bC << endl;
	cout << "Мой класс complex: a*b = " << a << " + " << b << " = " << a * b << endl;
	cout << "Встроенный класс complex: a*b = " << aC << " + " << bC << " = " << aC * bC << endl;
	cout << "Мой класс complex: a/b = " << a << " + " << b << " = " << a / b << endl;
	cout << "Встроенный класс complex: a/b = " << aC << " + " << bC << " = " << aC / bC << endl;
	cout << "Мой класс complex: a+c = " << a << " + " << c << " = " << a + c << endl;
	cout << "Встроенный класс complex: a+b = " << aC << " + " << c << " = " << aC + c << endl;
	//cout << "Мой класс Rational: a+b = " << num1 << " + " << num2 << " = " << num1 + num2 << endl;
	//cout << "Мой класс Rational: a-b = " << num1 << " - " << num2 << " = " << num1 - num2 << endl;
	//cout << "Мой класс Rational: a*b = " << num1 << " * " << num2 << " = " << num1 * num2 << endl;
	//cout << "Мой класс Rational: a/b = " << num1 << " / " << num2 << " = " << num1 / num2 << endl;
	//cout << "Мой класс Rational: a = " << num1 << " числитель = " << num1.Nominator() << " знаменатель = " << num1.Denominator() << " double = " << double(num1) << endl;
	/*cout << "f(" << a << ") = " << f(a) << endl;
	cout << "f(" << num2 << ") = " << f(num2) << endl;
	cout << "f(" << c << ") = " << f(c) << endl;*/
	return 0;
}