#pragma once

#include <iostream>
#include <cmath>
using namespace std;
export module Math;
export class Complex;
export class Rational;


namespace Math {
	class Complex {
	private:
		double m_mod;
		double m_arg;

	public:
		double getComplexArg(double re, double im) {
			const double PI = acos(-1.0);
			if (re > 0) {
				return atan(im / re);
			}
			else if (re < 0) {
				if (im >= 0) {
					return PI + atan(im / re);
				}
				else {
					return -PI + atan(im / re);
				}
			}
			else {
				if (im > 0) {
					return PI / 2;
				}
				else {
					return -PI / 2;
				}
			}
		}

		Complex()
			: m_mod(0.0), m_arg(0.0)
		{ }

		Complex(double re, double im) {
			m_mod = sqrt(re * re + im * im);
			m_arg = getComplexArg(re, im);
		}

		Complex(double mod)
			: m_mod(mod), m_arg(0.0)
		{ }

		static Complex FromExponentialForm(double mod, double arg) {
			double re = mod * cos(arg);
			double im = mod * sin(arg);
			Complex cur = Complex(re, im);
			return cur;
		}

		static Complex FromAlgebraicForm(double re, double im) {
			Complex cur = Complex(re, im);
			return cur;
		}

		double Re() {
			double re = m_mod * cos(m_arg);
			return re;
		}

		double Im() {
			double im = m_mod * sin(m_arg);
			return im;
		}

		double Mod() {
			return m_mod;
		}

		double Arg() {
			return m_arg;
		}

		explicit operator double() {
			return this->Re();
		}

		Complex operator- () {
			return Complex(-this->Re(), -this->Im());
		}

		Complex* operator++ () {
			double cur_re = this->Re() + 1;
			double im = this->Im();
			this->m_arg = sqrt(cur_re * cur_re + im * im);
			this->m_mod = this->getComplexArg(cur_re, im);
			return this;
		}

		Complex operator++ (int x) {
			double cur_re = this->Re() + 1;
			double im = this->Im();
			Complex cur = Complex(cur_re, im);
			return cur;
		}

		Complex* operator-- () {
			double cur_re = this->Re() - 1;
			double im = this->Im();
			this->m_arg = sqrt(cur_re * cur_re + im * im);
			this->m_mod = this->getComplexArg(cur_re, im);
			return this;
		}

		Complex operator-- (int x) {
			double cur_re = this->Re() - 1;
			double im = this->Im();
			Complex cur = Complex(cur_re, im);
			return cur;
		}

		Complex* operator+= (Complex x) {
			this->Re() + x.Re();
			this->Im() + x.Im();
			return this;
		}

		Complex* operator-= (Complex x) {
			this->Re() - x.Re();
			this->Im() - x.Im();
			return this;
		}

		Complex* operator*= (Complex x) {
			double cur_re = this->Re() * x.Re() - this->Im() * x.Im();
			double cur_im = this->Re() * x.Im() + this->Im() * x.Re();
			this->m_mod = sqrt(cur_re * cur_re + cur_im * cur_im);
			this->m_arg = this->getComplexArg(cur_re, cur_im);
			return this;
		}

		Complex* operator/= (Complex x) {
			double cur_re = (this->Re() * x.Re() + this->Im() * x.Im()) / (x.Re() * x.Re() + x.Im() * x.Im());
			double cur_im = (this->Im() * x.Re() - this->Re() * x.Im()) / (x.Re() * x.Re() + x.Im() * x.Im());
			this->m_mod = sqrt(cur_re * cur_re + cur_im * cur_im);
			this->m_arg = this->getComplexArg(cur_re, cur_im);
			return this;
		}

		friend Complex operator+(Complex lhs, Complex rhs);
		friend Complex operator-(Complex lhs, Complex rhs);
		friend Complex operator*(Complex lhs, Complex rhs);
		friend Complex operator/(Complex lhs, Complex rhs);
		friend Complex operator""i(long double x);
		friend Complex operator""i(unsigned long long x);
		friend ostream& operator<< (ostream& os, const Complex complex);
	};

	Complex operator+ (Complex lhs, Complex rhs) {
		double re = lhs.Re() + rhs.Re();
		double im = lhs.Im() + rhs.Im();
		return Complex(re, im);
	}

	Complex operator- (Complex lhs, Complex rhs) {
		double re = lhs.Re() - rhs.Re();
		double im = lhs.Im() - rhs.Im();
		return Complex(re, im);
	}

	Complex operator*(Complex lhs, Complex rhs) {
		double re = lhs.Re() * rhs.Re() - lhs.Im() * rhs.Im();
		double im = lhs.Re() * rhs.Im() + lhs.Im() * rhs.Re();
		return Complex(re, im);
	}

	Complex operator/ (Complex lhs, Complex rhs) {
		double re = (lhs.Re() * rhs.Re() + lhs.Im() * rhs.Im()) / (rhs.Re() * rhs.Re() + rhs.Im() * rhs.Im());
		double im = (lhs.Im() * rhs.Re() - lhs.Re() * rhs.Im()) / (rhs.Re() * rhs.Re() + rhs.Im() * rhs.Im());
		return Complex(re, im);
	}

	Complex operator""i(long double x) {
		return Complex(0, x);
	}

	Complex operator""i(unsigned long long x) {
		return Complex(0, x);
	}

	int FindGretestCommonDivisor(int a, int b) {
		if (a < b)
		{
			int temp = a;
			a = b;
			b = temp;
		}
		int r;
		while (b != 0)
		{
			r = a % b;
			a = b;
			b = r;
		}
		return a;
	}

	int FindLeastCommonMultiple(int x, int y) {
		return abs(x * y) / FindGretestCommonDivisor(x, y);
	}

	ostream& operator<< (ostream& os, Complex complex) {
		if (complex.Im() < 0) {
			return os << complex.Re() << complex.Im() << "i";
		}
		else
			return os << complex.Re() << "+" << complex.Im() << "i";
	}


	class Rational {
	private:
		int m_nominator;
		int m_denominator;

	public:
		void reduce() {
			if (m_nominator < 0 && m_denominator < 0) {
				m_nominator = -m_nominator;
				m_denominator = -m_denominator;
			}
			else if (m_nominator > 0 && m_denominator < 0) {
				m_nominator = -m_nominator;
				m_denominator = -m_denominator;
			}
			int nod = FindGretestCommonDivisor(m_nominator, m_denominator);
			m_nominator = m_nominator / nod;
			m_denominator = m_denominator / nod;
		}

		Rational() {
			m_nominator = 0;
			m_denominator = 1;
		}

		Rational(int x, int y) {
			m_nominator = x;
			m_denominator = y;
			this->reduce();
		}

		Rational(int x) {
			m_nominator = x;
			m_denominator = 1;
		}

		int Nominator() {
			return m_nominator;
		}

		int Denominator() {
			return m_denominator;
		}

		explicit operator double() {
			return m_nominator / m_denominator;
		}

		Rational operator- () {
			return Rational(-this->m_nominator, this->m_denominator);
		}

		Rational* operator++ () {
			this->m_nominator + 1;
			this->reduce();
			return this;
		}

		Rational operator++ (int x) {
			Rational cur = *this;
			m_nominator + 1; this->reduce();
			return cur;
		}

		Rational* operator-- () {
			this->m_nominator - 1;
			this->reduce();
			return this;
		}

		Rational operator-- (int x) {
			Rational cur = *this;
			m_nominator - 1;
			this->reduce();
			return cur;
		}

		Rational* operator+= (Rational x) {
			int no = m_nominator * x.m_denominator + m_denominator * x.m_nominator;
			int de = m_denominator * x.m_denominator;
			m_nominator = no;
			m_denominator = de;
			reduce();
			return this;
		}

		Rational* operator-= (Rational x) {
			int no = m_nominator * x.m_denominator - m_denominator * x.m_nominator;
			int de = m_denominator * x.m_denominator;
			m_nominator = no;
			m_denominator = de;
			reduce();
			return this;
		}

		Rational* operator*= (Rational x) {
			int no = m_nominator * x.m_nominator;
			int de = m_denominator * x.m_denominator;
			m_nominator = no;
			m_denominator = de;
			reduce();
			return this;
		}

		Rational* operator/= (Rational x) {
			int no = m_nominator * x.m_denominator;
			int de = m_denominator * x.m_nominator;
			m_nominator = no;
			m_denominator = de;
			reduce();
			return this;
		}

		friend Rational operator+ (Rational lhs, Rational rhs);
		friend Rational operator- (Rational lhs, Rational rhs);
		friend Rational operator* (Rational lhs, Rational rhs);
		friend Rational operator/ (Rational lhs, Rational rhs);
		friend bool operator== (Rational lhs, Rational rhs);
		friend ostream& operator<< (ostream& os, const Rational& rational);
	};

	Rational operator+ (Rational lhs, Rational rhs) {
		int no = lhs.m_nominator * rhs.m_denominator + lhs.m_denominator * rhs.m_nominator;
		int de = lhs.m_denominator * rhs.m_denominator;
		Rational now = Rational(no, de);
		now.reduce();
		return now;
	}

	Rational operator- (Rational lhs, Rational rhs) {
		int no = lhs.m_nominator * rhs.m_denominator - lhs.m_denominator * rhs.m_nominator;
		int de = lhs.m_denominator * rhs.m_denominator;
		Rational now = Rational(no, de);
		now.reduce();
		return now;
	}

	Rational operator* (Rational lhs, Rational rhs) {
		int no = lhs.m_nominator * rhs.m_nominator;
		int de = lhs.m_denominator * rhs.m_denominator;
		Rational now = Rational(no, de);
		now.reduce();
		return now;
	}

	Rational operator/ (Rational lhs, Rational rhs) {
		int no = lhs.m_nominator * rhs.m_denominator;
		int de = lhs.m_denominator * rhs.m_nominator;
		Rational now = Rational(no, de);
		now.reduce();
		return now;
	}

	bool operator== (Rational lhs, Rational rhs) {
		if (lhs.m_nominator == rhs.m_nominator && lhs.m_denominator == rhs.m_denominator) {
			return true;
		}
		else {
			return false;
		}
	}

	ostream& operator<< (ostream& os, const Rational& rational) {
		return os << rational.m_nominator << "/" << rational.m_denominator;
	}
}

