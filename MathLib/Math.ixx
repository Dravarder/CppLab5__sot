#include <iostream>
#include <math.h>


export module Math;
export namespace Math
//123

{
	export class Complex
	{
	private:

		double m_mod;
		double m_arg;

	public:
		Complex()
		{
			m_mod = 0;
			m_arg = 0;
		}
		Complex(double real, double image)
		{
			m_mod = sqrt(real * real + image * image);
			m_arg = atan2(image, real);
			if (m_mod < 6.6613381477509392e-15) {
				m_mod = 0;
			}
			if ((m_arg > 3.1415926535897) && (m_arg < 3.1415926535898)) {
				m_arg = 0;
			}
		}
		Complex(double a)
		{
			m_mod = a;
			m_arg = 0;
		}
		static Complex FromExponentialForm(double mod, double arg)
		{

			return Complex(mod * cos(arg), mod * sin(arg));
		}
		static Complex FromAlgebraicForm(double real, double img)
		{
			Complex a;
			a.m_mod = sqrt(real * real + img * img);
			a.m_arg = atan2(img, real);
			return a;
		}
		double Re() const
		{
			return m_mod * cos(m_arg);
		}
		double Im() const
		{
			return m_mod * sin(m_arg);
		}
		double Mod() const
		{
			return m_mod;
		}
		double Arg() const
		{
			return m_arg;
		}
		explicit operator double() const {
			return this->Re();
		}
		Complex operator-() const
		{
			return Complex(-1 * this->Re(), -1 * this->Im());
		}
		Complex& operator++()
		{
			double real = this->Re() + 1;
			double img = this->Im();
			this->m_mod = sqrt(real * real + img * img);
			this->m_arg = atan2(img, real);
			return *this;
		}
		Complex operator++ (int) {
			Complex copy(*this);
			double real = this->Re() + 1;
			double img = this->Im();
			this->m_mod = sqrt(real * real + img * img);
			this->m_arg = atan2(img, real);
			return copy;
		}
		Complex& operator--()
		{
			double real = this->Re() - 1;
			double img = this->Im();
			this->m_mod = sqrt(real * real + img * img);
			this->m_arg = atan2(img, real);
			return *this;
		}
		Complex operator-- (int) {
			Complex copy(*this);
			double real = this->Re() - 1;
			double img = this->Im();
			this->m_mod = sqrt(real * real + img * img);
			this->m_arg = atan2(img, real);
			return copy;
		}

		Complex& operator+=(Complex a)
		{
			double real = Re();
			double img = Im();
			real += a.Re();
			img += a.Im();
			m_mod = sqrt(real * real + img * img);
			m_arg = atan2(img, real);
			return *this;
		}
		Complex& operator-=(Complex a)
		{
			double real = Re();
			double img = Im();
			real -= a.Re();
			img -= a.Im();
			m_mod = sqrt(real * real + img * img);
			m_arg = atan2(img, real);
			return *this;
		}
		Complex& operator*=(Complex a)
		{
			m_mod *= a.m_mod;
			m_arg += a.m_arg;
			return *this;
		}
		Complex& operator/=(Complex a)
		{
			m_mod /= a.m_mod;
			m_arg -= a.m_arg;
			return *this;
		}
		friend Complex operator+ (const Complex& left, const Complex& right);
		friend Complex operator- (const Complex& left, const Complex& right);
		friend Complex operator* (const Complex& left, const Complex& right);
		friend Complex operator/ (const Complex& left, const Complex& right);
		friend Complex operator ""i(long double im);
		friend Complex operator ""i(unsigned long long im);
		friend std::ostream& operator<<(std::ostream& stream, const Complex& a);
		friend Complex Poww(const Complex& r, double B);
	};

	export Complex operator+(const Complex& left, const Complex& right)
	{
		double real = left.Re() + right.Re();
		double img = left.Im() + right.Im();
		return Complex(real, img);
	}

	export Complex operator-(const Complex& left, const Complex& right)
	{
		double real = left.Re() - right.Re();
		double img = left.Im() - right.Im();
		return Complex(real, img);
	}

	export Complex operator*(const Complex& left, const Complex& right)
	{
		double real = left.Mod() * right.Mod();
		double img = left.Arg() + right.Arg();
		double real1 = real * cos(img);
		double img1 = real * sin(img);
		return Complex(real1, img1);
	}

	export Complex operator/(const Complex& left, const Complex& right)
	{
		double real = left.Mod() / right.Mod();
		double img = left.Arg() - right.Arg();
		double real1 = real * cos(img);
		double img1 = real * sin(img);
		return Complex(real1, img1);
	}

	export Complex operator ""i(unsigned long long im) {
		return Complex(0.0, static_cast<double>(im));
	}

	export Complex operator ""i(long double im) {
		return Complex(0.0, static_cast<double>(im));
	}
	export std::ostream& operator<<(std::ostream& stream, const Complex& a) {
		if (a.Im() < 0) {
			stream << a.Re() << " " << a.Im() << "i";
		}
		else {
			stream << a.Re() << " + " << a.Im() << "i";
		}
		return stream;
	}

	export int FindGreatestCommonDivisor(int a, int b) {
		if (a < 0) {
			a *= -1;
		}
		if (b < 0) {
			b *= -1;
		}
		while (a != 0 && b != 0) {
			if (a > b) {
				a %= b;
			}
			else {
				b %= a;
			}
		}
		return a + b;
	}

	export int FindLeastCommonMultiple(int a, int b) {
		return abs(a * b) / FindGreatestCommonDivisor(a, b);
	}

	export Complex Sin(const Complex& z) {
		Complex a = 2 * z;
		return Complex(sin((double)a));
	}

	export class Rational {
	private:
		int m_nominator;
		int m_denominator;
		void normalizator() {
			if (m_denominator < 0) {
				m_denominator *= -1;
				m_nominator *= -1;
			}
			int nod = FindGreatestCommonDivisor(m_nominator, m_denominator);
			m_nominator /= nod;
			m_denominator /= nod;
		}
	public:
		Rational() {
			m_nominator = 0;
			m_denominator = 1;
		}

		Rational(int x, int y) {
			m_nominator = x;
			m_denominator = y;
			this->normalizator();
		}

		Rational(int x) {
			m_nominator = x;
			m_denominator = 1;
		}
		int Nominator() const {
			return m_nominator;
		}
		int Denominator() const {
			return m_denominator;
		}
		explicit operator double() const {
			return static_cast<double>(m_nominator) / m_denominator;
		}
		Rational operator-() {
			Rational a(*this);
			a.m_nominator *= -1;
			return a;
		}
		Rational& operator++ () {
			Rational a(*this);
			m_nominator += a.m_denominator;
			normalizator();
			return *this;
		}
		Rational operator++ (int) {
			Rational a(*this);
			m_nominator += a.m_denominator;
			normalizator();
			return a;
		}
		Rational& operator-- () {
			Rational a(*this);
			m_nominator -= a.m_denominator;
			normalizator();
			return *this;
		}
		Rational operator-- (int) {
			Rational a(*this);
			m_nominator -= a.m_denominator;
			normalizator();
			return a;
		}
		Rational& operator+=(Rational a) {
			int Nominator = m_nominator * a.m_denominator + m_denominator * a.m_nominator;
			int Denominator = m_denominator * a.m_denominator;
			m_nominator = Nominator;
			m_denominator = Denominator;
			normalizator();
			return *this;
		}
		Rational& operator-=(Rational a) {
			int Nominator = m_nominator * a.m_denominator - m_denominator * a.m_nominator;
			int Denominator = m_denominator * a.m_denominator;
			m_nominator = Nominator;
			m_denominator = Denominator;
			normalizator();
			return *this;
		}
		Rational& operator*=(Rational a) {
			m_denominator *= a.m_denominator;
			m_nominator *= a.m_nominator;
			normalizator();
			return *this;
		}
		Rational& operator/=(Rational a) {
			m_denominator *= a.m_nominator;
			m_nominator *= a.m_denominator;
			normalizator();
			return *this;
		}
		friend Rational operator+ (const Rational& left, const Rational& right);
		friend Rational operator- (const Rational& left, const Rational& right);
		friend Rational operator* (const Rational& left, const Rational& right);
		friend Rational operator/(const Rational& left, const Rational& right);
		friend bool operator==(const Rational& left, const Rational& right);
		friend bool operator>(const Rational& left, const Rational& right);
		friend bool operator<(const Rational& left, const Rational& right);
		friend bool operator>=(const Rational& left, const Rational& right);
		friend bool operator<=(const Rational& left, const Rational& right);
		friend std::ostream& operator<< (std::ostream& stream, const Rational& rational);
		friend Rational Poww(const Rational& r, double B);
	};

	export Rational operator+ (const Rational& left, const Rational& right) {
		int denominator = FindLeastCommonMultiple(left.m_denominator, right.m_denominator);
		int nominator = denominator / left.m_denominator * left.m_nominator;
		nominator += denominator / right.m_denominator * right.m_nominator;
		return Rational(nominator, denominator);
	}

	export Rational operator-(const Rational& left, const Rational& right)
	{
		int denominator = FindLeastCommonMultiple(left.m_denominator, right.m_denominator);
		int nominator = denominator / left.m_denominator * left.m_nominator;
		nominator -= denominator / right.m_denominator * right.m_nominator;
		return Rational{ nominator, denominator };
	}

	export Rational operator*(const Rational& left, const Rational& right)
	{
		return Rational(left.m_nominator * right.m_nominator, right.m_denominator * left.m_denominator);
	}

	export Rational operator/(const Rational& left, const Rational& right)
	{
		return Rational(left.m_nominator * right.m_denominator, left.m_denominator * right.m_nominator);
	}

	export bool operator==(const Rational& left, const Rational& right)
	{
		return left.m_nominator == right.m_nominator && left.m_denominator == right.m_denominator;
	}

	export bool operator>(const Rational& left, const Rational& right)
	{
		int lcm = FindLeastCommonMultiple(left.m_denominator, right.m_denominator);
		return lcm / left.m_denominator * left.m_nominator > lcm / right.m_denominator * right.m_nominator;
	}

	export bool operator<(const Rational& left, const Rational& right)
	{
		int lcm = FindLeastCommonMultiple(left.m_denominator, right.m_denominator);
		return lcm / left.m_denominator * left.m_nominator < lcm / right.m_denominator * right.m_nominator;
	}

	export bool operator>=(const Rational& left, const Rational& right)
	{
		int lcm = FindLeastCommonMultiple(left.m_denominator, right.m_denominator);
		return lcm / left.m_denominator * left.m_nominator >= lcm / right.m_denominator * right.m_nominator;
	}

	export bool operator<=(const Rational& left, const Rational& right)
	{
		int lcm = FindLeastCommonMultiple(left.m_denominator, right.m_denominator);
		return lcm / left.m_denominator * left.m_nominator <= lcm / right.m_denominator * right.m_nominator;
	}

	export std::ostream& operator<<(std::ostream& stream, const Rational& rational) {
		stream << rational.m_nominator << "/" << rational.m_denominator;
		return stream;
	}

	export Rational Poww(const Rational& r, double B) {
		Rational res = r;
		if (B == 0) {
			return 1;
		}
		else {
			if (B < 0) {
				B = B * (-1);
				while (B - 1 > 0) {
					res = res * r;
					B--;
				}
				res = 1 / res;
			}
			else {
				while (B - 1 > 0) {
					res = res * r;
					B--;
				}
			}
		}
		return res;
	}

	export Complex Poww(const Complex& r, double B) {
		Complex res = r;

		if (B == 0) {
			return 1;
		}
		else {
			if (B < 0) {
				B = B * (-1);
				while (B - 1 > 0) {
					res = res * r;
					B--;
				}
				res = 1 / res;
			}
			else {
				while (B - 1 > 0) {
					res = res * r;
					B--;
				}
			}
		}
		return res;
	}


	export Complex f(const Complex& c) {
		Complex a = Complex(0.0, 1.0);
		Complex result = Poww(c, 3) + (1 + 2 * a) * Poww(c, 2) + (1 - 2 * a) * Poww(c, -5 * (double)a);
		return result;
	}

	export Rational f(const Rational& r) {
		Rational a(1, 2);
		Rational result = Poww(r, 3) + (1 + 2 * a) * Poww(r, 2) + (1 - 2 * a) * Poww(r, -5 * (double)a);
		return result;
	}

	export double f(const double& b) {
		double a = 0.5;
		double result = pow(b, 3) + (1 + 2 * a) * pow(b, 2) + (1 - 2 * a) * pow(b, -5 * a);
		return result;
	}
}

