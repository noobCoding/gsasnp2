#pragma  once
//Polynomial Fit
#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<map>
using namespace std;

#define	maxN	50000
#define STEPS	12

typedef vector<double> vec;

struct myPoint {
	double x;
	double y;
};

struct SplineSet {
	double a;
	double b;
	double c;
	double d;
	double x;
};

class vcRegression
{
private:
	void estimateTangents(vector<double>& x, vector<double>& y, vector<double>& m)
	{
		int n = x.size();
		vector<double> delta(n);
		//slopes;
		for (int k = 0; k < n - 1; k++) {
			delta[k] = (y[k + 1] - y[k]) / (x[k + 1] - x[k]);
		}
		//average tangent at data point
		m[0] = delta[0];
		m[n - 1] = delta[n - 2];
		for (int k = 1; k < n - 1; k++) {
			m[k] = (delta[k - 1] + delta[k]) / 2; //error corrected.
		}
		for (int k = 0; k < n - 1; k++) {
			if (delta[k] == 0.)
				m[k] = m[k + 1] = 0.;
			else {
				double alpha = m[k] / delta[k];
				double beta = m[k + 1] / delta[k];
				double dist = alpha * alpha + beta * beta;
				if (dist > 9.) {
					double tau = 3. * delta[k] / sqrt(dist);
					m[k] = tau * alpha;
					m[k + 1] = tau * beta;
				}
			}
		}
	}

	// monotonic cubic splines variables
	vector<double> xx;
	vector<double> yy;
	vector<double> m;

public:
	vcRegression() {};
	~vcRegression() {};

	void monotoneCubicSplineInit(vector<myPoint>& fittingpoints)
	{
		int n = (int)fittingpoints.size();
		if (n < 2) return;

		//if (!isMonotonicPolyline(fittingpoints)) return;
		map<double, double> monotonicPoints;
		for (int i = 0; i < n; i++) {
			monotonicPoints[fittingpoints[i].x] = fittingpoints[i].y;
		}

		xx.resize(n), yy.resize(n);

		int k = 0;
		for (auto mi : monotonicPoints) {
			xx[k] = mi.first;
			yy[k++] = mi.second;
		}

		//tangent table		
		m.resize(n);
		for (auto vi : m) vi = 0;
		estimateTangents(xx, yy, m);
	}

	double monotoneCubicSplineInterp(double x)
	{
		int n = (int)xx.size();
		myPoint tmpPoint;

		//Hermite spline	
		for (int i = 0; i < n - 1; i++)
			if ((double)xx[i] < x && x <= (double)xx[i + 1])
			{
				double h = xx[i + 1] - xx[i];
				//double t = (x - xx[i]) / (xx[i+1] - xx[i]); // 0 < t <= 1
				double t = (x - 0.5 - xx[i]) / h;
				double t2 = t*t;
				double t3 = t2*t;
				double h00 = 2.*t3 - 3.*t2 + 1.;
				double h10 = t3 - 2.*t2 + t;
				double h01 = -2.*t3 + 3.*t2;
				double h11 = t3 - t2;
				double yf = h00*yy[i] + h10*h*m[i] + h01*yy[i + 1] + h11*h*m[i + 1];
				
				//tmpPoint.x = xx[i] + t*h + 0.5;
				tmpPoint.x = x;
				tmpPoint.y = yf;
				break;
			}

		//for (int i = 0; i < n - 1; i++)
		//	//if ((double)xx[i] < x && x <= (double)xx[i + 1])
		//{
		//	double h = xx[i + 1] - xx[i];
		//	//double t = (x - xx[i]) / (xx[i+1] - xx[i]); // 0 < t <= 1
		//	//double t = (x - 0.5 - xx[i]) / h;
		//	for (int k = 0; k < STEPS; k++) {
		//		double t = (double)k / STEPS;
		//		double t2 = t*t;
		//		double t3 = t2*t;
		//		double h00 = 2.*t3 - 3.*t2 + 1.;
		//		double h10 = t3 - 2.*t2 + t;
		//		double h01 = -2.*t3 + 3.*t2;
		//		double h11 = t3 - t2;
		//		double yf = h00*yy[i] + h10*h*m[i] + h01*yy[i + 1] + h11*h*m[i + 1];
		//		tmpPoint.x = xx[i] + t*h + 0.5;
		//		//tmpPoint.x = x;
		//		tmpPoint.y = yf + 0.5;
		//		curve.push_back(tmpPoint);
		//	}
		//}

		return tmpPoint.y;
	}

	vector<SplineSet> spline(vec &x, vec &y)
	{
		map<double, double> xx; // <x-value, sum(y) value>
		map<double, int> xs; // <x-value, no of y>
		int n = x.size() - 1;

		for (int i = 0; i < n + 1; i++) {
			xx[x[i]] += y[i];
			xs[x[i]] ++;
		}

		x.clear(); x.resize(0);
		y.clear(); y.resize(0);
		double maxY = 0.;

		for (map<double, double>::iterator mi = xx.begin(); mi != xx.end(); mi++) {
			x.push_back(mi->first);
			y.push_back(mi->second / xs[mi->first]); // average of all y(s) in the same pin x
			if (maxY < mi->second / xs[mi->first]) maxY = mi->second / xs[mi->first];	// best y value	
		}

		x.push_back(maxN);
		y.push_back((maxY + y.back()) / 2); // last boundary = 1/2 * (max y value + last y value)

											// original
		n = (int)x.size() - 1;
		vec a;
		a.insert(a.begin(), y.begin(), y.end());
		vec b(n);
		vec d(n);
		vec h;

		for (int i = 0; i < n; ++i)
			h.push_back(x[i + 1] - x[i]);

		vec alpha;
		for (int i = 0; i < n; ++i)
			alpha.push_back(3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1]);

		vec c(n + 1);
		vec l(n + 1);
		vec mu(n + 1);
		vec z(n + 1);
		l[0] = 1;
		mu[0] = 0;
		z[0] = 0;

		for (int i = 1; i < n; ++i)
		{
			l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
			mu[i] = h[i] / l[i];
			z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
		}

		l[n] = 1;
		z[n] = 0;
		c[n] = 0;

		for (int j = n - 1; j >= 0; --j)
		{
			c[j] = z[j] - mu[j] * c[j + 1];
			b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
			d[j] = (c[j + 1] - c[j]) / 3 / h[j];
		}

		vector<SplineSet> output_set(n);
		for (int i = 0; i < n; ++i)
		{
			output_set[i].a = a[i];
			output_set[i].b = b[i];
			output_set[i].c = c[i];
			output_set[i].d = d[i];
			output_set[i].x = x[i];
		}
		return output_set;
	}

	double Nspline(vector<SplineSet> nsp, double p) {
		int n = (int)nsp.size() - 1;
		for (int i = 0; i < n; i++)
			if (nsp[i].x <= p && p <= nsp[i + 1].x)
				return (nsp[i].d*pow((p - nsp[i].x), 3) + nsp[i].c*pow((p - nsp[i].x), 2) + nsp[i].b*(p - nsp[i].x) + nsp[i].a);
		return 0.;
	}

	//void Bspline(vector<double> xx, vector<double> yy, vector<double> &esCoef) { // cubic/Bezier curve
	//	char choice = 'y';
	//	int n, i, j;
	//	float h[10], a, b, c, d, sum, s[10] = { 0 }, x[10], F[10], f[10], p, m[10][10] = { 0 }, temp;
	//	cout << "No of samples ? "; cin >> n;
	//	cout << "enter all sample points: " << endl;
	//	for (i = 0; i < n; i++)
	//		cin >> x[i] >> f[i];
	//	for (i = n - 1; i > 0; i--)
	//	{
	//		F[i] = (f[i] - f[i - 1]) / (x[i] - x[i - 1]);
	//		h[i - 1] = x[i] - x[i - 1];
	//	}
	//	//*********** formation of h, s , f matrix **************//
	//	for (i = 1; i < n - 1; i++)
	//	{
	//		m[i][i] = 2 * (h[i - 1] + h[i]);
	//		if (i != 1)
	//		{
	//			m[i][i - 1] = h[i - 1];
	//			m[i - 1][i] = h[i - 1];
	//		}
	//		m[i][n - 1] = 6 * (F[i + 1] - F[i]);
	//	}
	//	//***********  forward elimination **************//
	//	for (i = 1; i < n - 2; i++)
	//	{
	//		temp = (m[i + 1][i] / m[i][i]);
	//		for (j = 1; j <= n - 1; j++)
	//			m[i + 1][j] -= temp*m[i][j];
	//	}
	//	//*********** back ward substitution *********//
	//	for (i = n - 2; i > 0; i--)
	//	{
	//		sum = 0;
	//		for (j = i; j <= n - 2; j++)
	//			sum += m[i][j] * s[j];
	//		s[i] = (m[i][n - 1] - sum) / m[i][i];
	//	}
	//	while (choice == 'y')
	//	{
	//		cout << endl << "Enter x  : "; cin >> p;
	//		for (i = 0; i < n - 1; i++)
	//			if (x[i] <= p&&p <= x[i + 1])
	//			{
	//				a = (s[i + 1] - s[i]) / (6 * h[i]);
	//				b = s[i] / 2;
	//				c = (f[i + 1] - f[i]) / h[i] - (2 * h[i] * s[i] + s[i + 1] * h[i]) / 6;
	//				d = f[i];
	//				sum = a*pow((p - x[i]), 3) + b*pow((p - x[i]), 2) + c*(p - x[i]) + d;
	//			}
	//		cout << "coefficients of sub interval : " << endl << a << endl << b << endl << c << endl << d;
	//		cout << endl << "Functional value is: " << endl << sum;
	//		cout << "wanna continue (y/n) ? "; cin >> choice;
	//	}
	//	cin.get();
	//}

	void Logarithmic(vector<double> x, vector<double> y, vector<double> &esCoef) {
		int Ndata = (int)x.size();
		double* newX = new double[Ndata];
		double Sxx = 0., Sxy = 0., Syy = 0.;
		double meanNewX = 0., meanY = 0.;

		for (int i = 0; i < Ndata; i++) {
			newX[i] = log(x[i]);
			meanNewX += newX[i];
		}
		meanNewX /= Ndata;

		for (int i = 0; i < Ndata; i++) {
			meanY += y[i];
		}
		meanY /= Ndata;

		for (int i = 0; i < Ndata; i++) {
			Sxx += (newX[i] - meanNewX)*(newX[i] - meanNewX);
		}

		for (int i = 0; i < Ndata; i++) {
			Sxy += (newX[i] - meanNewX)*(y[i] - meanY);
		}

		for (int i = 0; i < Ndata; i++) {
			Syy += (y[i] - meanY)*(y[i] - meanY);
		}

		double b1 = Sxy / Sxx;
		double b0 = meanY - b1 * meanNewX;
		esCoef.push_back(b1);
		esCoef.push_back(b0);

		delete[] newX, newX = 0;
	}

	void Polynomial(vector<double> x, vector<double> y, int nDegree, vector<double> &esCoef) { // n is the degree of Polynomial 
		int Ndata = (int)x.size();
		int n = nDegree;
		double* X = new double[2 * n + 1];				//Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
		for (int i = 0; i < 2 * n + 1; i++)
		{
			X[i] = 0;
			for (int j = 0; j < Ndata; j++)
				X[i] += pow(x[j], i);					//consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
		}

		double* a = new double[n + 1];					//B is the Normal matrix(augmented) that will store the equations, 'a' is for value of the final coefficients
		double** B = new double*[n + 1];
		for (int i = 0; i < n + 1; i++) {
			B[i] = new double[n + 2];
		}

		for (int i = 0; i <= n; i++)
			for (int j = 0; j <= n; j++)
				B[i][j] = X[i + j];						//Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix

		double* Y = new double[n + 1];					//Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
		for (int i = 0; i < n + 1; i++)
		{
			Y[i] = 0;
			for (int j = 0; j < Ndata; j++)
				Y[i] += pow(x[j], i)*y[j];				//consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
		}

		for (int i = 0; i <= n; i++)
			B[i][n + 1] = Y[i];							//load the values of Y as the last column of B(Normal Matrix but augmented)

		n = n + 1;										//n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations

		for (int i = 0; i < n; i++)                    //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
			for (int k = i + 1; k < n; k++)
				if (B[i][i] < B[k][i])
					for (int j = 0; j <= n; j++)
					{
						double temp = B[i][j];
						B[i][j] = B[k][j];
						B[k][j] = temp;
					}

		for (int i = 0; i < n - 1; i++)					//loop to perform the Gauss elimination
			for (int k = i + 1; k < n; k++)
			{
				double t = B[k][i] / B[i][i];
				for (int j = 0; j <= n; j++)
					B[k][j] = B[k][j] - t*B[i][j];		//make the elements below the pivot elements equal to zero or eliminate the variables
			}

		for (int i = n - 1; i >= 0; i--)                //back-substitution
		{												//x is an array whose values correspond to the values of x,y,z..
			a[i] = B[i][n];								//make the variable to be calculated equal to the rhs of the last equation
			for (int j = 0; j < n; j++)
				if (j != i)								//then subtract all the lhs values except the coefficient of the variable whose value  is being calculated
					a[i] = a[i] - B[i][j] * a[j];
			a[i] = a[i] / B[i][i];						//now finally divide the rhs by the coefficient of the variable to be calculated
			esCoef.push_back(a[i]);
		}

		//cout << "\nThe values of the coefficients are as follows:\n";
		//for (int i = 0; i < n; i++)
		//	cout << "x^" << i << "=" << a[i] << endl;            // Print the values of x^0,x^1,x^2,x^3,....    
		//cout << "\nHence the fitted Polynomial is given by:\ny=";
		//for (int i = 0; i < n; i++)
		//	cout << " + (" << a[i] << ")" << "x^" << i;
		//cout << "\n";

		delete[] a, a = 0;
		delete[] X, X = 0;
		delete[] Y, Y = 0;
		for (int i = 0; i < nDegree + 1; i++) {
			delete[] B[i];
		}
		delete[] B, B = 0;
	}
};
