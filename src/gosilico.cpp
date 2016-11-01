#include <iostream>
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>

using namespace std;

//#define DEBUG
//#define TEST
#ifdef TEST
#include <armadillo> // for testing
#endif

#define EQ_NUM 2
#define EPSILON 0.000000000001

typedef array<array<double, EQ_NUM>, EQ_NUM> Mat;
typedef array<double, EQ_NUM> Vec;

// lVec <= lVec + rVec
void operator+=(Vec& lVec, Vec& rVec) {
	for (size_t i = 0; i < lVec.size(); i++) {
		lVec[i] += rVec[i];
	}
}

// is a almost equal to b?
bool inline equal(double &a, double &b) {
	return a == b || fabs(a - b) < EPSILON;
}

// true if some element is not equal to val
bool operator!=(Vec& lVec, double val) {
	for (size_t i = 0; i < lVec.size(); i++) {
		if (lVec[i] != val)
			return true;
	}
	return false;
}

const double Lambda = 0.5, k_eq1 = 0.14, k_eq2 = 2.7, nu1 = 6.6, nu2 = 6.2;
const double delta1 = 38, delta2 = 58, cSalt = 0.15;

#define COMMON_EXPR(nu) pow((Lambda - (nu1 + delta1)*q1 - (nu2 + delta2)*q2) / cSalt, nu)
// the equations
#define EQ1 k_eq1 * COMMON_EXPR(nu1) * c_eq1 - q1
#define EQ2 k_eq2 * COMMON_EXPR(nu2) * c_eq2 - q2

// partial derivatives d(eq_i)/d(q_i)
#define dCOMMON_EXPR_dq1(nu) k_eq1 * nu1 * c_eq1 * (-nu1-delta1) * COMMON_EXPR(nu-1) / cSalt
#define dCOMMON_EXPR_dq2(nu) k_eq2 * nu2 * c_eq2 * (-nu2-delta2) * COMMON_EXPR(nu-1) / cSalt

#define dEQ1_dq1 dCOMMON_EXPR_dq1(nu1) - 1
#define dEQ1_dq2 dCOMMON_EXPR_dq1(nu1)
#define dEQ2_dq1 dCOMMON_EXPR_dq2(nu2)
#define dEQ2_dq2 dCOMMON_EXPR_dq2(nu2) - 1

#define PARALLEL_LINES 1
#define SAME_LINES 2
/*
 * solve a 2x2 linear equation: AX=B
 * A00.X0 + A01.X1 = B0
 * A10.X0 + A11.X1 = B1
 * returns non-zero when no solution
 */
int solve(Mat A, Vec B, Vec &X) {
	if (A[1][0]) {
		double ratio = -A[0][0] / A[1][0];
		A[0][0] = 0;
		A[0][1] += ratio * A[1][1];
		B[0] += ratio * B[1];
		if (!B[0] && !A[0][1])
			return SAME_LINES;
		else if (!A[0][1])
			return PARALLEL_LINES;
		X[1] = B[0] / A[0][1];
		if (A[1][0]) {
			X[0] = (B[1] - A[1][1] * X[1]) / A[1][0];
		}
	} else {
		if (!A[1][1]) { // eq2 is bogus
			return -1;
		}
		X[1] = B[1] / A[1][1];
		if (!A[0][0]) {
			if (A[0][1] * X[1] == B[0])
				return SAME_LINES; // i.e. same points
			else
				return PARALLEL_LINES; // i.e. two distinct points
		}
		X[0] = (B[0] - A[0][1] * X[1]) / A[0][0];
	}
	return 0;
}

#ifdef TEST
// test the simple solver against a well known library
void testMySolver() {
	srand(time(NULL));
	int rounds = 1000;

	while (rounds--) {
		arma::mat A = arma::randu<arma::mat>(EQ_NUM, EQ_NUM);
		arma::vec b = arma::randu<arma::vec>(EQ_NUM);
		arma::vec x;
		bool found = true;
		try {
			x = arma::solve(A, b, arma::solve_opts::no_approx);
		} catch (runtime_error e) {
			found = false;
		}

		Vec _X, _B = {b[0], b[1]};
		Mat _A;
		_A[0]= {A(0,0),A(0,1)};
		_A[1]= {A(1,0),A(1,1)};
		int err = solve(_A, _B, _X);
		if (!found && !err) {
			cout << "wrong result";
			return;
		}
		if (found && !(equal(x[0], _X[0]) && equal(x[1], _X[1]))) {
			cout.precision(10);
			cout << "A=\n" << A(0, 0) << " " << A(0, 1) << " =" << b[0] << '\n'
			<< A(1, 0) << " " << A(1, 1) << " =" << b[1] << '\n';

			cout << "x[0]=" << x[0] << " x[1]=" << x[1] << '\n';
			cout << "_X[0]=" << _X[0] << " _X[1]=" << _X[1] << '\n';
			return;
		}
	}
	cout << "test passed!";
}
#endif

// L1 distance to zero vector
double norm(Vec X) {
	double s = 0;
	for (size_t i = 0; i < X.size(); i++) {
		s = max(s, abs(X[i]));
	}
	return s;
}

#define RANDX 8*((rand()) % 100) / (double) 100000

bool solveEQ(double c_eq1, double c_eq2, Vec &Q) {
	Mat J; // Jacobian matrix
	Vec X = { 0, 0 }, F, Y;

	int limit = 500; // max iterations
	double q1, q2;
	do {
		q1 = X[0];
		q2 = X[1];
		// update F(X) and J(X)
		F = {-(EQ1), -(EQ2)};
		J[0]= {dEQ1_dq1, dEQ1_dq2};
		J[1]= {dEQ2_dq1, dEQ2_dq2};

#ifdef DEBUG
		cout << "X0=" << X[0] << " X1=" << X[1] << '\n';
		cout << "F0=" << F[0] << " F1=" << F[1] << '\n';
		cout << "J00=" << J[0][0] << " J01=" << J[0][1] << " J10=" << J[1][0]
		<< " J11=" << J[1][1] << '\n';
		cout << "Y0=" << Y[0] << " Y1=" << Y[1] << " ||Y||=" << norm(Y)
		<< " ||F||=" << norm(F) << '\n';
		cout << "----------------" << '\n';
#endif

		// calculate the equation system JY=-F.
		if (solve(J, F, Y) != 0 || isinf(F[0]) || isinf(F[1])) {
			X[0] = RANDX;
			X[1] = RANDX;
			//cout << "reset X0=" << X[0] << " X1=" << X[1] << '\n';
			continue;
		}
		// next step
		X += Y;

	} while (norm(F) > EPSILON && limit--);

	if (limit < 0) {
		return false;
	}
	Q[0] = X[0];
	Q[1] = X[1];
	return true;
}

int main() {
	srand(time(NULL));
	cout.precision(10);
	double c_eq1 = 0.1, c_eq2 = 0.02;
	Vec Q;
	std::ofstream file1("q1.txt"), file2("q2.txt");
	clock_t begin = clock();

	double max = 10, step1 = 0.0001, step2 = 0.00001;
	int count = 0;
	for (double i = 0; i <= max; i += 0.1) {
		for (double j = 0; j <= max; j += 0.1) {
			c_eq1 = i * step1;
			c_eq2 = j * step2;

			if (solveEQ(c_eq1, c_eq2, Q)) {
				count++;
				file1 << c_eq1 << ' ' << c_eq2 << ' ' << Q[0] << '\n';
				file2 << c_eq1 << ' ' << c_eq2 << ' ' << Q[1] << '\n';
			}
		}
	}

	clock_t end = clock();
	double elapsed = double(end - begin) / (CLOCKS_PER_SEC / 1000);
	cout << "elapsed: " << elapsed << "ms" << " wrote " << count
			<< " point pairs (q1,q2)" << '\n';

	return 0;
}
