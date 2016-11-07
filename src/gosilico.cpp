#include <iostream>
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <functional>

using namespace std;

//#define DEBUG
//#define TEST
#ifdef TEST
#include <armadillo> // for testing
#endif

#define EQ_NUM 2
#define EPSILON 0.000000000001
#define RANDX 8*((rand()) % 100) / (double) 100000

// construct the equation system
#define COMMON_EXPR(nu) pow((Lambda - (nu1 + delta1)*q1 - (nu2 + delta2)*q2) / cSalt, nu)
#define EQ1 k_eq1 * COMMON_EXPR(nu1) * c_eq1 - q1
#define EQ2 k_eq2 * COMMON_EXPR(nu2) * c_eq2 - q2
// partial derivatives d(eq_i)/d(q_i)
#define dCOMMON_EXPR_dq1(nu) k_eq1 * nu1 * c_eq1 * (-nu1-delta1) * COMMON_EXPR(nu-1) / cSalt
#define dCOMMON_EXPR_dq2(nu) k_eq2 * nu2 * c_eq2 * (-nu2-delta2) * COMMON_EXPR(nu-1) / cSalt
#define dEQ1_dq1 dCOMMON_EXPR_dq1(nu1) - 1
#define dEQ1_dq2 dCOMMON_EXPR_dq1(nu1)
#define dEQ2_dq1 dCOMMON_EXPR_dq2(nu2)
#define dEQ2_dq2 dCOMMON_EXPR_dq2(nu2) - 1

typedef array<array<double, EQ_NUM>, EQ_NUM> Mat;
typedef array<double, EQ_NUM> Vec;
typedef std::function<void(Vec &X, Mat &J, Vec &F)> EQ_functor;

// lVec <= lVec + rVec
void operator+=(Vec& lVec, Vec& rVec) {
	for (size_t i = 0; i < lVec.size(); i++) {
		lVec[i] += rVec[i];
	}
}

/**
 * Newton-Raphson 2D non-linear solver
 */
class NR_2DSolver {
	const int max_iterations = 50; // approximation steps
	/*
	 * solve a 2x2 linear equation: AX=B
	 * A00.X0 + A01.X1 = B0
	 * A10.X0 + A11.X1 = B1
	 * returns non-zero when there is no unique solution
	 */
	int solveLinear(Mat A, Vec B, Vec &X) {
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
	// L1 distance to zero vector
	double norm(Vec X) {
		double s = 0;
		for (size_t i = 0; i < X.size(); i++) {
			s = max(s, abs(X[i]));
		}
		return s;
	}

public:
	enum results {
		PARALLEL_LINES, SAME_LINES
	};
	// current is the start point, and will be the resulting approximation
	bool solve(EQ_functor &evaluate,
			Vec &current) {
		Mat J; // Jacobian matrix
		Vec F, Y;
		int limit = max_iterations;

		do {
			// update F(X) and J(X) with X=current
			evaluate(current, J, F);
			// calculate the equation system JY=-F.
			if (solveLinear(J, F, Y) != 0 || isinf(F[0]) || isinf(F[1])) {
				// safety measure in case of bad start point
				current[0] = RANDX;
				current[1] = RANDX;
				cout << "bad start point! " << "reset X0=" << current[0]
						<< " X1=" << current[1]
						<< '\n';
				continue;
			}

			// next approximation step
			current += Y;
		} while (norm(F) > EPSILON && limit--);

		if (limit < 0) { // the current point is not close enough
			return false;
		}
		return true;
	}

};

/**
 * models the parameterized equation system
 */
class MyEQ {
	const double Lambda = 0.5, k_eq1 = 0.14, k_eq2 = 2.7, nu1 = 6.6, nu2 = 6.2;
	const double delta1 = 38, delta2 = 58, cSalt = 0.15;

	const double max_ceq1 = 0.001, max_ceq2 = 0.0001; // domain intervals
	const double delta_c1, delta_c2;

public:
	double c_eq1, c_eq2;
	MyEQ(int n) :
			delta_c1(max_ceq1 / n), delta_c2(max_ceq2 / n) {
		set_c1c2_index(0, 0);
	}
	// update equation parameters
	void set_c1c2(double c1, double c2) {
		c_eq1 = c1;
		c_eq2 = c2;
	}
	// generate grid-like pairs (c1,c2) useful for plotting
	void set_c1c2_index(int i, int j) {
		set_c1c2(delta_c1 * i, delta_c2 * j);
	}
	// X <= F(X), J <= dF(X)/dX
	void evaluate(Vec &X, Mat &J, Vec &F) {
		double q1, q2;
		q1 = X[0];
		q2 = X[1];
		F = {-(EQ1), -(EQ2)};
		J[0]= {dEQ1_dq1, dEQ1_dq2};
		J[1]= {dEQ2_dq1, dEQ2_dq2};
	}
	void operator()(Vec &X, Mat &J, Vec &F) {
		evaluate(X, J, F);
	}
};

int main() {
	Vec Q;
	std::ofstream file1("q1.txt"), file2("q2.txt");
	clock_t begin = clock();

	const int dim = 60;
	MyEQ myEq(dim);
	NR_2DSolver solver;
	int count = 0;
	for (int i = 0; i <= dim; ++i) {
		for (int j = 0; j <= dim; ++j) {
			myEq.set_c1c2_index(i, j);
			EQ_functor evaluator(myEq);
			if (solver.solve(evaluator, Q)) {
				++count;
				file1 << myEq.c_eq1 << ' ' << myEq.c_eq2 << ' ' << Q[0] << '\n';
				file2 << myEq.c_eq1 << ' ' << myEq.c_eq2 << ' ' << Q[1] << '\n';
			}
		}
	}

	clock_t end = clock();
	double elapsed = double(end - begin) / (CLOCKS_PER_SEC / 1000);
	cout << "execution time: " << elapsed << "ms" << ", for " << count
			<< " pairs (q1,q2)" << '\n';
	return 0;
}

#ifdef TEST
// is a almost equal to b?
bool inline equal(double &a, double &b) {
	return a == b || fabs(a - b) < EPSILON;
}
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
