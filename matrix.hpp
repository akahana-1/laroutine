#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <limits>
#include <vector>
#include <utility>
#include <iomanip>
#include <algorithm>
#include <cmath>

typedef std::vector<std::vector<double>> M;
using ld=long double;

constexpr ld epsilon = std::numeric_limits<ld>::epsilon();

const ld tol = 1e-10;

class Matrix {

	public:
		// Constructor
		Matrix();
		Matrix(int row, int column);
		Matrix(int row, int column, double init);
		Matrix(M& m);
		Matrix(const Matrix& src);

		// operator overload
		inline std::vector<double>& operator[](std::size_t index){
			return dat[index];
		}
		inline const std::vector<double>& operator[](std::size_t index) const{
			return dat[index];
		}
		Matrix operator*(Matrix const& rvalue);
		Matrix operator+(Matrix const& rvalue);
		Matrix operator-(Matrix const& rvalue);
		bool operator==(Matrix const& rvalue);

		int row, column;

		ld determinant();
		std::vector<double> solve(std::vector<double>& rvalue);
		Matrix inverse(), triu(), tril(), diag(), eig(), transpose();
		void transform();

		static Matrix eye(int size){
			auto m = Matrix(size, size);
			for(int i = 0;i < size;++i){
				m[i][i] = 1.;
			}
			return m;
		}


	private:
		M dat;

		int pivotsearch(M& mat, int col);
		Matrix gauss();
		int lu_decomp(Matrix& lu, Matrix& p);
		void qr_decomp(Matrix& q, Matrix& r);
		Matrix householder();
		std::vector<double> inverseiterate(double sigma, std::vector<double> vec = std::vector<double>(0));
		std::vector<double> helpersolve(Matrix const& lu, std::vector<double> const& b);
};

std::vector<double> operator*(Matrix& lvalue, std::vector<double>& rvalue);
Matrix operator*(double const& k, Matrix const& A);

#endif
