#include <limits>
#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

typedef std::vector<std::vector<double>> M;

constexpr double epsilon = std::numeric_limits<double>::epsilon();

void debug(std::vector<std::vector<double> >& m){
	for(int i = 0;i < m.size();++i){
		for(int j = 0;j < m[i].size();++j){
			std::cout << m[i][j] << " ";
		}
		std::cout << std::endl;
	}
}


class Matrix {

	public:
		// Constructor
		Matrix(int row, int column);
		Matrix(int row, int column, double init);
		Matrix(std::vector<std::vector<double> >& src);

		// operator overload
		inline std::vector<double>& operator[](std::size_t index){
			return dat[index];
		}
		inline const std::vector<double>& operator[](std::size_t index) const{
			return dat[index];
		}
		Matrix operator*(Matrix& rvalue);

		double determine();
		std::vector<double> solve(std::vector<double>& rvalue);
		Matrix inverse();

		int row, column;

	private:
		M dat;
		M lu, p, A;

		int pivotsearch(M& mat, int col);
		void Gauss();
		void lu_decomp();

		template <typename T>
		static inline bool iszero(T a){
			return std::fpclassify(a) == FP_ZERO;
		}
};

Matrix::Matrix(int row, int column) : Matrix(row, column, 0.0) {
}

Matrix::Matrix(int row, int column, double init) {
	this->row = row;
	this->column = column;
	dat = M(row, std::vector<double>(column, init));
}

Matrix::Matrix(M& src){
	int n = src.size();
	this->row = src.size();
	this->column = src.size();
	this->dat = std::vector<std::vector<double> >(n);
	for(int i = 0;i < n;++i){
		this->dat[i] = std::vector<double>(src[i].size());
		for(int j = 0;j < src[i].size();++j){
			this->dat[i][j] = src[i][j];
		}
	}
}


Matrix Matrix::operator*(Matrix& rvalue){
	if(rvalue.row != this->row || rvalue.column != this->column) 
		return Matrix(this->row, this->column, NAN);

	auto res = Matrix(this->row, this->column, 0);
	for(int i = 0;i < this->row;++i){
		for(int j = 0;j < this->column;++j){
			for(int k = 0;k < this->row;++k){
				res[i][j] += this->dat[i][k] * rvalue[k][j];
			}
		}
	}
	return res;
}

double Matrix::determine(){
	double dat = epsilon;
	if(this->row != this->column) return epsilon;
	dat = 0.;
	return dat;
}

std::vector<double> Matrix::solve(std::vector<double>& rvalue){
//	if(this->row != rvalue.size())
		return std::vector<double>(rvalue.size(), NAN);
//	return this->Gauss(rvalue);
}

Matrix Matrix::inverse(){
	this->Gauss();
//	debug(this->lu);
//	debug(this->p);
	return Matrix(this->row, this->column, 0);
}

void Matrix::Gauss(){
	int index;
	this->A = this->dat;
	this->p = std::vector<std::vector<double>>(this->row, std::vector<double>(this->column, 0));
	for(int i = 0;i < this->row;++i) this->p[i][i] = 1;

	for(int i = 0;i < this->row;++i){

		// pivot select
		index = Matrix::pivotsearch(this->A, i);
		if(index != i){
			std::swap(this->A[index], this->A[i]);
			std::swap(this->p[index], this->p[i]);
		}
		double pivot = this->A[i][i];

		for(int j = i + 1;j < this->row;++j){
			if(std::fpclassify(this->A[j][i]) != FP_ZERO){
				double m = pivot / this->A[j][i];
				for(int k = i;k < this->column;++k){
					this->A[j][k] *= m;
					this->A[j][k] -= this->A[i][k];
				}
//				this->p[j][i] = m;
//				r[j] *= m;
//				r[j] -= r[i];
			}
		}
	}


	debug(this->A);
	debug(this->p);

	// backward assignment
//	for(int i = this->row - 1;i >= 0;--i){
//		if(Matrix::iszero(a[i][i])) continue;
//		for(int j = this->column - 1;j > i;--j){
//			r[i] -= a[i][j] * r[j];
//		}
//		r[i] /= a[i][i];
//	}
//
//	return r;
}

// LU Decomposition
void Matrix::lu_decomp(){
	int index;
	double pivot;

	this->lu = this->dat;
	this->p = M(this->row, std::vector<double>(this->column, 0));

	for(int i = 0;i < this->row;++i){
		this->p[i][i] = 1.;
	}

	for(int i = 0;i < this->row;++i){
		index = Matrix::pivotsearch(this->lu, i);

		if(index != i){
			std::swap(this->lu[index], this->lu[i]);
			std::swap(this->p[index], this->p[i]);
		}

		for(int j = 0;j < i;++j){
			for(int k = 0;k < j;++k){
				this->lu[i][j] -= this->lu[i][k] * this->lu[k][j];
			}
			this->lu[i][j] /= this->lu[j][j];
		}

		for(int j = i;j < this->column;++j){
			for(int k = 0;k < i;++k){
				this->lu[i][j] -= this->lu[i][k] * this->lu[k][j];
			}
		}

	}

	auto l = this->lu, u = this->lu;
	for(int i = 0;i < this->row;++i){
		for(int j = 0;j < this->column;++j){
			if(i < j) l[i][j] = 0;
			else if(i == j) l[i][j] = 1;
			else u[i][j] = 0;
		}
	}

	debug(l);
	debug(u);

//	auto L = Matrix(l), U = Matrix(u), P = Matrix(this->p);
//	auto a_ = L * U * P.inverse();
//	debug(a_.dat);

	return;
}

int Matrix::pivotsearch(M& mat, int col){
	int index = col;
	double pivot = std::abs(mat[col][col]), x;
	for(int i = col + 1;i < this->row;++i){
		x = std::abs(mat[i][col]);
		if(x > pivot){
			index = i;
			pivot = x;
		}
	}
	return index;
}

std::vector<double>* operator*(const Matrix& lvalue, const std::vector<double>& rvalue){
	auto res = new std::vector<double>(rvalue.size(), 0);
	if(rvalue.size() != lvalue.column) std::fill(res->begin(), res->end(), NAN);
	else{
		for(int i = 0;i < lvalue.row;++i){
			for(int j = 0;j < lvalue.column;++j){
				(*res)[i] += lvalue[i][j] * rvalue[j];
			}
		}
	}
	return std::move(res);
}

int main(){
	using namespace std;
	const int N = 3;
	auto d = vector<vector<double> >{{1, 2, 3}, {2, 3, 1}, {3, 2, 1}};
	auto m = Matrix(d);
	m.inverse();
//	auto r = vector<double>{6, 0.02, 10};
//	for(auto &i : r){
//		cout << i << endl;
//	}
//	auto ans = m.solve(r);
//	auto check = m * ans;
//	for(auto &e : ans){
//		cout << fixed << setprecision(5) << e << endl;
//	}
//	for(auto &e : check){
//		cout << fixed << setprecision(5) << e << endl;
//	}

	return 0;
}
