#include <limits>
#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <tuple>

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
		Matrix(M& m);
		Matrix(const Matrix& src);

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

		static Matrix eye(int size){
			auto m = Matrix(size, size);
			for(int i = 0;i < size;++i){
				m[i][i] = 1.;
			}
			return m;
		}

		int row, column;

	private:
		M dat;

		int pivotsearch(M& mat, int col);
		std::tuple<Matrix, Matrix> gauss(), lu_decomp();

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

Matrix::Matrix(const Matrix& src) : Matrix(src.row, src.column) {
	for(int i = 0;i < this->row;++i){
		for(int j = 0;j < this->column;++j){
			(*this)[i][j] = src[i][j];
		}
	}
}

Matrix::Matrix(M& m){
	this->row = m.size();
	this->column = m[0].size();
	this->dat = M(this->row, std::vector<double>(this->column));
	for(int i = 0;i < this->row;++i){
		for(int j = 0;j < this->column;++j){
			this->dat[i][j] = m[i][j];
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
			if(Matrix::iszero(res[i][j])) res[i][j] = 0;
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
//	return this->gauss(rvalue);
}

Matrix Matrix::inverse(){
	this->gauss();
//	debug(this->lu);
//	debug(this->p);
	return Matrix(this->row, this->column, 0);
}

std::tuple<Matrix, Matrix> Matrix::gauss(){
	int index;
	//this->A = this->dat;
	auto A = *this;
	auto P = Matrix::eye(this->row);

	auto indices = std::vector<int>(this->row);

	for(int i = 0;i < this->row - 1;++i){

		for(int j = 0;j < this->row;++j) indices[j] = j;

		// pivot select
		auto p = Matrix::eye(this->row);
		index = Matrix::pivotsearch(A.dat, i);
		if(index != i){
			std::swap(A.dat[index], A.dat[i]);
			std::swap(indices[index], indices[i]);
			std::swap(p.dat[index], p.dat[i]);
	//		std::cout << i << " " << index << std::endl;
		}
		double pivot = A[i][i];

		for(int j = i + 1;j < this->row;++j){
			if(!Matrix::iszero(A[j][i])){
				double m = pivot / A[j][i];
				for(int k = i;k < column;++k){
					A[j][k] *= m;
					A[j][k] -= A[i][k];
				}
				p[j][indices[i]] = -1;
				p[j][indices[j]] = m;
			}
		}

//		debug(P.dat);
//		debug(p.dat);

		auto t = i == 0 ? P * p : p * P;

		P = Matrix(t);
	}

	debug(A.dat);
//
	debug(P.dat);

	auto A_ = P * (*this);

	debug(A_.dat);

	return std::forward_as_tuple(A, P);

}

// LU Decomposition
std::tuple<Matrix, Matrix> Matrix::lu_decomp(){
	int index;
	double pivot;

	auto LU = *this;
	auto P = Matrix::eye(this->row);

	for(int i = 0;i < this->row;++i){
		index = Matrix::pivotsearch(LU.dat, i);

		if(index != i){
			std::swap(LU.dat[index], LU.dat[i]);
			std::swap(P.dat[index], P.dat[i]);
		}

		for(int j = 0;j < i;++j){
			for(int k = 0;k < j;++k){
				LU[i][j] -= LU[i][k] * LU[k][j];
			}
			LU[i][j] /= LU[j][j];
		}

		for(int j = i;j < this->column;++j){
			for(int k = 0;k < i;++k){
				LU[i][j] -= LU[i][k] * LU[k][j];
			}
		}

	}

	return std::forward_as_tuple(LU, P);
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
	auto d = vector<vector<double> >{{1, 1, 1}, {3, 2, 2}, {2, -1, 3}};
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
