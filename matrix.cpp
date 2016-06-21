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
		inline std::vector<double>& operator[](int index){
			return dat[index];
		}
		std::vector<double> operator*(std::vector<double>& rvalue);

		double determine();
		std::vector<double> solve(std::vector<double>& rvalue);
		Matrix inverse();

	private:
		int row, column;
		M dat;
		M lu, elem;

		int pivotsearch(M& mat, int col);
		std::vector<double> Gauss(std::vector<double>& rvalue);
		void LU();
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

std::vector<double> Matrix::operator*(std::vector<double>& rvalue){
	std::vector<double> res = std::vector<double>(rvalue.size(), 0);
	if(rvalue.size() != this->column) return res;
	for(int i = 0;i < this->row;++i){
		for(int j = 0;j < this->column;++j){
			res[i] += this->dat[i][j] * rvalue[j];
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
	if(this->row != rvalue.size())
		return std::vector<double>(rvalue.size(), 0);
	return this->Gauss(rvalue);
}

Matrix Matrix::inverse(){
	return Matrix(this->row, this->column, 0);
}

std::vector<double> Matrix::Gauss(std::vector<double>& rvalue){
	int index;
	auto a = this->dat;
	std::vector<double> r = rvalue;
	for(int i = 0;i < this->row;++i){

		// pivot select
		index = Matrix::pivotsearch(a, i);
		if(index != i){
			std::swap(a[index], a[i]);
			std::swap(r[index], r[i]);
		}
		double pivot = a[i][i];

		// forward delete
		for(int j = i + 1;j < this->row;++j){
			if(std::fpclassify(a[j][i]) != FP_ZERO){
				double m = pivot / a[j][i];
				for(int k = i;k < this->column;++k){
					a[j][k] *= m;
					a[j][k] -= a[i][k];
				}
				r[j] *= m;
				r[j] -= r[i];
			}
		}
	}

	// backward assignment
	for(int i = this->row - 1;i >= 0;--i){
		if(Matrix::iszero(a[i][i])) continue;
		for(int j = this->column - 1;j > i;--j){
			r[i] -= a[i][j] * r[j];
		}
		r[i] /= a[i][i];
	}

	return r;
}

// LU Decomposition
void Matrix::LU(){
	int index;
	double pivot;
	auto a = this->dat;

	this->lu = M(this->row, std::vector<double>(this->column, 0));
	this->elem = M(this->row, std::vector<double>(this->column, 0));

	for(int i = 0;i < this->row;++i){
		this->elem[i][i] = 1.;
	}

	for(int i = 0;i < this->row;++i){
	}
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

int main(){
	using namespace std;
	const int N = 3;
	auto d = vector<vector<double> >{{3, 5, 2}, {0.001, 0.002, 0.05}, {2, 100, 3}};
	auto r = vector<double>{6, 0.02, 10};
	for(auto &i : r){
		cout << i << endl;
	}
	auto m = Matrix(d);
	auto ans = m.solve(r);
	auto check = m * ans;
	for(auto &e : ans){
		cout << fixed << setprecision(5) << e << endl;
	}
	for(auto &e : check){
		cout << fixed << setprecision(5) << e << endl;
	}
}
