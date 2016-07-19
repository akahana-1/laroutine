#include "matrix.hpp"

void debug(std::vector<std::vector<double> >& m){
	for(int i = 0;i < m.size();++i){
		for(int j = 0;j < m[i].size();++j){
			std::cout << m[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

bool iszero(double a){
	return fabs(a - 0) < tol;
}

Matrix::Matrix() : Matrix(0, 0){
}

Matrix::Matrix(int row, int column) : Matrix(row, column, 0.0) {
}

Matrix::Matrix(int row, int column, double init) {
	this->row = row;
	this->column = column;
	this->dat = M(row, std::vector<double>(column, init));
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


Matrix Matrix::operator*(Matrix const& rvalue){
	if(rvalue.row != this->row || rvalue.column != this->column) 
		return Matrix(this->row, this->column, NAN);

	auto res = Matrix(this->row, this->column, 0);
	for(int i = 0;i < this->row;++i){
		for(int j = 0;j < this->column;++j){
			for(int k = 0;k < this->row;++k){
				res[i][j] += this->dat[i][k] * rvalue[k][j];
			}
			if(iszero(res[i][j])) res[i][j] = 0;
		}
	}
	return res;
}

bool Matrix::operator==(Matrix const& rvalue){
	bool res = true;
	if(this->row != rvalue.row || this->column != rvalue.column) return false;
	for(int i = 0;i < this->row;++i){
		for(int j = 0;j < this->column;++j){
			res = res && iszero(this->dat[i][j] - rvalue[i][j]);
		}
	}
	return res;
}

ld Matrix::determinant(){
	ld det = NAN;
	bool f = false;
	if(this->row != this->column) return det;
	Matrix lu = *this, p = Matrix::eye(this->row);
	this->lu_decomp(lu, p);
	det = 1;
	for(int i = 0;i < this->row;++i) det *= lu[i][i];
	for(int i = 0;i < this->row;++i) {
		f = true; det *= p[i][i] == 1 ? 1 : -1;
	}

	return f ? -1 * det : det;
}

std::vector<double> Matrix::solve(std::vector<double>& rvalue){
	if(this->row != rvalue.size())
		return std::vector<double>(rvalue.size(), NAN);
	Matrix lu, p;
	this->lu_decomp(lu, p);
	Matrix pinv = p.inverse();

	auto L = lu.tril(), U = lu.triu();
	for(int i = 0;i < this->row;++i) L[i][i] = 1.;
	debug(L.dat);
	debug(U.dat);
	debug(p.dat);
	auto b = pinv * rvalue;

	auto res = this->helpersolve(lu, b);
//	res = pinv * res;

	return res;
}

Matrix Matrix::inverse(){
	Matrix inv(this->row, this->column);
	std::vector<double> b;

	Matrix lu, p;
	this->lu_decomp(lu, p);


	for(int j = 0;j < this->column;++j){
		b = std::vector<double>(this->row, 0.);
		for(int i = 0;i < this->row;++i){
			if(p[i][j]) b[i] = 1;
		}
		auto ans = this->helpersolve(lu, b);
		for(int i = 0;i < this->row;++i){
			inv[i][j] = ans[i];
		}
	}

	return inv;
}

Matrix Matrix::tril(){
	auto m = Matrix(this->row, this->column);
	for(int i = 0;i < this->row;++i){
		for(int j = 0;j <= i;++j){
			m[i][j] = (*this)[i][j];
		}
	}
	return m;
}

Matrix Matrix::triu(){
	auto m = Matrix(this->row, this->column);
	for(int i = 0;i < this->row;++i){
		for(int j = i;j < this->column;++j){
			m[i][j] = (*this)[i][j];
		}
	}
	return m;
}

Matrix Matrix::gauss(){
	int index;

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
		}
		double pivot = A[i][i];

		for(int j = i + 1;j < this->row;++j){
			if(!iszero(A[j][i])){
				double m = pivot / A[j][i];
				for(int k = i;k < column;++k){
					A[j][k] *= m;
					A[j][k] -= A[i][k];
				}
				p[j][indices[i]] = -1;
				p[j][indices[j]] = m;
			}
		}

		auto t = i == 0 ? P * p : p * P;

		P = Matrix(t);
	}

	return A;

}

void Matrix::lu_decomp(Matrix& lu, Matrix& p){
	int index;
	double pivot, div;

	lu = *this;
	p = Matrix::eye(this->row);

	/*
	for(int j = 0;j < this->column;++j){
		index = Matrix::pivotsearch(lu.dat, j);

		if(index != j){
			std::swap(lu.dat[index], lu.dat[j]);
			std::swap(p.dat[index], p.dat[j]);
		}

		// 軸選択をしても改善されない時
		if(iszero(lu[j][j])) {
			lu = Matrix(this->row, this->column, NAN);
			break;
		}

		for(int i = j + 1;i < this->row;++i){
			lu[i][j] /= lu[j][j];
			for(int k = j + 1;k < this->column;++k){
				lu[i][k] -= lu[i][j] * lu[j][k];
			}
		}
	}
	*/

	for(int i = 0;i < this->row;++i){
		index = Matrix::pivotsearch(lu.dat, i);

		if(index != i){
			std::swap(lu.dat[index], lu.dat[i]);
			std::swap(p.dat[index], p.dat[i]);
		}

		// 軸選択をしても改善されない時
		if(iszero(lu[i][i])) {
			lu = Matrix(this->row, this->column, NAN);
			break;
		}

		for(int j = 0;j < i;++j){
			div = 1. / lu.dat[j][j];
			for(int k = 0;k < j;++k){
				lu.dat[i][j] -= lu.dat[i][k] * lu.dat[k][j];
			}
			lu.dat[i][j] *= div;
		}

		for(int j = i;j < this->column;++j){
			for(int k = 0;k < i;++k){
				lu.dat[i][j] -= lu.dat[i][k] * lu.dat[k][j];
			}
		}
	}

}

std::vector<double> Matrix::helpersolve(Matrix const& lu, std::vector<double> const& b){
	auto ans = b;

	for(int i = 1;i < ans.size();++i){
		for(int j = 0;j < i;++j){
			ans[i] -= lu[i][j] * ans[j];
		}
		ans[i] = iszero(ans[i]) ? 0. : ans[i];
	}

	for(int i = lu.row - 1;0 <= i;--i){
		for(int j = i + 1;j < lu.column;++j){
			ans[i] -= lu[i][j] * ans[j];
		}
		ans[i] /= lu[i][i];
		ans[i] = iszero(ans[i]) ? 0. : ans[i];
	}

	return ans;
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

std::vector<double> operator*(Matrix& lvalue, std::vector<double>& rvalue){
	auto res = std::vector<double>(rvalue.size(), 0);
	if(rvalue.size() != lvalue.column) std::fill(res.begin(), res.end(), NAN);
	else{
		for(int i = 0;i < lvalue.row;++i){
			for(int j = 0;j < lvalue.column;++j){
				res[i] += lvalue[i][j] * rvalue[j];
			}
		}
	}
	return res;
}
