#include <iostream>
#include <random>
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

ld norm(std::vector<double> const& vec, int dim = 2){
	ld x = 0.;
	for(auto&& e : vec){
		x += pow(e, dim);
	}
	return pow(x, 1. / (double)dim);
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

Matrix Matrix::operator+(Matrix const& rvalue){
	if(rvalue.row != this->row || rvalue.column != this->column) 
		return Matrix(this->row, this->column, NAN);

	auto res = *this;
	for(int i = 0;i < this->row;++i){
		for(int j = 0;j < this->column;++j){
			res[i][j] += rvalue[i][j];
			if(iszero(res[i][j])) res[i][j] = 0;
		}
	}
	return res;
}

Matrix Matrix::operator-(Matrix const& rvalue){
	if(rvalue.row != this->row || rvalue.column != this->column) 
		return Matrix(this->row, this->column, NAN);

	auto res = *this;
	for(int i = 0;i < this->row;++i){
		for(int j = 0;j < this->column;++j){
			res[i][j] -= rvalue[i][j];
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
	int detsign;
	bool f = false;
	if(this->row != this->column) return det;
	Matrix lu = *this, p = Matrix::eye(this->row);
	detsign = this->lu_decomp(lu, p);
	det = 1;
	for(int i = 0;i < this->row;++i) det *= lu[i][i];

	return detsign * det;
}

std::vector<double> Matrix::solve(std::vector<double>& rvalue){
	if(this->row != rvalue.size())
		return std::vector<double>(rvalue.size(), NAN);
	Matrix lu, p;
	ld det = 1.;
	this->lu_decomp(lu, p);

	// 正則判定
	for(int i = 0;i < this->row;++i){
		det *= lu[i][i];
	}
	if(iszero(det)) return std::vector<double>(rvalue.size(), NAN);

	auto b = p * rvalue;
	auto res = this->helpersolve(lu, b);

	return res;
}

Matrix Matrix::eig(){
	Matrix ak, q(this->row, this->column), r(this->row, this->column), result(this->row, 2 * this->column), hb, eigvec = Matrix::eye(this->row);
	int m = this->row;

	bool isdiag = *(this) == this->diag();

	hb = isdiag ? *(this) : this->householder();
	ak = hb;

	// Calc eigen value (QR method)
	while(m != 1){
		if(iszero(ak[m - 1][m - 2])) --m;
		ak.qr_decomp(q, r);
		ak = r * q;
		eigvec = eigvec * q;
	}

	for(int j = 0;j < this->column;++j){
		double nx = 0.;
		for(int i = 0;i < this->row;++i){
			nx += pow(eigvec[i][j], 2);
		}
		nx = std::sqrt(nx);
		for(int i = 0;i < this->row;++i){
			eigvec[i][j] /= nx;
		}
	}

	// Calc eigen vector
	for(int i = 0;i < this->column;++i){
		result[i][i] = ak[i][i];
		std::vector<double> est = std::vector<double>(this->row);
		for(int j = 0;j < this->row;++j){
			est[j] = eigvec[j][i];
		}
		auto v = this->inverseiterate(result[i][i] + 1.e-4, est);
		for(int j = 0;j < this->row;++j){
			result[j][this->column + i] = v[j];
		}
	}

	return result;
}

std::vector<double> Matrix::inverseiterate(double sigma, std::vector<double> vec){
	Matrix G = (sigma * Matrix::eye(this->row)) - *(this);
	G = G.inverse();
	double nx = 0., resid = 0.;

	if(!vec.size()){
		vec = std::vector<double>(this->row);
		std::random_device rnd;
		std::mt19937 mt(rnd());

		std::uniform_real_distribution<> dist(0.0, 1.0);

		for(int i = 0;i < this->row;++i){
			vec[i] = dist(mt);
		}
		nx = norm(vec);
		for(int i = 0;i < this->row;++i){
			vec[i] /= nx;
		}
	}

	for(int i = 0;;++i){
		auto x = G * vec;
		nx = norm(x);
		resid = 0.;
		for(int i = 0;i < this->row;++i){
			x[i] /= nx;
			resid += pow(fabs(x[i]) - fabs(vec[i]), 2);
		}
		resid = std::sqrt(resid);
		if(iszero(resid)) break;
		std::copy(x.begin(), x.end(), vec.begin());
	}

	return vec;
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
			m[i][j] = this->dat[i][j];
		}
	}
	return m;
}

Matrix Matrix::triu(){
	auto m = Matrix(this->row, this->column);
	for(int i = 0;i < this->row;++i){
		for(int j = i;j < this->column;++j){
			m[i][j] = this->dat[i][j];
		}
	}
	return m;
}

Matrix Matrix::diag(){
	auto m = Matrix(this->row, this->column);
	for(int i = 0;i < this->row;++i){
		m[i][i] = this->dat[i][i];
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

int Matrix::lu_decomp(Matrix& lu, Matrix& p){
	int index, detsign = 1;
	double pivot, div;

	lu = *this;
	p = Matrix::eye(this->row);

	for(int i = 0;i < this->row;++i){
		index = Matrix::pivotsearch(lu.dat, i);

		if(index != i){
			std::swap(lu.dat[index], lu.dat[i]);
			std::swap(p.dat[index], p.dat[i]);
			detsign *= -1;
		}

		// 軸選択をしても0となる=>特異行列
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

	return detsign;
}

void Matrix::qr_decomp(Matrix& q, Matrix& r){
	q = *this;
	r = Matrix::eye(this->row);
	double dotprod, norm;
	
	// Gram-Schmit
	for(int j = 0;j < this->column;++j){
		norm = 0.;
		for(int k = 0;k < j;++k){
			dotprod = 0.;
			for(int i = 0;i < this->row;++i){
				dotprod += this->dat[i][j] * q[i][k];
			}
			for(int i = 0;i < this->row;++i){
				q[i][j] -= dotprod * q[i][k];
			}
			r[k][j] = dotprod;
		}
		for(int i = 0;i < this->row;++i){
			norm += std::pow(q[i][j], 2);
		}
		norm = sqrt(norm);
		for(int i = 0;i < this->row;++i){
			q[i][j] /= norm;
		}
		r[j][j] = norm;
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

Matrix Matrix::householder(){
	Matrix res = *this, p = Matrix::eye(this->row);
	double s, c;
	std::vector<double> w(this->row), q(this->row), aw, waw;

	for(int i = 0;i < this->row - 1;++i){
		s = c = 0.;
		for(int k = i + 1;k < this->row;++k){
			s += pow(res[k][i], 2);
		}
		s = (std::signbit(res[i + 1][i]) ? -1 : 1) * sqrt(s);

		c = 1 / (pow(s, 2) + res[i + 1][i] * s);

		for(int k = 0;k < this->row;++k){
			if(k <= i) w[k] = 0.;
			else w[k] = res[k][i];
		}
		w[i + 1] += s;

		aw = res * w;

		p = Matrix::eye(this->row);
		for(int i = 0;i < this->row;++i){
			for(int j = 0;j < this->column;++j){
				p[i][j] -= c * w[i] * w[j];
			}
		}

		res = p * res * p;
	}

	return res;
}

void Matrix::transform(){
	Matrix res, q, r;
	res = this->eig();
	debug(res.dat);
}

Matrix Matrix::transpose(){
	Matrix m = {this->column, this->row};
	for(int i = 0;i < this->row;++i){
		for(int j = 0;j < this->column;++j){
			m[j][i] = this->dat[i][j];
		}
	}
	return m;
}

std::vector<double> operator*(Matrix& lvalue, std::vector<double>& rvalue){
	auto res = std::vector<double>(rvalue.size(), 0);
	if(rvalue.size() != lvalue.column) std::fill(res.begin(), res.end(), NAN);
	else{
		for(int i = 0;i < lvalue.row;++i){
			for(int j = 0;j < lvalue.column;++j){
				res[i] += lvalue[i][j] * rvalue[j];
			}
			res[i] = iszero(res[i]) ? 0. : res[i];
		}
	}
	return res;
}

Matrix operator*(double const& k, Matrix const& A){
	Matrix res = A;
	for(int i = 0;i < A.row;++i){
		for(int j = 0;j < A.column;++j){
			res[i][j] *= k;
		}
	}
	return res;
}
