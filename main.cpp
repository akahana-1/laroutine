#include <iostream>
#include "matrix.hpp"

int main(){
	int row, column;
	ld det;
	std::cin >> row >> column;

	auto m = Matrix(row, column);
	std::vector<double> b = std::vector<double>(row);
	for(int i = 0;i < row;++i){
		for(int j = 0;j < column;++j){
			std::cin >> m[i][j];
		}
	}
	for(int i = 0;i < row;++i){
		std::cin >> b[i];
	}

	std::cout << "determinant" << std::endl;
	std::cout << m.determinant() << std::endl;

	auto eigen = m.eig();
	std::cout << "eigen values" << std::endl;
	for(int i = 0;i < eigen.row;++i){
		std::cout << "t" << i << " : " << eigen[i][i] << std::endl;
	}

	std::cout << "eigen vectors" << std::endl;
	for(int j = eigen.row;j < eigen.column;++j){
		std::cout << "v" << j - eigen.row << " : ";
		for(int i = 0;i < eigen.row;++i){
			std::cout << eigen[i][j] << " ";
		}
		std::cout << std::endl;
	}

	auto x = m.solve(b);
	std::cout << "solve result" << std::endl;
	for(int i = 0;i < x.size();++i){
		std::cout << "x" << i << " : " << x[i] << std::endl;
	}

	return 0;
}
