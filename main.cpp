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

	auto res = m.solve(b);

	auto check = m * res;

	for(auto &&x : res){
		std::cout << x << std::endl;
	}

	for(auto &&x : check){
		std::cout << x << std::endl;
	}

	std::cout << m.determinant() << std::endl;

	return 0;
}
