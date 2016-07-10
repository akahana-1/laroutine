#include "matrix.hpp"

int main(){
	int row, column;
	ld det;
	std::cin >> row >> column;
	auto m = Matrix(row, column);
	for(int i = 0;i < row;++i){
		for(int j = 0;j < column;++j){
			std::cin >> m[i][j];
		}
	}
	std::cin >> det;

	std::cout << m.determinant() << std::endl;
	return 0;
}
