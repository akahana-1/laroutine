CXX=g++
CXXFLAGS= -std=c++11 -O2

all: main.cpp matrix.cpp
	$(CXX) $^ $(CXXFLAGS)
