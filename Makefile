# Makefile for Multigrid Solver (flat version, no subfolders)

CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
EXEC = multigrid_solver

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o $(EXEC) *.txt

run: all
	./$(EXEC)
