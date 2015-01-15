SDSL := ${HOME}/opt/sdsl
CC := g++ -std=c++11 -O2 -g -Wall -I${SDSL}/include -L${SDSL}/lib

lz_random_test: lz_index.hpp lz_random_test.cpp
	${CC} lz_random_test.cpp -o lz_random_test -lsdsl

lz_memory_visualization: lz_index.hpp lz_memory_visualization.cpp
	${CC} lz_memory_visualization.cpp -o lz_memory_visualization -lsdsl

.PHONY: clean

clean:
	rm -f lz_random_test *.o
