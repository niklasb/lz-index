SDSL := ${HOME}/opt/sdsl

lz_index: lz_index.cpp
	g++ -std=c++11 -O0 -g -Wall -I${SDSL}/include -L${SDSL}/lib lz_index.cpp -o lz_index -lsdsl

test: lz_index
	./lz_index
