all: main

main:		main.cpp
	g++ -O3 -DNDEBUG -std=c++17 -Wall -I./include/parlaylib/include -I./include -I/usr/include/boost1.78 -o main main.cpp -L/data/bhuan102/libspatialindex/lib64 -lspatialindex -pthread

clean:
	rm -f test main


# g++ -O3 -DNDEBUG -std=c++17 -Wall -Iinclude -o main main.cpp --L/data/bhuan102/libspatialindex/lib64 lspatialindex
# g++ -g -O3 -DDEBUG -std=c++17 -Wall -Iinclude -o main main.cpp -L/data/bhuan102/libspatialindex/lib64 -lspatialindex