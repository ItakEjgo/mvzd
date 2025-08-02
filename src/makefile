all: main

CXX := /usr/bin/g++
# CXX := clang++

main:		main.cpp geobase.h zdtree_node.hpp zdtree.hpp morton.hpp zdtest.hpp seq_zdtree.hpp hilbert.c hilbert.h binary_rtree.hpp
	$(CXX) -O3 -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../include -o main main.cpp hilbert.c -L/usr/local/lib -Wl,-R/usr/local/lib -ljemalloc
	# g++ -O3 -DNDEBUG -mcx16 -march=native -DPARLAY_SEQUENTIAL -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../include -o main main.cpp hilbert.c -L/usr/local/lib -Wl,-R/usr/local/lib -ljemalloc
	# g++ -O3 -DNDEBUG -mcx16 -march=native -DHOMEGROWN -DPARLAY_USE_STD_ALLOC -pthread -std=c++17 -Wall -I../parlaylib/include -I../include -o main main.cpp hilbert.c -L/usr/local/lib -Wl,-R/usr/local/lib -ljemalloc
	# g++ -O3 -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../include -o main main.cpp hilbert.c -L/usr/local/lib -Wl,-R/usr/local/lib -ljemalloc -DPARLAY_SEQUENTIAL
	# g++ -O3 -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../include -o main main.cpp -L/usr/local/lib -Wl,-R/usr/local/lib -ljemalloc
	# g++ -O3 -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../include -o main main.cpp -L/usr/local/lib -ljemalloc
	# g++ -O3 -DDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../include -o main main.cpp -L/usr/local/lib -ljemalloc

clean:
	rm -f test main
