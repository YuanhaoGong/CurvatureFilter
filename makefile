cf:main.cpp
	g++ -o cf -O3 main.cpp -march=native `pkg-config opencv --libs --cflags`