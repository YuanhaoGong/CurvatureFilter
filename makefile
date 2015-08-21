cf:main.cpp
	g++ -o cf -O3 main.cpp -march=native `pkg-config opencv --libs --cflags`
mcf:main_MultiScale.cpp
	g++ -o mcf -O3 main_MultiScale.cpp -march=native `pkg-config opencv --libs --cflags`