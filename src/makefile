Mutation_rate_estimator: Mutation_rate_estimator.o Newton_method.o sketch.o MurmurHash3.o
	g++ -std=c++17 Mutation_rate_estimator.o Newton_method.o sketch.o MurmurHash3.o -o ../Mutation_rate_estimator -lmpfr -lgmp

Mutation_rate_estimator.o: Mutation_rate_estimator.cpp 
	g++ -std=c++17 -c Mutation_rate_estimator.cpp

Newton_method: Newton_method.o
	g++ -std=c++17 Newton_method.o -o Newton_method -lmpfr -lgmp

Newton_method.o: Newton_method.cpp
	g++ -std=c++17 -c Newton_method.cpp 

sketch: sketch.o MurmurHash3.o
	g++ -std=c++17 sketch.o MurmurHash3.o -o sketch

sketch.o: sketch.cpp 
	g++ -std=c++17 -c sketch.cpp

MurmurHash3: MurmurHash3.o
	g++ -std=c++17 MurmurHash3.o -o MurmurHash3

MurmurHash3.o: MurmurHash3.cpp
	g++ -std=c++17 -c MurmurHash3.cpp

clean:
	rm -f *.o

