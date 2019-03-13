OBJS = val_test01_solved val_test02_solved mmult1 omp_solved2 omp_solved3 omp_solved4 omp_solved5 omp_solved6
all : $(OBJS)

val_test01_solved : val_test01_solved.cpp
	g++ val_test01_solved.cpp -o val_test01_solved
val_test02_solved : val_test02_solved.cpp
	g++ val_test02_solved.cpp -o val_test02_solved
mmult1 : MMult1.cpp
	g++ -std=c++11 -fopenmp -O3 -march=native MMult1.cpp -o MMult1
omp_solved2 : omp_solved2.c
	g++ -std=c++11 -fopenmp -O3 -march=native omp_solved2.c -o omp_solved2
omp_solved3 : omp_bug3_solved.c
	g++ -std=c++11 -fopenmp -O3 -march=native omp_solved3.c -o omp_solved3
omp_solved4 : omp_solved4.c
	g++ -std=c++11 -fopenmp -O3 -march=native omp_solved4.c -o omp_solved4
omp_solved5 : omp_solved5.c
	g++ -std=c++11 -fopenmp -O3 -march=native omp_solved5.c -o omp_solved5
omp_solved6 : omp_solved6.c
	g++ -std=c++11 -fopenmp -O3 -march=native omp_solved6.c -o omp_solved6
clean :
	rm $(OBJS)
