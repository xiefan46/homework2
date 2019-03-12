OBJS = val_test01_solved val_test02_solved mmult1 omp_bug2_solved omp_bug3_solved omp_bug4_solved omp_bug5_solved omp_bug6_solved
all : $(OBJS)

val_test01_solved : val_test01_solved.cpp
	g++ val_test01_solved.cpp -o val_test01_solved
val_test02_solved : val_test02_solved.cpp
	g++ val_test02_solved.cpp -o val_test02_solved
mmult1 : MMult1.cpp
	g++ -std=c++11 -fopenmp -O3 -march=native MMult1.cpp -o MMult1
omp_bug2_solved : omp_bug2_solved.c
	g++ -std=c++11 -fopenmp -O3 -march=native omp_bug2_solved.c -o omp_bug2_solved
omp_bug3_solved : omp_bug3_solved.c
	g++ -std=c++11 -fopenmp -O3 -march=native omp_bug3_solved.c -o omp_bug3_solved
omp_bug4_solved : omp_bug4_solved.c
	g++ -std=c++11 -fopenmp -O3 -march=native omp_bug4_solved.c -o omp_bug4_solved
omp_bug5_solved : omp_bug5_solved.c
	g++ -std=c++11 -fopenmp -O3 -march=native omp_bug5_solved.c -o omp_bug5_solved
omp_bug6_solved : omp_bug6_solved.c
	g++ -std=c++11 -fopenmp -O3 -march=native omp_bug6_solved.c -o omp_bug6_solved
clean :
	rm $(OBJS)