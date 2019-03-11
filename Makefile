OBJS = val_test01_solved val_test02_solved mmult1
all : $(OBJS)
	
val_test01_solved : val_test01_solved.cpp
	g++ val_test01_solved.cpp -o val_test01_solved 
val_test02_solved : val_test02_solved.cpp
	g++ val_test02_solved.cpp -o val_test02_solved
mmult1 : MMult1.cpp
	g++ -std=c++11 -march=native -O2  MMult1.cpp -o MMult1
clean : 
	rm $(OBJS)
