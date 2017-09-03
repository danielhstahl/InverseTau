INCLUDES= -I../FangOost -I../FunctionalUtilities -I../ODESolver 
GCCVAL=g++

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	GCCVAL=g++-7
endif
invertTau:main.o 
	$(GCCVAL) -std=c++14 -O3 -pthread --coverage -g  main.o $(INCLUDES) -o invertTau -fopenmp

main.o: main.cpp 
	$(GCCVAL) -std=c++14 -O3 -pthread --coverage   -g  -c main.cpp   $(INCLUDES) -fopenmp

clean:
	-rm *.o test *.out



