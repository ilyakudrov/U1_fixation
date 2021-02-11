INCLUDE_PATH=../include
SRC_PATH=../src
OBJ_PATH=../obj
GCC=g++
#FLAGS=-std=c++11 -O0 -D DOUBLE -g
FLAGS=-std=c++11 -O3 -D DOUBLE
CMINPACK_PATH=/home/ilya/code/cminpack

OBJFILES=$(CMINPACK_PATH)/hybrd1.o $(CMINPACK_PATH)/enorm.o $(CMINPACK_PATH)/dpmpar.o \
 $(CMINPACK_PATH)/hybrd.o $(CMINPACK_PATH)/qrfac.o $(CMINPACK_PATH)/qform.o \
  $(CMINPACK_PATH)/dogleg.o $(CMINPACK_PATH)/r1updt.o $(CMINPACK_PATH)/r1mpyq.o \
  $(CMINPACK_PATH)/fdjac1.o

minpack_root_finding: $(OBJ_PATH)/minpack_root_finding.o $(OBJFILES)
	$(GCC) $(OBJ_PATH)/minpack_root_finding.o $(OBJFILES) -lf2c -lm -o minpack_root_finding

$(OBJ_PATH)/minpack_root_finding.o: minpack_root_finding.cpp
	$(GCC) $(FLAGS) -c minpack_root_finding.cpp -o $(OBJ_PATH)/minpack_root_finding.o

gnu_root_finding: $(OBJ_PATH)/gnu_root_finding.o
	$(GCC) $(OBJ_PATH)/gnu_root_finding.o -lgsl -lgslcblas -lm -o gnu_root_finding

$(OBJ_PATH)/gnu_root_finding.o: gnu_root_finding.cpp
	$(GCC) $(FLAGS) -c gnu_root_finding.cpp -o $(OBJ_PATH)/gnu_root_finding.o

maximization_test: $(OBJ_PATH)/maximization_test.o $(OBJ_PATH)/link.o $(OBJ_PATH)/matrix.o $(OBJ_PATH)/roots.o $(OBJ_PATH)/maximization.o
	$(GCC) $(OBJ_PATH)/maximization_test.o $(OBJ_PATH)/link.o $(OBJ_PATH)/matrix.o $(OBJ_PATH)/roots.o $(OBJ_PATH)/maximization.o -lgsl -lgslcblas -lm -o maximization_test

$(OBJ_PATH)/maximization_test.o: maximization_test.cpp
	$(GCC) $(FLAGS) -c maximization_test.cpp -o $(OBJ_PATH)/maximization_test.o

$(OBJ_PATH)/link.o: $(SRC_PATH)/link.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/link.cpp -o $(OBJ_PATH)/link.o

$(OBJ_PATH)/matrix.o: $(SRC_PATH)/matrix.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/matrix.cpp -o $(OBJ_PATH)/matrix.o

$(OBJ_PATH)/roots.o: $(SRC_PATH)/roots.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/roots.cpp -o $(OBJ_PATH)/roots.o

$(OBJ_PATH)/maximization.o: $(SRC_PATH)/maximization.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/maximization.cpp -o $(OBJ_PATH)/maximization.o

clean:
	rm -rf $(OBJ_PATH)/*.o
