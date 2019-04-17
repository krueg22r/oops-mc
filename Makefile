CC=g++-4.9
CFLAGS=-Wall -std=c++11  -O3
INC=-I/home/eigen
ARCH = -arch x86_64
CPP_FILES=$(wildcard *.cpp)
SRC=$(CPP_FILES)
OBJ=$(SRC:.cpp=.o)
EXECUTABLE=exe

TARGET=mc

default: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) -o $@ $(OBJ)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

clean: 
	-rm $(OBJ)

