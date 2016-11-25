CC = g++
CPP_FILES = $(wildcard src/*.cpp)
OBJ_FILES = $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
CFLAGS    = -Wall -Wextra -Wunused -Werror -Wpedantic -g -std=c++11
O2 		  = -O2
O0	 	  = -O0
OPT 	  = $(O0)
NDEBUG	  = -DNDEBUG
CFLAGS    += $(OPT) # $(NDEBUG)

.PHONY: all clean 

all: test run

test:
	$(CC) $(CFLAGS) src/testCDW.cpp -o bin/testCDW

run:
	$(CC) $(CFLAGS) src/runCDW.cpp -o bin/runCDW

clean:
	rm -f bin/testCDW
	rm -f bin/runCDW

