.PHONY: clean

# Test for needed library files.
CPATH=


CC = gcc


CFLAGS =
CFLAGS += -Wall -Wextra
CFLAGS += -fopenmp 

# CFLAGS += -O0
# CFLAGS += -O2
CFLAGS += -O3 -flto


SOFLAGS = -fPIC -shared


LDFLAGS = -lm

LIBSRC =

# library = $(C_LIB_PATH)

LIBSRC =

CFLAGS += -D'PLOT'


all: artificial_matrix.exe artificial_matrix_generation_double.o artificial_matrix_generation_float.o ordered_set.o libartificial_matrix_generation_double.so libartificial_matrix_generation_float.so


artificial_matrix.exe: artificial_matrix.c artificial_matrix_generation.c ordered_set.c $(LIBSRC)
	$(CC) $(CFLAGS) -D'VERBOSE' $^ -o $@ $(LDFLAGS)


libartificial_matrix_generation_double.so: artificial_matrix_generation.c ordered_set.o
	$(CC) $(CFLAGS) -D'ValueType=double' $(SOFLAGS) $^ -o $@ $(LDFLAGS)
libartificial_matrix_generation_float.so: artificial_matrix_generation.c ordered_set.o
	$(CC) $(CFLAGS) -D'ValueType=float' $(SOFLAGS) $^ -o $@ $(LDFLAGS)


artificial_matrix_generation_double.o: artificial_matrix_generation.c
	$(CC) $(CFLAGS) -D'ValueType=double' $^ -c -o $@ $(LDFLAGS)
artificial_matrix_generation_float.o: artificial_matrix_generation.c
	$(CC) $(CFLAGS) -D'ValueType=float' $^ -c -o $@ $(LDFLAGS)

ordered_set.o: ordered_set.c
	$(CC) $(CFLAGS) $^ -c -o $@ $(LDFLAGS)




clean: 
	$(RM) *.exe a.out *.so *.o


