.PHONY: clean


CC = gcc


CFLAGS =
CFLAGS += -Wall -Wextra
CFLAGS += -fopenmp 

# CFLAGS += -O0
CFLAGS += -O2


SOFLAGS = -fPIC -shared


LDFLAGS = -lm


LIBSRC =

ifdef C_LIB_PATH
    library = ${C_LIB_PATH}

    LIBSRC = $(library)/plot/plot.c $(library)/parallel_io.c $(library)/string_util.c

    CFLAGS += -I$(library) -D'PLOT'
    # CFLAGS += -I$(library) -D'PLOT' -D'VERBOSE'
endif


all: artificial_matrix.exe artificial_matrix_generation_double.o artificial_matrix_generation_float.o ordered_set.o artificial_matrix_generation_double.so artificial_matrix_generation_float.so


artificial_matrix.exe: artificial_matrix.c artificial_matrix_generation.c ordered_set.c $(LIBSRC)
	$(CC) $(CFLAGS) -D'VERBOSE' $^ -o $@ $(LDFLAGS)


artificial_matrix_generation_double.so: artificial_matrix_generation.c ordered_set.o
	$(CC) $(CFLAGS) -D'ValueType=double' $(SOFLAGS) $^ -o $@ $(LDFLAGS)
artificial_matrix_generation_float.so: artificial_matrix_generation.c ordered_set.o
	$(CC) $(CFLAGS) -D'ValueType=float' $(SOFLAGS) $^ -o $@ $(LDFLAGS)


artificial_matrix_generation_double.o: artificial_matrix_generation.c
	$(CC) $(CFLAGS) -D'ValueType=double' $^ -c -o $@ $(LDFLAGS)
artificial_matrix_generation_float.o: artificial_matrix_generation.c
	$(CC) $(CFLAGS) -D'ValueType=float' $^ -c -o $@ $(LDFLAGS)

ordered_set.o: ordered_set.c
	$(CC) $(CFLAGS) $^ -c -o $@ $(LDFLAGS)




clean: 
	$(RM) *.exe a.out *.so *.o


