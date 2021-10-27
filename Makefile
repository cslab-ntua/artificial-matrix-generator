.PHONY: clean


CC = gcc


CFLAGS =
CFLAGS += -Wall -Wextra
CFLAGS += -fopenmp 

# CFLAGS += -O0
CFLAGS += -O2


LDFLAGS = -lm


LIBSRC =

ifdef C_LIB_PATH
    library = ${C_LIB_PATH}

    LIBSRC = $(library)/plot/plot.c $(library)/parallel_io.c

    CFLAGS += -I$(library) -D'PLOT'
endif


all: artificial_matrix.exe artificial_matrix_generation_double.o artificial_matrix_generation_float.o sorted_set.o


artificial_matrix.exe: artificial_matrix.c artificial_matrix_generation.c sorted_set.c $(LIBSRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)


artificial_matrix_generation_double.o: artificial_matrix_generation.c
	$(CC) $(CFLAGS) -D'ValueType=double' $^ -c -o $@ $(LDFLAGS)
artificial_matrix_generation_float.o: artificial_matrix_generation.c
	$(CC) $(CFLAGS) -D'ValueType=float' $^ -c -o $@ $(LDFLAGS)
sorted_set.o: sorted_set.c
	$(CC) $(CFLAGS) $^ -c -o $@ $(LDFLAGS)




clean: 
	$(RM) *.exe a.out *.so *.o


