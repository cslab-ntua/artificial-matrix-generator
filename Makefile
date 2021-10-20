.PHONY: clean


library = ${C_LIB_PATH}


CC = gcc


CFLAGS =
CFLAGS += -Wall -Wextra
CFLAGS += -fopenmp 
CFLAGS += -I$(library)

# CFLAGS += -O0
CFLAGS += -O2


LDFLAGS = -lm


# PLOT = 1

ifdef PLOT
    LIBSRC = $(library)/plot/plot.c $(library)/parallel_io.c
    CFLAGS += -D'PLOT'
else
    LIBSRC =
endif


all: artificial_matrix.exe


artificial_matrix.exe: artificial_matrix.c artificial_matrix_generation.c  sorted_set.c $(LIBSRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)


clean: 
	$(RM) *.exe a.out *.so


