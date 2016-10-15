OBJS = src/fast_corr.o

CC = gcc
ifdef notshared
ifeq ($(notshared),yes)
EXEC = main.exe
CFLAGS =
OFLAGS =
endif
else
EXEC = src/c_fast_corr.so
CFLAGS = -fPIC
OFLAGS = -shared
endif

INCL = -I/home/tom/code/gsl/include/ -fopenmp -O3
LIBS = -lgsl -lgslcblas -L/home/tom/code/gsl/lib -lm -fopenmp -O3
.SUFFIXES : .c .o
%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

#$(OBJS): $(INCL)

.PHONY : clean

clean:
	rm -f $(OBJS) main.exe src/c_fast_corr.so
	rm -f *~
