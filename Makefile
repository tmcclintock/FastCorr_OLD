OBJS = src/fastcorr.o

CC = gcc
ifdef notshared
ifeq ($(notshared),yes)
EXEC = main.exe
CFLAGS =
OFLAGS =
endif
else
EXEC = src/c_fastcorr.so
CFLAGS = -fPIC
OFLAGS = -shared
endif

INCL = -I/home/tmcclintock/code/gsl/include/ -O3
LIBS = -lgsl -lgslcblas -L/home/tmcclintock/code/gsl/lib -lm -O3
.SUFFIXES : .c .o
%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

#$(OBJS): $(INCL)

.PHONY : clean

clean:
	rm -f $(OBJS) main.exe src/c_fastcorr.so
	rm -f src/*~ src/*.pyc
