OBJS = fastcorr.o higher_orders.o

CC = gcc
ifdef notshared
ifeq ($(notshared),yes)
EXEC = main.exe
CFLAGS =
OFLAGS =
endif
else
EXEC = _fastcorr.so
CFLAGS = -fPIC
OFLAGS = -shared
endif

INCL = -I${GSLI} -O3
LIBS = -lgsl -lgslcblas -L${GSLL} -lm -O3
.SUFFIXES : .c .o
%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

#$(OBJS): $(INCL)

.PHONY : clean

clean:
	rm -f $(OBJS) *.pyc main.exe _fastcorr.so
