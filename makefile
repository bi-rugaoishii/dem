CC = g++
CFLAGS = -O3  -g -mtune=native -march=native
LIBS = -lm
OBJS = main.cpp
PROGRAM = myDEM3d

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(PROGRAM) -L./ -lmylib

clean:
	rm -f *.o *~ $(PROGRAM)
