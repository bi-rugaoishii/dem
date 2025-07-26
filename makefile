CC = g++
CFLAGS = -O3  -g -mtune=native -march=native -fopenmp -pg
LIBS = -lm 
OBJS = main.o particleGroup.o demCalc.o readFiles.o particle.o boundingBox.o triangles.o fileOutput.o timeAdvTest.o
PROGRAM = myDEM3d

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(PROGRAM)

%.o : %.cpp
	${CC} ${CFLAGS} -c $<

clean:
	rm -f $(PROGRAM)  $(OBJS)
