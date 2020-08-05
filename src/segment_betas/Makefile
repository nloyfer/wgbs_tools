CC=g++
CFLAGS=-std=c++11

ODIR=obj

DEPS = segmentor.h

OBJ = segmentor.o main.o

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ *.cpp $(CFLAGS)

all: segmentor

segmentor: $(OBJS)
	g++ -o $@ *.cpp $(CFLAGS) $(LIBS)


.PHONY: clean

clean:
	rm -f ./*.o ./segmentor
