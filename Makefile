INC_DIR = -I include
CC = g++

BDIR = build
SDIR = src

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	GL = -lglut -lGL
endif
ifeq ($(UNAME_S), Darwin)
	GL = -framework OpenGL -framework GLUT
	INC_DIR += -I/usr/local/include
endif

_OBJS = filehandler.o graph.o socialgraph.o drawing.o main.o
OBJS = $(patsubst %, $(BDIR)/%, $(_OBJS))

$(BDIR)/%.o : $(SDIR)/%.cpp
	$(CC) $(INC_DIR) -O2 -fopenmp -c -o $@ $<
	
main: $(OBJS)
	$(CC) -o main $(INC_DIR) $(OBJS) -fopenmp $(GL)

clean:
	rm $(OBJS) main
