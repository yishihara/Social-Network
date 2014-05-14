INC_DIR = -I include


globjects = filehandler.o graph.o socialgraph.o drawing.o main.o

vpath 

%.o : src/%.cpp
	g++ $(INC_DIR) -O2 -fopenmp -c src/$*.cpp
	
main: $(globjects)
	g++ -o main $(INC_DIR) $(globjects) -lglut -fopenmp

clean:
	rm $(globjects) main
