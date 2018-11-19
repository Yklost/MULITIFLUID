CC=mpicxx
CFLAGS=-O3 -std=gnu++11
SRC= memory.cpp apply_BC.cpp calc_RHS.cpp physics.cpp initial_condition.cpp main.cpp numerical_flux.cpp mpiio.cpp mpigeom.cpp collisions.cpp
OBJ=${SRC:.cpp=.o}

%.o: %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

code: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)
	rm data_*

clean:
	rm *.o rti_fv data_*
        
