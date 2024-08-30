CC = gcc -g
CFLAGS = -O2 -march=native -fopenmp -Wall -Wno-unused-variable -std=gnu99
LDFLAGS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -L/home/boye4060/projects/NRml_TOV/GRHayL/build/lib -lhdf5 -lgrhayl

INCLUDEDIRS= -I/usr/include/hdf5/serial -I/home/boye4060/projects/NRml_TOV/GRHayL/build/include
all: NRml_TOV NRml_TOV_Driver_main.o
NRml_TOV_Driver_main.o: NRml_TOV_Driver_main.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c NRml_TOV_Driver_main.c -o NRml_TOV_Driver_main.o $(LDFLAGS)

NRml_TOV:	NRml_TOV_Driver_main.o
	$(CC)	NRml_TOV_Driver_main.o -o NRml_TOV $(LDFLAGS) -lm -lgomp

clean:
	rm -f *.o */*.o *~ */*~ ./#* *.dat *.avi *.png NRml_TOV
