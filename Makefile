# 'make'	build executable file
# 'make clean'	removes all *.o and executalbe file

# define the C compiler
CC	= gcc

# define any compile-time flags
CFLAGS = -Wall -g -fPIC -O1 -std=c99

# define openmp flags
OPENMP  = -fopenmp
CUOPENMP  = -Xcompiler -fopenmp

# define the direction containing header file
INCLUDES= -I/usr/local/include -I./ -I./lattice

# define the library path
LFLAGS	= -L/usr/local/lib

# define any libraries to link to executable
LIBS	= -lm -lgsl -lgslcblas

# define the C object files
OBJS	= data_struct.o estimator.o  lattice/diluted_bilayer_heisenberg.o kernal/heisenberg_model.o

#define the directory for object
OBJSDIR = object

# define the executable file
MAIN	= exe

all: $(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(MAIN) *.o $(LIBS) $(LFLAGS) $(INCLUDES) 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $^

lib: $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $^



# clean the executable file and object files
clean:
	$(RM) *.o $(MAIN)
