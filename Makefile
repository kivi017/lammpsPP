# Declaration of variables
CC = g++
DEBUG_FLAG = #-g
OPTI_FLAGS = #-O3
WARN_FLAGS = 
CC_FLAGS = $(DEBUG_FLAG) $(OPTI_FLAGS) $(WARN_FLAGS) #-fsanitize=address,undefined  

# File names
EXEC = lammps_dump_read
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
HEADERS = $(wildcard *.hpp)
LIBFLAGS = #-lprofiler

all: $(EXEC)

# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(CC_FLAGS) $(OBJECTS) -o $(EXEC) $(LIBFLAGS)

# To obtain object files
%.o: %.cpp
	$(CC) -c $(CC_FLAGS) $< -o $@ 

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)
