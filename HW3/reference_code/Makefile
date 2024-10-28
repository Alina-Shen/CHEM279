CXX = g++

# Compilation flags: -Wall for compiler warnings, -g for debugging information
CXXFLAGS =

# Include directories for header files
INCLUDES = #-I/usr/include on Datahub

# Library paths, if Armadillo is not in a standard location
LFLAGS = #-L/usr/lib on Datahub

# Libraries to link against, including Armadillo
LIBS = -larmadillo

# Source files
SRCS = main.cpp EH.cpp AO.cpp util.cpp

# Object files, derived from SRCS
OBJS = $(SRCS:.cpp=.o)

# Executable name
MAIN = hw3

# Targets
.PHONY: depend clean

all:    $(MAIN)
	@echo  hw3 has been compiled

$(MAIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# Rule to convert .cpp to .o
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)
