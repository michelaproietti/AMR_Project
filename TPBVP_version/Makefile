CFLAGS =  -I./build -Wall -I/usr/include/eigen3 -g -fPIC -std=c++11 
LDFLAGS = -lifopt_core -lifopt_ipopt 

OBJS = build/v_repExtAMRProject.o build/v_repLib.o

OS = $(shell uname -s)
ECHO=@

ifeq ($(OS), Linux)
	CFLAGS += -D__linux
	OPTION = -shared
	EXT = so
else
	CFLAGS += -D__APPLE__
	OPTION = -dynamiclib -current_version 1.0
	EXT = dylib
endif

TARGET = ./bin/libv_repExtAMRProject.$(EXT)

# change the following line according to the path of your V-REP main folder 
VREP_HOME = /home/sofia/Scaricati/CoppeliaSim_Edu_V4_0_0_Ubuntu18_04

default: v_repExtAMRProjectLib
	cp ./bin/libv_repExtAMRProject.${EXT} ${VREP_HOME}

v_repExtAMRProjectLib: $(OBJS)
		@echo "Linking $(OBJS) to $(TARGET)"
		$(ECHO)$(CXX) $(CFLAGS) $(OBJS) $(OPTION) -o $(TARGET) $(LDFLAGS)

%.o: %.c
		@echo "Compiling $< to $@"
		$(ECHO)$(CXX) $(CFLAGS) -c $< -o $@

./build/%.o: ./%.cpp
		@echo "Compiling $< to $@"
		$(ECHO)$(CXX) $(CFLAGS) -c $< -o $@

clean:
		@echo "Cleaning $(OBJS) $(TARGET)"
		$(ECHO)rm -rf $(OBJS) $(TARGET)
