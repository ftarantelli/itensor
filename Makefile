# 1. Put this file in the same folder as your 'driver' code 
#    (the code containing the 'main' function).

# 2. Edit LIBRARY_DIR to point at the location of your ITensor Library
#    source folder (this is the folder that has options.mk in it).
#    Also, edit TDVP_DIR to point at the location of the TDVP source files.
LIBRARY_DIR=$(HOME)/Desktop/C++/itensor
TDVP_DIR=$(LIBRARY_DIR)/tdvp

# 3. If your 'main' function is in a file called 'myappname.cc', then
#    set APP to 'myappname'. Running 'make' will compile the app.
#    Running 'make debug' will make a program called 'myappname-g'
#    which includes debugging symbols and can be used in gdb (Gnu debugger);
APP=tmpo

# 4. Add any headers your program depends on here. The make program
#    will auto-detect if these headers have changed and recompile your app.
HEADERS=$(TDVP_DIR)/tdvp.h $(PWD)/$(APP).hpp

# 5. For any additional .cc (source) files making up your project,
#    add their full filenames here.
CCFILES=$(APP).cpp

#################################################################
#################################################################
#################################################################
#################################################################


include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/core.h

CCFLAGS+=-I$(TDVP_DIR)
CCGFLAGS+=-I$(TDVP_DIR)

#Mappings --------------
OBJECTS=$(patsubst %.cpp,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------

%.o: %.cpp $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cpp $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(APP) $(LIBFLAGS)
	#./$(APP)

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)

clean:
	rm -fr .debug_objs *.o $(APP) $(APP)-g

mkdebugdir:
	mkdir -p .debug_objs
	
	
CC = g++ -g -Wall

foqt:
	$(CC) foqt.cpp -o foqt.out -larmadillo -llapack -lblas
	#./foqt.out

qtraj:
	$(CC) qtraj.cpp -o qtraj.out -larmadillo -llapack -lblas
	./qtraj.out

cleanall:
	rm foqt.out

data:
	rm data/*.dat

