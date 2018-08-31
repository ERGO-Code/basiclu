include config.mk

#-------------------------------------------------------------------------------
# compile the package
#-------------------------------------------------------------------------------

# The default creates the static and shared library in the lib/ subdirectory
# and compiles the example program 'maxvolume'.

default: static shared example/maxvolume

#-------------------------------------------------------------------------------
# files
#-------------------------------------------------------------------------------

DEP_FILES = $(wildcard src/*.h) $(wildcard include/*.h) Makefile config.mk
SRC_FILES = $(wildcard src/*.c)
OBJ_FILES = $(patsubst src/%.c, build/%.o, $(SRC_FILES))

#-------------------------------------------------------------------------------
# create the static library
#-------------------------------------------------------------------------------

static: lib/$(AR_TARGET)

lib/$(AR_TARGET): $(OBJ_FILES)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $@

#-------------------------------------------------------------------------------
# create the shared library
#-------------------------------------------------------------------------------

shared: lib/$(SO_TARGET)

lib/$(SO_TARGET): $(OBJ_FILES)
	$(CC99) $(CF) $(SO_OPTS) $^ -o $@ $(LDLIBS)
	( cd lib; ln -sf $(SO_TARGET) $(SO_PLAIN) )
	( cd lib; ln -sf $(SO_TARGET) $(SO_MAIN) )

#-------------------------------------------------------------------------------
# compile example binary (use static linkage)
#-------------------------------------------------------------------------------

example/maxvolume: example/maxvolume.c example/mmio.c lib/$(AR_TARGET)
	$(CC99) $(CF) -I./include -o $@ $^ $(LDLIBS)

#-------------------------------------------------------------------------------
# compile each object file from its source file
#-------------------------------------------------------------------------------

build/%.o: src/%.c $(DEP_FILES)
	$(CC99) $(CF) -I./include -c $< -o $@

#-------------------------------------------------------------------------------
# clean and purge
#-------------------------------------------------------------------------------

clean:
	$(RM) $(OBJ_FILES)

purge: clean
	$(RM) lib/$(AR_TARGET) lib/$(SO_TARGET) lib/$(SO_PLAIN) lib/$(SO_MAIN)
	$(RM) example/maxvolume
