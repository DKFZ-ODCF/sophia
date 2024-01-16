# Compiler
CXX = x86_64-conda_cos6-linux-gnu-g++

INCLUDE_DIR = ./include
BUILD_DIR = ./build
SRC_DIR = ./src
TESTS_DIR = ./tests

# Compiler flags
LIBRARY_FLAGS := -lz -lm -lrt -lboost_system -lboost_iostreams -lboost_program_options -ldl
LDFLAGS := $(LDFLAGS) -flto=auto -rdynamic -no-pie
CXXFLAGS := -I$(INCLUDE_DIR) $(CXXFLAGS) -std=c++17 -flto=auto -Wall -Wextra -c -fmessage-length=0 -Wno-attributes

ifeq ($(static),true)
	LD_BEGIN_FLAGS := -L$(boost_lib_dir)
	LD_END_FLAGS := $(LDFLAGS) -static -static-libgcc -static-libstdc++ $(LIBRARY_FLAGS)
else
    LD_BEGIN_FLAGS :=
	LD_END_FLAGS := $(LDFLAGS) $(LIBRARY_FLAGS)
endif

ifeq ($(develop),true)
	# NOTE: Generally, it is a good idea to compile with -O0 during development, because it seems
	#       that thus the compiler actually catches some binary dependencies during linking that
	#       will otherwise be missed.
	CXXFLAGS := $(CXXFLAGS) -O0 -ggdb3 -DDEBUG -fno-inline
	LD_END_FLAGS := $(LD_END_FLAGS) -Wl,-O0 -ggdb3 -DDEBUG -fno-inline
else
	# Ignore some leftover unused variables from SvEvent::assessBreakpointClonalityStatus.
	CXXFLAGS := $(CXXFLAGS) -O3 -DNDEBUG -Wno-unused-variable
endif

# Source files
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)

# Test source files
TESTS = $(wildcard $(TESTS_DIR)/*.cpp)

# Object files should have .o instead of .cpp.
# Note, we put the objects for production and tests both into the build directory.
OBJECTS = $(SOURCES: $(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o) $(BUILD_DIR)/strtk.o

# Binaries
BINARIES = sophiaMref sophia sophiaAnnotate

# Default rule
all: test $(BINARIES)

vpath %.h $(INCLUDE_DIR)
vpath %.hpp $(INCLUDE_DIR)
vpath %.cpp $(SRC_DIR) $(TESTS_DIR)


$(BUILD_DIR):
	mkdir -p $@

# Retrieve and build strtk.
$(INCLUDE_DIR)/strtk.hpp:
	wget -c https://raw.githubusercontent.com/ArashPartow/strtk/d2b446bf1f7854e8b08f5295ec6f6852cae066a2/strtk.hpp -O $(INCLUDE_DIR)/strtk.hpp

# General compilation rule for object files that have matching .h files.
$(BUILD_DIR)/%.o: %.cpp %.h $(INCLUDE_DIR)/strtk.hpp Makefile | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rules for sophia
$(BUILD_DIR)/sophia.o: $(SRC_DIR)/sophia.cpp Makefile | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@
sophia: $(BUILD_DIR)/global.o \
		$(BUILD_DIR)/Alignment.o \
		$(BUILD_DIR)/Breakpoint.o \
		$(BUILD_DIR)/ChosenBp.o \
		$(BUILD_DIR)/ChrConverter.o \
		$(BUILD_DIR)/Hg37ChrConverter.o \
		$(BUILD_DIR)/Hg38ChrConverter.o \
		$(BUILD_DIR)/SamSegmentMapper.o \
		$(BUILD_DIR)/Sdust.o \
		$(BUILD_DIR)/SuppAlignment.o \
		$(BUILD_DIR)/HelperFunctions.o \
		$(BUILD_DIR)/GlobalAppConfig.o \
		$(BUILD_DIR)/sophia.o
	$(CXX) $(LD_BEGIN_FLAGS) -o $@ $^ $(LD_END_FLAGS)

# Rules for sophiaAnnotate
$(BUILD_DIR)/sophiaAnnotate.o: $(SRC_DIR)/sophiaAnnotate.cpp Makefile | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@
sophiaAnnotate: $(BUILD_DIR)/global.o \
				$(BUILD_DIR)/Alignment.o \
				$(BUILD_DIR)/AnnotationProcessor.o \
				$(BUILD_DIR)/Breakpoint.o \
				$(BUILD_DIR)/BreakpointReduced.o \
				$(BUILD_DIR)/ChrConverter.o \
				$(BUILD_DIR)/Hg37ChrConverter.o \
				$(BUILD_DIR)/Hg38ChrConverter.o \
				$(BUILD_DIR)/DeFuzzier.o \
				$(BUILD_DIR)/GermlineMatch.o \
				$(BUILD_DIR)/MrefEntry.o \
				$(BUILD_DIR)/MrefEntryAnno.o \
				$(BUILD_DIR)/MrefMatch.o \
				$(BUILD_DIR)/SuppAlignment.o \
				$(BUILD_DIR)/SuppAlignmentAnno.o \
				$(BUILD_DIR)/Sdust.o \
				$(BUILD_DIR)/ChosenBp.o \
				$(BUILD_DIR)/SvEvent.o \
				$(BUILD_DIR)/HelperFunctions.o \
				$(BUILD_DIR)/GlobalAppConfig.o \
				$(BUILD_DIR)/sophiaAnnotate.o
	$(CXX) $(LD_BEGIN_FLAGS) -o $@ $^ $(LD_END_FLAGS)

# Rules for sophiaMref
$(BUILD_DIR)/sophiaMref.o: $(SRC_DIR)/sophiaMref.cpp Makefile | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@
sophiaMref: $(BUILD_DIR)/global.o \
			$(BUILD_DIR)/Alignment.o \
			$(BUILD_DIR)/GlobalAppConfig.o \
			$(BUILD_DIR)/ChrConverter.o \
			$(BUILD_DIR)/Hg37ChrConverter.o \
			$(BUILD_DIR)/Hg38ChrConverter.o \
			$(BUILD_DIR)/HelperFunctions.o \
			$(BUILD_DIR)/SuppAlignment.o \
			$(BUILD_DIR)/SuppAlignmentAnno.o \
			$(BUILD_DIR)/Sdust.o \
			$(BUILD_DIR)/ChosenBp.o \
			$(BUILD_DIR)/MrefEntry.o \
			$(BUILD_DIR)/MrefEntryAnno.o \
			$(BUILD_DIR)/MrefMatch.o \
			$(BUILD_DIR)/MasterRefProcessor.o \
			$(BUILD_DIR)/Breakpoint.o \
			$(BUILD_DIR)/BreakpointReduced.o \
			$(BUILD_DIR)/GermlineMatch.o \
			$(BUILD_DIR)/DeFuzzier.o \
			$(BUILD_DIR)/sophiaMref.o
	$(CXX) $(LD_BEGIN_FLAGS) -o $@ $^ $(LD_END_FLAGS)

# There are usually no .h files for test files, so we need a separate rule for them.
$(BUILD_DIR)/%.o: $(TESTS_DIR)/%.cpp Makefile | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

testRunner: \
		$(BUILD_DIR)/global.o \
		$(BUILD_DIR)/ChrConverter.o \
		$(BUILD_DIR)/Hg38ChrConverter.o \
		$(BUILD_DIR)/Hg38ChrConverter_test.o \
		$(BUILD_DIR)/SuppAlignment.o \
		$(BUILD_DIR)/GlobalAppConfig.o \
		$(BUILD_DIR)/SuppAlignment_test.o
	$(CXX) $(LD_BEGIN_FLAGS) -o testRunner $^ $(LDFLAGS) $(LIBRARY_FLAGS) -Wl,-Bdynamic -lgtest -lgtest_main -pthread

# Rule for running the tests
test: testRunner
	./testRunner

# Rule for clean
.PHONY: clean clean-all
clean:
	rm -f $(BUILD_DIR)/*.o $(BINARIES)

clean-all: clean
	rm -f $(INCLUDE_DIR)/strtk.hpp
