# Compiler
CXX = x86_64-conda_cos6-linux-gnu-g++

INCLUDE_DIR = ./include
BUILD_DIR = ./build
SRC_DIR = ./src

# Compiler flags
LDFLAGS := $(LDFLAGS) -flto=auto -lz -lboost_system -lboost_iostreams -lboost_program_options
CXXFLAGS := -I$(INCLUDE_DIR) $(CXXFLAGS) -std=c++17 -O3 -flto=auto -Wall -Wextra -c -fmessage-length=0 -Wno-attributes
ifeq ($(STATIC),true)
	CXXFLAGS := $(CXXFLAGS) -static -static-libgcc -static-libstdc++
	LD_BEGIN_FLAGS := -L$(boost_lib_dir)
	# Note that without the additional -lrt flag, the static linking fails for sophiaMref.
	LD_END_FLAGS := $(LDFLAGS) -static -static-libgcc -static-libstdc++ -lrt
else
    LD_BEGIN_FLAGS :=
	LD_END_FLAGS := $(LDFLAGS)
endif


# Source files
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)

# Object files should have .o instead of .cpp.
OBJECTS = $(SOURCES: $(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o) $(BUILD_DIR)/strtk.o

# Binaries
BINARIES = sophia sophiaAnnotate sophiaMref

# Default rule
all: $(BINARIES)

vpath %.h $(INCLUDE_DIR)
vpath %.hpp $(INCLUDE_DIR)
vpath %.cpp $(SRC_DIR)


$(BUILD_DIR):
	mkdir -p $@

# Retrieve and build strtk.
$(INCLUDE_DIR)/strtk.hpp:
	wget -c https://raw.githubusercontent.com/ArashPartow/strtk/d2b446bf1f7854e8b08f5295ec6f6852cae066a2/strtk.hpp -O $(INCLUDE_DIR)/strtk.hpp

#$(BUILD_DIR)/strtk.o: $(INCLUDE_DIR)/strtk.hpp | $(BUILD_DIR)
#	$(CXX) $(LIBDIR_FLAGS) $(CXXFLAGS) -c $(INCLUDE_DIR)/strtk.hpp -o $@

# General compilation rule for object files.
$(BUILD_DIR)/%.o: %.cpp $(INCLUDE_DIR)/strtk.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule for sophia
sophia: $(BUILD_DIR)/Alignment.o \
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

# Rule for sophiaAnnotate
sophiaAnnotate: $(BUILD_DIR)/AnnotationProcessor.o \
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
				$(BUILD_DIR)/SvEvent.o \
				$(BUILD_DIR)/HelperFunctions.o \
				$(BUILD_DIR)/GlobalAppConfig.o \
				$(BUILD_DIR)/sophiaAnnotate.o
	$(CXX) $(LD_BEGIN_FLAGS) -o $@ $^ $(LD_END_FLAGS)

# Rule for sophiaMref
sophiaMref: $(BUILD_DIR)/GlobalAppConfig.o \
			$(BUILD_DIR)/ChrConverter.o \
			$(BUILD_DIR)/Hg37ChrConverter.o \
			$(BUILD_DIR)/Hg38ChrConverter.o \
			$(BUILD_DIR)/HelperFunctions.o \
			$(BUILD_DIR)/SuppAlignment.o \
			$(BUILD_DIR)/SuppAlignmentAnno.o \
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

# Rule for clean
.PHONY: clean clean-all
clean:
	rm -f $(BUILD_DIR)/*.o $(BINARIES)

clean-all: clean
	rm -f $(INCLUDE_DIR)/strtk.hpp
