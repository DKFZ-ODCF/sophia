# Compiler
CXX = x86_64-conda_cos6-linux-gnu-g++

INCLUDE_DIR = ./include
BUILD_DIR = ./build
SRC_DIR = ./src


INCLUDE_FLAGS= -Iinclude -I$$CONDA_PREFIX/include
LIBDIR_FLAGS = -L$$CONDA_PREFIX/lib

# Compiler flags
CXXFLAGS = $(LIBDIR_FLAGS) -std=c++17 -O3 -Wall -Wextra -flto  -c -fmessage-length=0 -Wno-attributes $(INCLUDE_FLAGS)
LDFLAGS = $(LIBDIR_FLAGS)
ifeq ($(STATIC), "true")
	# Note that static compilation needs static libraries for boost, which are not available in
	# conda!
	CXXFLAGS := $(CXXFLAGS) -static -static-libgcc -static-libstdc++
	LDFLAGS := $(LIBDIR_FLAGS) -static -static-libgcc -static-libstdc++
endif


# Source files
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)

# Object files should have .o instead of .cpp.
OBJECTS = $(SOURCES: $(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Binaries
BINARIES = sophia sophiaAnnotate sophiaMref

# Default rule
all: $(BINARIES)

$(BUILD_DIR):
	mkdir -p $@

# Rule for object files
$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(LIBDIR_FLAGS) $(CXXFLAGS) -c $< -o $@

download_strtk: include/strtk.hpp
	wget -c https://github.com/ArashPartow/strtk/raw/master/strtk.hpp -O include/strtk.hpp

vpath %.h includes
vpath %.cpp src

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
        $(BUILD_DIR)/sophia.o \
        download_strtk
	$(CXX) $(LDFLAGS) -lboost_program_options -o $@ $^

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
				$(BUILD_DIR)/sophiaAnnotate.o \
				download_strtk
	$(CXX) $(LDFLAGS) -lz -lboost_system -lboost_iostreams -o $@ $^

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
			$(BUILD_DIR)/sophiaMref.o \
			download_strtk
	$(CXX) $(LDFLAGS) -lz -lboost_system -lboost_iostreams -lboost_program_options -o $@ $^

# Rule for clean
.PHONY: clean
clean:
	rm -f $(BUILD_DIR)/*.o $(BINARIES)
