# Compiler
CXX = x86_64-conda_cos6-linux-gnu-g++

INCLUDE_DIR = ./include
BUILD_DIR = ./build
SRC_DIR = ./src


INCLUDE_FLAGS= -Iinclude -I$$CONDA_PREFIX/include
LIBDIR_FLAGS = -L$$CONDA_PREFIX/lib

# Compiler flags
CXXFLAGS = -O3 -Wall -Wextra -static -static-libgcc -static-libstdc++ -flto  -c -fmessage-length=0 -Wno-attributes
ifeq ($(STATIC), "true")
	CXXFLAGS := -static -static-libgcc -static-libstdc++ $(CXXFLAGS)
endif
CXXFLAGS := $(LIBDIR_FLAGS) -std=c++17 $(INCLUDE_FLAGS) $(CXXFLAGS)

# Source files
SOURCES = $(wildcard src/*.cpp)

# Object files should have .o instead of .cpp.
OBJECTS = $(SOURCES:.cpp=.o)

# Binaries
BINARIES = sophia sophiaAnnotate sophiaMref

# Default rule
all: $(BINARIES)

$(BUILD_DIR):
	mkdir -p $@

# Rule for object files
$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(INCLUDE_FLAGS) $(LIBDIR_FLAGS) $(CXXFLAGS) -c $< -o $@

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
	$(CXX) $(LIBDIR_FLAGS) -lboost_program_options -o $@ $^

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
	$(CXX) $(LIBDIR_FLAGS) -lz -lboost_system -lboost_iostreams $(CXXFLAGS) -o $@ $^

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
	$(CXX) $(LIBDIR_FLAGS) -lz -lboost_system -lboost_iostreams -lboost_program_options $(CXXFLAGS) -o $@ $^

# Rule for clean
.PHONY: clean
clean:
	rm -f $(OBJECTS) $(BINARIES)
