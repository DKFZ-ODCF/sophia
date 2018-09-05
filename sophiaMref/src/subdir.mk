################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Alignment.cpp \
../src/AnnotationProcessor.cpp \
../src/Breakpoint.cpp \
../src/BreakpointReduced.cpp \
../src/ChosenBp.cpp \
../src/ChrConverter.cpp \
../src/DeFuzzier.cpp \
../src/GermlineMatch.cpp \
../src/MasterRefProcessor.cpp \
../src/MrefEntry.cpp \
../src/MrefEntryAnno.cpp \
../src/MrefMatch.cpp \
../src/SamSegmentMapper.cpp \
../src/Sdust.cpp \
../src/SuppAlignment.cpp \
../src/SuppAlignmentAnno.cpp \
../src/SvEvent.cpp 

OBJS += \
./src/Alignment.o \
./src/AnnotationProcessor.o \
./src/Breakpoint.o \
./src/BreakpointReduced.o \
./src/ChosenBp.o \
./src/ChrConverter.o \
./src/DeFuzzier.o \
./src/GermlineMatch.o \
./src/MasterRefProcessor.o \
./src/MrefEntry.o \
./src/MrefEntryAnno.o \
./src/MrefMatch.o \
./src/SamSegmentMapper.o \
./src/Sdust.o \
./src/SuppAlignment.o \
./src/SuppAlignmentAnno.o \
./src/SvEvent.o 

CPP_DEPS += \
./src/Alignment.d \
./src/AnnotationProcessor.d \
./src/Breakpoint.d \
./src/BreakpointReduced.d \
./src/ChosenBp.d \
./src/ChrConverter.d \
./src/DeFuzzier.d \
./src/GermlineMatch.d \
./src/MasterRefProcessor.d \
./src/MrefEntry.d \
./src/MrefEntryAnno.d \
./src/MrefMatch.d \
./src/SamSegmentMapper.d \
./src/Sdust.d \
./src/SuppAlignment.d \
./src/SuppAlignmentAnno.d \
./src/SvEvent.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++1y -I"/home/umuttoprak/cppProjectsCevelop/sophia/include" -O3 -Wall -c -fmessage-length=0 -static -flto -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


