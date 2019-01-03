################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/bamlib.cpp 

OBJS += \
./src/bamlib.o 

CPP_DEPS += \
./src/bamlib.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -fPIC -std=c++14 -DSEQAN_HAS_ZLIB -DSEQAN_HAS_EXECINFO -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -I/usr/local/include -I/home/eliot/anaconda2/include -O3 -DNDEBUG -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


