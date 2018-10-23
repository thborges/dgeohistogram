################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../dataset_specs.c \
../ehist.c \
../ehist_intermed.c \
../histogram.c \
../main.c \
../minskew.c \
../sthist.c 

O_SRCS += \
../dataset_specs.o \
../ehist.o \
../ehist_intermed.o \
../histogram.o \
../main.o \
../minskew.o \
../sthist.o 

OBJS += \
./dataset_specs.o \
./ehist.o \
./ehist_intermed.o \
./histogram.o \
./main.o \
./minskew.o \
./sthist.o 

C_DEPS += \
./dataset_specs.d \
./ehist.d \
./ehist_intermed.d \
./histogram.d \
./main.d \
./minskew.d \
./sthist.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I/usr/include/gdal -I"/home/venom/Documents/DadosEspaciais/dgeohistogram/auxs" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


