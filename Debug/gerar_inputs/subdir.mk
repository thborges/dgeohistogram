################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../gerar_inputs/gerar_inputs.c 

OBJS += \
./gerar_inputs/gerar_inputs.o 

C_DEPS += \
./gerar_inputs/gerar_inputs.d 


# Each subdirectory must supply rules for building sources it contributes
gerar_inputs/%.o: ../gerar_inputs/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I/usr/include/gdal -I"/home/venom/Documents/Dados Espaciais/dgeohistogram/auxs" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


