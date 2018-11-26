################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../auxs/geosext.cpp \
../auxs/ogrext.cpp 

C_SRCS += \
../auxs/dataset.c \
../auxs/glibwrap.c \
../auxs/lru-buffer.c \
../auxs/rtree-distributor.c \
../auxs/rtree-lazy.c \
../auxs/rtree-reinsert.c \
../auxs/rtree-star.c \
../auxs/rtree-wquery.c \
../auxs/rtree.c \
../auxs/utils.c \
../auxs/wkbconvert.c 

O_SRCS += \
../auxs/dataset.o \
../auxs/geosext.o \
../auxs/glibwrap.o \
../auxs/lru-buffer.o \
../auxs/ogrext.o \
../auxs/rtree-distributor.o \
../auxs/rtree-lazy.o \
../auxs/rtree-reinsert.o \
../auxs/rtree-star.o \
../auxs/rtree-wquery.o \
../auxs/rtree.o \
../auxs/utils.o \
../auxs/wkbconvert.o 

OBJS += \
./auxs/dataset.o \
./auxs/geosext.o \
./auxs/glibwrap.o \
./auxs/lru-buffer.o \
./auxs/ogrext.o \
./auxs/rtree-distributor.o \
./auxs/rtree-lazy.o \
./auxs/rtree-reinsert.o \
./auxs/rtree-star.o \
./auxs/rtree-wquery.o \
./auxs/rtree.o \
./auxs/utils.o \
./auxs/wkbconvert.o 

CPP_DEPS += \
./auxs/geosext.d \
./auxs/ogrext.d 

C_DEPS += \
./auxs/dataset.d \
./auxs/glibwrap.d \
./auxs/lru-buffer.d \
./auxs/rtree-distributor.d \
./auxs/rtree-lazy.d \
./auxs/rtree-reinsert.d \
./auxs/rtree-star.d \
./auxs/rtree-wquery.d \
./auxs/rtree.d \
./auxs/utils.d \
./auxs/wkbconvert.d 


# Each subdirectory must supply rules for building sources it contributes
auxs/%.o: ../auxs/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -DHAVE_ISNAN=1 -I"/home/kronos/Documents/DadosEspaciais/dgeohistogram/auxs" -I/usr/include/gdal -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

auxs/%.o: ../auxs/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DHAVE_ISNAN=1 -I/usr/include/gdal -I"/home/kronos/Documents/DadosEspaciais/dgeohistogram/auxs" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


