TARGET:=threetubekin
OPTFLAGS:=-O2 -Wall
CFLAGS:=-I$(OH_SDK_BASE)/libsrc/include/HDU -I$(OH_SDK_BASE)/libsrc/include/HLU $(OPTFLAGS) -Dlinux -DNDEBUG
CXXFLAGS:=$(CFLAGS) -std=c++11
LDFLAGS:=-flto -L$(OH_SDK_BASE)/libsrc/HDU -L$(OH_SDK_BASE)/libsrc/HLU -lHDU -lHD -lrt -lm -lstdc++ -lpthread 
CC:=gcc
CXX:=g++
LD:=g++
HEADERS:=*.h
OBJS:=Rk8_coeffs.o

%.o: %.c $(HEADERS)
	$(CC) -c $(CFLAGS) $< -o $@

%.o: %.cpp $(HEADERS)
	$(CXX) -c $(CXXFLAGS) $< -o $@

default: $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $(TARGET)

.PHONY: clean
clean:
	rm -f $(OBJS) OmniDriver
