TARGET:=kintest
OPTFLAGS:=-O3
CFLAGS:=-I/usr/include/eigen3 $(OPTFLAGS) -Dlinux -DNDEBUG -g
CXXFLAGS:=$(CFLAGS) -std=c++11 -fsanitize=address
LDFLAGS:=-flto -fsanitize=address
CC:=gcc
CXX:=g++
LD:=g++
HEADERS:=*.h
OBJS:=Rk8_Coeffs.o\
		ReferenceCannulaKinematics.o\
		TaggedInterval.o

%.o: %.c $(HEADERS)
	$(CC) -c $(CFLAGS) $< -o $@

%.o: %.cpp $(HEADERS)
	$(CXX) -c $(CXXFLAGS) $< -o $@

default: $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $(TARGET)

.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)
