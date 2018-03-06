CC		= g++
CFLAGS  = -O3 -std=c++11 -I../third_party/eigen
LDFLAGS = 

TARGET = mdsim
SRCS	= main.cpp md_utils.cpp file_utils.cpp
OBJS	= main.o md_utils.o file_utils.o

RM		= rm -f


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

main.o: main.cpp file_utils.hpp md_utils.hpp
	$(CC) $(CFLAGS) -c $<

md_utils.o: md_utils.cpp md_utils.hpp file_utils.hpp
	$(CC) $(CFLAGS) -c $<

file_utils.o: file_utils.cpp file_utils.hpp
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) $(OBJS)
