#Makefile
CFLAGS += $(shell root-config --cflags)
LDFLAGS += $(shell root-config --libs --ldflags)

SRCS_H=$(wildcard *.cxx)
SRCS_M=$(wildcard *.c)
OBJS=$(SRCS_H:.cxx= .o)
EXEC=$(SRCS_M:.c=.exe)

all: $(EXEC) $(OBJS)

%.o: %.cxx
	g++ -c $(CFLAGS) $< -o $@
%.exe: %.c $(OBJS)
	g++ $(CFLAGS) $(LDFLAGS) $< -o $@ $(OBJS)

clean:
	rm $(EXEC)
