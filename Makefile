#Makefile
CFLAGS += $(shell root-config --cflags)
LDFLAGS += $(shell root-config --libs --ldflags)

SRCS_H=$(wildcard *.cxx)
SRCS_M=$(wildcard *.C)
OBJS=$(SRCS_H:.cxx= .o)
EXEC=$(SRCS_M:.C=.exe)

all: $(EXEC) $(OBJS)
	echo $(SRCS_M) 

%.o: %.cxx
	g++ -c $(CFLAGS) $< -o $@
%.exe: %.C $(OBJS)
	g++ $(CFLAGS) $(LDFLAGS) $< -o $@ $(OBJS)

clean:
	rm $(EXEC)

