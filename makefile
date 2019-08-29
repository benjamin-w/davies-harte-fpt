CC = gcc
OPTIM = -O3 
CFLAGS += -Wall 

OBJFILES = dh_main.o dh_functions.o

LDFLAGS = -lfftw3 -lm -lgsl

TARGET = dh

$(TARGET): $(OBJFILES)  
	$(CC) $(OPTIM) $(LIBRARYPATHS) $(LDFLAGS) -o $@ $^ 

.c.o:
	$(CC) $(OPTIM) $(INCLUDEPATHS) $(CFLAGS)   -c -o $@ $^

clean:
	rm -f $(OBJFILES) $(TARGET) *~
