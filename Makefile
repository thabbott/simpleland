SHELL=/bin/bash
CC=gcc -c
LD=gcc
INCLUDE = -I./include
LINK = -lm
CFLAGS = -fPIC -g
LDFLAGS = -shared

TARGETS = lib/libland.so py/land.so
SOURCES = $(shell echo src/*.c)
HEADERS = $(shell echo include/*.h)
OBJECTS = $(SOURCES:.c=.o)

all : $(TARGETS)

lib/libland.so: $(OBJECTS)
	$(LD) $(CFLAGS) $(INCLUDE) $(OBJECTS) -o $@ $(LDFLAGS) $(LINK)

py/land.so:
	python setup.py build_ext
	python setup.py install --prefix=./py
	rm -rf build py/src/*.c

%.o : %.c
	$(CC) $(INCLUDE) $(CFLAGS) $< -o $@

clean :
	rm -rf src/*.o lib/* py/lib/*
