#
# Makefile for core code
# Created by Dietrich Geisler
#

CC = gcc
CFLAGS = -I

all: core.c
	$(CC) -o a.out core.c -lm

clean:
	rm -rf *.o
	rm -f a.out