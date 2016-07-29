CC=gcc
CFLAGS= -std=c99 -pedantic -Wall -g3

#####
# Instructions to make hw2
#####

hw2: hw2.o
	${CC} ${CFLAGS} -o hw2 hw2.o
