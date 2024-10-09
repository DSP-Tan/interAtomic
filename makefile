CC = gcc
CFLAGS = -Wall
objects = main.o makebonds.o AtomTypes.o SkipToCommand.o angles.o parseSupercell.o
all: mai

mai: $(objects)
	$(CC) $(CFLAGS) -O2 -o mai $(objects) -lm

main.o: main.c funcs.h
	$(CC) $(CFLAGS) -c main.c

AtomTypes.o: AtomTypes.c funcs.h
	$(CC) $(CFLAGS) -c AtomTypes.c

makebonds.o: makebonds.c funcs.h
	$(CC) $(CFLAGS) -c makebonds.c

angles.o: angles.c funcs.h
	$(CC) $(CFLAGS) -c angles.c

parseSupercell.o: parseSupercell.c funcs.h
	$(CC) $(CFLAGS) -c parseSupercell.c

SkipToCommand.o: SkipToCommand.c
	$(CC) $(CFLAGS) -c SkipToCommand.c

clean:
	rm *.o 



# gcc tests/testGetCellSize.c parseSupercell.c SkipToCommand.c -lm
