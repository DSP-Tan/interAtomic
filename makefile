CC = gcc
CFLAGS = -Wall
objects = BondsTwoChecks.o makebonds.o AtomTypes.o SkipToCommand.o angles.o
all: mai

mai: $(objects)
	$(CC) $(CFLAGS) -O2 -o mai $(objects) -lm

BondsTwoChecks.o: BondsTwoChecks.c funcs.h
	$(CC) $(CFLAGS) -c BondsTwoChecks.c

AtomTypes.o: AtomTypes.c funcs.h
	$(CC) $(CFLAGS) -c AtomTypes.c

makebonds.o: makebonds.c funcs.h
	$(CC) $(CFLAGS) -c makebonds.c

angles.o: angles.c funcs.h
	$(CC) $(CFLAGS) -c angles.c


SkipToCommand.o: SkipToCommand.c
	$(CC) $(CFLAGS) -c SkipToCommand.c

clean:
	rm *.o *~


	

