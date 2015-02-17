CC1=gcc -Wall
CC2=g++ -Wall

all: clients

clients: test.c laserobstacleavoid.cc
	$(CC1) test.c -o test -I/usr/local/include/player-3.0 -L/usr/local/lib -lplayerc -lm
	$(CC2) laserobstacleavoid.cc -o laserobstacleavoid -I/usr/local/include/player-3.0 -L/usr/local/lib -lplayerc++ -lm

clean:
	@rm test laserobstacleavoid
