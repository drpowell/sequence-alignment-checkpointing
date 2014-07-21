
CC = g++
CFLAGS = -ggdb -Wall
LIBS = 

all: ukk_checkp ukk_2str dpa_checkp dpa_2str makestrings

makestrings: makestrings.o
	$(CC) makestrings.o -o makestrings $(LIBS) -lm

ukk_checkp: ukk_checkp.o
	$(CC) ukk_checkp.o -o ukk_checkp $(LIBS)

ukk_2str: ukk_2str.o
	$(CC) ukk_2str.o -o ukk_2str $(LIBS)

dpa_checkp: dpa_checkp.o
	$(CC) dpa_checkp.o -o dpa_checkp $(LIBS)

dpa_2str: dpa_2str.o
	$(CC) dpa_2str.o -o dpa_2str $(LIBS)

tarball:
	./tar.pl align2str_checkp.tar.gz align2str_checkp README COPYRIGHT Makefile dpa_2str.cc dpa_checkp.cc ukk_2str.cc ukk_checkp.cc ukk_noalign.h


clean:
	rm -f *.o core ukk_checkp ukk_2str dpa_checkp dpa_2str makestrings align2str_checkp.tar.gz

.cc.o:
	@echo ......................................Compiling $< to $@
	$(CC) -c $(CFLAGS) -o $@ $<



