
CC = g++
#CFLAGS = -ggdb -Wall -I/usr/include/g++-include
CFLAGS = -ggdb -Wall
LIBS = 

all: dpa_linear dpa_lcheckp ukk_linear ukk_lcheckp

ukk_lcheckp: ukk_lcheckp.o
	$(CC) ukk_lcheckp.o -o ukk_lcheckp $(LIBS)

ukk_linear: ukk_linear.o
	$(CC) ukk_linear.o -o ukk_linear $(LIBS)

dpa_linear: dpa_linear.o
	$(CC) dpa_linear.o -o dpa_linear $(LIBS)

dpa_lcheckp: dpa_lcheckp.o
	$(CC) dpa_lcheckp.o -o dpa_lcheckp $(LIBS)

tarball:
	./tar.pl align2str_linear_checkp.tar.gz align2str_linear_checkp README COPYRIGHT Makefile dpa_linear.cc dpa_lcheckp.cc ukk_linear.cc ukk_lcheckp.cc ukk_linear.h common.h

clean:
	rm -f *.o core dpa_linear dpa_lcheckp ukk_linear ukk_lcheckp align2str_linear_checkp.tar.gz

.cc.o:
	@echo ......................................Compiling $< to $@
	$(CC) -c $(CFLAGS) -o $@ $<

ukk_lcheckp.o: ukk_lcheckp.cc ukk_linear.h
ukk_linear.o: ukk_linear.cc ukk_linear.h

