
HOME=/home/rodrigo/mycodes
CC=clang -Wall
CFLAGS=-march=native -O4
LIBS=-lm

DEPS = contact_map.h
OBJ = chemical_map.o pdb_map.o protein_map.o contact_map.o

contact_map: $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean install

clean:
	rm -f *.o
	rm contact_map
	rm $(HOME)/bin/contact_map

install: contact_map
	ln -s $(PWD)/contact_map $(HOME)/bin/contact_map
