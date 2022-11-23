all: rng.o api.o sha256.o dme-ref.o dme-opt.o kat-ref kat-opt kat.txt

rng.o: rng.c rng.h
	gcc -c rng.c -Wall -Wextra -Werror -pedantic -std=gnu99 -O3

api.o: api.c api.h sha256.h dme.h
	gcc -c api.c -Wall -Wextra -Werror -pedantic -std=gnu99 -O3

sha256.o: sha256.c sha256.h
	gcc -c sha256.c -Wall -Wextra -Werror -pedantic -std=gnu99 -O3

dme-ref.o: dme.c dme.h rng.h
	gcc -o dme-ref.o -c dme.c -Wall -Wextra -Werror -pedantic -std=gnu99 -O3 -DREFERENCE

dme-opt.o: dme.c dme.h rng.h
	gcc -o dme-opt.o -c dme.c -Wall -Wextra -Werror -pedantic -std=gnu99 -O3

kat-ref: kat.c dme-ref.o rng.o api.o sha256.o
	gcc -o kat-ref kat.c dme-ref.o rng.o api.o sha256.o -Wall -Wextra -Werror -pedantic -std=c99 -O3

kat-opt: kat.c dme-opt.o rng.o api.o sha256.o
	gcc -o kat-opt kat.c dme-opt.o rng.o api.o sha256.o -Wall -Wextra -Werror -pedantic -std=c99 -O3

kat.txt: kat-opt
	./kat-opt > kat.txt

clean:
	rm -f *.o
	rm -f kat-ref
	rm -f kat-opt
	rm -f kat.txt

.PHONY: all clean
