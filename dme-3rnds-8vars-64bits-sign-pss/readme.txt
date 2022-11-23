An implementation of the "dme-3rnds-8vars-64bits-sign-pss" signature
scheme developed by Ignacio Luengo (iluengo@ucm.es).

Files:

   * dme.h, dme.c: encryption, decryption and key-generation core functions,
       including a generic C99 implementation and a AVX2 optimized version.
   * api.h, api.c: API functions required by NIST.
   * kat.c: a simple test program that outputs the kat.txt file.
   * rng.h, rng.c: pseudorandom number generator (currently it is just a
       wrapper of the standard rand() function in libc, so it should never
       be used in production).
   * sha256.h, sha256.c: the reference implementation of the sha2-256 hash
       function.
   * makefile: the usual recipes to build everything.
   * readme.txt: this file.

To compile everything, you just need to run "make" on the folder that contains
all the files. This produces:

   * rng.o
   * dme-ref.o
   * dme-opt.o
   * sha256.o
   * api.o

   * kat-ref
   * kat-opt

   * kat.txt

The object files dme-ref.o and dme-opt.o correspond to the generic C99 and
AVX2 optimized version, respectively. The same for the kat-ref and kat-opt
executables, which can be used to generate the kat.txt file (in the makefile,
the optimized version is run to create kat.txt). The last lines of kat.txt
contain average timings for the crypto_sign(), crypto_sign_open() and
crypto_sign_keypair() API functions. Use "tail kat.txt -n 4" to see these
timings (in microseconds, as measured by the standard time() libc function).

To clean the folder, just run "make clean". This will remove the object
and executable files created by "make" and also the kat.txt file.

In case you just want the vanilla C99 version (maybe because you don't have
a CPU that supports the AVX2 instruction set), you can create kat-ref with
"make kat-ref" and the generate the kat.txt with "./kat-ref > kat.txt".
