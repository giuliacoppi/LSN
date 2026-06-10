CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS = -larmadillo

# L'eseguibile dipende solo da main.o e random.o
main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe $(LIBS)
	
# Se modifichi functions.h o random.h, main.cpp viene ricompilato automaticamente
main.o : main.cpp functions.h random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm -f *.o main.exe seed.out