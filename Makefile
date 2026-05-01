all: pocketfft_lib

pocketfft_lib:
# 	g++ -O3 -march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize -fopenmp-simd -DNDEBUG -std=c++11 -fPIC -pthread -c pocketfft_c.cc -o pocketfft_c.o
	g++ -O3 -march=native -mtune=native -funroll-loops -ftree-vectorize \
        -fopenmp-simd -DNDEBUG -std=c++11 -fPIC -pthread -c ./pocketfft/pocketfft_c.cc \
        -o ./bin/pocketfft_c.o
	ar rcs bin/libpocketfft_c.a bin/pocketfft_c.o

test:
	odin build . -o:speed -extra-linker-flags:"-lstdc++ -lpthread -lm" -out:pocketfft_test.exe

run:
	./pocketfft_test.exe

clean:
	rm -f bin/libpocketfft_c.a bin/pocketfft_c.o pocketfft_test.exe
