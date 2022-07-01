TARGET:=diffusion.c
OUTPUT:=diffusion
OPTIONS:=-O3

all: gcc clang bench

bench:
	./bench.sh

gcc:
	gcc $(OPTIONS) $(TARGET) -o $(OUTPUT)_gcc

clang:
	clang $(OPTIONS) $(TARGET) -o $(OUTPUT)_clang

clean:
	-rm $(OUTPUT)_gcc
	-rm $(OUTPUT)_clang
	-rm *.csv
	-rm *.png
	-rm *.txt
