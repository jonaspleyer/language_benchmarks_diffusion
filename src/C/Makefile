TARGET:=diffusion.c
OUTPUT:=diffusion
OPTIONS:=-O3

gcc:
	gcc $(OPTIONS) $(TARGET) -o $(OUTPUT)_gcc

clang:
	clang $(OPTIONS) $(TARGET) -o $(OUTPUT)_clang

clean:
	-rm $(OUTPUT)_gcc
	-rm $(OUTPUT)_clang
