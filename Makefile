#	Makefile
CC = g++
TARGETS = clean

OBJSDP = geometry_common.o geometry_sample.o integralGL.o general_func.o core_matrix.o solve.o main.o
CPPSDP = geometry_common.cpp geometry_sample.cpp integralGL.cpp general_func.cpp core_matrix.cpp solve.cpp main.cpp
Appname = app
LIBS = -I/usr/include/eigen3

MAIN: $(Appname) clean
	@printf "#\n# the App has been built!!! \n#\n"

$(Appname): CPP
	$(CC) -o $@ $(OBJSDP)

CPP:
	$(CC) -c $(CPPSDP) -std=c++17 $(LIBS)

clean: clean_msg
	rm -f *.o

clean_msg:
	@printf "#\n# cleaning obj files: \n#\n"
