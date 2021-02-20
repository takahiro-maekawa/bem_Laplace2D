#	Makefile
CC = g++
TARGETS = clean

FILES := geometry_common geometry_sample integralGL general_func core_matrix solve main
Appname = app
LIBS = -I/usr/include/eigen3

MAIN: $(Appname) clean
	@printf "#\n# the App has been built!!! \n#\n"

$(Appname): CPP
	$(CC) -o $@ $(addsuffix .o,$(FILES))

CPP:
	$(CC) -c $(addsuffix .cpp,$(FILES)) -std=c++17 $(LIBS)

clean: clean_msg
	rm -f *.o

clean_msg:
	@printf "#\n# cleaning obj files: \n#\n"
