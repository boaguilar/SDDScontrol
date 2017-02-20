main: SDDS_Main.o mt64.o read_files.o 
	g++ -std=c++11 SDDS_Main.o mt64.o read_files.o -o SDDS
mt64.o: mt64.c mt64.h
	g++ -std=c++11 -o mt64.o -c mt64.c
read_files.o: read_files.cpp read_files.h
	g++ -std=c++11 -o read_files.o -c read_files.cpp
SDDS_Main.o: SDDS_Main.cpp
	g++ -std=c++11 -o SDDS_Main.o -c SDDS_Main.cpp 
clean:
	rm -f *.o SDDS
