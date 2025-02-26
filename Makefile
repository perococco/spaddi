all:
	cd src && ./config.sh && make all

clean:
	cd src && test -e Makefile && make clean

clean_all:
	cd src && test -e Makefile && make clean_all
