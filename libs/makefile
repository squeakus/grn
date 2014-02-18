all:
	gcc -fPIC -std=c99 -c grnlib.c
	gcc -shared grnlib.o -o libgrn.so
	sudo mv libgrn.so /usr/lib
