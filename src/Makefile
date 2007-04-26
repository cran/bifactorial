bifactorial.so : ave.o krit.o pv.o Rcpp.o misc.o
	g++ -shared -o $@ $^ -L/usr/lib/R/lib -lR

ave.o : ave.cpp
	g++ -I/usr/share/R/include -I/usr/share/R/include -fpic -g -O2 -c $^ -o $@

krit.o : krit.cpp
	g++ -I/usr/share/R/include -I/usr/share/R/include -fpic -g -O2 -c $^ -o $@

misc.o : misc.cpp
	g++ -I/usr/share/R/include -I/usr/share/R/include -fpic -g -O2 -c $^ -o $@

pv.o : pv.cpp
	g++ -I/usr/share/R/include -I/usr/share/R/include -fpic -g -O2 -c $^ -o $@

Rcpp.o : Rcpp.cpp
	g++ -I/usr/share/R/include -I/usr/share/R/include -fpic -g -O2 -c $^ -o $@

clean:
	rm -f bifactorial.so ave.o krit.o misc.o pv.o Rcpp.o