CXXFLAGS=-O3

nbody: nbody.cpp
	g++ -O3 nbody.cpp -o nbody

# Test case 1
# Solar simulation.
solar_seq.out: nbody
	date
	./nbody "planet" 200 5000000 10000 > solar_seq.out # maybe a minutes
	date

solar_seq.pdf: solar.out
	python plot.py solar_seq.out solar_seq.pdf 1000

# Test case 2
# 100 particles where time step size = 1, and number of steps = 10000.
hundred_seq.out: nbody
	date
	./nbody 100 1 10000 100 > hundred_seq.out
	date

hundred_seq.pdf: 100p_seq.out
	python3 plot.py hundred_seq.out hundred_seq.pdf 1000

# Test case 3
# 1000 particles where time step size = 1, and number of steps = 10000.
thousand_seq.out: nbody
	date
	./nbody 1000 1 10000 100 > thousand_seq.out
	date

thousand_seq.pdf: 1000p_seq.out
	python3 plot.py thousand_seq.out thousand_seq.pdf 1000
