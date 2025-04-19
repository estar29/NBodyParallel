nbody_par: nbody-parallel.cpp
	   g++ nbody-parallel.cpp -o nbody_par


# Test case 1
# Same solar case as the sequential code
solar_par.out: nbody_par
	date
	./nbody_par "planet" 200 5000000 10000 > solar_par.out # maybe a minutes
	date

solar_par.pdf: solar_par.out
	date
	python plot.py solar_par.out solar_par.pdf 1000
	date


# Test case 2
# 100 particles where time step size = 1, and number of steps = 10000.
hundred_par.out: nbody_par
	date
	./nbody_par 100 1 10000 100 > hundred_par.out
	date

hundred_par.pdf: hundred_par.out
	python3 plot.py hundred_par.out hundred_par.pdf 1000


# Test case 3
# 1000 particles where time step size = 1, and number of steps = 10000.
thousand_par.out: nbody_par
	date
	./nbody_par 1000 1 10000 100 > thousand_par.out
	date

thousand_par.pdf: thousand_par.out
	python3 plot.py thousand_par.out thousand_par.pdf 1000
