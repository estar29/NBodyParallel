nbody_par: nbody-parallel.cpp
	   g++ nbody-parallel.cpp -o nbody_par

# Test case with a simple 3 particle system.
test.out: nbody_par
	       ./nbody_par 3 1 100 5 > test.out

test.pdf: test.out
	  python3 plot.py test.out test.pdf 1000

# Test case 1
# Same solar case as the sequential code to check if parallel implementation
# is correct.
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
hun_part_par.out: nbody_par
	date
	./nbody_par 100 1 10000 100 > hun_part_par.out
	date

# Test case 3
# 1000 particles where time step size = 1, and number of steps = 10000.
thou_part_par.out: nbody_par
	date
	./nbody_par 1000 1 10000 100 > thou_part_par.out
	date
