all :
	cd submodules && make all
	
test : all
	python tests/tests.py

clean :
	cd submodules && make clean
	
