all :
	cd submodules && make all
	
test : all
	scripts/executionScript.sh ./tests/tests.py --logInfo

clean :
	cd submodules && make clean
	
