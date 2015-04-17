all :
	cd submodules && make all
	
test : all
	src/executionScript.sh ./tests/tests.py --logInfo

testLong : all
	src/executionScript.sh ./tests/tests.py --logInfo --longTests

clean :
	cd submodules && make clean
	
