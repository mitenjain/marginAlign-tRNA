all :
	cd submodules && make all
	
test : all
	scripts/executionScript.sh ./tests/tests.py --logInfo

testLong : all
	scripts/executionScript.sh ./tests/tests.py --logInfo --longTests

clean :
	cd submodules && make clean
	
