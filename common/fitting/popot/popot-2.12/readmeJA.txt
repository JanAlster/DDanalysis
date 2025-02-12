po 
cmake .. -G"Unix Makefiles"
je taky potreba zmenit 
build/src/CMakeFiles/PyPopot.dir/link.txt
-lboost_python
na
-lboost_python-py34
