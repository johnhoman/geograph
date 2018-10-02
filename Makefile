PYTHON_BASE=/usr/local
PYTHON_INCLUDE_DIR=$(PYTHON_BASE)/include/python3.6m
PYTHON_LIB_DIR=$(PYTHON_BASE)/lib

geograph:
	  gcc -pthread -g -fPIC -Wall -DNDEBUG \
		    -I$(PYTHON_INCLUDE_DIR) -L$(PYTHON_LIB_DIR) \
				geograph.c -o test -lm -lpython3.6m \
				-Wl,--no-as-needed -ldl -lutil
