MKDIR_P = mkdir -p

OUT_DIR = plots

directories: ${OUT_DIR}

.PHONY: directories thesis.pdf all # clean

all: directories thesis.pdf all # clean


${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}

thesis.pdf: 
	cd thesis && $(MAKE) 

#clean:
#
