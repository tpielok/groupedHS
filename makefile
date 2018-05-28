MKDIR_P = mkdir -p

OUT_DIR = plots

.PHONY: all

all: directories thesis.pdf all # clean

directories: ${OUT_DIR}

${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}

thesis.pdf: 
	cd thesis && $(MAKE) 

#clean:
#
