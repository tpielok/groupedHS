MKDIR_P = mkdir -p

OUT_DIR = plots

.PHONY: all

all: directories thesis.pdf all

directories: ${OUT_DIR}

${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}

thesis.pdf:
	cd src && $(MAKE)
	cd thesis && $(MAKE) 

clean:
	cd src && $(MAKE) clean
	cd thesis && $(MAKE) clean

