RDIR= .

RFILES := $(wildcard $(RDIR)/*.R)
OUT_FILES:= $(RFILES:.R=.Rout)

all: $(OUT_FILES)

$(RDIR)/%.Rout: $(RDIR)/%.R 
	R CMD BATCH $<

R: $(OUT_FILES)

clean:
	rm -fv $(OUT_FILES) 

.PHONY: all 



