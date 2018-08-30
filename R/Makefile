

all: install

clean:
	rm -rf ./partitionstability.Rcheck
	rm -f partitionstability*.tar.gz

build: clean
	R CMD build partitionstability

install: clean
	R CMD INSTALL partitionstability

check: build
	R CMD check partitionstability_*.tar.gz

cran: build
	R CMD check --as-cran partitionstability_*.tar.gz

doc: 
	rm -f partitionstability.pdf
	R -e 'library("devtools"); devtools::document("partitionstability")'
	R CMD Rd2pdf partitionstability
