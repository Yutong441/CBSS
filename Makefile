# Build CBSS environment
install:
	mkdir software -p; \
	# install niimath
	git clone https://github.com/rordenlab/niimath.git; \
	mv niimath software; cd software/niimath; mkdir build -p; cd build; cmake ..; make; cd ../../../; \
	# python environment
	cd software; python3 -m venv cbss_venv; cd ..; \
	. software/cbss_venv/bin/activate; pip install pip --upgrade; pip install -r requirements.txt; deactivate;\
	# compile fortran
	cd CBSS/module; gfortran -c -free -fPIC bobyqa.f90; cd ../..; \
	cd CBSS; mv module/bobyqa.o src/; mv module/bobyqa_module.mod src/; cd ..; \
	. software/cbss_venv/bin/activate; cd CBSS/src; python -m numpy.f2py -llapack -c -I. bobyqa.o -m fort *.f90; cd -
