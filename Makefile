# all: fluffy
# fluffy: interp.o fluffy.o
# 	gfortran interp.o fluffy.o -o fluffy

FLAGS= -g -pg
all: fluffy
fluffy: clouddata.o fluffycloud.o wrapper.o
	gfortran $(FLAGS) clouddata.o fluffycloud.o wrapper.o -o fluffy -pg

clouddata.o: clouddata.f
	gfortran $(FLAGS) -c clouddata.f

fluffycloud.o: fluffycloud.f
	gfortran $(FLAGS) -c fluffycloud.f

# sizes.o: sizes.f
# 	gfortran -c sizes.f

# fluffymod.o: fluffymod.f
# 	gfortran -c fluffymod.f
#
wrapper.o: wrapper.f
	gfortran $(FLAGS) -c wrapper.f

# all:
# 	@echo "awesome"
# PHONY: all
# all: ; $(info $$var is [${var}])echo Hello world

clean:
	rm -r *.o *.mod fluffy
