LIBDIR=lib
MEKDDIR=..
INCDIR1=$(MEKDDIR)/src
INCDIR2=$(MEKDDIR)/src/Extra_code
INCDIR3=$(MEKDDIR)/src/higgs_properties
INCDIR4=$(MEKDDIR)/src/MadGraphSrc
INCDIR5=$(MEKDDIR)/interface
PDF_DIR=$(MEKDDIR)/src/PDFTables
CXXFLAGS= -pedantic -Wunused-variable -Werror=sign-compare -Wunused-but-set-variable -fPIC -O -I $(INCDIR1) -I $(INCDIR2) -I $(INCDIR3) -I $(INCDIR4) -I $(INCDIR5) -I.
CFLAGS= -pedantic -Wunused-variable -Werror=sign-compare -Wunused-but-set-variable -fPIC -std=c99
LIBFLAGS= -L$(LIBDIR)
UNAME := $(shell uname)

########## ROOT PART ##########
ROOTLIB=$(ROOTSYS)/lib
ROOTINCLUDES=$(ROOTSYS)/include
###############################


objects1_Z4l=$(INCDIR4)/Z4l/qq_Z4l_DN_OF.o $(INCDIR4)/Z4l/qq_Z4l_DN_OFpA.o $(INCDIR4)/Z4l/qq_Z4l_DN_SF.o $(INCDIR4)/Z4l/qq_Z4l_DN_SFpA.o $(INCDIR4)/Z4l/qq_Z4l_UP_OF.o $(INCDIR4)/Z4l/qq_Z4l_UP_OFpA.o $(INCDIR4)/Z4l/qq_Z4l_UP_SF.o $(INCDIR4)/Z4l/qq_Z4l_UP_SFpA.o

objects1_ZZ=$(INCDIR4)/ZZ/qq_ZZ_DN_OF.o $(INCDIR4)/ZZ/qq_ZZ_UP_OF.o $(INCDIR4)/ZZ/qq_ZZ_DN_SF.o $(INCDIR4)/ZZ/qq_ZZ_UP_SF.o $(INCDIR4)/ZZ/qq_ZZ_DN_OFpA.o $(INCDIR4)/ZZ/qq_ZZ_UP_OFpA.o $(INCDIR4)/ZZ/qq_ZZ_DN_SFpA.o $(INCDIR4)/ZZ/qq_ZZ_UP_SFpA.o

objects1_spin0=$(INCDIR4)/Spin0/gg_Spin0_OF.o $(INCDIR4)/Spin0/gg_Spin0_SF.o $(INCDIR4)/Spin0/gg_Spin0_OFpA.o $(INCDIR4)/Spin0/gg_Spin0_SFpA.o $(INCDIR4)/Spin0/qq_Spin0_DN_OF.o $(INCDIR4)/Spin0/qq_Spin0_UP_OF.o $(INCDIR4)/Spin0/qq_Spin0_DN_SF.o $(INCDIR4)/Spin0/qq_Spin0_UP_SF.o $(INCDIR4)/Spin0/qq_Spin0_DN_OFpA.o $(INCDIR4)/Spin0/qq_Spin0_UP_OFpA.o $(INCDIR4)/Spin0/qq_Spin0_DN_SFpA.o $(INCDIR4)/Spin0/qq_Spin0_UP_SFpA.o $(INCDIR4)/Spin0/Spin0_OF.o $(INCDIR4)/Spin0/Spin0_SF.o $(INCDIR4)/Spin0/Spin0_OFpA.o $(INCDIR4)/Spin0/Spin0_SFpA.o

objects1_spin1=$(INCDIR4)/Spin1/qq_Spin1_DN_OF.o $(INCDIR4)/Spin1/qq_Spin1_UP_OF.o $(INCDIR4)/Spin1/qq_Spin1_DN_SF.o $(INCDIR4)/Spin1/qq_Spin1_UP_SF.o $(INCDIR4)/Spin1/qq_Spin1_DN_OFpA.o $(INCDIR4)/Spin1/qq_Spin1_UP_OFpA.o $(INCDIR4)/Spin1/qq_Spin1_DN_SFpA.o $(INCDIR4)/Spin1/qq_Spin1_UP_SFpA.o $(INCDIR4)/Spin1/Spin1_OF.o $(INCDIR4)/Spin1/Spin1_SF.o $(INCDIR4)/Spin1/Spin1_OFpA.o $(INCDIR4)/Spin1/Spin1_SFpA.o

objects1_spin2=$(INCDIR4)/Spin2/gg_Spin2_OF.o $(INCDIR4)/Spin2/gg_Spin2_SF.o $(INCDIR4)/Spin2/gg_Spin2_OFpA.o $(INCDIR4)/Spin2/gg_Spin2_SFpA.o $(INCDIR4)/Spin2/qq_Spin2_DN_OF.o $(INCDIR4)/Spin2/qq_Spin2_UP_OF.o $(INCDIR4)/Spin2/qq_Spin2_DN_SF.o $(INCDIR4)/Spin2/qq_Spin2_UP_SF.o $(INCDIR4)/Spin2/qq_Spin2_DN_OFpA.o $(INCDIR4)/Spin2/qq_Spin2_UP_OFpA.o $(INCDIR4)/Spin2/qq_Spin2_DN_SFpA.o $(INCDIR4)/Spin2/qq_Spin2_UP_SFpA.o $(INCDIR4)/Spin2/Spin2_OF.o $(INCDIR4)/Spin2/Spin2_SF.o $(INCDIR4)/Spin2/Spin2_OFpA.o $(INCDIR4)/Spin2/Spin2_SFpA.o

objects1_DY=$(INCDIR4)/DY/qq_DY_DN_2l.o $(INCDIR4)/DY/qq_DY_UP_2l.o $(INCDIR4)/DY/qq_DY_DN_2lpA.o $(INCDIR4)/DY/qq_DY_UP_2lpA.o

objects1_spin0_2=$(INCDIR4)/Spin0/gg_Spin0_2l.o $(INCDIR4)/Spin0/qq_Spin0_DN_2l.o $(INCDIR4)/Spin0/qq_Spin0_UP_2l.o $(INCDIR4)/Spin0/gg_Spin0_2lpA.o $(INCDIR4)/Spin0/qq_Spin0_DN_2lpA.o $(INCDIR4)/Spin0/qq_Spin0_UP_2lpA.o

objects1_spin1_2=$(INCDIR4)/Spin1/qq_Spin1_DN_2l.o $(INCDIR4)/Spin1/qq_Spin1_UP_2l.o $(INCDIR4)/Spin1/qq_Spin1_DN_2lpA.o $(INCDIR4)/Spin1/qq_Spin1_UP_2lpA.o

objects1_spin2_2=$(INCDIR4)/Spin2/gg_Spin2_2l.o $(INCDIR4)/Spin2/qq_Spin2_DN_2l.o $(INCDIR4)/Spin2/qq_Spin2_UP_2l.o $(INCDIR4)/Spin2/gg_Spin2_2lpA.o $(INCDIR4)/Spin2/qq_Spin2_DN_2lpA.o $(INCDIR4)/Spin2/qq_Spin2_UP_2lpA.o


objects2=$(INCDIR4)/HelAmps_HEF_UFO_bkg.o $(INCDIR4)/HelAmps_HEF_MEKD_spin0.o $(INCDIR4)/HelAmps_HEF_MEKD_spin1.o $(INCDIR4)/HelAmps_HEF_MEKD_spin2.o $(INCDIR4)/HelAmps_HZZ_Unitary_bkgpA.o $(INCDIR4)/HelAmps_HEF_MEKD_spin0pA.o $(INCDIR4)/HelAmps_HEF_MEKD_spin1pA.o $(INCDIR4)/HelAmps_HEF_MEKD_spin2pA.o $(INCDIR4)/HelAmps_HEF_MEKD_spinX.o $(INCDIR4)/Parameters_MEKD.o $(INCDIR4)/read_slha.o

objects2_2=$(INCDIR4)/HelAmps_HEF_MEKD2_BKG_DY_2l.o $(INCDIR4)/HelAmps_HEF_MEKD2_BKG_DY_2lpA.o $(INCDIR4)/HelAmps_HEF_MEKD2_spin0_2l.o $(INCDIR4)/HelAmps_HEF_MEKD2_spin0_2lpA.o $(INCDIR4)/HelAmps_HEF_MEKD2_spin1_2l.o $(INCDIR4)/HelAmps_HEF_MEKD2_spin1_2lpA.o $(INCDIR4)/HelAmps_HEF_MEKD2_spin2_2l.o


objects3=$(INCDIR2)/MEKD_CalcHEP_PDF.o $(INCDIR2)/MEKD_CalcHEP_Extra_functions.o $(INCDIR2)/MEKD_MG_Boosts.o $(INCDIR3)/hggeffective.o $(PDF_DIR)/pdt.o

objects4=$(INCDIR1)/MEKD_MG.o


objects1=$(objects1_ZZ) $(objects1_Z4l) $(objects1_spin0) $(objects1_spin1) $(objects1_spin2) $(objects1_DY) $(objects1_spin0_2) $(objects1_spin1_2) $(objects1_spin2_2)

objects=$(objects1) $(objects2) $(objects2_2) $(objects3) $(objects4)


targetlib=$(LIBDIR)/libMEKD_MG.so
target=runKD_MAD
test=MEKD_Test



# test: $(test)

all: $(target)




############################
### Linux building block ###
############################

ifeq ($(UNAME), Linux)
# do something Linux-y
ifeq ($(MEKD_STANDALONE), Yes)
CXXFLAGS+= -D MEKD_STANDALONE
endif
$(targetlib): $(objects)
	$(CXX) $(CXXFLAGS) -fPIC -shared -Wl,-soname,libMEKD_MG.so -o $(targetlib) $(objects) $(LIBFLAGS)


ifeq ($(MEKD_COMPILE_WITH_ROOT), Yes)
CXXFLAGS+= -D MEKD_with_ROOT -std=c++0x

$(target): $(targetlib)
	$(CXX) $(CXXFLAGS) `$(ROOTSYS)/bin/root-config --ldflags` -I $(ROOTINCLUDES) $(target).cpp -o $(target) -L. $(LIBFLAGS) -L $(ROOTLIB) `$(ROOTSYS)/bin/root-config --libs --glibs` $(targetlib)

$(test): $(targetlib)
	$(CXX) $(CXXFLAGS) `$(ROOTSYS)/bin/root-config --ldflags` -I $(ROOTINCLUDES) $(test).cpp -o $(test) -L. $(LIBFLAGS) -L $(ROOTLIB) `$(ROOTSYS)/bin/root-config --libs --glibs` $(targetlib)
else
$(target): $(targetlib)
	$(CXX) $(CXXFLAGS) $(target).cpp -o $(target) -L. $(LIBFLAGS) $(targetlib)

$(test): $(targetlib)
	$(CXX) $(CXXFLAGS) $(test).cpp -o $(test) -L. $(LIBFLAGS) $(targetlib)
endif

endif




#############################
### Mac OS building block ###
#############################

ifeq ($(UNAME), $(filter $(UNAME),OSX Darwin))
# do something OSX-y
ifeq ($(MEKD_STANDALONE), Yes)
CXXFLAGS+= -D MEKD_STANDALONE
endif
$(targetlib): $(objects)
	$(CXX) $(CXXFLAGS) -shared -Wl,-install_name,libMEKD_MG.so -o $(targetlib) $(objects) $(LIBFLAGS)


ifeq ($(MEKD_COMPILE_WITH_ROOT), Yes)
CXXFLAGS+= -D MEKD_with_ROOT -std=c++0x
$(target): $(targetlib)
	$(CXX) $(CXXFLAGS) `$(ROOTSYS)/bin/root-config --ldflags` -I $(ROOTINCLUDES) $(target).cpp -o $(target) -L. $(LIBFLAGS) -L $(ROOTLIB) `$(ROOTSYS)/bin/root-config --libs --glibs` $(targetlib)

$(test): $(targetlib)
	$(CXX) $(CXXFLAGS) `$(ROOTSYS)/bin/root-config --ldflags` -I $(ROOTINCLUDES) $(test).cpp -o $(test) -L. $(LIBFLAGS) -L $(ROOTLIB) `$(ROOTSYS)/bin/root-config --libs --glibs` $(targetlib)
else
$(target): $(targetlib)
	$(CXX) $(CXXFLAGS) $(target).cpp -o $(target) -L. $(LIBFLAGS) $(targetlib)

$(test): $(targetlib)
	$(CXX) $(CXXFLAGS) $(test).cpp -o $(test) -L. $(LIBFLAGS) $(targetlib)
endif

endif



.PHONY: clean

clean:
	rm -f $(targetlib); \
	rm -f $(target); \
	rm -f $(test); \
	rm -f $(objects)
