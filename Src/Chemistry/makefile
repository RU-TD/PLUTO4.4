# Makefile for Chemistry sources   

.SUFFIXES: .f .f90

VPATH        += $(SRC)/Chemistry
INCLUDE_DIRS += -I$(SRC)/Chemistry

CHEM_OBJ = chemistry.o
F_OBJ    = opkda2.o opkda1.o opkdmain.o lkrm_commons.o lkrm_utils.o lkrm_rates.o lkrm_ode.o lkrm_rates_photo.o lkrm_xsecs.o lkrm_user_photo.o lkrm.o 
OBJ     += $(CHEM_OBJ) $(F_OBJ)
HEADERS += chemistry.h mod_defs.h

LDFLAGS += -lifcore -limf

$(CHEM_OBJ):  $(HEADERS)

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $<

.f90.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $<

%.o: %.mod

