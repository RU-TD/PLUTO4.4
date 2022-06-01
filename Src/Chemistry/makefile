# Makefile for Chemistry sources   

.SUFFIXES: .f .f90

VPATH        += $(SRC)/Chemistry
INCLUDE_DIRS += -I$(SRC)/Chemistry

CHEM_OBJ = chemistry.o
#prizmo objects
objs = opkda1.o
objs += opkda2.o
objs += opkdmain.o
objs += prizmo_fit.o
#objs += prizmo_spline.o
objs += prizmo_commons.o
objs += prizmo_utils.o
objs += prizmo_loaders.o
objs += prizmo_shielding.o
objs += prizmo_rates.o
objs += prizmo_rates_photo.o
objs += prizmo_rates_heating.o
objs += prizmo_flux.o
objs += prizmo_tdust.o
objs += prizmo_heating_photo.o
objs += prizmo_heating_photoelectric.o
objs += prizmo_heating_H2diss.o
objs += prizmo_heating_CR.o
objs += prizmo_heating.o
objs += prizmo_cooling_H2.o
objs += prizmo_cooling_CO.o
objs += prizmo_cooling_atomic.o
objs += prizmo_cooling_dust.o
objs += prizmo_cooling_chemical.o
objs += prizmo_cooling.o
objs += prizmo_ode.o
objs += prizmo_attenuate.o
objs += prizmo_core.o
objs += prizmo.o
OBJ     += $(CHEM_OBJ) $(objs)
HEADERS += chemistry.h mod_defs.h

LDFLAGS += -lm

$(CHEM_OBJ):  $(HEADERS)

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $<

.f90.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $<

%.o: %.mod

