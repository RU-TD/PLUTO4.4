# Add objects file for Particles module

VPATH        += $(SETUP_DIR)/particles_picogna
INCLUDE_DIRS += -I$(SETUP_DIR)/particles_picogna

OBJ +=  particles_set.o particles_boundary.o particles_output.o particles_init.o
OBJ +=  particles_tools.o particles_dtype.o
HEADERS += particles.h

particles_set.o:  particles.h
particles_init.o: particles.h
particles_boundary.o:  particles.h
particles_output.o: particles.h
particles_tools.o: particles.h
particles_dtype: particles.h
