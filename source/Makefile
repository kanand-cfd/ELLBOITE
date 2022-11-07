FC = gfortran

#Options of compilation
OPT = -fdefault-real-8 -fdefault-double-8 -O3 -C -g -fbacktrace -funroll-loops -fno-range-check -ffree-line-length-none -Wall -fbounds-check
#( -C -g -fbounds-check -Wall -finit-real=snan )-fdefault-integer-8

BLAS = /usr/local/lib/libblas.a

LAPACK = /usr/local/lib/liblapack.a

# Name of executable
EXE = ellbox

# Defining the objects variables
MODOBJS = \
mod_quaternion.o \
mod_param_phys.o \
mod_ellipsoid_particle.o \
minpack_reduc.o 


FOBJS = \
main.o \
readparams.o \
init_param.o \
init_particle.o \
adv_particle.o \
euler_integration.o \
quaternion_integration.o \
box_boundary_condition.o \
closest_point_plane.o \
max_encroachment_ellipsoid.o \
wall_ellipsoid_rebound.o \
collision.o \
ellipsoid_contact_detection_new.o \
gradient_descent_ellipsoid.o \
ellipsoid_contact_optimization.o \
ellipsoid_pair_depth.o \
eval_func.o \
eval_func2.o \
calculate_ellipsoid_matrix.o \
collision_response.o \
boundary_particle.o \
stat_particle.o \
stat_particle_new.o \
transform_basis.o \
invert_ndim3_matrix.o \
save_particle.o

MOD = *.mod

LIBS = $(BLAS) $(LAPACK)

OBJS = $(MODOBJS) $(FOBJS) $(EOBJS)

$(EXE):	$(OBJS)
	$(FC) $(OPT) $(OBJS) $(LIBS) -o $@

clean:
	/bin/rm -f $(OBJS) $(MOD) $(EXE)

.SUFFIXES: .o .f90 .f


.f90.o: 
	$(FC) -c $(OPT) $<

.f:
	$(FC) -c $(OPT) $<
