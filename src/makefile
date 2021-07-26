.SUFFIXES: .f90
.SUFFIXES: .o
FC = mpif90
FCOPT = -O2 -fpp # -CB
EXE = MD.out
SRCS = subroutine_function.f90 MD_mod.f90 fix_compute_mod.f90\
	  fix_freeze.f90 fix_langevin.f90 fix_langevin2.f90 fix_move_linear.f90 fix_set_pos.f90 fix_addforce.f90\
	  fix_rigid.f90\
	  fix.f90\
	  compute_temp.f90 compute_chunk_1d.f90 compute_compv.f90\
	  compute_group_group.f90 compute_group_group_hf.f90\
	  compute_correlate.f90 compute_ave.f90\
	  compute.f90\
	  output.f90\
	  MD.f90
OBJS = $(SRCS:.f90=.o)

all: $(EXE)

$(EXE): $(OBJS)
	$(FC) $^ -o $@ $(FCOPT)

.f90.o:
	$(FC) -c $< $(FCOPT)

.o.out:
	$(FC) $< -o $@ $(FCOPT)

clean:
	rm -rf ./*.o ./*.mod ./*.out $(EXE)

