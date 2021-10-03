CC 	= g++

CFLAGS 	= -O3 -fomit-frame-pointer -fopenmp -std=c++11
#CFLAGS 	= -O3 -march=i586 -fomit-frame-pointer -I/home/maruyama/libs/include
#LIBS = -L/home/maruyama/libs/lib -lfftw3 
LIBS = -lfftw3_threads -lfftw3 -lm

OBJS   	= main.o\
	  add_tuv.o\
	  anderson.o\
	  cal_Coulomb.o\
	  cal_exchem.o\
	  cal_exnew.o\
	  cal_grad.o\
	  cal_LJ.o\
	  cal_pmv.o\
	  cal_potential.o\
	  cal_rms.o\
	  calculate.o\
	  cell.o\
	  excp_fd.o\
	  fft3d.o\
	  initialize.o\
	  initialize_g.o\
	  initialize_tuv.o\
	  iterate.o\
	  output.o\
	  output_cuv.o\
	  output_grad.o\
	  output_guv.o\
	  output_huv.o\
	  output_xmu.o\
	  read_input.o\
	  read_tuv.o\
	  set_fname.o\
	  set_solvent.o\
	  solute.o\
	  solvent_read.o\
	  solvent_spline.o\
	  solvent_spline2.o\
	  spline.o\
	  splint.o\
	  write_tuv.o\
	  
PROGRAM = 3drism
.cc.o:
	$(CC) $(CFLAGS) -c $< -o $@

all:	$(PROGRAM)

$(PROGRAM) : $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) $(LIBS) -o $(PROGRAM)

clean:;	@rm -f $(OBJS) core

$(OBJS) : 




