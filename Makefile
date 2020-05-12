all:
	nvcc -Xcompiler -fopenmp -lm -arch sm_35 -rdc=true -default-stream per-thread cubo_rec2.cu -o cuborec2
