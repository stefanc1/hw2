/**
 *  \file mandelbrot_joe.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include "render.hh"

#define WIDTH 1000
#define HEIGHT 1000

int
mandelbrot(double x, double y) {
	int maxit = 511;
	double cx = x;
	double cy = y;
	double newx, newy;

	int it = 0;
	for (it = 0; it < maxit && (x*x + y*y) < 4; ++it) {
		newx = x*x - y*y + cx;
		newy = 2*x*y + cy;
		x = newx;
		y = newy;
	}
	return it;
}

int
main (int argc, char* argv[])
{
	double minX = -2.1;
	double maxX = 0.7;
	double minY = -1.25;
	double maxY = 1.25;
	double t_start, t_elapsed;

	int height, width;
	if (argc == 3) {
		height = atoi (argv[1]);
		width = atoi (argv[2]);
		assert (height > 0 && width > 0);
	} else {
		fprintf (stderr, "usage: %s <height> <width>\n", argv[0]);
		fprintf (stderr, "where <height> and <width> are the dimensions of the image.\n");
		return -1;
	}
	double it = (maxY - minY)/height;
	double jt = (maxX - minX)/width;
	double x, y;

	int rank, np;
	MPI_Status status;

	//MPI Start
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* Get process id */
	MPI_Comm_size (MPI_COMM_WORLD, &np);	/*Get number of processes*/

	MPI_Barrier (MPI_COMM_WORLD);

	if(rank == 0){
		t_start = MPI_Wtime();
	}
	int N = height / np;
	int blockSize = N * width;
	int sendBuffer[blockSize];

	y = minY + rank*N*it;
	for (int i = 0; i < N; ++i) {
		x = minX;
		for (int j = 0; j < width; ++j) {
			sendBuffer[(i*width) + j] = mandelbrot(x,y);
			x += jt;
		}
		y += it;
	}

	if (rank == 0) {
		int leftOver = height - (N*np);
		int 

		int receiveBuffer = (int *)malloc(np * blockSize 8 sizeof(int));
	}





}

/* eof */
