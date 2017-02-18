/**
 *  \file mandelbrot_ms.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include "render.hh"

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

	// master
	if(rank == 0){
		t_start = MPI_Wtime();
		printf("Number of MPI processes: %d\r\n", np);


		MPI_Status status;
		int terminationFlag = -1;


		int finalArray[width*height];

		// store the currentRow at the end needed when master received it
		int receiveBuff[width + 1];

		int currRow = 0;

		// set currProcessor to 1 instead of 0 because we are excluding master tread
		int currProcessor = 1;


		while (currRow < height) {

			// init the send to the slaves
			// need to exclude the master processor
			if (currProcessor < np) {
				MPI_Send(&currRow, 1, MPI_INT, currProcessor, 0, MPI_COMM_WORLD);
				currProcessor++;
			}
			else {
				MPI_Recv(receiveBuff, width+1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Send(&currRow, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
				memcpy(finalArray + receiveBuff[width]*width, receiveBuff, width*sizeof(int));
			}

			currRow++;
		}

		// receive the last rows from the slaves and indicate termination
		for (int i = 1; i < np; i++) {
			MPI_Recv(receiveBuff, width+1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Send(&terminationFlag, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			memcpy(finalArray + receiveBuff[width]*width, receiveBuff, width*sizeof(int));
		}


		// Generate Image
		gil::rgb8_image_t img(height, width);
		auto img_view = gil::view(img);

		for (int k = 0; k < height; ++k) {
			for (int p = 0; p < width; ++p) {
				img_view(p, k) = render(receiveBuff[ (k*width) + p] / 512.0);
			}
		}

		t_elapsed = MPI_Wtime () - t_start; // compute the overall time taken
	    printf("Total time: %f\r\n", t_elapsed);
		gil::png_write_view("mandelbrot-ms.png", const_view(img));

	}

	// slaves
	else {
		MPI_Status status;
		int currRow;
		int sendBuff[width + 1];

		while(1) {
			MPI_Recv(&currRow, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			if (currRow == -1) {
				break;
			}	

			y = minY + currRow*it;
			x = minX;
			for (int i = 0; i < width; ++i) {
				sendBuff[i] = mandelbrot(x,y);
				x += jt;
			}
			sendBuff[width] = currRow;
			MPI_Send(sendBuff, width+1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		
	}


	MPI_Finalize();

	return 0;
}

/* eof */
