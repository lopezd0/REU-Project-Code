
//#include "stdafx.h"
#define _CRT_SECURE_NO_DEPRECATE
#define M	512      // height of image
#define N	512     // width of image

#include "pgm.h"
#include <complex.h>
#include "../fftw3/fftw3.h"

float Estimate_ReferenceLight_fast(float* I1, float* I2, float step);
float Estimate_ReferenceLight_original(float* I1, float* I2, float step);
void correlation_fftw(int width, int height, float* realData1, float* imagData1, float* realData2, float* imagData2, double& peak);
int find_maximum_index(float *a, int n);
int find_minimum_index(float *a, int n);

void main(int argc, char *argv[])
{

	clock_t start, end;
	//  struct timeval beginTime, endTime, elapsedTime;

	FILE *fpin, *fpinI1, *fpinI2, *fpout;

	float* img = new float[N*M];
	float* I1 = new float[N*M];
	float* I2 = new float[N*M];

	//read the orignal image lena
	fpin = fopen("lena.dat", "rb");
	if (fpin == NULL)
	{
		puts("Can't find input file: lena.dat");
		exit(1);
	}
	fread(img, N*M, sizeof(float), fpin);
	fclose(fpin);
	printf("Done with reading lena.dat\n");

	//read the first image I1
	fpinI1 = fopen("I1.dat", "rb");
	if (fpinI1 == NULL)
	{
		puts("Can't find input file: I1.dat");
		exit(1);
	}
	fread(I1, N*M, sizeof(float), fpinI1);
	fclose(fpinI1);
	printf("Done with reading I1.dat\n");

	//read the second image I2
	fpinI2 = fopen("I2.dat", "rb");
	if (fpinI2 == NULL)
	{
		puts("Can't find input file: I2.dat");
		exit(1);
	}
	fread(I2, N*M, sizeof(float), fpinI2);
	fclose(fpinI2);
	printf("Done with reading I2.dat\n");

	////verify correctness of data
	//printf("%f  %f \n", I1[0], I1[1]);
	//printf("%f  %f \n", I2[0], I2[1]);
	//printf("%f  %f \n", 255*img[256*M+200], 255*img[200*M+200]);
	  
	//floatWritePGM("out.pgm", img, M, N);


	start = clock();

	float step = 0.1;
	float R = Estimate_ReferenceLight_fast(I1, I2, step);
	//float R = Estimate_ReferenceLight_original(I1, I2, step);

	end = clock();

	printf("The estimated reference light intensity is %f\n\n", R);
	printf("The image size is %d and %d\n", M, N);
	printf("Time required for execution: %f seconds\n\n", (double)(end - start) / CLOCKS_PER_SEC);


	delete[] I1;
	delete[] I2;
	delete[] img;

	getchar();

}

//function to estimate the reference light using the fast algorithm based on minimizing an error function
float Estimate_ReferenceLight_fast(float* I1, float* I2, float step)
{
	float Rtest;
	// compute the range of search
	//you code here to calculate Rmin and Rmax and replace the following numbers
	float Rmin = 1.0;
	float Rmax = 2.5;

	//I0 is the sum of Reference intensity and object intensity
	float* I0 = new float[N*M];
	
	//y used to store all the errors or correlation peaks
	float* y = new float[1000];

	//calculate the number of elements of y based on the step size
	int ysize = (int)((Rmax - Rmin) / step);

	int idx = 0;

	for (int r = (int)(Rmin / step); r < (int)(Rmax / step); r++)
	{
		//	Rtest = r/100.0; //step size of 0.01
		Rtest = r*step;
		printf("Rtest= %f\n\n", Rtest);

		float sum = 0.0;

		for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
		{
			
			float i1 = I1[i*N + j]; //get a pixel of I1
			float i2 = I2[i*N + j]; //get a pixel of I2
			//To see the result image, call the function floatWritePGM("out.pgm", I1, M, N);

			
			float b = (2 * Rtest *Rtest + i1 + i2);
			float ac = 2 * (i1*i1 + i2*i2 + 4.0 * powf(Rtest, 4));
			float i0 = 0.5*b - 0.5*sqrt(b*b - ac);
			float o2 = abs(i0 - Rtest*Rtest);

			sum += o2;
		}

		y[idx] = sum;
		idx++;
	}

	idx = find_minimum_index(y, ysize);

//	printf("correlation peak at the index of %d\n", idx);

	Rtest = Rmin + idx*step;

	// Cleanup.
	delete[] I0;
	delete[] y;

	return Rtest;
}

////function to estimate the reference light using the original algorithm basedon maximizing correlation peak
float Estimate_ReferenceLight_original(float* I1, float* I2, float step)
{
	float Rtest;
	// compute the range of search
	//you code here to calculate Rmin and Rmax and replace the following numbers
	float Rmin = 1.0;
	float Rmax = 2.5;

	//I0 is the sum of Reference intensity and object intensity
	float* I0 = new float[N*M];
	//Hr is a complex image, need two arrays, Hr1 for real part and Hr2 for imaginary part 
	float* Hr1 = new float[N*M];
	float* Hr2 = new float[N*M];
	//Ht is a complex image, need two arrays, Ht1 for real part and Ht2 for imaginary part 
	float* Ht1 = new float[N*M];
	float* Ht2 = new float[N*M];
	//y used to store all the errors or correlation peaks
	float* y = new float[1000];
	
	//calculate the number of elements of y based on the step size
	int ysize = (int) ((Rmax - Rmin) / step);


	int idx = 0;
	for (int r = (int)(Rmin/step); r < (int)(Rmax/step); r++)
	{    
	
		Rtest = r*step;
	//	printf("Rtest= %f\n\n", Rtest);

		for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
		{
			
			float i1 = I1[i*N + j]; //get a pixel of I1
			float i2 = I2[i*N + j]; //get a pixel of I2
			//To see the result image, call the function floatWritePGM("out.pgm", I1, M, N);

			float b = (2 * Rtest *Rtest + i1 + i2);
			float ac = 2 * (i1*i1 + i2*i2 + 4.0 * powf(Rtest, 4));
			float i0 = 0.5*b - 0.5*sqrt(b*b - ac);
			I0[i*N + j] = i0;

			
			Ht1[i*N + j] = 0.5*(i1 - i2);
			Ht2[i*N + j] = -0.5*(i1 - i2);
			Hr1[i*N + j] = (i1 - i0) / 2.0;
			Hr2[i*N + j] = (i2 - i0) / 2.0;
			
		}
		//printf("The Rtest is : %f%f\n", Rtest);
		double peakRT = 0, peakRR=0;
		correlation_fftw(N, M, Ht1, Ht2, Hr1, Hr2, peakRT);
		correlation_fftw(N, M, Hr1, Hr2, Hr1, Hr2, peakRR);
		y[idx] = (float)peakRT/peakRR;
		idx++;
	}

	idx = find_maximum_index(y, ysize);
	
//	printf("correlation peak at the index of %d\n", idx);

	Rtest = Rmin  + idx*step;

	// Cleanup.
	delete[] I0;
	delete[] Hr1;
	delete[] Hr2;
	delete[] Ht1;
	delete[] Ht2;
	delete[] y;

	return Rtest;
}

// Perform correlation between two complex images and return the largest correlation peak
//	first image: realData1 + i*imagData1,  second image:  realData2 + i*imagData2
void correlation_fftw(int width, int height, float* realData1, float* imagData1, float* realData2, float* imagData2, double& peak)
{
	// The size of the whole image.
	int size = height * width;

	// Allocate input and output buffers for the first image
	fftw_complex* buf1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* size);
	fftw_complex* buf2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* size);
	// Copy in image data as real values for the transform.
	for (int i = 0; i < size; i++) {
		buf1[i][0] = (double)realData1[i];
		buf1[i][1] = (double)imagData1[i]; 
	}
	// Transform to frequency space.
	fftw_plan pFwd = fftw_plan_dft_2d(width, height, buf1, buf2, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(pFwd);
	fftw_destroy_plan(pFwd);

	// Allocate input and output buffers for the second image
	fftw_complex* buf3 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* size);
	fftw_complex* buf4 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* size);
	// Copy in image data as real values for the transform.
	for (int i = 0; i < size; i++) {
		buf3[i][0] = (double)realData2[i];
		buf3[i][1] = (double)imagData2[i];
	}
	// Transform to frequency space.
	fftw_plan pFwd2 = fftw_plan_dft_2d(width, height, buf3, buf4, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(pFwd2);
	fftw_destroy_plan(pFwd2);

	fftw_complex* buf5 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* size);
	for (int i = 0; i < size; i++) {
		buf5[i][0] = buf2[i][0] * buf4[i][0] + buf2[i][1] * buf4[i][1];
		buf5[i][1] = buf2[i][0] * buf4[i][1] - buf2[i][1] * buf4[i][0];
	}

	// Transform back to image space.
	fftw_plan pBack = fftw_plan_dft_2d(width, height, buf5, buf1, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(pBack);
	fftw_destroy_plan(pBack);

	// Have to scale the output values to get back to the original.
	for (int i = 0; i < size; i++) {
		buf1[i][0] = buf1[i][0] / size;
		buf1[i][1] = buf1[i][1] / size;
	}

	double max=0.0; 
	// Copy the magnitude of the result back into the image.
	for (int i = 0; i < size; i++) {
		double re = buf1[i][0];
		double im = buf1[i][1];
		double mag = sqrt(re*re + im*im);
	//	magnitude[i] = (float)mag;
		if (i == 0) 
			max = mag;
		else
		{
			if (mag > max)
				max = mag;
		}
		
	}

	peak = max;

	// Cleanup.
	fftw_free(buf1);
	fftw_free(buf2);
	fftw_free(buf3);
	fftw_free(buf4);
	fftw_free(buf5);
}


int find_maximum_index(float *a, int n) {
	int c, max, index;

	max = a[0];
	index = 0;

	for (c = 1; c < n; c++) {
		if (a[c] > max) {
			index = c;
			max = a[c];
		}
	}

	return index;
}

int find_minimum_index(float *a, int n) {
	int c, min, index;

	min = a[0];
	index = 0;

	for (c = 1; c < n; c++) {
		if (a[c] < min) {
			index = c;
			min = a[c];
		}
	}

	return index;
}
