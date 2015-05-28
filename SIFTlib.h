/* Program for Scale Invariant Feature Transform */
/* Written by:  R R Iyer
 * Institution: ISI Kolkata
 * Date:        27 May 2015
*/
#include "PGMlib.h"
#include <math.h>

/* Generates a nxn Gaussian filter with sigma = sigma */
double **generate_gaussian_filter(double sigma, int n)
{
	double **gKernel = allocate_dynamic_matrix_double(n,n);
    double r, s = 2.0 * sigma * sigma;

    // sum is for normalization
    double sum = 0.0;

    // generate 5x5 kernel
    int x,y;
    for (x = -(n-1)/2; x <= (n-1)/2; x++)
    {
        for(y = -(n-1)/2; y <= (n-1)/2; y++)
        {
            r = sqrt(x*x + y*y);
            gKernel[x + (n-1)/2][y + (n-1)/2] = (exp(-(r*r)/s))/(3.14159265 * s);
            sum += gKernel[x + (n-1)/2][y + (n-1)/2];
        }
    }


	return gKernel;
}

/* Generates a nxn normalised Gaussian filter with sigma = sigma */
double **generate_normalised_gaussian_filter(double sigma, int n)
{
	double **gKernel = allocate_dynamic_matrix_double(n,n);
    double r, s = 2.0 * sigma * sigma;

    // sum is for normalization
    double sum = 0.0;

    // generate 5x5 kernel
    int x,y;
    for (x = -(n-1)/2; x <= (n-1)/2; x++)
    {
        for(y = -(n-1)/2; y <= (n-1)/2; y++)
        {
            r = sqrt(x*x + y*y);
            gKernel[x + (n-1)/2][y + (n-1)/2] = (exp(-(r*r)/s))/(3.14159265 * s);
            sum += gKernel[x + (n-1)/2][y + (n-1)/2];
        }
    }

    // normalize the Kernel
    int i,j;
    for(i = 0; i < 5; ++i)
        for(j = 0; j < 5; ++j)
            gKernel[i][j] /= sum;

	return gKernel;
}

/* Function for convolution mask
* Inputs: in: matrix of size mxn
*         coeffs: matrix of size kxk (The filter mask)
*/
double **convolve(int **in, int m, int n, double **coeffs, int K)
{
	double **out = allocate_dynamic_matrix_double(m-K+1,n-K+1);
	int i, j, ii, jj;

	for (i = K / 2; i < m - K / 2; ++i) // iterate through image
	{
	    for (j = K / 2; j < n - K / 2; ++j)
	    {
	        double sum = 0; // sum will be the sum of input data * coeff terms

	        for (ii = - K / 2; ii <= K / 2; ++ii) // iterate over kernel
	        {
	            for (jj = - K / 2; jj <= K / 2; ++jj)
	            {
	                sum += in[i + ii][j +jj] * coeffs[ii + K / 2][jj + K / 2];
   	            }
	        }
	        out[i-K/2][j-K/2] = sum;
	    }

	}
	return out;
}

/* Function to downsample the picture to the given octave 
*  Arguments: data: The pic to be downsampled
*                o: Specify octave
*             down: Empty object to store result into
*/
void down_sample_pic(PGMData *data, int o, PGMData *down)
{
    int i,j;
    down->width= (int) ((double)data->width / pow(2,o));
    down->height= (int) ((double)data->height / pow(2,o));
    down->max_gray = data->max_gray;
    for(i=0; i<down->width; i++)
    {
    	for(j=0;j<down->height; j++)
    	{
    		int x = (int) pow(2,o);
    		down->pixels[i][j] = data->pixels[x*i][x*j];
    	}
    }
}

/* Function to calculate key_points based on local maxima
*  Arguments: data: The picture
*                o: The octave 
*          key_loc: Empty array to store the location of key_points
*                   Even terms contain x coordinate, odd terms cotain y coordinates
*          key_mag: Empty array to be passed as argument to store magnitude of key_points
*  Returns: number of key points
*/
int calculate_keypoints(PGMData *data, int o, int *key_loc, double *key_mag)
{
	printf("Calculating keypoints \n");
	int i,j,k,l, loc_count=0, mag_count=0;
	double **gauss, **h0, **h1, **h2, **h3;
	double sigma0 = 1.6;
	if(o!=0)
    {
		PGMData temp;
		temp.width = data->width;
		temp.height = data->height;
		temp.max_gray = data->max_gray;
		temp.pixels = allocate_dynamic_matrix(temp.width, temp.height);
		for(i=0; i<temp.width; i++)
		{
			for(j=0; j<temp.width; j++)
			{
				temp.pixels[i][j] = data->pixels[i][j];
			}
		}
		down_sample_pic(&temp,o,data);
    }

	for(i=0; i<4; i++)
	{
		double temp = 4*o + i;
		temp = temp/4;
		double sigma = sigma0 * pow(2,temp);
		printf("%f\n",sigma);
		gauss = generate_gaussian_filter(sigma,5);
	    if(i==0)
		{
	    	h0 = convolve(data->pixels, data->width, data->height, gauss, 5);
		}
		if(i==1)
		{
			h1 = convolve(data->pixels, data->width, data->height, gauss, 5);
		}
		if(i==2)
		{
			h2 = convolve(data->pixels, data->width, data->height, gauss, 5);
		}
		if(i==3)
		{
			h3 = convolve(data->pixels, data->width, data->height, gauss, 5);
		}
	}
	printf("...Applied filters \n");
	for(i=0; i<data->width-4; i++)
    {
        for(j=0; j<data->height-4; j++)
        {
    	    h0[i][j] = h1[i][j] - h0[i][j];
    	    h1[i][j] = h2[i][j] - h1[i][j];
    	    h2[i][j] = h3[i][j] - h2[i][j];
        }
    }
    PGMData result = make_PGM(h0, data->width-4, data->height-4);
    writePGM("h0.pgm",&result,1);
    result = make_PGM(h1, data->width-4, data->height-4);
    writePGM("h1.pgm",&result,1);
    result = make_PGM(h2, data->width-4, data->height-4);
    writePGM("h2.pgm",&result,1);
    printf("... Wrote images \n");

    for(i=1; i<data->width-5; i++)
    {
        for(j=1; j<data->height-5; j++)
        {
        	double array[26];
            int count=0;
        	for(k=-1; k<2; k++)
        	{
        		for(l=-1; l<2; l++)
        		{
                    array[count] = h0[i-k][j-l];
                    count++;
        		}
        	}
        	for(k=-1; k<2; k++)
        	{
        		for(l=-1; l<2; l++)
        		{
                    if(!(k==0 && l==0))
                        array[count] = h1[i-k][j-l];
                        count++;
        		}
        	}
        	for(k=-1; k<2; k++)
        	{
        		for(l=-1; l<2; l++)
        		{
                    array[count] = h2[i-k][j-l];
                    count++;
        		}
        	}
        	//printf("count = %d\n",count);
        	double max = array[0], min = array[0];
        	for(k=0; k<26; k++)
        	{
        		if(array[k]>max) max = array[k];
        		if(array[k]<min) min = array[k];
        	}
        	if(h1[i][j]>max || h1[i][j]<min)
        	{
        		key_loc[loc_count] = i;
        		loc_count++;
        		key_loc[loc_count] = j;
        		loc_count++;
        		key_mag[mag_count] = h1[i][j];
        		mag_count++;
        	}
        }
    }
    printf("Generated %d key points\n",mag_count);
    return mag_count;
}

/* To find the inverse of a 3x3 matrix 
* a is the input matrix
* inverse is the empty matrix to store the inverse into
*/
void find_inverse(double a[3][3], double inverse[3][3])
{
	int i,j;
	double determinant;
	for(i=0;i<3;i++)
	      determinant = determinant + (a[0][i]*(a[1][(i+1)%3]*a[2][(i+2)%3] - a[1][(i+2)%3]*a[2][(i+1)%3]));

	for(i=0;i<3;i++)
	    for(j=0;j<3;j++)
	        inverse[i][j] = ((a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3]) - (a[(i+1)%3][(j+2)%3]*a[(i+2)%3][(j+1)%3]))/ determinant;
}


/* Function to calculate key_points based on curve fitting with contrast setting and edge rejection
*  Arguments: data: The picture
*                o: The octave 
*          key_loc: Empty array to store the location of key_points
*                   Even terms contain x coordinate, odd terms cotain y coordinates
*          key_mag: Empty array to be passed as argument to store magnitude of key_points
*  Returns: number of key points
*/
int calculate_keypoints_mod(PGMData *data, int o, int *key_loc, double *key_mag)
{
	printf("Calculating keypoints \n");
	int i,j,loc_count=0, mag_count=0;
	double **gauss, **h0, **h1, **h2, **h3;
	double sigma0 = 1.6;

	if(o!=0)
    {
		PGMData temp;
		temp.width = data->width;
		temp.height = data->height;
		temp.max_gray = data->max_gray;
		temp.pixels = allocate_dynamic_matrix(temp.width, temp.height);
		for(i=0; i<temp.width; i++)
		{
			for(j=0; j<temp.width; j++)
			{
				temp.pixels[i][j] = data->pixels[i][j];
			}
		}
		down_sample_pic(&temp,o,data);
    }

	for(i=0; i<4; i++)
	{
		double temp = 4*o + i;
		temp = temp/4;
		double sigma = sigma0 * pow(2,temp);
		printf("%f\n",sigma);
		gauss = generate_normalised_gaussian_filter(sigma,5);
	    if(i==0)
		{
	    	h0 = convolve(data->pixels, data->width, data->height, gauss, 5);
		}
		if(i==1)
		{
			h1 = convolve(data->pixels, data->width, data->height, gauss, 5);
		}
		if(i==2)
		{
			h2 = convolve(data->pixels, data->width, data->height, gauss, 5);
		}
		if(i==3)
		{
			h3 = convolve(data->pixels, data->width, data->height, gauss, 5);
		}
	}
	printf("...Applied filters \n");
	for(i=0; i<data->width-4; i++)
    {
        for(j=0; j<data->height-4; j++)
        {
    	    h0[i][j] = h1[i][j] - h0[i][j];
    	    h1[i][j] = h2[i][j] - h1[i][j];
    	    h2[i][j] = h3[i][j] - h2[i][j];
        }
    }
/*    double max0=h0[0][0], max1=h1[0][0], max2=h2[0][0], min0=h0[0][0], min1=h1[0][0], min2=h2[0][0];
	for(i=0; i<data->width-4; i++)
    {
        for(j=0; j<data->height-4; j++)
        {
            if(h0[i][j]>max0) max0 = h0[i][j];
            if(h0[i][j]<min0) min0 = h0[i][j];
            if(h1[i][j]>max1) max1 = h1[i][j];
            if(h1[i][j]<min1) min1 = h1[i][j];
            if(h2[i][j]>max2) max2 = h2[i][j];
            if(h2[i][j]<min2) min2 = h2[i][j];
        }
    }
	for(i=0; i<data->width-4; i++)
    {
        for(j=0; j<data->height-4; j++)
        {
            h0[i][j] = (h0[i][j]-min0)/(max0-min0);
            h1[i][j] = (h1[i][j]-min1)/(max1-min1);
            h2[i][j] = (h2[i][j]-min2)/(max2-min2);
        }
    }*/
	PGMData result = make_PGM(h0, data->width-4, data->height-4);
    writePGM("h0.pgm",&result,1);
    result = make_PGM(h1, data->width-4, data->height-4);
    writePGM("h1.pgm",&result,1);
    result = make_PGM(h2, data->width-4, data->height-4);
    writePGM("h2.pgm",&result,1);
    printf("... Wrote images \n");
    double min = -67;
    double max = 100;
    double count[11];
    for(i=0; i<11; i++) count[i]=0;
    double **D_cap = allocate_dynamic_matrix_double(data->width-5,data->height-5);
    for(i=1; i<data->width-5; i++)
    {
        for(j=1; j<data->height-5; j++)
        {
        	double diff2_D[3][3];
        	double diff_D[3];
        	diff2_D[0][0] = (h1[i-1][j]+h1[i+1][j]-2*h1[i][j]);
        	diff2_D[1][1] = (h1[i][j-1]+h1[i][j+1]-2*h1[i][j]);
        	diff2_D[2][2] = 0.7071*(h2[i][j]+h0[i][j]-2*h1[i][j]);
            diff2_D[0][1] = diff2_D[1][0] = 0.25*(h1[i+1][j+1] - h1[i+1][j-1] - h1[i-1][j+1] + h1[i-1][j-1]);
            diff2_D[0][2] = diff2_D[2][0] = 0.25*0.8409*(h2[i+1][j] - h0[i+1][j] - h2[i-1][j] + h0[i-1][j]);
            diff2_D[1][2] = diff2_D[2][1] = 0.25*0.8409*(h2[i][j+1] - h2[i][j-1] - h0[i][j+1] + h0[i][j-1]);
            diff_D[0] = 0.5*(h1[i+1][j] - h1[i-1][j]);
            diff_D[1] = 0.5*(h1[i][j+1] - h1[i][j-1]);
            diff_D[2] = 0.5*0.8409*(h2[i][j] - h0[i][j]);

            double inverse[3][3];
            find_inverse(diff2_D,inverse);
            double x_dash =-(inverse[0][0]*diff_D[0] + inverse[0][1]*diff_D[1] + inverse[0][2]*diff_D[2]);
            double y_dash = -(inverse[1][0]*diff_D[0] + inverse[1][1]*diff_D[1] + inverse[1][2]*diff_D[2]);
            double s_dash = -(inverse[2][0]*diff_D[0] + inverse[2][1]*diff_D[1] + inverse[2][2]*diff_D[2]);

            D_cap[i][j] = h1[i][j] + 0.5 * (diff_D[0]*x_dash + diff_D[1]*y_dash + diff_D[2]*s_dash);
            double tr = diff2_D[0][0] + diff2_D[1][1];
            double det = diff2_D[0][0]*diff2_D[1][1] - diff2_D[1][0]*diff2_D[0][1];
            double ratio = tr*tr/(det+0.000001);
            if(ratio<20) D_cap[i][j] = 100;
        }
    }
    printf("Cmae najbjc\n");
    for(i=1; i<data->width-5; i++)
    {
    	for(j=1; j<data->height-5; j++)
    	{
            if(D_cap[i][j]>min) min=D_cap[i][j];
            if(D_cap[i][j]<max) max=D_cap[i][j];
            if(D_cap[i][j]>0.0 && D_cap[i][j]<=0.1)  count[0]++;
            if(D_cap[i][j]>0.1 && D_cap[i][j]<=0.2)  count[1]++;
            if(D_cap[i][j]>0.2 && D_cap[i][j]<=0.3)  count[2]++;
            if(D_cap[i][j]>0.3 && D_cap[i][j]<=0.4)  count[3]++;
            if(D_cap[i][j]>0.4 && D_cap[i][j]<=0.5)  count[4]++;
            if(D_cap[i][j]>0.54 && D_cap[i][j]<=0.56)  count[5]++;
            if(D_cap[i][j]>0.6 && D_cap[i][j]<=0.7)  count[6]++;
            if(D_cap[i][j]>0.7 && D_cap[i][j]<=0.8)  count[7]++;
            if(D_cap[i][j]>0.8 && D_cap[i][j]<=0.9)  count[8]++;
            if(D_cap[i][j]>0.9 && D_cap[i][j]<=1.0)  count[9]++;

            if(D_cap[i][j]==100) count[10]++;
            if(fabs(D_cap[i][j])>0.1 && fabs(D_cap[i][j])<0.7)
            {
            	key_loc[loc_count]=i;
            	loc_count++;
            	key_loc[loc_count]=j;
            	loc_count++;
            	key_mag[mag_count]=h1[i][j];
            	mag_count++;
            }
    	}
    }
    for(i=0; i<11; i++)
    {
    	printf("%f ",count[i]);
    }

    printf("\nGenerated %d key points between %f and %f \n",mag_count,max,min);
    return mag_count;
}
