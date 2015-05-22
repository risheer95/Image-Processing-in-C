/* Written by:  R R Iyer
 * Institution: ISI Kolkata
 * Date:        22 May 2015
*/


int find_least_combination(int number, int size)
{
	int min = number;
	int i;
	for(i=0; i<size; i++)
	{
		int temp = number&0x01;
		number = number>>1;
		number = number + temp*pow(2,size-1);
		if(number<min) min=number;
	}
	return min;
}
/*Function that returns the LBP matrix of a given image
 * Arguments: data: The PGM data
 *            radius: The radius of window to be considered
 *            no_of_points: No of points in the window*/
PGMData calculate_LBP(PGMData *data, int radius, int no_of_points)
{
	printf("Finding the LBP Matrix\n");
	int i,j,k, jdash, kdash;
	double del_theta = 2*3.14159265/no_of_points, delx, dely;
	int **lbp;
	lbp = allocate_dynamic_matrix(data->width, data->height);
	for(j=0+radius; j<data->width-radius; j++)
	{
	    for(k=0+radius; k<data->height-radius; k++)
	    {
	    	lbp[j-radius][k-radius]=0;

		    for(i=0; i<no_of_points; i++)
	        {
		        dely = ceil(sin(del_theta*i)*radius);
		        delx = ceil(cos(del_theta*i)*radius);
		        lbp[j-radius][k-radius]=lbp[j-radius][k-radius]<<1;
		        if(j+delx>=0 && j+delx<data->width)
		        	jdash=j+delx;
		        else
		        	jdash = j;

		        if(k+dely>=0 && k+dely<data->height)
		        	kdash=k+dely;
		        else
		        	kdash = k;

		        if(data->pixels[jdash][kdash]>data->pixels[j][k])
		        {
		        	lbp[j-radius][k-radius]=lbp[j-radius][k-radius]|0x01;
		        }
//		        printf("%d %d %d\n",j,k,lbp[j-1][k-1]);
	        }

		    lbp[j][k] = find_least_combination(lbp[j][k], no_of_points);
	    }
	}

	PGMData result;
	result.width = data->width-2*radius;
	result.height = data->height-2*radius;
	result.max_gray = pow(2,no_of_points);

	result.pixels = allocate_dynamic_matrix(data->width, data->height);
	for(i=0; i<result.width; i++)
	{
		for(j=0; j<result.height; j++)
		{
			result.pixels[i][j] = lbp[i][j];
//			printf("%d\t",lbp[i][j]);
		}
//		printf("\n");
	}
	char ver=1;
	writePGM("LBP.pgm",&result,ver);
	deallocate_dynamic_matrix(lbp, result.height);

	return result;
}
