/* Written by:  R R Iyer
 * Institution: ISI Kolkata
 * Date:        22 May 2015
*/

/*To find the least combination of a combination of bits
 * after rotating it right and comparing with the rest.
 */

#define S(X1,X2) ((X1>=X2)?1:0)
#define f(x,y) ((x==y)?1:0)
#define t(x,c) ((x>=c)?1:0)
/* Function to find least combo of numbers by ROR
 * Arguments: number and size
 */
int find_least_combination(int number, int size)
{
	int min = number;
	int i;
	for(i=0; i<size; i++)
	{
		int temp = number&0x0001;
		number = number>>1;
		number = number + temp*pow(2,size-1);
		if(number<min) min=number;
	}
	return min;
}

/*Function to calculate uniformity
 * Arguments: LBP: Vector of size P+1 where the first element is the
 *                 central element
 *            size: P
 */
int uniformity(int *LBP, int size)
{
	int i, U=abs(S(LBP[size-1],LBP[0]) - S(LBP[1],LBP[0]));
	for(i=1; i<size; i++)
	{
		U = U + abs(S(LBP[i],LBP[0]) - S(LBP[i-1],LBP[0]));
	}
return U;
}

/* Function to find the distribution of the LBP
 * Argument: The PGM image containing the LBP of the original image
 *           k
 */
int h_k(PGMData *data, int k)
{
	int i,j, h=0;
	for(i=0; i<data->width; i++)
	{
		for(j=0; j<data->height; j++)
		{
			h = h + f(data->pixels[i][j],k);
		}
	}
	return h;
}

/* Function to find the histogram of the given image
 * Arguments: data: The image whose histogram is to be evaluated
 *            h: Empty array to store the result in*/
void calculate_histogram(PGMData *data, int *h)
{
	int i;
	for(i=0; i<data->max_gray; i++)
	{
		h[i] = h_k(data,i);
//		printf("%d %d\n",i,h[i]);
	}
//printf("\n");
}
/*Function that returns the LBP matrix of a given image
 * Arguments: data: The PGM data
 *            radius: The radius of window to be considered
 *            no_of_points: No of points in the window
 *
 */
PGMData calculate_LBP(PGMData *data, int radius, int no_of_points)
{
	printf("Finding the LBP Matrix\n");
	int i,j,k, jdash, kdash;
	double del_theta = 2*3.14159265/no_of_points, delx, dely;
	int **lbp;
	lbp = allocate_dynamic_matrix(data->width, data->height);
	int max=0;
	for(j=0+radius; j<(data->width-radius); j=j+2*radius+1)
	{
	    for(k=0+radius; k<(data->height-radius); k=k+2*radius+1)
	    {
	    	lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=0;

		    for(i=0; i<no_of_points; i++)
	        {
		        dely = ceil(sin(del_theta*i)*radius);
		        delx = ceil(cos(del_theta*i)*radius);

		        lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]<<1;
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
		        	lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]|0x01;
		        }
	        }

//		    lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)] = find_least_combination(lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)], no_of_points);

		    if(lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]>max)
		    	max = lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)];
	    }
	}

	PGMData result;
	result.width = (data->width-2*radius)/(2*radius+1);
	result.height = (data->height-2*radius)/(2*radius+1);
	result.max_gray = max;

	result.pixels = allocate_dynamic_matrix(data->width, data->height);
	for(i=0; i<result.width; i++)
	{
		for(j=0; j<result.height; j++)
		{
			result.pixels[i][j] = lbp[i][j];
		}
	}
	char ver=1;
	deallocate_dynamic_matrix(lbp, result.height);

	writePGM("LBP.pgm",&result,ver);

	return result;
}

/*Function that returns the Rotational Invariant Uniform LBP matrix of a given image
 * Arguments: data: The PGM data
 *            radius: The radius of window to be considered
 *            no_of_points: No of points in the window
 */

PGMData calculate_LBPriu2(PGMData *data, int radius, int no_of_points)
{
	printf("Finding the Uniform Rotational Invariant LBP Matrix\n");
	int i,j,k, jdash, kdash;
	double del_theta = 2*3.14159265/no_of_points, delx, dely;
	int **lbp;
	lbp = allocate_dynamic_matrix(data->width, data->height);
	int temp[16];
	for(j=0+radius; j<(data->width-radius); j=j+2*radius+1)
	{
	    for(k=0+radius; k<(data->height-radius); k=k+2*radius+1)
	    {
	    	lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=0;

	    	temp[0] = data->pixels[j][k];
	    	for(i=0; i<no_of_points; i++)
	        {
		        dely = ceil(sin(del_theta*i)*radius);
		        delx = ceil(cos(del_theta*i)*radius);

		        if(j+delx>=0 && j+delx<data->width)
		        	jdash=j+delx;
		        else
		        	jdash = j;

		        if(k+dely>=0 && k+dely<data->height)
		        	kdash=k+dely;
		        else
		        	kdash = k;

		        temp[i+1] = data->pixels[jdash][kdash];
		    }

	    	int U = uniformity(temp, no_of_points);

	    	if(U<=2)
	    	{
	    		for(i=0; i<no_of_points; i++)
	    			lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)] = lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)] + S(temp[i+1],temp[0]);
	    	}

	    	else
	    		lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)] = no_of_points+1;
	    }
	}

	PGMData result;
	result.width = (data->width-2*radius)/(2*radius+1);
	result.height = (data->height-2*radius)/(2*radius+1);
	result.max_gray = no_of_points+2;

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
	deallocate_dynamic_matrix(lbp, result.height);

	writePGM("LBPriu2.pgm",&result,ver);
	return result;

}

/*Function that returns the LBP matrix of a given image
 * Arguments: data: The PGM data
 *            radius: The radius of window to be considered
 *            no_of_points: No of points in the window
 *                       f: An array to store the Haralick parameters
 */
double *calculate_completed_LBP(PGMData *data, int radius, int no_of_points, double f[39])
{
	printf("\nFinding the Completed LBP Matrix\n");
	int i,j,k, jdash, kdash;
	double average=0;
	for(i=0; i<data->width; i++)
	{
		for(j=0; j<data->height; j++)
		{
			average = average + data->pixels[i][j];
		}
	}
	average = average/(data->width * data->height);
	double del_theta = 2*3.14159265/no_of_points, delx, dely;
	int **clbp_s, **clbp_m, **clbp_c;
	clbp_s = allocate_dynamic_matrix((data->width-2*radius)/(2*radius+1), (data->height-2*radius)/(2*radius+1));
	clbp_m = allocate_dynamic_matrix((data->width-2*radius)/(2*radius+1), (data->height-2*radius)/(2*radius+1));
	clbp_c = allocate_dynamic_matrix((data->width-2*radius)/(2*radius+1), (data->height-2*radius)/(2*radius+1));

	int Sp[16], Mp[16];
	for(j=0+radius; j<(data->width-radius); j=j+2*radius+1)
	{
	    for(k=0+radius; k<(data->height-radius); k=k+2*radius+1)
	    {
	    	clbp_s[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=0;
	    	clbp_m[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=0;
	    }
	}
	for(j=0+radius; j<(data->width-radius); j=j+2*radius+1)
	{
	    for(k=0+radius; k<(data->height-radius); k=k+2*radius+1)
	    {
	    	for(i=0; i<no_of_points; i++)
	        {
		        dely = ceil(sin(del_theta*i)*radius);
		        delx = ceil(cos(del_theta*i)*radius);

		        if(j+delx>=0 && j+delx<data->width)
		        	jdash=j+delx;
		        else
		        	jdash = j;

		        if(k+dely>=0 && k+dely<data->height)
		        	kdash=k+dely;
		        else
		        	kdash = k;

		        Sp[i] = (data->pixels[jdash][kdash]>=data->pixels[j][k])?1:0;
		        Mp[i] = (data->pixels[jdash][kdash]>=average)?1:0;

		        clbp_s[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]+=Sp[i];
		        clbp_m[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]+=Mp[i];
		    }
	        clbp_c[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=t(data->pixels[j][k],average);
	    }
	}

	PGMData result;
	result.width = (data->width-2*radius)/(2*radius+1);
	result.height = (data->height-2*radius)/(2*radius+1);
	result.max_gray = data->max_gray;

	result.pixels = allocate_dynamic_matrix(result.width, result.height);

	for(i=0; i<result.width; i++)
	{
		for(j=0; j<result.height; j++)
		{
			result.pixels[i][j] = clbp_c[i][j];
		}
	}
    double f1[13],f2[13],f3[13];

    create_cooccurance_matrix(&result,1,0,f1);
//	calculate_histogram(&result,h_c);
    deallocate_dynamic_matrix(result.pixels,result.height);
	result.width = (data->width-2*radius)/(2*radius+1);
	result.height = (data->height-2*radius)/(2*radius+1);
	result.max_gray = data->max_gray;

	result.pixels = allocate_dynamic_matrix(result.width, result.height);

	for(i=0; i<result.width; i++)
	{
		for(j=0; j<result.height; j++)
		{
			result.pixels[i][j] = clbp_m[i][j];
//			printf("%d ",clbp_m[i][j]);
		}
	}
    create_cooccurance_matrix(&result,1,0,f2);

//	calculate_histogram(&result,h_m);
//	writePGM("clbpm.pgm",&result,ver);
    deallocate_dynamic_matrix(result.pixels,result.height);
	result.width = (data->width-2*radius)/(2*radius+1);
	result.height = (data->height-2*radius)/(2*radius+1);
	result.max_gray = data->max_gray;

	result.pixels = allocate_dynamic_matrix(result.width, result.height);

	for(i=0; i<result.width; i++)
	{
		for(j=0; j<result.height; j++)
		{
			result.pixels[i][j] = clbp_s[i][j];
		}
	}
    create_cooccurance_matrix(&result,1,0,f3);

//	calculate_histogram(&result,h_s);
    for(i=0; i<13; i++)
    {
    	f[i] = f1[i];
    	f[i+13] = f2[i];
    	f[i+26] = f3[i];
    }

	deallocate_dynamic_matrix(clbp_s, result.height);
	deallocate_dynamic_matrix(clbp_c, result.height);
	deallocate_dynamic_matrix(clbp_m, result.height);

	deallocate_dynamic_matrix(result.pixels, result.height);

//	writePGM("Histogram.pgm",&result,ver);
	return f;
}


/*Function that returns the Co-occurrence of Adjacent LBP matrix of a given image
 * Arguments: data: The PGM data
 *            radius: The radius of window to be considered
 *            no_of_points: No of points in the window
 *            f: To store the Haralick parameters
 */
double *calculate_CoALBP(PGMData *data, int radius, int no_of_points, double f[13])
{

	printf("\nFinding the LBP Matrix\n");
	int i,j,k, jdash, kdash;
	double del_theta = 2*3.14159265/no_of_points, delx, dely;
	int **lbp;
	lbp = allocate_dynamic_matrix((data->width-2*radius)/(2*radius+1), (data->height-2*radius)/(2*radius+1));

	for(j=0+radius; j<(data->width-radius); j=j+2*radius+1)
	{
	    for(k=0+radius; k<(data->height-radius); k=k+2*radius+1)
	    {
	    	lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=0;

		    for(i=0; i<no_of_points; i++)
	        {
		        dely = ceil(sin(del_theta*i)*radius);
		        delx = ceil(cos(del_theta*i)*radius);

		        lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]<<1;
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
		        	lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]|0x01;
		        }
	        }

	    }
	}

    PGMData result;
    result.width = (data->width - 2*radius)/(2*radius+1);
    result.height = (data->height - 2*radius)/(2*radius+1);
    int max=lbp[0][0];
	for(j=0; j<(data->width - 2*radius)/(2*radius+1)-1; j++)
	{
		for(k=0; k<(data->height - 2*radius)/(2*radius+1)-1; k++)
		{
			if(max<lbp[j][k]) max = lbp[j][k];
		}
	}

    result.max_gray = max;
//    printf("\nMax = %d\n",max);
    result.pixels = allocate_dynamic_matrix(result.width, result.height);
	for(j=0; j<(data->width - 2*radius)/(2*radius+1)-1; j++)
	{
		for(k=0; k<(data->height - 2*radius)/(2*radius+1)-1; k++)
		{
			result.pixels[j][k] = (int) lbp[j][k];
		}
	}
    create_cooccurance_matrix(&result,1,0,f);
    deallocate_dynamic_matrix(result.pixels, result.height);
    deallocate_dynamic_matrix(lbp, (data->height-2*radius)/(2*radius+1));

    return f;
}

/*  Function to generate the Map M such that
 *  P_theta_ri(r, del_r) = M*P_theta(r,delr);
 */
int **generate_M(int N)
{
	int i,j, idash, jdash, id=1;
	int **M = allocate_dynamic_matrix((int) pow(2,N),(int) pow(2,N));

	for(i=0; i<(int) pow(2,N); i++)
	{
		for(j=0; j<(int) pow(2,N); j++)
		{
              M[i][j] = (-1);
		}
	}

	for(i=0; i<(int) pow(2,N); i++)
	{
		for(j=0; j<(int)pow(2,N); j++)
		{
			if(M[i][j]== -1)
			{
				idash = i >> (N/2);
				jdash = j >> (N/2);
				M[i][j] = id;
				M[jdash][idash] = id;
				id++;
			}
		}
	}
    return M;
}

/*Function that returns the Co-occurrence of Adjacent LBP matrix of a given image
 * Arguments: data: The PGM data
 *            radius: The radius of window to be considered
 *            no_of_points: No of points in the window
 *            theta: The angle to be specified
 *            f: To store the Haralick parameters
 */
double* calculate_RIV_LBP(PGMData *data, int radius, int no_of_points, double theta, double f[13])
{

	printf("Finding the Rot Invariant LBP Matrix\n");
	int i,j,k, jdash, kdash;
	double **P = allocate_dynamic_matrix_double((data->width - 2*radius)/(2*radius+1),(data->height - 2*radius)/(2*radius+1));
	double del_theta = 2*3.14159265/no_of_points, delx, dely;
	int **lbp;
	lbp = allocate_dynamic_matrix((data->width - 2*radius)/(2*radius+1), (data->height - 2*radius)/(2*radius+1));

	for(j=0+radius; j<(data->width-radius); j=j+2*radius+1)
	{
	    for(k=0+radius; k<(data->height-radius); k=k+2*radius+1)
	    {
	    	lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=0;

		    for(i=0; i<no_of_points; i++)
	        {
		        dely = ceil(sin(del_theta*i + theta)*radius);
		        delx = ceil(cos(del_theta*i + theta)*radius);

		        lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]<<1;
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
		        	lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]=lbp[(j-radius)/(2*radius + 1)][(k-radius)/(2*radius + 1)]|0x01;
		        }
	        }
	    }
	}
	int **M = generate_M(no_of_points);
	for(j=0; j<(data->width - 2*radius)/(2*radius+1); j++)
		{
			for(k=0; k<(data->height - 2*radius)/(2*radius+1); k++)
			{
				int dx = ceil(radius*cos(theta));
				int dy = ceil(radius*sin(theta));

				int jdash,kdash;
				if(j+dx<0 || j+dx>=(data->width - 2*radius)/(2*radius+1))
					jdash=j;
				else jdash=j+dx;

				if(k+dy<0 || k+dy>=(data->height - 2*radius)/(2*radius+1))
					kdash=k;
				else kdash=k+dy;

				P[j][k] = (double) M[lbp[j][k]][lbp[jdash][kdash]];
			}
		}

    calculate_haralick_parameters_RIVLBP(P,f,(data->width - 2*radius)/(2*radius+1),(data->height - 2*radius)/(2*radius+1));
    deallocate_dynamic_matrix(lbp,(data->height - 2*radius)/(2*radius+1) );
    deallocate_dynamic_matrix_double(P, (data->height - 2*radius)/(2*radius+1));
    deallocate_dynamic_matrix(M, (int) pow(2,no_of_points));

    return f;
}
