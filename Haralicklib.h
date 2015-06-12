/* Written by:  R R Iyer
 * Institution: ISI Kolkata
 * Date:        22 May 2015
*/



#include <math.h>
/* Function to calculate Haralick parameters from the concurrence matrix
 * Argurments: P: The cooccurance matrix
 *             data_max_gray: max gray value
 *             fx: The array to store the result*/
double *calculate_haralick_parameters(double **P, int data_max_gray, double fx[13])
{

//	printf("Calculating Haralick parameters\n");
    /*Normalise the matrix*/
	double N=0;
    int i,j,k;
    for(i=0;i<data_max_gray;i++)
    {
    	for(j=0;j<data_max_gray;j++)
    	{
    		N = N + P[i][j];
    	}
    }
    for(i=0;i<data_max_gray;i++)
    {
    	for(j=0;j<data_max_gray;j++)
    	{
    		P[i][j] = P[i][j]/N;

    	}

    }
//    printf("N = %f\tmax gray = %d\n",N,data_max_gray);
    double Px[17000], Py[17000], ux, uy,u, sx, sy, Pxplusy[33000], Pxminusy[33000];

    for(i=0; i<data_max_gray; i++)
    {
    	Px[i]=0;
       	Py[i]=0;
       	for(j=0; j<data_max_gray; j++)
       	{
       		Px[i] = Px[i]+P[i][j];
       		Py[i] = Py[i]+P[j][i];
       	}
    }

    ux=0; uy=0, u=0;
    for(i=0;i<data_max_gray;i++)
    {
      	ux = ux + i*Px[i];
       	uy = uy + i*Py[i];
    }

    for(k=0; k<2*data_max_gray-1; k++)
    {
       	Pxplusy[k]=0;
       	Pxminusy[k]=0;
    }



    for(i=0; i<data_max_gray; i++)
    {
       	for(j=0; j<data_max_gray; j++)
       	{
       		Pxplusy[i+j] = Pxplusy[i+j]+P[i][j];
       		if(i>=j)
       			Pxminusy[i-j] = Pxminusy[i-j]+P[i][j];
       		else
       			Pxminusy[j-i] = Pxminusy[j-i]+P[i][j];

      	}
    }

    u = 0.5* ux+ 0.5*uy;

    sx = 0; sy = 0;
    for(i=0; i<data_max_gray; i++)
    {
      	sx = sx+Px[i]*(i-ux)*(i-ux);
       	sy = sy+Py[i]*(i-uy)*(i-uy);
    }


    sx = sqrt(sx);
    sy = sqrt(sy);


    double HXY1, HXY2, HX=0, HY=0;
    for(i=0; i<data_max_gray; i++)
    {
        HX = HX - Px[i]*log(Px[i]+0.00000001);
  	    HY = HY - Py[i]*log(Py[i]+0.00000001);
    }


    HXY1 = 0; HXY2 = 0;
    for(i=0; i<data_max_gray; i++)
    {
     	for(j=0; j<data_max_gray; j++)
     	{
   		    HXY1 = HXY1 - P[i][j]*log(Px[i]*Py[j] +0.00000001);
            HXY2 = HXY2 - Px[i]*Py[j]*log(Px[i]*Py[j] +0.00000001);

       	}
   	}



//    printf("ux = %g\t\tuy = %g\t\tsx = %g\t\tsy=%g\t\tu = %f\n",ux,uy,sx,sy,u);
//    printf("HXY1 = %g\t\tHXY2 = %g\t\tHX = %g\t\tHY = %g\n",HXY1,HXY2,HX,HY);


    for(i=0; i<14; i++)
        fx[i]=0;
    for(i=0; i<data_max_gray; i++)
    {
    	for(j=0;j<data_max_gray; j++)
    	{
    		fx[0] = fx[0] + P[i][j]*P[i][j];
    		fx[2] = fx[2] + (i-ux)*(i-uy)*P[i][j]/(sx*sy +0.00000001);
    		fx[3] = fx[3] + (i-u)*(i-u)*P[i][j];
    		fx[4] = fx[4] + P[i][j]/(1+(i-j)*(i-j));
    		fx[8] = fx[8] - P[i][j]*log(P[i][j]+0.00000001);
    		fx[10] = fx[10] - Pxminusy[i]*log(Pxminusy[i]+0.00000001);
    	}
    }

    for(k=0; k<data_max_gray; k++)
    {
    	for(i=0; i<data_max_gray; i++)
    	{
    		for(j=0; j<data_max_gray; j++)
    		{
    			if(i-j == k)
    				fx[1]=fx[1]+ P[i][j]*k*k;
    		}
    	}
    }


    for(i=0; i<2*data_max_gray-1; i++)
    {
    	fx[5] = fx[5] + i*Pxplusy[i];
    	fx[7] = fx[7] - Pxplusy[i]*log(Pxplusy[i]+0.00000001);
    }

    for(i=0; i<2*data_max_gray-1; i++)
    {
    	fx[6] = fx[6] + (i-fx[5])*(i-fx[5])*Pxplusy[i];
    }

    for(i=0; i<data_max_gray; i++)
    {
    	double temp=0;
    	for(j=0; j<data_max_gray; j++)
    	{
    		temp = temp + j*Pxminusy[j];
    	}
    	fx[9] = fx[9] + (i-temp)*(i-temp)*Pxminusy[i];
    }

    double t;
    if(HX>HY) t=HX;
    else t=HY;
    fx[11] = (fx[8]-HXY1)/(t+0.00000001);
    fx[12] = sqrt(1-exp(-2*(HXY2 - fx[8])));

    return fx;

}

/* To create co-occurance matrix
 * Arguements: data: The PGM image
 *             delta: The distance
 *             angle: The config in degrees
 */
double *create_cooccurance_matrix(PGMData *data, int delta, int angle, double f[13])
{
//    printf("Creating Co-occurance matrix\n");
	int delx, dely;
    if(angle == 0)
    {
    	delx = 0;
    	dely = delta;
    }
    if(angle == 45)
    {
    	delx = -delta;
    	dely = delta;
    }
    if(angle == 90)
    {
    	delx=-delta;
        dely = 0;
    }
    if(angle == 135)
    {
    	delx = -delta;
    	dely = -delta;
    }


    double **P;
    P = allocate_dynamic_matrix_double((data->max_gray), (data->max_gray));

    int i,j,x,y;

    for(i=0;i<data->max_gray;i++)
    	for(j=0;j<data->max_gray;j++)
    		P[i][j]=0;

    int temp1, temp2;
    for(i=0;i<data->max_gray;i++)
    {
    	for(j=0;j<data->max_gray;j++)
    	{
    		for(x=0;x<data->width;x++)
    		{
    			for(y=0;y<data->height;y++)
    			{
                    temp1 = (x+delx >= 0)? x+delx:-1;
                    temp2 = (y+dely >= 0)? y+dely:-1;
                    if(temp1>=0 && temp2 >=0)
    				if((data->pixels[x][y]==i && data->pixels[temp1][temp2]==j))// || (data->pixels[x][y]==j && data->pixels[temp1][temp2]==i))
    				    P[i][j]=P[i][j]+1;
    			}
    		}
    	}
    }


//    printf("Calculating Haralick parameters\n");
    calculate_haralick_parameters(P,data->max_gray,f);
/*
    printf("Haralicks PARAMETERS\n");
    for(i=0; i<13; i++)
    {
    	printf("f%d = %g\t\t",i+1, f[i]);
    	if(i%3==2) printf("\n");
    }
*/
    deallocate_dynamic_matrix_double(P, data->max_gray);

//    printf("\n");
    return f;
}


void write_haralick(int no, double *f, char *name, int n)
{
    FILE *train_file;
    int i;
    train_file = fopen(name, "ab");
    if (train_file == NULL) {
        perror("Cannot open file to write\n");
        exit(EXIT_FAILURE);
    }
    fprintf(train_file,"%d ",no);
    for(i=0; i<n; i++)
    {
    	fprintf(train_file,"%f ",f[i]);
    }
    fprintf(train_file,"\n");
}


double *calculate_haralick_parameters_RIVLBP(double **P, double fx[13], int m, int n)
{

//	printf("Calculating Haralick parameters\n");
    /*Normalise the matrix*/
	double N=0;
    int i,j,k;
    for(i=0;i<m;i++)
    {
    	for(j=0;j<n;j++)
    	{
    		N = N + P[i][j];
    	}
    }
    for(i=0;i<m;i++)
    {
    	for(j=0;j<n;j++)
    	{
    		P[i][j] = P[i][j]/N;

    	}

    }
//    printf("N = %f\tmax gray = %d\n",N,data_max_gray);
    double Px[17000], Py[17000], ux, uy,u, sx, sy, Pxplusy[33000], Pxminusy[33000];

    for(i=0; i<m; i++)
    {
    	Px[i]=0;
       	Py[i]=0;
       	for(j=0; j<n; j++)
       	{
       		Px[i] = Px[i]+P[i][j];
       		Py[i] = Py[i]+P[j][i];
       	}
    }

    ux=0; uy=0, u=0;
    for(i=0;i<m;i++)
    {
      	ux = ux + i*Px[i];
    }
    for(i=0; i<n; i++)
    {
       	uy = uy + i*Py[i];
    }

    for(k=0; k<m+n-1; k++)
    {
       	Pxplusy[k]=0;
       	Pxminusy[k]=0;
    }



    for(i=0; i<m; i++)
    {
       	for(j=0; j<n; j++)
       	{
       		Pxplusy[i+j] = Pxplusy[i+j]+P[i][j];
       		if(i>=j)
       			Pxminusy[i-j] = Pxminusy[i-j]+P[i][j];
       		else
       			Pxminusy[j-i] = Pxminusy[j-i]+P[i][j];

      	}
    }

    u = 0.5* ux+ 0.5*uy;

    sx = 0; sy = 0;
    for(i=0; i<m; i++)
    {
      	sx = sx+Px[i]*(i-ux)*(i-ux);
    }
    for(i=0; i<n; i++)
    {
       	sy = sy+Py[i]*(i-uy)*(i-uy);
    }


    sx = sqrt(sx);
    sy = sqrt(sy);


    double HXY1, HXY2, HX=0, HY=0;
    for(i=0; i<m; i++)
    {
        HX = HX - Px[i]*log(Px[i]+0.00000001);
    }
    for(i=0; i<n; i++)
    {
  	    HY = HY - Py[i]*log(Py[i]+0.00000001);
    }


    HXY1 = 0; HXY2 = 0;
    for(i=0; i<m; i++)
    {
     	for(j=0; j<n; j++)
     	{
   		    HXY1 = HXY1 - P[i][j]*log(Px[i]*Py[j] +0.00000001);
            HXY2 = HXY2 - Px[i]*Py[j]*log(Px[i]*Py[j] +0.00000001);
       	}
   	}



//    printf("ux = %g\t\tuy = %g\t\tsx = %g\t\tsy=%g\t\tu = %f\n",ux,uy,sx,sy,u);
//    printf("HXY1 = %g\t\tHXY2 = %g\t\tHX = %g\t\tHY = %g\n",HXY1,HXY2,HX,HY);


    for(i=0; i<14; i++)
        fx[i]=0;
    for(i=0; i<m; i++)
    {
    	for(j=0;j<n; j++)
    	{
    		fx[0] = fx[0] + P[i][j]*P[i][j];
    		fx[2] = fx[2] + (i-ux)*(i-uy)*P[i][j]/(sx*sy +0.00000001);
    		fx[3] = fx[3] + (i-u)*(i-u)*P[i][j];
    		fx[4] = fx[4] + P[i][j]/(1+(i-j)*(i-j));
    		fx[8] = fx[8] - P[i][j]*log(P[i][j]+0.00000001);
    		fx[10] = fx[10] - Pxminusy[i]*log(Pxminusy[i]+0.00000001);
    	}
    }


    int x;
    if(m>n) x=m;
    else x=n;
    for(k=0; k<x; k++)
    {
    	for(i=0; i<m; i++)
    	{
    		for(j=0; j<n; j++)
    		{
    			if(i-j == k)
    				fx[1]=fx[1]+ P[i][j]*k*k;
    		}
    	}
    }

    for(i=0; i<m+n-1; i++)
    {
    	fx[5] = fx[5] + i*Pxplusy[i];
    	fx[7] = fx[7] - Pxplusy[i]*log(Pxplusy[i]+0.00000001);
    }

    for(i=0; i<m+n-1; i++)
    {
    	fx[6] = fx[6] + (i-fx[5])*(i-fx[5])*Pxplusy[i];
    }

    for(i=0; i<m; i++)
    {
    	double temp=0;
    	for(j=0; j<n; j++)
    	{
    		temp = temp + j*Pxminusy[j];
    	}
    	fx[9] = fx[9] + (i-temp)*(i-temp)*Pxminusy[i];
    }

    double t;
    if(HX>HY) t=HX;
    else t=HY;
    fx[11] = (fx[8]-HXY1)/(t+0.00000001);
    fx[12] = sqrt(1-exp(-2*(HXY2 - fx[8])));

    return fx;

}

void write_SVM_format(double f[], char name[], int clas, int n)
{
    FILE *train_file;
    int i;
    train_file = fopen(name, "ab");
    if (train_file == NULL) {
        perror("Cannot open file to write\n");
        exit(EXIT_FAILURE);
    }
    fprintf(train_file,"%d \t",clas);
    for(i=0; i<n-1; i++)
    {
    	fprintf(train_file,"%d:%f \t",i+1,f[i]);
    }
    fprintf(train_file,"%d:%f\n",n,f[n-1]);
}
