/* Program to read a PGM file and find its Co-occurance matrix and Harlick parameters*/
/* Written by:  R R Iyer
 * Institution: ISI Kolkata
 * Date:        18 May 2015
*/

/*
 * Example of a P2 coded PGM file
 * a P5 file would have the pixel values in binary
 * Row value, Col Value, max pixel size are in ASCII decimal
 *
 * 1) Magic number (P2 or P5)
 * 2) Empty space
 * 3) Row numbers (width)
 * 4) Empty space
 * 5) Column numbers (height)
 * 6) Empty space
 * 7) Max pixel value
 * 8) Single white space
 * 9) Row wise values
 *
 * P2
 * # feep.pgm
 * 24 7
 * 15
 * 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 * 0  3  3  3  3  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15 15 15 15  0
 * 0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0 15  0
 * 0  3  3  3  0  0  0  7  7  7  0  0  0 11 11 11  0  0  0 15 15 15 15  0
 * 0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0  0  0
 * 0  3  0  0  0  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15  0  0  0  0
 * 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 *
 * For finding the co-occurance matrix P, we use the following formula
 *
 * P(i,j) = (sigma x=1 to M)(sigma y=1 to N){ 1 if I(x,y) = i and I(x+delx,y+dely)=j
 *                                          {  0 otherwise
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/* To get the upper 8 bits of a number */
#define UP8(num) (((num) & 0x0000FF00) >> 8)
/* To get the lower 8 bits of a number */
#define LO8(num) ((num) & 0x000000FF)

/* Structure of PGMdata */
typedef struct _PGMData
{
    int width;      // no of rows (width)
    int height;      // no of columns (height)
    int max_gray; // max pixel value
    int **pixels; // The 2D array containing all pixel values
}PGMData;

/*Dynamically allocate memory size for the matrix to store pixel values
 * Total size = row * col
 */
int **allocate_dynamic_matrix(int width, int height)
{
    int **ret_val;
    int i;

    ret_val = (int **)malloc(sizeof(int *) * width);
    if (ret_val == NULL)
    {
        perror("Memory allocation failure");
        exit(1);
    }

    for (i = 0; i < width; ++i)
    {
        ret_val[i] = (int *)malloc(sizeof(int) * height);
        if (ret_val[i] == NULL)
        {
            perror("Memory allocation failure");
            exit(1);
        }
    }

    return ret_val;
}


double **allocate_dynamic_matrix_d(int width, int height)
{
    double **ret_val;
    int i;

    ret_val = (double **)malloc(sizeof(double *) * width);
    if (ret_val == NULL)
    {
        perror("Memory allocation failure");
        exit(1);
    }

    for (i = 0; i < width; ++i)
    {
        ret_val[i] = (double *)malloc(sizeof(double) * height);
        if (ret_val[i] == NULL)
        {
            perror("Memory allocation failure");
            exit(1);
        }
    }

    return ret_val;
}

double *allocate_dynamic_matrix_a(int width)
{
    double *ret_val;
    //int i;

    ret_val = (double *)malloc(sizeof(double *) * width);
    if (ret_val == NULL)
    {
        perror("Memory allocation failure");
        exit(1);
    }

    return ret_val;
}
/* Deallocate memory for the matrix
 * Total size = row * col
 */
void deallocate_dynamic_matrix(int **mat, int row)
{
    int i;

    for (i = 0; i < row; ++i)
        free(mat[i]);
    free(mat);
}

/* As seen in the example, a PGM file may contain comments.
 * Comments start with # and continue for a single line.
 * If the character read is not a comment, go back to the
 * character and read it again
 */
void skip_comments(FILE *fp)
{
    int ch;
    char line[100];

    while ((ch = fgetc(fp)) != EOF && isspace(ch));

    if (ch == '#')
    {
        fgets(line, sizeof(line), fp);
        skip_comments(fp);
    }
    else
        fseek(fp, -1, SEEK_CUR);
}

/* To read a PGM file given a filename and storage space*/
PGMData* readPGM(const char *file_name, PGMData *data)
{
    FILE *pgm_file;     // File pointer
    char version[3];    // To get the version P2 or P5
    char ver;           // Flag for version
    int i, j;           // Counters/ Iteration variables
    int lower, upper;   // To read the upper 8 bits and lower 8 bits

    /* Opening file*/
    pgm_file = fopen(file_name, "rb");
    if (pgm_file == NULL)
    {
        perror("Cannot open file");
        exit(1);
    }
    //skip_comments(pgm_file);
    /*Reading the File version and assigning appropriate flags */
    fgets(version, sizeof(version), pgm_file);
    if (!strcmp(version, "P5"))
        ver = 1;
    else if(!strcmp(version,"P2"))
       	ver = 0;
    else
    {
    	fprintf(stderr, "Wrong file type!\n");
        exit(1);
    }
    printf("Opened file\n");

    /* Reading width, height, max gray values */
    skip_comments(pgm_file);
    fscanf(pgm_file, "%d", &data->height);
    skip_comments(pgm_file);
    fscanf(pgm_file, "%d", &data->width);
    skip_comments(pgm_file);
    fscanf(pgm_file, "%d", &data->max_gray);
    fgetc(pgm_file);


    /* Allocating memory dynamically */
    data->pixels = allocate_dynamic_matrix(data->width, data->height);

    /* If version is P5, reading the values as binary */
    if(ver)
    {
        if (data->max_gray > 255)
            for (i = 0; i < data->width; ++i)
                for (j = 0; j < data->height; ++j)
                {
                    upper = fgetc(pgm_file);
                    lower = fgetc(pgm_file);
                    data->pixels[i][j] = (upper << 8) + lower;
                }
        else
            for (i = 0; i < data->width; ++i)
                for (j = 0; j < data->height; ++j)
                {
                    lower = fgetc(pgm_file);
                    data->pixels[i][j] = lower;
                }
    }

    /* If version is P2 reading values as ASCII */
    else
    {
        for (i = 0; i < data->width; ++i)
            for (j = 0; j < data->height; ++j)
                fscanf(pgm_file, "%d", &data->pixels[i][j]);
    }

    fclose(pgm_file);
    return data;
}

/* Normalize all values to given value */
PGMData* normalisePGM(PGMData *data, int value)
{
	int i, j;
	if(value>data->max_gray)
	    for(i=0; i<data->width; i++)
		    for(j=0; j<data->height; j++)
			    data->pixels[i][j]*=(value/data->max_gray);

	else
		for(i=0; i<data->width; i++)
		    for(j=0; j<data->height; j++)
			    data->pixels[i][j]/=(data->max_gray/value);

	data->max_gray = value;
	printf("Successfully read file\n");
	return data;
}

/* To write a structure into a PGM file */
void writePGM(const char *filename, PGMData *data, char ver)
{
    FILE *pgm_file;
    int i, j;
    int hi, lo;

    pgm_file = fopen(filename, "wb");
    if (pgm_file == NULL) {
        perror("Cannot open file to write");
        exit(EXIT_FAILURE);
    }

    if (ver)
        fprintf(pgm_file, "P5 ");
    else
    	fprintf(pgm_file, "P2 ");
    fprintf(pgm_file, "%d %d ", data->height, data->width);
    fprintf(pgm_file, "%d ", data->max_gray);

    if(ver)
    {
        if (data->max_gray > 255)
        {
            for (i = 0; i < data->width; ++i)
            {
                for (j = 0; j < data->height; ++j)
                {
                    hi = UP8(data->pixels[i][j]);
                    lo = LO8(data->pixels[i][j]);
                    fputc(hi, pgm_file);
                    fputc(lo, pgm_file);
                }
            }
        }
        else
        {
            for (i = 0; i < data->width; ++i)
                for (j = 0; j < data->height; ++j)
                {
                    lo = LO8(data->pixels[i][j]);
                    fputc(lo, pgm_file);
                }
        }
    }

    else
    {
    	for (i = 0; i < data->width; ++i)
    	{
    	    for (j = 0; j < data->height; ++j)
    	        {
    	    	    fprintf(pgm_file,"%d ",data->pixels[i][j]);
    	        }
    	}

    }
    fclose(pgm_file);
    //deallocate_dynamic_matrix(data->pixels, data->row);
}

/* To print the structure's data */
void printPGM(PGMData picture)
{
    int i;
    int j;
    printf("Printing the PGM image:\n");

	printf("Width = %d\tHeight=%d\tMax gray value = %d\n",picture.width,picture.height,picture.max_gray);

	    for(i=0; i<picture.width; i++)
	    {
	    	for(j=0; j<picture.height; j++)
	    	{
	    		printf("\t%d",picture.pixels[i][j]);
	    	}
	    printf("\n");
	    }
}

double *calculate_harlick_parameters(double **P, int data_max_gray, double *fx)
{

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

    printf("N = %f\tmax gray = %d\n",N,data_max_gray);
    double Px[1024], Py[1024], ux, uy,u, sx, sy, Pxplusy[2048], Pxminusy[2048];
/*
    Px = allocate_dynamic_matrix_a(data_max_gray);
    Py = allocate_dynamic_matrix_a(data_max_gray);
    Pxplusy = allocate_dynamic_matrix_a(2*data_max_gray-1);
    Pxminusy = allocate_dynamic_matrix_a(2*data_max_gray-1);
*/
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

       		u = u+ i*j*P[i][j];

       	}
    }



    sx = 0; sy = 0;
    for(i=0; i<data_max_gray; i++)
    {
      	sx = sx+Px[i]*(i-ux)*(i-ux);
       	sy = sy+Py[i]*(i-uy)*(i-uy);
    }


    sx = sqrt(sx);
    sy = sqrt(sy);


    double HXY1, HXY2, HX=0, HY=0;
    double **Q;
    Q = allocate_dynamic_matrix_d((data_max_gray), (data_max_gray));
    for(i=0; i<data_max_gray; i++)
    {
        if(Px[i]>0)
            HX = HX - Px[i]*log(Px[i]);
            //printf("%f\t%f\n",Px[i],HX);
       	if(Py[i]>0)
    	    HY = HY - Py[i]*log(Py[i]);
    }


    HXY1 = 0; HXY2 = 0;
    for(i=0; i<data_max_gray; i++)
    {
     	for(j=0; j<data_max_gray; j++)
     	{
            if(Px[i]*Py[j]!=0)
       		    HXY1 = HXY1 - P[i][j]*log(Px[i]*Py[j]);
         	if(Px[i]*Py[j]!=0)
         		HXY2 = HXY2 - Px[i]*Py[j]*log(Px[i]*Py[j]);

         	Q[i][j]=0;
         	for(k=0;k<data_max_gray;k++)
         		if(Px[i]*Py[j]!=0)
         		Q[i][j] = Q[i][j]+ P[i][k]*P[j][k]/(Px[i]*Py[j]);
         }
     }

       printf("ux = %g\tuy = %g\tsx = %g\tsy=%g\n",ux,uy,sx,sy);
       printf("HXY1 = %g\tHXY2 = %g\tHX = %g\tHY = %g\n",HXY1,HXY2,HX,HY);


    for(i=0; i<14; i++)
        fx[i]=0;
    for(i=0; i<data_max_gray; i++)
    {
    	for(j=0;j<data_max_gray; j++)
    	{
    		fx[0] = fx[0] + P[i][j]*P[i][j];
    		fx[2] = fx[2] + (i-ux)*(i-uy)*P[i][j]/(sx*sy);
    		fx[3] = fx[3] + (i-u)*(i-u)*P[i][j];
    		fx[4] = fx[4] + P[i][j]/(1+(i-j)*(i-j));
    		if(P[i][j]>0)
    		    fx[8] = fx[8] - P[i][j]*log(P[i][j]);
    		if(Pxminusy[i]>0)
    		    fx[10] = fx[10] - Pxminusy[i]*log(Pxminusy[i]);
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
    	if(Pxplusy[i]>0)
    	    fx[7] = fx[7] - Pxplusy[i]*log(Pxplusy[i]);
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
    fx[11] = (fx[8]-HXY1)/t;
    fx[12] = sqrt(1-exp(-2*(HXY2 - fx[8])));


    return fx;

}


/* To create co-occurance matrix */
void create_cooccurance_matrix(PGMData *data, int delta, int angle)
{
	printf("Come on I don't crash\n");
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
    P = allocate_dynamic_matrix_d((data->max_gray), (data->max_gray));

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
    deallocate_dynamic_matrix(data->pixels, data->width);
    /*
    printf("\nPRINTING COOCCURANCE MATRIX with theta = %d\n",angle);

    for(i=0;i<data->max_gray;i++)
    {
    	for(j=0;j<data->max_gray;j++)
    	{
    		printf("%g\t",P[i][j]);
    	}
    printf("\n");
    }*/

    double f[14];
    calculate_harlick_parameters(P,data->max_gray,f);

    printf("HARLICKs PARAMETERS\n");
    for(i=0; i<13; i++)
    {
    	printf("f%d = %g\t\t",i+1, f[i]);
    	if(i%3==2) printf("\n");
    }

    printf("\n");


}


int main(int argc, char *argv[])
{
    if(argc!=2)
    {
    	perror("File name not specified");
    	exit(1);
    }

    PGMData picture;

    readPGM(argv[1], &picture);
 //   printPGM(picture);
    create_cooccurance_matrix(&picture,1,45);
    return 0;
}
