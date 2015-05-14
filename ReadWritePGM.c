/* Program to read a PGM file */
/* Written by:  R R Iyer
 * Institution: ISI Kolkata
 * Date:        14 May 2015
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
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

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

int main(int argc, char *argv[])
{
    if(argc!=4)
    {
    	perror("File name not specified");
    	exit(1);
    }

    PGMData picture;

    readPGM(argv[1], &picture);
    printPGM(picture);
    normalisePGM(&picture,63);
    printPGM(picture);
    writePGM(argv[2],&picture,1);
    writePGM(argv[3],&picture,0);
    return 0;
}
