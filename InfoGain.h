/* To calculate the InfoGain of each parameter stored in SVM data format */
/* Author: RR Iyer
 * Institution: ISI Kolkata
 * Date: 11-Jun-2015
 */


int is_present(double S, double Tf[], int n)
{
	int i, res=-1;
	for(i=0; i<n; i++)
	{
		if(S - Tf[i]< 0.000009)
		{
			res=i;
			break;
		}
	}
	return res;
}

void form_matrix(int clas, double f[], int n, double **Sx, int ii)
{
	Sx[ii][0] = (double) clas;
	int j;
	for(j=1; j<=n; j++)
	{
		Sx[ii][j] = trunc(f[j-1]*10000)/10000;
	}
}

void find_info_gain(double **S, double *x, int n_cl, int n_sam, int n_feat)
{
    int i,j,k, count_t=0;
    double pC[1000], Tf[6000], pTf[6000];
    for(i=0; i<1000; i++) pC[i]=0;

    for(i=0; i<n_sam; i++)
    {
    	int temp = (int) S[i][0];
    	pC[temp-1]++;
    }
    double sum_c=0;
    for(i=0; i<n_cl; i++)
    {
    	sum_c = sum_c + pC[i];
    }

    for(i=0; i<n_cl; i++)
    {
    	pC[i] = pC[i]/sum_c;
    }

    for(k=1; k<=n_feat; k++)
    {
    	for(i=0; i<6000; i++) pTf[i]=0;
    	for(i=0; i<6000; i++) Tf[i]=0;
    	count_t=0;

    	for(i=0; i<n_sam; i++)
    	{
    		int pres = is_present(S[i][k],Tf,count_t);
    		if(pres==-1)
    		{
    			Tf[count_t]=trunc(S[i][k]*10000)/10000;
    			pTf[count_t]++;
    			count_t++;
            }
    		else
    		{
    			pTf[pres]++;
    		}
    	}
    	sum_c=0;
    	for(i=0; i<count_t; i++)
    	{
    		sum_c = sum_c + pTf[i];
    	}
    	for(i=0; i<count_t; i++)
    	{
    		pTf[i] = pTf[i]/sum_c;
    	}
    	double **M = allocate_dynamic_matrix_double(count_t, n_cl);
        for(i=0; i<count_t; i++)
        {
        	for(j=0; j<n_cl; j++)
        	{
        		M[i][j] = 0;
        	}
        }

        for(i=0; i<n_sam; i++)
        {
        	int pos = is_present(S[i][k], Tf, count_t);
		    int temp = (int)S[i][0];
		    M[pos][temp-1]++;
        }
        sum_c=0;
        for(i=0; i<count_t; i++)
        {
        	for(j=0; j<n_cl; j++)
        	{
        		sum_c = sum_c+M[i][j];
        	}
        }
        for(i=0; i<count_t; i++)
        {
        	for(j=0; j<n_cl; j++)
        	{
        		M[i][j] = M[i][j]/sum_c;
        	}
        }
        double term1 = 0;
        for(i=0; i<n_cl; i++)
        {
        	term1 = term1 - pC[i]*log(pC[i]+0.00001);
        }
        double term2 = 0;
        for(i=0; i<count_t; i++)
        {
        	for(j=0; j<n_cl; j++)
        	{
        		term2 = term2 + M[i][j] * log (M[i][j]/(pC[j]*pTf[i] + 0.00001) + 0.00001);
        	}
        }
        x[k-1] = term1 + term2;
        deallocate_dynamic_matrix_double(M,count_t);
    }
}
