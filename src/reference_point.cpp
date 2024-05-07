#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "rand.h"

int nchoosek(int n, int k)
{
	int i;
	double prod;
	prod = 1.0;

	if (n == 0 && k == 0)
		return 1;
	else
	{
		for (i = 1; i <= k; i++)
			prod = prod * (double)((double)(n + 1 - i) / (double)i);

		return (int)(prod + 0.5);
	}
}

void add_ref_points(int p)
{

	int i, j, l, m, e, no, tp;
	double delta, limit, temp, beta;
	delta = (double)(1.0 / (double)p);
	m = nchoosek(nobj + p - 1, p);

	for (i = 0; i < nobj - 1; i++)
	{
		e = nref;
		while (e <= nref + m - 1)
		{
			if (i == 0)
				limit = 0;
			else
			{
				limit = 0.0;
				for (j = 0; j < i; j++)
					limit = limit + ref_pt[e][j];
			}
			for (j = 0; j <= (int)(((1 - limit) / delta) + 0.5); j++)
			{
				beta = delta * (double)j;

				tp = (int)(((1 - limit - beta) / delta) + 0.5);
				no = nchoosek(nobj - i - 2 + tp, tp);
				for (l = e; l < e + no; l++)
				{
					ref_pt[l][i] = beta;
				}

				e = e + no;
			}
		}
	}

	for (i = nref; i < nref + m; i++)
	{
		temp = 0.0;
		for (j = 0; j < nobj - 1; j++)
		{
			temp = temp + ref_pt[i][j];
		}
		ref_pt[i][nobj - 1] = 1 - temp;
		created_around[i] = 0;
	}
	nref = nref + m;
	/*for(i=0;i<m;i++)
	{
		for (j=0;j<nobj;j++)
		{
			ref_pt[i][j]=ref_pt[i][j]*0.25;
			if (j==nobj-1)
				ref_pt[i][j] += 0.75;
		}
	}*/
}

void create_ref_points(int p)
{

	int i, j, l, m, e, no, tp;
	double delta, limit, temp, beta;
	float value;
	FILE *fp;

	if (p == 0)
	{
		fp = fopen("pref_ref.txt", "r");
		for (i = 0; i < nref; i++)
		{
			temp = 0.0;
			for (j = 0; j < nobj; j++)
			{
				fscanf(fp, "%f\t", &value);
				ref_pt[i][j] = (double)value;
				temp += ref_pt[i][j];
			}
			/* Normalizing the reference points */
			for (j = 0; j < nobj; j++)
				ref_pt[i][j] = ref_pt[i][j] / temp;

			fscanf(fp, "\n");
		}
		fclose(fp);
	}
	else
	{
		/* Das and Dennis's approach to create structred reference points */

		delta = (double)(1.0 / (double)p);
		m = nchoosek(nobj + p - 1, p);
		nref = m;

		for (i = 0; i < nobj - 1; i++)
		{
			e = 0;
			while (e <= m - 1)
			{
				if (i == 0)
					limit = 0;
				else
				{
					limit = 0.0;
					for (j = 0; j < i; j++)
						limit = limit + ref_pt[e][j];
				}
				for (j = 0; j <= (int)(((1 - limit) / delta) + 0.5); j++)
				{
					beta = delta * (double)j;

					tp = (int)(((1 - limit - beta) / delta) + 0.5);
					no = nchoosek(nobj - i - 2 + tp, tp);

					for (l = e; l < e + no; l++)
					{
						ref_pt[l][i] = beta;
					}

					e = e + no;
				}
			}
		}
		for (i = 0; i < m; i++)
		{
			temp = 0.0;
			for (j = 0; j < nobj - 1; j++)
			{
				temp = temp + ref_pt[i][j];
			}
			ref_pt[i][nobj - 1] = 1 - temp;
			created_around[i] = 0;
		}

		/* if nobj is high use multilayered reference points */
		if (nobj > 5)
		{
			for (i = 0; i < m; i++)
			{
				for (j = 0; j < nobj; j++)
				{
					ref_pt[i][j] = ref_pt[i][j] * 0.5 + 0.5 / (double)nobj;
				}
			}
			add_ref_points(p + 1);
			/*scaling = 2.0 * scaling;*/
		}
	}
	onref = nref;
}

int generate_ref_points(int p, double **pts)
{

	int i, j, l, m, e, no, tp;
	double delta, limit, temp, beta;

	delta = (double)(1.0 / (double)p);
	m = nchoosek(nobj + p - 1, p);

	for (i = 0; i < nobj - 1; i++)
	{
		e = 0;
		while (e <= m - 1)
		{
			if (i == 0)
				limit = 0;
			else
			{
				limit = 0.0;
				for (j = 0; j < i; j++)
					limit = limit + pts[e][j];
			}
			for (j = 0; j <= (int)(((1 - limit) / delta) + 0.5); j++)
			{
				beta = delta * (double)j;

				tp = (int)(((1 - limit - beta) / delta) + 0.5);
				no = nchoosek(nobj - i - 2 + tp, tp);
				/*fprintf(fp,"beta = %f\t tp = %i\t no = %i\t e = %i\n",beta,tp,no,e);*/
				for (l = e; l < e + no; l++)
				{
					pts[l][i] = beta;
				}

				e = e + no;
			}
		}
	}

	for (i = 0; i < m; i++)
	{
		temp = 0.0;
		for (j = 0; j < nobj - 1; j++)
		{
			temp = temp + pts[i][j];
		}
		pts[i][nobj - 1] = 1 - temp;
	}
	return m;
}
