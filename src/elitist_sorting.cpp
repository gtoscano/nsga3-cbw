#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <fmt/core.h>
#include "global.h"
#include "rand.h"

int remaining;		   /* Number of points that are to be selected from last front */
int limiting_front;	   /* Last front */
int selected_pop_size; /* size of pop S_t*/
int final_pop_size;
int nholes;						 /* Number of reference points having no pop member associated to them ie holes */
int *rank;						 /* Non-domination rank of each pop member */
int *selected_pop_indices;		 /* Indices of pop members upto last front */
int *final_pop_indices;			 /* Indices of pop members in array "selected_pop_indices" which are surely included in P_{t+1}*/
int *nclstr;					 /* Number of pop members associated with each reference point including last front ie "niche count"*/
int *nclstr_sure;				 /* Number of pop members associated with each reference point that are to be surely selected */
int *idx;						 /* contains the index of reference point associated with each population point */
int **clstr_info;				 /* clstr_info[i][j] contains indices of population members belongining to jth reference point */
int **clstr_info_limiting_front; /* clstr_info_limiting_front[i][j] contains indices of population members belongining to last front and to jth reference point */
int *within_clstr_rank;
int *hole; /* Indices of reference points with 0 niche count */
double *intercept;
double **trans_pop;
double **D; /* Perpendicular distance of each pop member from each ref point */
int incr_exhausted;

double extreme[2];
bool first = false;
int extreme_idx;

/* Routine to determine the Last Front */
void find_limiting_front(population *pop)
{
	int i, maxfront, total, max_feasible_front;
	int *front_count;

	if (curr_gen > 1)
	{
		front_count = (int *)malloc(2 * popsize * sizeof(int));
		maxfront = 0;
		max_feasible_front = 0;
		for (i = 0; i < 2 * popsize; i++)
		{
			front_count[i] = 0;
			rank[i] = pop->ind[i].rank;
			if ((rank[i] > max_feasible_front) && (pop->ind[i].constr_violation >= 0.0))
				max_feasible_front = rank[i];
			if (rank[i] > maxfront)
				maxfront = rank[i];

		}
		for (i = 0; i < 2 * popsize; i++)
			front_count[rank[i] - 1] += 1;
		total = 0;
		for (i = 0; i < maxfront; i++)
		{
			total += front_count[i];
			if (total > popsize)
			{
				limiting_front = i + 1;
				remaining = popsize - total + front_count[i];
				break;
			}
		}
		if (remaining > 0 && max_feasible_front < limiting_front)
		{
			remaining = 0;
			limiting_front += 1;
		}
		/*for(i=0;i<maxfront;i++)
			printf("front =%i count = %i\n",i,front_count[i]);*/

		free(front_count);
	}
	else
	{

		maxfront = 0;
		for (i = 0; i < popsize; i++)
		{

			rank[i] = pop->ind[i].rank;
			if (rank[i] > maxfront)
				maxfront = rank[i];
		}
		limiting_front = maxfront + 1;
		remaining = 0;
	}

	/*printf("lim : %i\n",limiting_front);
	printf("rem : %i\n",remaining);
	printf("maxf : %i\n",maxfront);*/
}

/* Routine to form population S_t and P_{t+1} for an unconstrained problem*/
void define_selected_pop()
{
	int i;
	final_pop_size = 0;
	selected_pop_size = 0;
	if (remaining > 0)
	{
		for (i = 0; i < 2 * popsize; i++)
		{
			if (rank[i] < limiting_front + 1)
			{
				selected_pop_indices[selected_pop_size] = i;
				selected_pop_size += 1;
				if (rank[i] < limiting_front)
				{
					final_pop_indices[final_pop_size] = selected_pop_size - 1;
					;
					final_pop_size += 1;
				}

			}
		}
	}
	else
	{
		if (curr_gen > 1)
		{
			for (i = 0; i < 2 * popsize; i++)
			{
				if (rank[i] < limiting_front)
				{
					selected_pop_indices[selected_pop_size] = i;
					selected_pop_size += 1;
					final_pop_indices[final_pop_size] = i;
					final_pop_size += 1;
				}
			}
		}
		else
		{
			for (i = 0; i < popsize; i++)
				selected_pop_indices[i] = i;
			selected_pop_size = popsize;
		}
	}
	/*printf("spop : %i\n",selected_pop_size);*/
}

/* Routine to form population S_t and P_{t+1} for an constrained problem*/
void define_feasible_selected_pop(population *pop)
{
	int i;
	final_pop_size = 0;
	selected_pop_size = 0;
	/* Here remaining > 0 implies that no. of feasible solutions is greater than popsize */
	if (remaining > 0)
	{
		for (i = 0; i < 2 * popsize; i++)
		{
			if (rank[i] < limiting_front + 1)
			{
				selected_pop_indices[selected_pop_size] = i;
				selected_pop_size += 1;

				if (rank[i] < limiting_front && final_pop_size < popsize)
				{
					final_pop_indices[final_pop_size] = selected_pop_size - 1;
					final_pop_size += 1;
				}
			}
		}
	}
	else
	{
		if (curr_gen > 1)
		{
			selected_pop_size = 0;
			for (i = 0; i < 2 * popsize; i++)
			{
				if (rank[i] < limiting_front)
				{
					if (pop->ind[i].constr_violation >= 0.0)
					{
						selected_pop_indices[selected_pop_size] = i;
						selected_pop_size += 1;
					}
					if (final_pop_size < popsize)
					{
						final_pop_indices[final_pop_size] = i;
						final_pop_size += 1;
					}
				}
			}
		}
		else
		{
			for (i = 0; i < popsize; i++)
			{
				if (pop->ind[i].constr_violation >= 0.0)
				{
					selected_pop_indices[selected_pop_size] = i;
					selected_pop_size += 1;
				}
				final_pop_indices[final_pop_size] = i;
				final_pop_size += 1;
			}
		}
	}
	/*	printf("spop : %i\n",selected_pop_size);*/
}

void translate_pop(population *pop)
{
	int i, j;
	trans_pop = (double **)malloc(selected_pop_size * sizeof(double *));
	for (i = 0; i < selected_pop_size; i++)
		trans_pop[i] = (double *)malloc(nobj * sizeof(double));
	for (i = 0; i < selected_pop_size; i++)
	{
		for (j = 0; j < nobj; j++)
			trans_pop[i][j] = pop->ind[selected_pop_indices[i]].obj[j] - ideal_point[j];
	}
	/*for(j=0;j<nobj;j++)
		printf("ideal pt : %f\t",ideal_point[j]);*/
}

void find_intercepts()
{
	int i, k, l, minind, flag, flag2;
	double *A;
	double asf, minasf, temp;
	double *nadir;
	double *fmax;

	A = (double *)malloc(nobj * nobj * (sizeof(double)));
	nadir = (double *)malloc(nobj * sizeof(double));
	fmax = (double *)malloc(nobj * sizeof(double));
	/* Determining Nadir Point */

	for (i = 0; i < nobj; i++)
	{
		nadir[i] = -1 * INF;
		fmax[i] = -1 * INF;
		for (k = 0; k < selected_pop_size; k++)
		{
			if (rank[selected_pop_indices[k]] == 1)
			{
				if (nadir[i] < trans_pop[k][i])
					nadir[i] = trans_pop[k][i];
			}
			if (fmax[i] < trans_pop[k][i])
				fmax[i] = trans_pop[k][i];
		}
		/*printf("\nnadir  = %f\n",nadir[i]);*/
	}

	/* Finding plane points to form Hyperplane */
	for (i = 0; i < nobj; i++)
	{
		intercept[i] = 1.0;
		asf = INF;
		/* Getting asf value of previous generation plane points */
		if (curr_gen > 1)
		{
			asf = -1.0;
			for (k = 0; k < nobj; k++)
			{
				plane_point[i][k] = plane_point[i][k] - ideal_point[k];
				temp = plane_point[i][k];
				if (temp < 1.0e-3)
					temp = 0;
				if (k != i)
					temp = temp * 1.0e6;
				if (temp > asf)
					asf = temp;
			}
			minasf = asf;
			minind = -1;
		}
		else
		{
			minind = 0;
			minasf = INF;
		}
		/* Identifying the minimum asf points */
		for (k = 0; k < selected_pop_size; k++)
		{
			if (rank[selected_pop_indices[k]] == 1)
			{
				asf = -1.0;
				for (l = 0; l < nobj; l++)
				{
					temp = trans_pop[k][l];
					if (temp < 1.0e-3)
						temp = 0.0;
					if (l != i)
						temp = temp * 1.0e6;
					if (temp > asf)
						asf = temp;
				}
			}
			/*printf("asf =  %f\t minasf = %f\n",asf,minasf);*/
			if (asf < minasf)
			{
				minasf = asf;
				minind = k;
			}
		}
		if (minind >= 0)
		{
			/*printf("minasf =  %f\n",minasf);*/
			for (k = 0; k < nobj; k++)
			{
				plane_point[i][k] = trans_pop[minind][k];
			}
		}
	}
	for (i = 0; i < nobj; i++)
	{
		for (k = 0; k < nobj; k++)
			A[i * nobj + k] = plane_point[i][k];
	}
	/* extending plane to axes and determining intercepts */

	flag = Gaussian_Elimination(A, nobj, intercept);
	if (flag == -1)
	{
		for (i = 0; i < nobj; i++)
		{
			intercept[i] = nadir[i];
		}
	}
	else
	{
		flag2 = 0;
		for (i = 0; i < nobj; i++)
		{
			if (intercept[i] <= 1e-6)
			{
				flag2 = 1;
				break;
			}
			intercept[i] = 1.0 / intercept[i];
			if (intercept[i] != intercept[i])
			{
				flag2 = 1;
				break;
			}
		}
		if (flag2 == 1)
		{
			for (i = 0; i < nobj; i++)
				intercept[i] = nadir[i];
		}
	}

	for (i = 0; i < nobj; i++)
	{
		if (intercept[i] <= 1e-6)
			intercept[i] = fmax[i];

		/*printf("intercept =  %f ideal pt = %f  nadir = %f\n",intercept[i],ideal_point[i],nadir[i]);*/

		for (k = 0; k < nobj; k++)
		{
			plane_point[i][k] = plane_point[i][k] + ideal_point[k];
			/*printf("plane pt[%i][%i] =  %f\n",i,k,plane_point[i][k]);*/
		}
	}

	free(nadir);
	free(fmax);
	free(A);
}

void associate()
{
	int i, j, k;
	double temp, d1, d2, lambda;
	/* Perpendicular distance calculation */
	for (i = 0; i < nref; i++)
	{
		for (j = 0; j < selected_pop_size; j++)
		{

			d1 = 0.0;
			lambda = 0.0;
			for (k = 0; k < nobj; k++)
			{
				d1 += trans_pop[j][k] * ref_pt[i][k] / intercept[k];
				lambda += ref_pt[i][k] * ref_pt[i][k];
			}
			lambda = sqrt(lambda);
			d1 = d1 / lambda;
			d2 = 0.0;
			for (k = 0; k < nobj; k++)
			{
				d2 += pow((trans_pop[j][k] / intercept[k] - d1 * ref_pt[i][k] / lambda), 2.0);
			}
			D[j][i] = sqrt(d2);
		}
		nclstr[i] = 0;
		nclstr_sure[i] = 0;
		clstr_info_limiting_front[0][i] = 0;
	}

	/* Associating each population member with the closest reference point */
	for (i = 0; i < selected_pop_size; i++)
	{
		temp = INF;
		for (j = 0; j < nref; j++)
		{
			if (temp > D[i][j])
			{
				temp = D[i][j];
				idx[i] = j;
			}
		}

		clstr_info[nclstr[idx[i]]][idx[i]] = i;
		nclstr[idx[i]] += 1;

		if (rank[selected_pop_indices[i]] < limiting_front)
			nclstr_sure[idx[i]] += 1;
	}
	/*for(i=0;i<nref;i++)
		printf("nclstr[%i] =  %i\n",i,nclstr[i]);
	for(i=0;i<selected_pop_size;i++)
		printf("id[%i] =  %i\n",i,idx[i]);*/
}

/* Routine to identify the closest pop member in a cluster and giving it a rank one while assigning ranks randomly to the remaining members */
void calculate_within_cluster_ranks([[maybe_unused]] population *pop)
{
	int i, j, best=0, tempi, flag;
	double temp;

	if (remaining > 0)
	{
		for (i = 0; i < nref; i++)
		{
			if (nclstr[i] > 0)
			{
				if (nclstr[i] > 1)
				{
					temp = INF;
					flag = 0;
					for (j = 0; j < nclstr[i]; j++)
					{
						if (rank[selected_pop_indices[clstr_info[j][i]]] == limiting_front && nclstr_sure[i] == 0)
						{
							if (temp > D[clstr_info[j][i]][i])
							{
								temp = D[clstr_info[j][i]][i];
								best = j;
							}
							flag = 1;
						}
					}
					if (flag == 1)
					{

						tempi = clstr_info[0][i];
						clstr_info[0][i] = clstr_info[best][i];
						clstr_info[best][i] = tempi;
					}
				}

				for (j = 0; j < nclstr[i]; j++)
				{
					if (rank[selected_pop_indices[clstr_info[j][i]]] == limiting_front)
					{

						clstr_info_limiting_front[clstr_info_limiting_front[0][i] + 1][i] = clstr_info[j][i];
						clstr_info_limiting_front[0][i] += 1;
					}
				}
			}
		}
	}
}

void remove_points()
{
	int i, j, k, t, n, flag;
	int *remove;
	remove = (int *)malloc(nref * sizeof(int));
	flag = 0;
	for (i = 0; i < nref; i++)
	{
		remove[i] = 0;
		if (nclstr[i] == 0)
		{
			/* remove reference points except the original ones */
			if (i >= onref)
				flag = 1;
			remove[i] = 1;
		}
	}
	if (flag == 1)
	{
		k = 0;
		i = onref;
		n = onref;
		while (n < nref)
		{
			t = 0;
			while (remove[i + k] == 1)
			{
				k += 1;
				t += 1;
			}
			n += t + 1;

			if (k > 0)
			{
				for (j = 0; j < nobj; j++)
					ref_pt[i][j] = ref_pt[i + k][j];
			}
			i += 1;
		}
		nref -= k;
		/*printf("reference points reduced from %i to %i\n",nref+k,nref);
		getchar();*/
	}
	for (i = 0; i < nref; i++)
		created_around[i] = 0;
	free(remove);
}

void define_final_pop()
{
	int i, j, nh, idef=0, nmax, nmin;
	int *select_count, *nothole;
	int temp;
	select_count = (int *)malloc(nref * sizeof(int)); /* contains the number of points selected from last front belonging to jth reference point */
	nothole = (int *)malloc(nref * sizeof(int));	  /* indices of reference points having atleast one pop member */

	nh = 0;
	for (i = 0; i < nref; i++)
	{
		select_count[i] = 0;
		if (nclstr[i] > 0)
		{
			nothole[nh] = i;
			nh += 1;
		}
	}
	if (remaining > 0)
	{
		i = 0;
		while (i < remaining)
		{
			temp = 100000;
			/* Identify minimum niche count clstr */
			for (j = 0; j < nh; j++)
			{
				if (temp > nclstr_sure[nothole[j]])
				{
					temp = nclstr_sure[nothole[j]];
					idef = nothole[j];
					if (temp == 0)
						break;
				}
			}

			if ((select_count[idef] + 1) > clstr_info_limiting_front[0][idef])
				nclstr_sure[idef] = 1000000;
			else
			{
				final_pop_indices[final_pop_size + i] = clstr_info_limiting_front[select_count[idef] + 1][idef];
				nclstr_sure[idef] += 1;
				select_count[idef] += 1;
				i += 1;
			}
		}
		for (i = 0; i < nref; i++)
		{
			nclstr[i] = 0;
		}
		for (i = 0; i < popsize; i++)
			nclstr[idx[final_pop_indices[i]]] += 1;
		nholes = 0;
		incr_exhausted = 1;
		for (i = 0; i < nref; i++)
		{
			if (nclstr[i] == 0)
				nholes += 1;
			if (nclstr[i] >= 1)
			{
				/*printf("\n%d",created_around[i]);*/
				if (created_around[i] < (nobj))
					incr_exhausted = 0;
			}
		}
		if (start_incr == 0)
		{
			nmax = nref - nholes;
			nmin = nref - nholes;
			for (i = 9; i > 0; i--)
			{
				if (nmax < active_ref_pts[i - 1])
					nmax = active_ref_pts[i - 1];
				if (nmin > active_ref_pts[i - 1])
					nmin = active_ref_pts[i - 1];
				active_ref_pts[i] = active_ref_pts[i - 1];
			}

			active_ref_pts[0] = nref - nholes;
			if ((nmax - nmin) < ((int)(0.01 * (double)onref) + 1))
			{
				printf("\nnmax = %d nmin = %d", nmax, nmin);
				start_incr = 1;
			}
		}
	}

	free(select_count);
	free(nothole);
}
void create_points_around(int clstrid)
{
	int i, j, k, flag, tnref, flag2, level, n, exist;
	double **temp, dist;
	flag = 0;
	temp = (double **)malloc(nref * sizeof(double *));
	for (i = 0; i < nref; i++)
		temp[i] = (double *)malloc(nobj * sizeof(double));
	tnref = nref;
	level = 1;
	while (level < 2 && flag == 0)
	{
		n = generate_ref_points(level, temp);

		for (i = 0; i < n; i++)
		{
			exist = 0;
			flag2 = 1;
			/* check if the introduced point is out of boundaries */
			for (k = 0; k < nobj; k++)
			{
				temp[i][k] = temp[i][k] * (1.0 / scaling) + ref_pt[clstrid][k] - (1.0 / ((double)nobj * scaling));
				if (temp[i][k] < 0)
					flag2 = 0;
			}
			if (flag2 == 1)
			{
				for (j = 0; j < tnref; j++)
				{
					/* Check if the introduced reference point already exists */
					dist = 0.0;
					for (k = 0; k < nobj; k++)
						dist += (ref_pt[j][k] - temp[i][k]) * (ref_pt[j][k] - temp[i][k]);
					dist = sqrt(dist);
					if (dist < 1e-3)
					{
						exist = 1;
						break;
					}
				}

				if (exist == 0)
				{
					/* Add the reference point */
					created_around[nref] = 0;
					for (k = 0; k < nobj; k++)
						ref_pt[nref][k] = temp[i][k];
					nref += 1;
					/*printf("reference points increased to %i\n",nref);
					getchar();*/
				}
			}
		}
		if (nref == tnref)
		{
			flag = 0;
			level += 1;
		}
		else
			flag = 1;
	}
	if (nref > tnref)
		created_around[clstrid] = 1;
	for (i = 0; i < tnref; i++)
		free(temp[i]);
	free(temp);
}
void create_points_around2(int clstrid)
{
	int i, j, k, flag, tnref, flag2, type, n, exist;
	double **temp, dist;

	flag = 0;
	temp = (double **)malloc(nref * sizeof(double *));
	for (i = 0; i < nref; i++)
		temp[i] = (double *)malloc(nobj * sizeof(double));
	tnref = nref;
	type = created_around[clstrid];
	while (created_around[clstrid] < nobj && flag == 0)
	{
		n = generate_ref_points(1, temp);
		type = created_around[clstrid];

		for (i = 0; i < n; i++)
		{
			if (i != type)
			{
				exist = 0;
				flag2 = 1;
				/* check if the introduced point is out of boundaries */
				for (k = 0; k < nobj; k++)
				{
					temp[i][k] = temp[i][k] / (2.0 * scaling) + ref_pt[clstrid][k] - (temp[type][k] / (2.0 * scaling));
					if (temp[i][k] < 0)
					{
						/*printf("\there1");*/
						flag2 = 0;
					}
				}
				if (flag2 == 1)
				{
					for (j = 0; j < tnref; j++)
					{
						/* Check if the introduced reference point already exists */
						dist = 0.0;
						for (k = 0; k < nobj; k++)
							dist += (ref_pt[j][k] - temp[i][k]) * (ref_pt[j][k] - temp[i][k]);
						dist = sqrt(dist);
						if (dist < 1e-3)
						{
							exist = 1;
							/*printf("\nexists");*/
							break;
						}
					}

					if (exist == 0)
					{
						/* Add the reference point */
						created_around[nref] = 0;
						for (k = 0; k < nobj; k++)
							ref_pt[nref][k] = temp[i][k];
						nref += 1;
						/*printf("reference points increased to %i\n",nref);
					getchar();*/
					}
				}
			}
		}
		created_around[clstrid] += 1;
		if (nref == tnref)
		{
			flag = 0;
		}
		else
			flag = 1;
	}
	/*if (flag==0 && created_around[clstrid]==nobj)
	{
		printf("\n");
		for (i=0;i<nobj;i++)
			printf("%f\t",ref_pt[clstrid][i]);
	}*/
	for (i = 0; i < tnref; i++)
		free(temp[i]);
	free(temp);
}

void check_for_increment()
{
	int i, tnref, approach;
	if (adaptive_increment == 1)
		approach = 0;
	else
		approach = 1;
	tnref = nref;
	for (i = 0; i < tnref; i++)
	{
		if (approach == 0)
		{
			if (nclstr[i] > 1 && created_around[i] == 0)
			{
				if (steps == 0)
				{
					if (i >= nobj)
						create_points_around(i);
				}
				else
					create_points_around(i);
			}
		}
		else
		{
			if (nclstr[i] > 1 && created_around[i] < nobj)
			{
				if (steps == 0)
				{
					if (i >= nobj)
						create_points_around2(i);
				}
				else
					create_points_around2(i);
			}
		}
	}
	/*if (nref>tnref)
	{
		printf("reference points increased from %i to %i\n",tnref,nref);
		getchar();
	}*/
}

void elitist_sorting(population *mixed_pop, population *new_pop)
{
	int i, temp_nref;
	double param;

	rank = (int *)malloc(2 * popsize * sizeof(int));
	final_pop_indices = (int *)malloc(popsize * sizeof(int));
	selected_pop_indices = (int *)malloc(2 * popsize * sizeof(int));

  static double extreme_global[2];
  double extreme_local[2];
  static bool label_global = false;
  bool label_local = false;
  static int extreme_global_idx;
  int extreme_local_idx;
  static double extreme_global_constr_violation = std::numeric_limits<double>::lowest();
  double extreme_local_constr_violation = std::numeric_limits<double>::lowest();
  auto which_popsize = (curr_gen <= 1)? popsize: popsize*2;

	find_limiting_front(mixed_pop);
	if (ncon > 0)
		define_feasible_selected_pop(mixed_pop);
	else
		define_selected_pop();

  for(int _idx=0; _idx < which_popsize; ++_idx) {
    if (label_global == false
        ||
        (extreme_global_constr_violation < 0 &&
         extreme_global_constr_violation < mixed_pop->ind[_idx].constr_violation)

        ||
        (mixed_pop->ind[_idx].obj[1] < extreme_global[1]
         && mixed_pop->ind[_idx].constr_violation == 0
        )
        )
    {
      label_global = true;
      extreme_global[0] = mixed_pop->ind[_idx].obj[0];
      extreme_global[1] = mixed_pop->ind[_idx].obj[1];
      extreme_global_idx = _idx;
      extreme_global_constr_violation = mixed_pop->ind[_idx].constr_violation;
    }
    if (label_local == false
        ||
        (extreme_local_constr_violation < 0 &&
         extreme_local_constr_violation < mixed_pop->ind[_idx].constr_violation)
        ||
        (mixed_pop->ind[_idx].obj[1] < extreme_local[1]
         && mixed_pop->ind[_idx].constr_violation == 0
        )
        )
    {
      label_local = true;
      extreme_local[0] = mixed_pop->ind[_idx].obj[0];
      extreme_local[1] = mixed_pop->ind[_idx].obj[1];
      extreme_local_idx = _idx;
      extreme_local_constr_violation = mixed_pop->ind[_idx].constr_violation;
    }
  }

	if (selected_pop_size > nobj)
	{
		intercept = (double *)malloc(nobj * sizeof(double));
		hole = (int *)malloc(nref * sizeof(int));
		D = (double **)malloc(selected_pop_size * sizeof(double *));
		clstr_info = (int **)malloc((selected_pop_size + 1) * sizeof(int *));
		within_clstr_rank = (int *)malloc(selected_pop_size * sizeof(int));
		clstr_info_limiting_front = (int **)malloc((selected_pop_size + 1) * sizeof(int *));
		for (i = 0; i < selected_pop_size + 1; i++)
		{
			if (i < selected_pop_size)
				D[i] = (double *)malloc(nref * sizeof(double));
			clstr_info[i] = (int *)malloc(nref * sizeof(int));
			clstr_info_limiting_front[i] = (int *)malloc(nref * sizeof(int));
		}
		nclstr = (int *)malloc(nref * sizeof(int));
		nclstr_sure = (int *)malloc(nref * sizeof(int));
		idx = (int *)malloc(selected_pop_size * sizeof(int));

		/* Normalization */
		translate_pop(mixed_pop);
		find_intercepts();

		/* Association */
    	associate();

		if (remaining > 0)
			calculate_within_cluster_ranks(mixed_pop);

		if (curr_gen > 1)
		{
			if (remaining > 0)
			{
				/* Niching */
				define_final_pop();

				/* Adaptive reference points */
				printf("No. of active reference points = %d\n", (nref - nholes));

				if (adaptive_increment >= 1 && start_incr == 1)
				{
					temp_nref = nref;
					param = (double)popsize / (double)(nref - nholes);
					if (curr_gen <= ngen)
						printf("\nparam = %f\n", param);
					if (nref < 10 * onref)
						check_for_increment();
					else
						remove_points();
					if (temp_nref < nref)
						printf("\nRef points increased from %d to %d", temp_nref, nref);
					if ((temp_nref == nref) && (param == 1.0) && (nholes > 0))
						remove_points();
					if ((temp_nref == nref) && param > 1.0 && incr_exhausted == 1)
					{
						scaling = 2.0 * scaling;
						remove_points();
						printf("Scaling Doubled >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
					}
					if (temp_nref > nref)
						printf("\nRef points reduced from %d to %d", temp_nref, nref);
				}
			}
		}

		/*fp=fopen("test.txt","w");
		fp1=fopen("preference.txt","r");
		fp2=fopen("crowding.txt","r");*/
	}
	if (curr_gen > 1)
	{

		if (remaining > 0)
		{
                          //selected_pop_indices [final_pop_indices[popsize-1]] = extreme_local_idx;
                          for (i = 0; i < popsize; i++) {
                          copy_ind(
                              &mixed_pop
                                   ->ind[selected_pop_indices[final_pop_indices[i]]],
                              &new_pop->ind[i]);
                        }
/*
                  copy_ind(
                      &mixed_pop ->ind[extreme_local_idx],
                      &new_pop->ind[popsize-1]);
                      */


		}
		else
		{
			for (i = 0; i < popsize; i++)
			{
				/*printf("fp[%i] = %i\n",i,final_pop_indices[i]);*/
				copy_ind(&mixed_pop->ind[final_pop_indices[i]], &new_pop->ind[i]);
			}
		}
	}
    for (i = 0; i < popsize; i++) {
        auto selected_i =  selected_pop_indices[final_pop_indices[i]];
        auto dst_uuid_idx = fmt::format("{}_{}", curr_gen, i);
        if (selected_i < popsize) {
          // Parent population.
          auto src_uuid_idx = fmt::format("{}_{}", curr_gen-1, selected_i);
          uuid_surviving_pop_registry[dst_uuid_idx]  = uuid_surviving_pop_registry[src_uuid_idx];
        }
        else {
          auto src_uuid_idx = fmt::format("{}_{}", curr_gen, selected_i-popsize);
           uuid_surviving_pop_registry[dst_uuid_idx]  = uuid_registry[src_uuid_idx];
        }
    }

	free(rank);
	free(selected_pop_indices);
	free(final_pop_indices);
	if (selected_pop_size > nobj)
	{
		free(intercept);
		free(idx);
		free(nclstr);
		free(nclstr_sure);
		free(within_clstr_rank);
		free(hole);
		for (i = 0; i < selected_pop_size + 1; i++)
		{
			if (i < selected_pop_size)
			{
				free(D[i]);
				free(trans_pop[i]);
			}
			free(clstr_info[i]);
			free(clstr_info_limiting_front[i]);
		}
		free(D);
		free(clstr_info);
		free(clstr_info_limiting_front);
		free(trans_pop);
	}
}
