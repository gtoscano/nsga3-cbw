#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#include "global.h"
#include "rand.h"


void find_ideal_point(population *pop)
{

	int i, j;
	double temp;

	for (i = 0; i < nobj; i++)
	{
		ideal_point[i] = INF;
		for (j = 0; j < popsize; j++)
		{
			if (pop->ind[j].constr_violation >= 0.0)
			{
				temp = pop->ind[j].obj[i];
				if (temp < ideal_point[i])
					ideal_point[i] = temp;
			}
		}
	}
}

void update_ideal_point(population *pop)
{
	int i, j;
	double temp;
	for (i = 0; i < nobj; i++)
	{
		for (j = 0; j < popsize; j++)
		{
			if (pop->ind[j].constr_violation >= 0.0)
			{
				temp = pop->ind[j].obj[i];
				if (temp < ideal_point[i])
					ideal_point[i] = temp;
			}
		}
	}
}
