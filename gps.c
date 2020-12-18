#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define N	10
#define M	3
#define NON_ZERO	30	
#define TMAX 100
#define EPS  (1.0e-6)
#define EPS_D (2.0e-6)
#define OMEGA  1.5
#define LEVEL_END_FLAG = -1

// CRS方式
typedef double crsdata[NON_ZERO];
typedef int crscol[NON_ZERO];
typedef int crsrow[N+1];
typedef int level_index_structure[N+1];

typedef double matrix[N][N];    // matrix
typedef double vector[N];       // vector
typedef int ivector[N];

int print_int_vector(ivector x)
{
    int i;
    for(i = 0; i < 10; i++){
        printf("vec[%d] = %d\n", i, x[i]);
    }
    return 0;
}

void minimal_digree_all(crsrow row, int *md, int *mdi){
	int i, rn;
	for (int i=0; i < N; i++){
		rn = row[i+1] - row[i];
		if (*md > rn){
			
			*md = rn;
			*mdi = i;
		}
	}
	return;
}

int minimal_digree(crsrow row, ivector level, const int si, const int ei){
	int i, rn;
	int md = N; 
	for (int i=0; i < N; i++){
		rn = row[i+1] - row[i];
		if (md > rn)
			md = rn;
	}
	return md;
}

int level_width(level_index_structure level_index, int level_depth){
	int i, width;
	int max_width = 0;
	for (i=1; i<level_depth; i++){
		width = level_index[i+1] - level_index[i];
		if (max_width < width)
			max_width = width;
	}
	return max_width;
}

int generate_level_structure(ivector level, level_index_structure level_index, ivector x_level, crscol col, crsrow row, int v){
	// level: index sorted by level increasing order, level_index: the index where the level is changed
	// x_level: the level of each index
	// Return: Level Number
	memset(x_level, 0, N * sizeof(int));
	int i, j, s, t, current_level, leveled_x;
	int loop_end_flag;
	level[0] = v;
	x_level[v] = 1;
	current_level = 2;
	level_index[0] = 1;
	level_index[1] = 1;
	leveled_x = 1;
	
	while(leveled_x < N){
		for (i = level_index[current_level - 2]; i < level_index[current_level - 1]; i++){
			s = level[i];
			for (j = row[s]; j < row[s+1]; j++){
				t = col[j];
				if (x_level[t] == 0){
					x_level[t] = current_level;
					level[leveled_x] = t;
					leveled_x++;
					if (leveled_x == N){
						level_index[current_level] = N;
						return current_level;
					}
				}
			}
			level_index[current_level] = leveled_x;
			current_level++;
		}
	}
}

void degree_increasing_sort(ivector i_arr, const int len, crsrow row){
	//TODO i_arr のindexをdegree increasing の順にソートする。
}

void pseudo_diameter(int *level_index_v, int *level_index_u, 
		int *x_level_v, int *x_level_u, crsrow row, crscol col, int *u, int *v, int *level_depth){
	int i, j, new_level_depth, u_width, min_u_width, u_index, u_depth; 
	level_index_structure level_index_0, level_index_1;
	ivector level_0, level_1, x_level_0, x_level_1;
	int *po_0, *po_1;
	po_0 = &level_0[0];
	po_1 = &level_1[0];
	min_u_width = NON_ZERO;
	minimal_digree_all(row, level_depth, v);
	new_level_depth = generate_level_structure(level_1, level_index_1, x_level_1, col, row, *v);
	*level_depth = new_level_depth;
	for (i = 0; i < TMAX; i++){
		if (i % 2 == 0){
			degree_increasing_sort(po_1 + level_index_1[new_level_depth - 1], N - level_index_1[new_level_depth - 1], row);
			for (j = level_index_1[new_level_depth - 1]; j < N; j++){
				new_level_depth = generate_level_structure(level_0, level_index_0, x_level_0, col, row, level_1[j]);
				if (*level_depth < new_level_depth){
					*v = level_1[j];
					break;
				}
				u_width = level_width(level_index_0, new_level_depth);
				if (min_u_width > u_width){
					min_u_width = u_width;
					u_index = j;
				}
				if (j == N-1){
					*u = level_1[u_index];
					u_depth = generate_level_structure(level_0, level_index_u, x_level_u, col, row, *u);
					level_index_v = level_index_1;
					return;
				}
			}
		}else{
			degree_increasing_sort(po_0 + level_index_0[new_level_depth - 1], N - level_index_0[new_level_depth - 1], row);
			for (j = level_index_0[new_level_depth - 1]; j < N; j++){
				new_level_depth = generate_level_structure(level_1, level_index_1, x_level_1, col, row, level_0[j]);
				if (*level_depth < new_level_depth){
					*v = level_0[j];
					break;
				}
				u_width = level_width(level_index_1, new_level_depth);
				if (min_u_width > u_width){
					min_u_width = u_width;
					u_index = j;
				}
				if (j == N-1){
					*u = level_0[u_index];
					u_depth = generate_level_structure(level_1, level_index_u, x_level_u, col, row, *u);
					level_index_v = level_index_0;
					return;
				}
			}
		}
	}
}

void minimizing_level_width(crsrow row, crscol col, ivector x_level, level_index_structure level_index, int *u, int *v, int *level_depth, int *chosed_element){
	int i, j, k, l, m, cen, c,  max_v_width, max_u_width;
	cen = 0;
	k = 0;
	l = 0;
	level_index_structure level_index_v, level_index_u;
	ivector x_level_v, x_level_u;
	ivector component_index = {0};
	ivector component_elements ={0};
	ivector component_sep_i = {0};
	ivector sorted_m = {0};
	pseudo_diameter(level_index_v, level_index_u, x_level_v, x_level_u, row, col, u, v, level_depth);
	for (i = 0; i < *level_depth; i++){
		level_index_v[i] = level_index_v[i + 1] - level_index_v[i];
		level_index_u[i] = level_index_u[i + 1] - level_index_u[i];
	}
	for (i = 0; i < N; i++){
		x_level_u[i] = *level_depth + 1 - x_level_u[i];
		if (x_level_v[i] == x_level_u[i]){
			x_level[i] = x_level_u[i];
			component_index[i] = -1;
			level_index[x_level[i] - 1]++;
		}else{
			cen++;
		}
	}
	m = 1;
	i = 0;
	for (i = 0;i < N;i++){
		if (component_index[i] == 0){
			component_elements[k] == i;
			l = k;
			k++;
			while(l < k){
				for (j = row[component_elements[l]]; j < row[component_elements[l] + 1]; j++){
					c = col[j];
					if (component_index[c] == 0){
						component_elements[k] == c;
						component_index[c] == m;
						k++;
					}
				}
				l++;
			}
			component_sep_i[m] = l;
			m++;
			cen -= l;
			if (cen == 0)
				break;
		}
	}
	// TODO component_size argsort >>>> sorted_m

	for(i = 0; i < m; i++){
		max_v_width = 0;
		max_u_width = 0;
		for(j = 0; j > *level_depth; j++){
			k = level_index_v[j] - level_index[j];	
			if (k > max_v_width)
				max_v_width = k;
			k = level_index_u[j] - level_index[j];	
			if (k > max_u_width)
				max_u_width = k;
		}
		if (max_v_width < max_u_width){
			if (i == 0)
				*chosed_element = 0;
			for(j = component_sep_i[sorted_m[i]]; j < component_sep_i[sorted_m[i] + 1]; j++){
				x_level[component_elements[j]] = x_level_v[component_elements[j]];
				level_index[x_level[component_elements[j]]]++;
			}
		}else if(max_v_width > max_u_width){
			if (i == 0)
				*chosed_element = 1;
			for(j = component_sep_i[sorted_m[i]]; j < component_sep_i[sorted_m[i] + 1]; j++){
				x_level[component_elements[j]] = x_level_u[component_elements[j]];
				level_index[x_level[component_elements[j]]]++;
			}
		}else{
			max_v_width = 0;
			max_u_width = 0;
			for(j = 0; j > *level_depth; j++){
				k = level_index_v[j];	
				if (k > max_v_width)
					max_v_width = k;
				k = level_index_u[j];	
				if (k > max_u_width)
					max_u_width = k;
			}
			if (max_v_width <= max_u_width){
				if (i == 0)
					*chosed_element = 0;

				for(j = component_sep_i[sorted_m[i]]; j < component_sep_i[sorted_m[i] + 1]; j++){
					x_level[component_elements[j]] = x_level_v[component_elements[j]];
					level_index[x_level[component_elements[j]]]++;
				}
			}else{
				if (i == 0)
					*chosed_element = 1;
				for(j = component_sep_i[sorted_m[i]]; j < component_sep_i[sorted_m[i] + 1]; j++){
					x_level[component_elements[j]] = x_level_u[component_elements[j]];
					level_index[x_level[component_elements[j]]]++;
				}
			}
		}
	}
	
	return;
}

void numbering_in_level(){}



void gps(crsrow row, crscol col, ivector numbered){
	int i, j, j0, j0_init, j1, k, l, l0, l1, d, cl, lvl,  u,v,w,level_depth, tmp, choosed_element, uv_switched_flag;
	level_index_structure level_width = {0};
	j = 0;
	l1 = 0;
	uv_switched_flag = 0;
	ivector x_level;
	ivector adj0 = {0};
	ivector adj1 = {0};
	minimizing_level_width(row, col, x_level, level_width,  &u, &v, &level_depth, &choosed_element);
	if (row[u + 1] - row[u] < row[v+1] - row[v]){
		tmp = u;
		u = v;
		v = tmp;
		for (i = 0; i < N; i++){
			x_level[i] = level_depth + 1 - x_level[i];
		}
		uv_switched_flag = 1;
	}
	numbered[0] = v;
	for (lvl = 0; lvl < level_depth; lvl++){
		l0 = l1;
		j0 = j1;
		j0_init = j0;
		j1 = l0;
		l1 = l0 + level_width[lvl];
		for (k = l0; k < l1; k++){
			w = numbered[k];
			for (i = row[w]; row[w+1]; i++){
				if (x_level[col[i]] == lvl + 1){
					for (l = l0; l < j0; l++){
						if (col[i] == numbered[l])
							break;
						if (l = j0 - 1){
							adj0[j0 - j0_init] = col[i];	
							j0++;
							break;
						}
					}
				}else if (x_level[col[i]] == lvl + 2){
					for (l = l1; l < j1; l++){
						if (col[i] == numbered[l])
							break;
						if (l = j1-1){
							adj1[j1 - l1] = col[i];	
							j1++;
							break;
						}
					}
				}
			}
			degree_increasing_sort(adj0, j0 - j0_init, row);
			degree_increasing_sort(adj1, j1 - l1, row);
			for (i = j0_init; i < j0; i++)
				numbered[i] = adj0[i - j0_init];
			for (i = l1; i < j1; i++){
				numbered[i] = adj1[i - l1];
			}
		}
		if (j0 != l1){
			// search another level lvl col
		}
	}
}

void main(){
	ivector arr = {4,5,6,7,8,9,10,1,2,3};
	int* po = &arr[0];
	po+=5;
	printf("== %d ==", po[2]);
	return;
}

