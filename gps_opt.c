#include "gps_opt.h"

#include <stdio.h>
#include <string.h>

#include "print_util.h"

#define PDTMAX 10

void minimal_digree_all(crsrow row, int *mdi) {
    int i, rn;
    int md = N;
    for (i = 0; i < N; i++) {
        rn = row[i + 1] - row[i];
        if (md > rn) {
            md = rn;
            *mdi = i;
        }
    }
    return;
}

int minimal_digree(crsrow row, ivector level, const int si, const int ei) {
    int rn;
    int md = N;
    for (int i = 0; i < N; i++) {
        rn = row[i + 1] - row[i];
        if (md > rn) md = rn;
    }
    return md;
}

int level_width(level_index_structure level_index, int level_depth) {
    int i, width;
    int max_width = 0;
    for (i = 1; i < level_depth; i++) {
        width = level_index[i + 1] - level_index[i];
        if (max_width < width) max_width = width;
    }
    return max_width;
}

int generate_level_structure(ivector level, level_index_structure level_index,
                             ivector x_level, crscol col, crsrow row, int v) {
    // level: index sorted by level increasing order, level_index: the index
    // where the level is changed x_level: the level of each index Return: Level
    // Number
    memset(x_level, 0, N * sizeof(int));
    int i, j, s, t, current_level, leveled_x;
    //	int pre_leveled_x=0;
    level[0] = v;
    x_level[v] = 1;
    current_level = 2;
    level_index[0] = 0;
    level_index[1] = 1;
    leveled_x = 1;
    // TODO どこにも隣接していない頂点が存在した場合どうするか
    while (leveled_x < N) {
        for (i = level_index[current_level - 2];
             i < level_index[current_level - 1]; i++) {
            s = level[i];
            for (j = row[s]; j < row[s + 1]; j++) {
                t = col[j];
                if (x_level[t] == 0) {
                    x_level[t] = current_level;
                    level[leveled_x] = t;
                    leveled_x++;
                    if (leveled_x == N) {
                        level_index[current_level] = N;
                        return current_level;
                    }
                }
            }
            level_index[current_level] = leveled_x;
        }
        //		if (pre_leveled_x==leveled_x)
        //		    break;
        //        pre_leveled_x = leveled_x;
        current_level++;
    }
    return -1;
}

void degree_increasing_sort(ivector i_arr, const int len, crsrow row) {
    // i_arr のindexをdegree increasing の順にソートする。
    int i, j, tmp;
    for (i = 0; i < len; i++) {
        for (j = 1; j < len - i; j++) {
            if (row[i_arr[j] + 1] - row[i_arr[j]] <
                row[i_arr[j - 1] + 1] - row[i_arr[j - 1]]) {
                tmp = i_arr[j];
                i_arr[j] = i_arr[j - 1];
                i_arr[j - 1] = tmp;
            }
        }
    }
}

void pseudo_diameter(int *level_index_v, int *level_index_u, int *x_level_v,
                     int *x_level_u, int *level_v, int *level_u, crsrow row,
                     crscol col, int *u, int *v, int *level_depth) {
    int i, j, k, new_level_depth, u_width, min_u_width, u_index;
    level_index_structure level_index_0, level_index_1;
    ivector level_0, level_1, x_level_0, x_level_1;
    int *po_0, *po_1;
    po_0 = &level_0[0];
    po_1 = &level_1[0];
    minimal_digree_all(row, v);
    new_level_depth = generate_level_structure(level_1, level_index_1,
                                               x_level_1, col, row, *v);
    //	for (i = 0; i<N; i++){
    //        new_level_depth = generate_level_structure(level_1, level_index_1,
    //        x_level_1, col, row, i); printf("level=%d\n", new_level_depth);
    //	}
    *level_depth = new_level_depth;
    for (i = 0; i < PDTMAX; i++) {
        min_u_width = NON_ZERO;
        //        printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        //	    printf("LOOP%d:v = %d, level_depth = %d\n", i, *v,
        //*level_depth);
        if (i % 2 == 0) {
            degree_increasing_sort(po_1 + level_index_1[new_level_depth - 1],
                                   N - level_index_1[new_level_depth - 1], row);

            for (j = level_index_1[new_level_depth - 1]; j < N; j++) {
                new_level_depth = generate_level_structure(
                    level_0, level_index_0, x_level_0, col, row, level_1[j]);
                //                printf("LOOP%d: new_level_depth = %d\n",
                //                level_1[j], new_level_depth);
                if (*level_depth < new_level_depth) {
                    *level_depth = new_level_depth;
                    *v = level_1[j];
                    break;
                }
                u_width = level_width(level_index_0, new_level_depth);
                if (min_u_width > u_width) {
                    min_u_width = u_width;
                    u_index = j;
                }
                if (j == N - 1) {
                    *u = level_1[u_index];
                    generate_level_structure(level_u, level_index_u, x_level_u,
                                             col, row, *u);
                    // TODO stop copying like this. use poiter.
                    for (k = 0; k < N; k++) {
                        level_index_v[k] = level_index_1[k];
                        level_v[k] = level_1[k];
                        x_level_v[k] = x_level_1[k];
                    }
                    return;
                }
            }
        } else {
            degree_increasing_sort(po_0 + level_index_0[new_level_depth - 1],
                                   N - level_index_0[new_level_depth - 1], row);
            for (j = level_index_0[new_level_depth - 1]; j < N; j++) {
                new_level_depth = generate_level_structure(
                    level_1, level_index_1, x_level_1, col, row, level_0[j]);
                //                printf("LOOP%d: new_level_depth = %d\n",
                //                level_0[j], new_level_depth);
                if (*level_depth < new_level_depth) {
                    *level_depth = new_level_depth;
                    *v = level_0[j];
                    break;
                }
                u_width = level_width(level_index_1, new_level_depth);
                if (min_u_width > u_width) {
                    min_u_width = u_width;
                    u_index = j;
                }
                if (j == N - 1) {
                    *u = level_0[u_index];
                    generate_level_structure(level_u, level_index_u, x_level_u,
                                             col, row, *u);
                    // TODO stop copying like this. use poiter.
                    for (k = 0; k < N; k++) {
                        level_index_v[k] = level_index_0[k];
                        level_v[k] = level_0[k];
                        x_level_v[k] = x_level_0[k];
                    }
                    return;
                }
            }
        }
    }
}

void minimizing_level_width(crsrow row, crscol col, ivector x_level,
                            level_index_structure level_index,
                            level_index_structure level_index_v,
                            level_index_structure level_index_u,
                            ivector level_v, ivector level_u, int *u, int *v,
                            int *level_depth, int *chosed_element) {
    int i, j, k, l, m, cen, c, max_v_width, max_u_width, tmp;
    cen = 0;
    memset(level_index, 0, N * sizeof(int));
    ivector x_level_v, x_level_u;
    ivector component_index = {0};
    ivector component_elements = {0};
    ivector component_sep_i = {0};
    ivector component_number = {0};
    ivector sorted_m = {0};
    pseudo_diameter(level_index_v, level_index_u, x_level_v, x_level_u, level_v,
                    level_u, row, col, u, v, level_depth);
    for (i = 0; i < *level_depth; i++) {
        level_index_v[i] = level_index_v[i + 1] - level_index_v[i];
        level_index_u[i] = level_index_u[i + 1] - level_index_u[i];
    }

    for (i = 0; i < N; i++) {
        x_level_u[i] = *level_depth + 1 - x_level_u[i];
        if (x_level_v[i] == x_level_u[i]) {
            x_level[i] = x_level_u[i];
            component_index[i] = -1;
            level_index[x_level[i] - 1]++;
        } else {
            cen++;
        }
    }

    m = 0;
    k = 0;
    for (i = 0; i < N; i++) {
        if (component_index[i] == 0) {
            m++;
            component_elements[k] = i;
            component_index[i] = m;
            l = k;
            k++;
            while (l < k) {
                for (j = row[component_elements[l]];
                     j < row[component_elements[l] + 1]; j++) {
                    c = col[j];
                    if (component_index[c] == 0) {
                        component_elements[k] = c;
                        component_index[c] = m;
                        k++;
                    }
                }
                l++;
            }
            component_sep_i[m] = l;
            cen -= l;
            if (cen == 0) break;
        }
    }
    for (i = 0; i < m; i++) {
        component_number[i] = component_sep_i[i + 1] - component_sep_i[i];
        sorted_m[i] = i;
    }
    for (i = 0; i < m; i++) {
        for (j = 1; j < m - i; j++) {
            if (component_number[j] > component_number[j - 1]) {
                tmp = component_number[j];
                component_number[j] = component_number[j - 1];
                component_number[j - 1] = tmp;
                tmp = sorted_m[j];
                sorted_m[j] = sorted_m[j - 1];
                sorted_m[j - 1] = tmp;
            }
        }
    }
    for (i = 0; i < m; i++) {
        max_v_width = 0;
        max_u_width = 0;
        for (j = 0; j < *level_depth; j++) {
            k = level_index_v[j] - level_index[j];
            if (k >= 1 && k + level_index[j] > max_v_width)
                max_v_width = k + level_index[j];
            k = level_index_u[*level_depth - j - 1] - level_index[j];
            if (k >= 1 && k + level_index[j] > max_u_width) max_u_width = k;
        }
        if (max_v_width < max_u_width) {
            if (i == 0) *chosed_element = 0;
            for (j = component_sep_i[sorted_m[i]];
                 j < component_sep_i[sorted_m[i] + 1]; j++) {
                x_level[component_elements[j]] =
                    x_level_v[component_elements[j]];
                level_index[x_level[component_elements[j]] - 1]++;
            }
        } else if (max_v_width > max_u_width) {
            if (i == 0) *chosed_element = 1;
            for (j = component_sep_i[sorted_m[i]];
                 j < component_sep_i[sorted_m[i] + 1]; j++) {
                x_level[component_elements[j]] =
                    x_level_u[component_elements[j]];
                level_index[x_level[component_elements[j]] - 1]++;
            }
        } else {
            //			max_v_width = 0;
            //			max_u_width = 0;
            //			for(j = 0; j < *level_depth; j++){
            //				k = level_index_v[j];
            //				if (k > max_v_width)
            //					max_v_width = k;
            //				k = level_index_u[*level_depth - j - 1];
            //				if (k > max_u_width)
            //					max_u_width = k;
            //			}
            //		    if (max_v_width <= max_u_width){
            if (i == 0) *chosed_element = 0;
            for (j = component_sep_i[sorted_m[i]];
                 j < component_sep_i[sorted_m[i] + 1]; j++) {
                x_level[component_elements[j]] =
                    x_level_v[component_elements[j]];
                level_index[x_level[component_elements[j]] - 1]++;
            }
            //			}else{
            //				if (i == 0)
            //					*chosed_element = 1;
            //				for(j = component_sep_i[sorted_m[i]]; j <
            //component_sep_i[sorted_m[i] + 1]; j++){
            //					x_level[component_elements[j]] =
            //x_level_u[component_elements[j]];
            //					level_index[x_level[component_elements[j]] -
            //1]++;
            //				}

            //			}
        }
    }
    return;
}

void gps(crsrow row, crscol col, ivector numbered) {
    int i, j, j0, j0_init, j1, k, l, l0, l1, n0, n1, l_first_flg, wcand, lvl,
        v_s, v_e, u_s, u_e, v, u, w, level_depth, tmp, choosed_element,
        uv_switched_flag;
    int min_degree, min_degree_n, md;
    level_index_structure level_index_v, level_index_u;
    ivector level_v, level_u;
    level_index_structure level_width = {0};
    l1 = 0;
    l_first_flg = 0;
    uv_switched_flag = 0;
    ivector x_level;
    ivector adj0 = {0};
    ivector adj1 = {0};
    minimizing_level_width(row, col, x_level, level_width, level_index_v,
                           level_index_u, level_v, level_u, &u, &v,
                           &level_depth, &choosed_element);
    if (row[u + 1] - row[u] < row[v + 1] - row[v]) {
        tmp = u;
        u = v;
        v = tmp;
        for (i = 0; i < N; i++) {
            x_level[i] = level_depth + 1 - x_level[i];
        }
        uv_switched_flag = 1;
    }
    numbered[0] = v;
    j1 = 1;
    if (uv_switched_flag) {
        v_s = N;
        u_s = N;
    } else {
        v_e = 0;
        u_e = 0;
    }
    for (lvl = 0; lvl < level_depth; lvl++) {
        l0 = l1;
        j0 = j1;
        l1 = l0 + level_width[lvl];
        j1 = l1;
        for (k = l0; k < j0; k++) {
            w = numbered[k];
            n0 = 0;
            n1 = 0;
            for (i = row[w]; i < row[w + 1]; i++) {
                if (x_level[col[i]] == lvl + 1) {
                    for (l = l0; l < j0; l++) {
                        if (col[i] == numbered[l]) {
                            l_first_flg = 1;
                            break;
                        }
                    }
                    if (!l_first_flg) {
                        adj0[n0] = col[i];
                        n0++;
                    } else {
                        l_first_flg = 0;
                    }

                } else if (x_level[col[i]] == lvl + 2) {
                    for (l = l1; l < j1; l++) {
                        if (col[i] == numbered[l]) {
                            l_first_flg = 1;
                            break;
                        }
                    }
                    if (!l_first_flg) {
                        adj1[n1] = col[i];
                        n1++;
                    } else {
                        l_first_flg = 0;
                    }
                }
            }
            degree_increasing_sort(adj0, n0, row);
            degree_increasing_sort(adj1, n1, row);
            for (i = 0; i < n0; i++) numbered[j0 + i] = adj0[i];
            for (i = 0; i < n1; i++) {
                numbered[j1 + i] = adj1[i];
            }
            j0 += n0;
            j1 += n1;
        }
        if (uv_switched_flag) {
            v_e = v_s;
            v_s = v_s - level_index_v[level_depth - lvl - 1];
            u_e = u_s;
            u_s = u_s - level_index_u[lvl];
        } else {
            v_s = v_e;
            v_e = v_e + level_index_v[lvl];
            u_s = u_e;
            u_e = u_e + level_index_u[level_depth - lvl - 1];
        }
        while (j0 < l1) {
            // search another level lvl col
            min_degree_n = -1;
            min_degree = N;
            for (i = v_s; i < v_e; i++) {
                wcand = level_v[i];
                if (x_level[wcand] == lvl + 1) {
                    for (j = l0; j < j0; j++) {
                        if (numbered[j] == wcand) {
                            l_first_flg = 1;
                            break;
                        }
                    }
                    if (!l_first_flg) {
                        md = row[wcand + 1] - row[wcand];
                        if (min_degree > md) {
                            min_degree = md;
                            min_degree_n = wcand;
                        }
                    } else {
                        l_first_flg = 0;
                    }
                }
            }
            for (i = u_s; i < u_e; i++) {
                wcand = level_u[N - 1 - i];
                if (x_level[wcand] == lvl + 1) {
                    for (j = l0; j < j0; j++) {
                        if (numbered[j] == wcand) {
                            l_first_flg = 1;
                            break;
                        }
                    }
                    if (!l_first_flg) {
                        md = row[wcand + 1] - row[wcand];
                        if (min_degree > md) {
                            min_degree = md;
                            min_degree_n = wcand;
                        }
                    } else {
                        l_first_flg = 0;
                    }
                }
            }
            numbered[j0] = min_degree_n;
            //			printf("numbered[%d]=%d\n",j0, numbered[j0]);
            j0_init = j0;
            j0++;
            for (k = j0_init; k < j0; k++) {
                w = numbered[k];
                n0 = 0;
                n1 = 0;
                for (i = row[w]; i < row[w + 1]; i++) {
                    if (x_level[col[i]] == lvl + 1) {
                        for (l = l0; l < j0; l++) {
                            if (col[i] == numbered[l]) break;
                            if (l == j0 - 1) {
                                adj0[n0] = col[i];
                                n0++;
                                break;
                            }
                        }
                    } else if (x_level[col[i]] == lvl + 2) {
                        for (l = l1; l < j1; l++) {
                            if (col[i] == numbered[l]) break;
                            if (l == j1 - 1) {
                                adj1[n1] = col[i];
                                n1++;
                                break;
                            }
                        }
                    }
                }
                degree_increasing_sort(adj0, n0, row);
                degree_increasing_sort(adj1, n1, row);
                for (i = 0; i < n0; i++) {
                    numbered[j0 + i] = adj0[i];
                }

                for (i = 0; i < n1; i++) {
                    numbered[j1 + i] = adj1[i];
                }
                j0 += n0;
                j1 += n1;
            }
        }
    }
    if (!(uv_switched_flag + choosed_element) % 2) {
        for (i = 0; i < N / 2; i++) {
            j = numbered[i];
            numbered[i] = numbered[N - 1 - i];
            numbered[N - 1 - i] = j;
        }
    }
    return;
}
void opt_mat(ivector numbered, crsdata data, crsdata optdata, crscol col,
             crscol optcol, crsrow row, crsrow optrow) {
    ivector optrow_n, i_numbered;
    int i, ni, j, k, k_init, l, tmp;
    double tmpdata;

    for (i = 0; i < N; i++) {
        ni = numbered[i];
        i_numbered[ni] = i;
        optrow_n[i] = row[ni + 1] - row[ni];
        optrow[i + 1] = optrow[i] + optrow_n[i];
    }
    k = 0;
    for (i = 0; i < N; i++) {
        k_init = k;
        for (j = row[numbered[i]]; j < row[numbered[i] + 1]; j++) {
            optcol[k] = i_numbered[col[j]];
            optdata[k] = data[j];
            for (l = k; l > k_init; l--) {
                if (optcol[l] < optcol[l - 1]) {
                    tmp = optcol[l];
                    optcol[l] = optcol[l - 1];
                    optcol[l - 1] = tmp;
                    tmpdata = optdata[l];
                    optdata[l] = optdata[l - 1];
                    optdata[l - 1] = tmpdata;
                }
            }
            k++;
        }
    }
}
