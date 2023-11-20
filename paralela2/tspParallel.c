/* WSCAD - 9th Marathon of Parallel Programming 
 * Simple Brute Force Algorithm for the 
 * Traveling-Salesman Problem
 * Author: Emilio Francesquini - francesquini@ic.unicamp.br
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
typedef struct {
    int to_town;
    int dist;
} d_info;

d_info **d_matrix;
int *dist_to_origin;
int nb_towns;
int myrank,NumProcs;
int present (int town, int *path) {
    if (path[town]) 
        return 1;
    return 0;
}


void tsp (int depth, int current_length, int *path,int *pathPresent,int *min_distance) {
    int i;
    if (current_length >= *min_distance) return;
    if (depth == nb_towns) {
        current_length += dist_to_origin[path[nb_towns - 1]];
        if (current_length < *min_distance)
            *min_distance = current_length;
    } 
    if(depth==1){
        int town, me, dist;
        me = path[depth - 1];
        for (i = myrank; i < nb_towns ; i+=NumProcs) {
            town = d_matrix[me][i].to_town;
            if (!present (town,pathPresent)) {
                path[depth] = town;
                pathPresent[town]=1;
                dist = d_matrix[me][i].dist;
                tsp (depth + 1, current_length + dist, path,pathPresent,min_distance);
                pathPresent[town]=0;
            }
        }
    }
    else{
        int town, me, dist;
        me = path[depth - 1];
        for (i = 0; i < nb_towns ; i++) {
            town = d_matrix[me][i].to_town;
            if (!present (town,pathPresent)) {
                path[depth] = town;
                pathPresent[town]=1;
                dist = d_matrix[me][i].dist;
                tsp (depth + 1, current_length + dist, path,pathPresent,min_distance);
                pathPresent[town]=0;
            }
        }
    }
}


void greedy_shortest_first_heuristic(int *x, int *y) {
    int i, j, k, dist;
    int *tempdist;

    tempdist = (int*) malloc(sizeof(int) * nb_towns);
    //Could be faster, albeit not as didactic.
    //Anyway, for tractable sizes of the problem it
    //runs almost instantaneously.
    for (i = 0; i < nb_towns; i++) {
        for (j = 0; j < nb_towns; j++) {
            int dx = x[i] - x[j];
            int dy = y[i] - y[j];
            tempdist [j] = dx * dx + dy * dy;
        }
        for (j = 0; j < nb_towns; j++) {
            int tmp = INT_MAX;
            int town = 0;
            for (k = 0; k < nb_towns; k++) {
                if (tempdist [k] < tmp) {
                    tmp = tempdist [k];
                    town = k;
                }
            }
            tempdist [town] = INT_MAX;
            d_matrix[i][j].to_town = town;
            dist = (int) sqrt (tmp);
            d_matrix[i][j].dist = dist;
            if (i == 0)
                dist_to_origin[town] = dist;
        }
    }

    free(tempdist);
}

void init_tsp(int *min_distance) {
    int st,*x, *y;

    *min_distance = INT_MAX;
    if(myrank==0){  
        st=scanf("%u", &nb_towns);
        if(st!=1)exit(1);
        MPI_Bcast(&nb_towns,1,MPI_INT,0,MPI_COMM_WORLD);
    }else{
        MPI_Bcast(&nb_towns,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    d_matrix = (d_info**) malloc (sizeof(d_info*) * nb_towns);
    for (int i = 0; i < nb_towns; i++)
        d_matrix[i] = (d_info*) malloc (sizeof(d_info) * nb_towns);
    dist_to_origin = (int*) malloc(sizeof(int) * nb_towns);
   
    x = (int*) malloc(sizeof(int) * nb_towns);
    y = (int*) malloc(sizeof(int) * nb_towns);
    
    if(myrank==0){
        for (int i = 0; i < nb_towns; i++) {
            st = scanf("%u %u", x + i, y + i);
            if (st != 2) exit(1);
        }
        MPI_Bcast(x,nb_towns,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(y,nb_towns,MPI_INT,0,MPI_COMM_WORLD);
    }else{
        MPI_Bcast(x,nb_towns,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(y,nb_towns,MPI_INT,0,MPI_COMM_WORLD);
    }
    greedy_shortest_first_heuristic(x, y);
    free(x);
    free(y);
}

int run_tsp(int *min_distance) {
    int i, *path,*pathPresent;

    path = (int*) malloc(sizeof(int) * nb_towns);
    pathPresent = (int*) malloc(sizeof(int) * nb_towns);
    memset(path,0,sizeof(int) * nb_towns);
    memset(pathPresent,0,sizeof(int) * nb_towns);
    pathPresent[0]=1;
    tsp(1, 0, path, pathPresent, min_distance);
    free(path);
    free(pathPresent);
    for (i = 0; i < nb_towns; i++)
        free(d_matrix[i]);
    free(d_matrix);
    free(dist_to_origin);
    return *min_distance;
}

int main (int argc, char **argv) {
    int num_instances, st,min_distance,answer;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
    if(myrank==0){
        st = scanf("%u", &num_instances);
        if (st != 1) exit(1);
        MPI_Bcast(&num_instances,1,MPI_INT,0,MPI_COMM_WORLD);
    }else{
        MPI_Bcast(&num_instances,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    while(num_instances-- > 0){
        init_tsp(&min_distance);
        run_tsp(&min_distance);
        MPI_Reduce(&min_distance,&answer,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
        if(myrank==0){
            printf("Meu rank= %d, min dist = %d\n",myrank,answer);
        }
    }
    MPI_Finalize();

    return 0;
}