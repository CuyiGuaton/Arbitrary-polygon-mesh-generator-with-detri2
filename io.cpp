/* Funciones para manejar entrada y salida de datos. */

#define _POSIX_C_SOURCE 200112L /* Para popen, pclose y posix_memalign. */

#include <iostream>
#include <vector> 
#include <algorithm>


#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "consts.h"
#include "triang.h"

#define filespath "input/"
#define filespathoutput "output/"

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define debug_block(fmt) do { if (DEBUG_TEST){ fmt }} while (0)
#define debug_print(fmt, ...) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)
#define debug_msg(fmt) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__,  __LINE__, __func__); } while (0)

/*geomview output*/
void write_geomview(double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, int print_triangles){

    int i,j;
    char cmd[1024] = "\0";
    strcat(cmd, filespathoutput);
    //strcat(cmd, ppath);
    strcat(cmd,".off");
    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    
    //imprimir puntos
    fprintf(fptr, "{ appearance  {+edge +face linewidth 2} LIST\n");
    fprintf(fptr, "OFF\n");
    fprintf(fptr,"%d %d 0\n", pnumber, num_region);
    for(i = 0; i < pnumber; i++)
        fprintf(fptr,"%.3f %.3f 0\n", r[2*i + 0], r[2*i + 1]);

  //imprimir polginos
    i = 0;
    while(i < i_mesh){
        int length_poly = mesh[i];
        i++;
        fprintf(fptr, "%d ", length_poly);
        for(j=0; j < length_poly;j++){
            fprintf(fptr, "%d ", mesh[i]);
            i++;
        }
        fprintf(fptr, "\n");
    }

    if(print_triangles){
        
        fprintf(fptr, "{ appearance  {+edge -face linewidth 2} LIST\n");
        int p0, p1,p2;
        for(i = 0; i < tnumber; i++){
            p0 = 3*i + 0;
            p1 = 3*i + 1;
            p2 = 3*i + 2;
            fprintf(fptr,"# %d %d\n", p0, p1);
            fprintf(fptr,"VECT 1 2 1 2 1\n");
            fprintf(fptr,"%.3f %.3f 0\n%.3f %.3f 0\n",	r[2*triangles[p0]+0], r[2*triangles[p0]+1], 
                                r[2*triangles[p1]+0], r[2*triangles[p1]+1]);
            fprintf(fptr,"0 1 1 1\n");

            fprintf(fptr,"# %d %d\n", p1, p2);
            fprintf(fptr,"VECT 1 2 1 2 1\n");
            fprintf(fptr,"%.3f %.3f 0\n%.3f %.3f 0\n",	r[2*triangles[p1]+0], r[2*triangles[p1]+1], 
                                r[2*triangles[p2]+0], r[2*triangles[p2]+1]);
            fprintf(fptr,"0 1 1 1\n");

            fprintf(fptr,"# %d %d\n", p0, p2);
            fprintf(fptr,"VECT 1 2 1 2 1\n");
            fprintf(fptr,"%.3f %.3f 0\n%.3f %.3f 0\n",	r[2*triangles[p0]+0], r[2*triangles[p0]+1], 
                                r[2*triangles[p2]+0], r[2*triangles[p2]+1]);
            fprintf(fptr,"0 1 1 1\n");
        }
    }

    fprintf(fptr," }\n");
    if(print_triangles){
        fprintf(fptr," }\n");
        fprintf(fptr," }\n");
    }
    /*
    std::sort( border_point.begin(), border_point.end() );
    border_point.erase( std::unique( border_point.begin(), border_point.end() ), border_point.end() );
    fprintf(fptr,"#Border vertices\n#");
    for ( i=0; i<border_point.size(); i++)
    {
        fprintf(fptr,"%d ", border_point[i]);
    }
    */
    fprintf(fptr,"\n");
    fclose(fptr);
}
