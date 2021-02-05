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
void write_geomview(double *r, int *p,  int pnumber,
									int tnumber, int ind_poly, 
									int *polygon, int id_start_poly,int *start_polygon, int print_triangles, char *ppath){

	int i,j;
	char cmd[1024] = "\0";
	strcat(cmd, filespathoutput);
	//strcat(cmd, ppath);
	strcat(cmd,".off");
	debug_print("Imprimiendo archivo en %s\n", cmd);
	FILE *fptr;
    fptr = fopen(cmd, "w");
	if(fptr == NULL){
		printf("Hmnito, no abrio el archivo :C\n");
		perror("fopen");
		exit(0);
	}
	
	/*
	if (print_triangles)
	{
		debug_print("Imprimiendo %d triangulos semilla", id_chose_seed_triangle);

		fprintf(fptr, "{ appearance  {-edge +face linewidth 2} LIST\n");
		for(j = 0; j < id_chose_seed_triangle; j++){

			i = chose_seed_triangle[j];
			fprintf(fptr,"{ OFF 3 1 1\n");
			fprintf(fptr,"%.3f %.3f 0\n%.3f %.3f 0\n%.3f %.3f 0\n",   r[2*p[3*i + 0] + 0], r[2*p[3*i + 0] + 1], 
																r[2*p[3*i + 1] + 0], r[2*p[3*i + 1] + 1], 
																r[2*p[3*i + 2] + 0], r[2*p[3*i + 2] + 1] );
			
			fprintf(fptr,"3 0 1 2  1.0 0.0 1.0 1 }\n");
		}
	}*/
	debug_print("Imprimiendo %d poligonos\n", id_start_poly);
	fprintf(fptr, "{ appearance  {+edge +face linewidth 2} LIST\n");
	fprintf(fptr, "OFF\n");
	fprintf(fptr,"%d %d 0\n", pnumber, id_start_poly);
	for(i = 0; i < pnumber; i++)
		fprintf(fptr,"%.3f %.3f 0\n", r[2*i + 0], r[2*i + 1]);
	
	int anterior = 0;
	for(j = 0,i = 0; i< id_start_poly; i ++){
		fprintf(fptr,"%d ", start_polygon[i] - anterior);
		for(; j <  start_polygon[i] ; j = j + 1 ){
				fprintf(fptr,"%d ", polygon[j]);
				
		}
		anterior = start_polygon[i] ;
		fprintf(fptr,"\n");
	}

	
	if(print_triangles){
		debug_msg("Imprimiendo malla delaunay\n");
		fprintf(fptr, "{ appearance  {+edge -face linewidth 2} LIST\n");
		int p0, p1,p2;
		for(i = 0; i < tnumber; i++){
			p0 = 3*i + 0;
			p1 = 3*i + 1;
			p2 = 3*i + 2;
			fprintf(fptr,"# %d %d\n", p0, p1);
			fprintf(fptr,"VECT 1 2 1 2 1\n");
			fprintf(fptr,"%.3f %.3f 0\n%.3f %.3f 0\n",	r[2*p[p0]+0], r[2*p[p0]+1], 
													r[2*p[p1]+0], r[2*p[p1]+1]);
			fprintf(fptr,"0 1 1 1\n");
			
			fprintf(fptr,"# %d %d\n", p1, p2);
			fprintf(fptr,"VECT 1 2 1 2 1\n");
			fprintf(fptr,"%.3f %.3f 0\n%.3f %.3f 0\n",	r[2*p[p1]+0], r[2*p[p1]+1], 
													r[2*p[p2]+0], r[2*p[p2]+1]);
			fprintf(fptr,"0 1 1 1\n");

			fprintf(fptr,"# %d %d\n", p0, p2);
			fprintf(fptr,"VECT 1 2 1 2 1\n");
			fprintf(fptr,"%.3f %.3f 0\n%.3f %.3f 0\n",	r[2*p[p0]+0], r[2*p[p0]+1], 
													r[2*p[p2]+0], r[2*p[p2]+1]);
			fprintf(fptr,"0 1 1 1\n");
		}
	}

	fprintf(fptr," }\n");
	if(print_triangles){
		fprintf(fptr," }\n");
		fprintf(fptr," }\n");
	}
	
	fprintf(fptr,"\n");
    fclose(fptr);
//printf("%d %d %d", 6.928203230275509 == 6.928203230275509, 6.928203230275509 > 6.928203230275509, 6.928203230275509 >= 6.928203230275509);
//printf("%d %d %d", 6.928203230275519 == 6.928203230275509, 6.928203230275519 > 6.928203230275509, 6.928203230275519 >= 6.928203230275509);
}

