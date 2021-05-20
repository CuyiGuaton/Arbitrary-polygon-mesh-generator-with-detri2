/* Funciones para manejar entrada y salida de datos. */

#define _POSIX_C_SOURCE 200112L /* Para popen, pclose y posix_memalign. */

#include <iostream>
#include <vector> 
#include <algorithm>
#include <set>
#include <list>


#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "consts.h"
#include "triang.h"
#include "io.h"

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


/*geomview output*/
void write_VEM(double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, int print_triangles){

    int i,j;
    char cmd[1024] = "\0";
    strcat(cmd, filespathoutput);
    //strcat(cmd, ppath);
    strcat(cmd,"vertice.txt");
    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    
    //imprimir puntos
    fprintf(fptr,"%d %d\n", pnumber, num_region);
    for(i = 0; i < pnumber; i++)
        fprintf(fptr,"%.3f %.3f\n", r[2*i + 0], r[2*i + 1]);

  //imprimir polginos
    i = 0;
    while(i < i_mesh){
        int length_poly = mesh[i];
        i++;
        fprintf(fptr, "%d ", length_poly);
        for(j=0; j < length_poly;j++){
            fprintf(fptr, "%d ", mesh[i] + 1);
            i++;
        }
        fprintf(fptr, "\n");
    }


    fprintf(fptr,"\n");
    fclose(fptr);
}

/*geomview output*/
void write_VEM_triangles(double *r, int *triangles, int *adj, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, std::list <int> &seed_bet){

    int i,j, lenght_poly;
    char cmd[1024] = "\0";
    strcat(cmd, filespathoutput);
    //strcat(cmd, ppath);
    strcat(cmd,"triangulos.txt");
    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    
    //imprimir puntos
    fprintf(fptr,"%d %d %d 0\n", pnumber, tnumber, num_region);
    for(i = 0; i < pnumber; i++)
        fprintf(fptr,"%.3f %.3f\n", r[2*i + 0], r[2*i + 1]);

        //imprimir triangulos
    for(i = 0; i < tnumber; i++)
        fprintf(fptr,"%d %d %d\n", triangles[3*i + 0] +1, triangles[3*i + 1] +1, triangles[3*i + 2] +1);

  //imprimir polginos
    std::set<int> s;
    j = 1;
    for( i = 0; i < tnumber; i++)
	{
		if(seed[i] == TRUE ){
			lenght_poly = look_triangles(i, s, triangles, adj, r);
            fprintf(fptr, "%d %ld ", j, s.size());
            
            for ( auto it = s.begin(); it != s.end(); it++ )
                fprintf(fptr, "%d ", *it + 1);
            fprintf(fptr, "\n");
            j++;
            s.clear();
        }
    }
    std::list <int> :: iterator ite;
    std::cout<<"iterando "<<seed_bet.size()<<" de "<<j<<std::endl;
    for(ite = seed_bet.begin(); ite != seed_bet.end(); ++ite){
            std::cout<<"seed: "<<*ite<<std::endl;
			lenght_poly = look_triangles(*ite, s, triangles, adj, r);
            fprintf(fptr, "%d %ld ", j, s.size());
            
            for ( auto it = s.begin(); it != s.end(); it++ )
                fprintf(fptr, "%d ", *it + 1);
            fprintf(fptr, "\n");
            j++;
            s.clear();
    }

    fprintf(fptr,"\n");
    fclose(fptr);
}

int look_triangles(int i, std::set<int> &s, int * triangles, int * adj, double *r) {

    s.insert(i);

    int ind_poly = 0;
	
	int initial_point = 0;
	int end_point = 0;
	
	int t0;
	int t1;	
	int t2;
    int ind0;
    int ind1;
    int ind2;
	int continuous;
	int k, j, aux;
	int origen;

    int num_FrontierEdges = count_FrontierEdges(i, adj);
    debug_print("Generando polinomio con triangulo %d FE %d\n", i, num_FrontierEdges);
    /*si tiene 3 se agregan y se corta el ciclo*/
    if (num_FrontierEdges == 3) {
        debug_print("T %d Tiene 3 Frontier edge, se guardan así\n", i);
        //poly[ind_poly] = triangles[3 * i + 0];
        ind_poly++;
        ////poly[ind_poly] = triangles[3 * i + 1];
        ind_poly++;
        ////poly[ind_poly] = triangles[3 * i + 2];
        ind_poly++;

        //visited[i] = TRUE;
        return ind_poly;
    } else if(num_FrontierEdges == 2) {
        debug_print("T %d Tiene 2 Frontier edge, es oreja, se usa como semilla para generar el poly\n", i);
        /*si tiene dos FE se agregan y se empieza el ciclo*/
        for(j = 0; j<3; j++){
            ind0 = 3*i + j;
            ind1 = 3*i + (j+1)%3;
            ind2 = 3*i + (j+2)%3;
            if(adj[ind0] == NO_ADJ && adj[ind1] == NO_ADJ){
                ////poly[ind_poly] = triangles[ind1];
                ind_poly++;
                ////poly[ind_poly] = triangles[ind2];
                ind_poly++;

                initial_point = triangles[ind1];
                end_point = triangles[ind0];  
            }
        }
    }else if (num_FrontierEdges == 1){
        debug_print("T %d Tiene 1 Frontier edge,se usa como FE initial\n", i);
        /*si tiene dos FE se agregan y se empieza el ciclo*/
        for(j = 0; j<3; j++){
            if(adj[3*i + j] == NO_ADJ){
                ////poly[ind_poly] = triangles[3*i + (j+1)%3];
                ind_poly++;
                initial_point = triangles[3*i + (j+1)%3];

                end_point = triangles[3*i + (j+2)%3];  
            }
        }
    }else {
        end_point = triangles[3*i + 0];
        initial_point = triangles[3*i + 0];
    }
    
    /*se marca como visitado */
    //visited[i] = TRUE;
    num_FrontierEdges = 0;
    k = i;
    aux = k;
    k = get_adjacent_triangle_share_endpoint(k, k, end_point, triangles, adj); /* cambia el indice */
    continuous = is_continuous(k, end_point, triangles);
    origen = aux;
//        debug_print("k %d origen %d, conti %d\n", k, origen, continuous);
    debug_print("T_inicial %d | Triangles %d %d %d | ADJ  %d %d %d\n", i, triangles[3*i + 0], triangles[3*i + 1], triangles[3*i + 2], adj[3*i + 0], adj[3*i + 1], adj[3*i + 2]);
    debug_print("initial_point %d endpoint %d | T_sig %d\n", initial_point, end_point, k);

    int triangugulo_initial = i;
    while (initial_point != end_point || triangugulo_initial != k) {
        s.insert(k);

        /*se marca el triangulo visto como visitado y se suma al area del poligono */
        
      //  visited[k] = TRUE;
        t0 = adj[3 * k + 0];
        t1 = adj[3 * k + 1];
        t2 = adj[3 * k + 2];

        num_FrontierEdges = count_FrontierEdges(k, adj);
        debug_print("FE %d | origen %d t %d | Triangles %d %d %d | ADJ  %d %d %d\n", num_FrontierEdges, origen, k, triangles[3*k + 0], triangles[3*k + 1], triangles[3*k + 2], adj[3*k + 0], adj[3*k + 1], adj[3*k + 2]);
        if(origen == -2)
            exit(0);
        if (num_FrontierEdges == 2 && continuous != -1) {
            /* ///////////////////si tiene 2 frontier edge se agregan a poly //////////////////////////////////// */

            if (t0 == NO_ADJ && t1 == NO_ADJ) {
                /*si endpoint es continua a t0  y t0-t1 son fe*/
                if (continuous == 1) {
                    ////poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;
                    ////poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                } else if (continuous == 0) {
                    ////poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;
                    ////poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];
                }
            } else if (t2 == NO_ADJ && t0 == NO_ADJ) {
                /*si endpoint es continua a t2  y t2-t0 son fe*/
                if (continuous == 0) {
                    ////poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;
                    ////poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                } else if (continuous == 2) {
                    ////poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;
                    ////poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                }
            } else if (t1 == NO_ADJ && t2 == NO_ADJ) {
                /*si endpoint es continua a t1 y t1-t2 son fe*/
                if (continuous == 2) {
                    ////poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;
                    ////poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                } else if (continuous == 1) {
                    ////poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;
                    //poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                }
            } else {
                fprintf(stderr, "** ERROR ** Adding 2 fronter edges\n");
                fprintf(stderr, "** ERROR ** k %d t %d %d %d ini %d end %d \n", k, t0, t1, t2, initial_point, end_point);
            }

            aux = k;
            k = get_adjacent_triangle_share_endpoint(k, -1, end_point, triangles, adj); /* se le permite volver al triangulo anterior */
            continuous = is_continuous(k, end_point, triangles);
            origen = aux;

        } else if (num_FrontierEdges == 1 && continuous != -1) {
            /* ///////////////////si solo se tiene 1 frontier edge //////////////////////////////////// */
            if (t0 == NO_ADJ) {
                /*si endpoint es continua a t0  y t0 es fe*/
                if (continuous == 1) {
                    //poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                } else if (continuous == 2) {
                    //poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                }
            } else if (t2 == NO_ADJ) {
                /*si endpoint es continua a t2  y t2 es fe*/
                if (continuous == 0) {
                    //poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                } else if (continuous == 1) {
                    //poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                }
            } else if (t1 == NO_ADJ) {
                /*si endpoint es continua a t1  y t1 es fe*/
                if (continuous == 2) {
                    //poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                } else if (continuous == 0) {
                    //poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                }
            } else {
                fprintf(stderr, "** ERROR ** Adding 1 fronter edges\n");
                fprintf(stderr, "** ERROR ** k %d t %d %d %d ini %d end %d conti %d\n", k, t0, t1, t2, initial_point, end_point, continuous);
            }
            /*si es continuo y tiene 1 fe no puede volver, ind si se guarda  o no*/
            aux = k;
            k = get_adjacent_triangle_share_endpoint(k, origen, end_point, triangles, adj); /* cambia el indice */
            continuous = is_continuous(k, end_point, triangles);
            origen = aux;
        } else {
            /*si no es continuo no puede regresar de donde venía*/
            aux = k;
            k = get_adjacent_triangle_share_endpoint(k, origen, end_point, triangles, adj); /* cambia el indice */
            continuous = is_continuous(k, end_point, triangles);
            origen = aux;
        }

    }
    
    return ind_poly;
}