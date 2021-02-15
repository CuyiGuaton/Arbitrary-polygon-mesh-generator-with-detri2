#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "consts.h"
#include "triang.h"
#include "polygon.h"
#include "mesh.h"
#include "BET_elimitation.h"

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define debug_block(fmt) do { if (DEBUG_TEST){ fmt }} while (0)
#define debug_print(fmt, ...) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)
#define debug_msg(fmt) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__,  __LINE__, __func__); } while (0)


//Remove_BE
int Remove_BE(int option, int *poly, int length_poly, int num_BE, int *triangles, int *adj, double *r, int tnumber, int *mesh, int i_mesh, int* trivertex){

    debug_msg("Removiendo barrier edge de "); debug_block(print_poly(poly, length_poly); );

    int *poly1 = (int *)malloc(length_poly*sizeof(int));
	int *poly2 = (int *)malloc(length_poly*sizeof(int));
	int length_poly1;
	int length_poly2;
    int num_BE_poly1;
    int num_BE_poly2;
    int v_be, v_other;
    int t1, t2;


    v_be = get_vertex_BarrierEdge(poly, length_poly);
    
    if(option == 0){
        //t1 = search_triangle_by_vertex_with_FrontierEdge(v_be, triangles, adj, tnumber);
        t1 = search_triangle_by_vertex_with_FrontierEdge_from_trivertex(v_be, triangles, adj, tnumber, trivertex);
        v_other = optimice1_max_area_criteria(&t1, v_be, poly, length_poly, poly1, &length_poly1, poly2, &length_poly2, num_BE, triangles, adj, r, tnumber);
        t2 = get_adjacent_triangle(t1, v_other, v_be, triangles, adj);
    }else if(option == 1){
        //t1 = search_triangle_by_vertex_with_FrontierEdge(v_be, triangles, adj, tnumber);
        t1 = search_triangle_by_vertex_with_FrontierEdge_from_trivertex(v_be, triangles, adj, tnumber, trivertex);
        v_other = optimice2_middle_edge(&t1, v_be, triangles, adj);
        //printf("v %d, t %d  | Triangles %d %d %d | ADJ  %d %d %d\n", v_be, t1, triangles[3*t1 + 0], triangles[3*t1 + 1], triangles[3*t1 + 2], adj[3*t1 + 0], adj[3*t1 + 1], adj[3*t1 + 2]);
        t2 = get_adjacent_triangle(t1, v_other, v_be, triangles, adj);
    }

    debug_print("Eliminando arista de %d - %d de los triangulos  %d y %d", v_be, v_other, t1, t2);
    adj[3*t1 + get_edge(t1, v_be, v_other, triangles)] = NO_ADJ;
    adj[3*t2 + get_edge(t2, v_be, v_other, triangles)] = NO_ADJ;

    //posible bug, si se parte un BE puede omitir un vertice
    //solución, verificar si el triangulo que particiona es válido, sino se cambia.


    length_poly1 = generate_polygon(t1, poly1, triangles, adj, r);
    length_poly2 = generate_polygon(t2, poly2, triangles, adj, r);

    
    num_BE_poly1 = count_BarrierEdges(poly1, length_poly1);
    num_BE_poly2 = count_BarrierEdges(poly2, length_poly2);

    debug_print("num_BE_poly1 %d, num_BE_poly2 %d\n", num_BE_poly1, num_BE_poly2);

    if(num_BE_poly1 > 0 && num_BE_poly2 == 0){						
        debug_msg("Guardando poly2 y enviando recursivamente poly1\n");
        i_mesh = save_to_mesh(mesh, poly2, i_mesh, length_poly2);	
        i_mesh = Remove_BE(option, poly1, length_poly1, num_BE_poly1, triangles, adj, r, tnumber, mesh, i_mesh, trivertex);
    }else if(num_BE_poly2 > 0 && num_BE_poly1 == 0){
        debug_msg("Guardando poly1 y enviando recursivamente poly2\n");
        i_mesh = save_to_mesh(mesh, poly1, i_mesh, length_poly1);	
        i_mesh = Remove_BE(option, poly2, length_poly2, num_BE_poly2, triangles, adj, r, tnumber, mesh, i_mesh, trivertex);
    }else if(num_BE_poly1 > 0 && num_BE_poly2 > 0){
        debug_msg("Enviando recursivamente poly1 y poly2\n");
        i_mesh = Remove_BE(option, poly1, length_poly1, num_BE_poly1, triangles, adj, r, tnumber, mesh, i_mesh, trivertex);
        i_mesh = Remove_BE(option, poly2, length_poly2, num_BE_poly2, triangles, adj, r, tnumber, mesh, i_mesh, trivertex);
    }else{
        debug_msg("Guardando poly1 y poly2\n");
        i_mesh = save_to_mesh(mesh, poly1, i_mesh, length_poly1);	
        i_mesh = save_to_mesh(mesh, poly2, i_mesh, length_poly2);	
    }
    free(poly1);
    free(poly2);
    return i_mesh;
}



/* Method 1 to split polygon, max area criteria */

/* Dado un poly con barrier edges
Optimiza la división de este y devuelve poly1 y poly2*/
int optimice1_max_area_criteria(int *t_original, int v_be, int *poly, int length_poly, int *poly1, int *length_poly1, int *poly2, int *length_poly2, int num_BE, int *triangles, int *adj, double *r, int tnumber){
    double A_poly, A1, A2, opt, r_prev, r_act;
    int v_other, aux, origen, t;
    t = *t_original;

    /*se calculca el valor optimo para el poligono */
    A_poly = get_signed_area_poly(poly, length_poly,r);
    opt = fabs(A_poly/(num_BE+1));

    debug_print("Area poly: %.2lf, opt = %.2lf\n", A_poly, opt);

    /*se calcula el otro vertice para partir poly*/
    v_other = search_next_vertex_to_split(t, v_be, -2, triangles, adj);

    debug_print("Agregar edge %d - %d del Triangulo %d \n", v_be, v_other,t);
    
    if(v_other == -1 || v_other == -2){
        debug_print("No se puede agregar el edge %d - %d del Triangulo %d\n", v_be, v_other, t);
        if(v_other == -2)
            fprintf(stderr, "Caso critico especial, no encuentra vertices para avanzar en la busqueda de eliminación de barries edge, pero es la primera iteración");
        exit(0);
    }

    debug_msg("Dividiendo poligono\n");
    /* Se divide el polygono en dos */
    split_poly(poly, length_poly, poly1, &(*length_poly1), poly2, &(*length_poly2), v_be, v_other);

    debug_msg("poly1: "); debug_block(print_poly(poly1, *length_poly1););
	debug_msg("poly2: "); debug_block( print_poly(poly2, *length_poly2););


    A1 = get_signed_area_poly(poly1, *length_poly1,r);
    A2 = get_signed_area_poly(poly2, *length_poly2,r);

    /* se calcula el r */
    r_prev = fabs(fmin(fabs(A1), fabs(A2)) - opt);
    r_act = 0.0;

    debug_print("A: %.2lf, A1: %.2lf, A2:  %.2lf, A1/A = %.2lf, A2/A = %.2lf, r_prev = %.2lf, r_act = %.2lf\n", A_poly, A1 , A2,  A1/A_poly, A2/A_poly, r_prev, r_act);
    origen = t;
    while (1){
        
        aux = t;
        t = get_adjacent_triangle_share_endpoint(t, origen, v_be, triangles, adj);
        origen = aux;
        v_other = search_next_vertex_to_split(t, v_be, origen, triangles, adj);
        

        debug_print("Agregar edge %d - %d del nuevo triangulo %d | origen = %d \n", v_be, v_other,t, origen); 
        
        if(v_other != -2){
            debug_msg("Dividiendo poligono de nuevo\n");

            split_poly(poly, length_poly, poly1, &(*length_poly1), poly2, &(*length_poly2), v_be, v_other);
            A1 =get_signed_area_poly(poly1, *length_poly1,r);
            A2 = get_signed_area_poly(poly2,*length_poly2,r);
            debug_msg("poly1: "); debug_block(print_poly(poly1, *length_poly1););
            debug_msg("poly2: "); debug_block( print_poly(poly2, *length_poly2););
            
            r_act = fabs(fmin(fabs(A1), fabs(A2)) - opt);
            debug_print("A: %.2lf, A1: %.2lf, A2:  %.2lf, A1/A = %.2lf, A2/A = %.2lf, r_prev = %.2lf, r_act = %.2lf\n", A_poly, A1 , A2,  A1/A_poly, A2/A_poly, r_prev, r_act);
        }
        if (r_act <= r_prev && v_other != -2){
            r_prev = r_act;
            debug_msg("Solución optima no encontrada, repitiendo\n");
        }else{
            debug_print("Se encontro la optimización con r_act %.2lf\n", r_act);
            v_other = search_prev_vertex_to_split(t, v_be, origen, triangles, adj);
            *t_original = t;
            return v_other;
            /*
            debug_print("Agregar edge %d - %d del nuevo triangulo %d | origen = %d \n", v_be, v_other,t, origen);
            split_poly(poly, length_poly, poly1, &(*length_poly1), poly2, &(*length_poly2), v_be, v_other);
            debug_msg("Poligonos generados\n");
			debug_msg("poly1: "); debug_block(print_poly(poly1, *length_poly1););
			debug_msg("poly2: "); debug_block( print_poly(poly2, *length_poly2); );
            debug_block(
            A1 =get_signed_area_poly(poly1, *length_poly1,r);
            A2 = get_signed_area_poly(poly2,*length_poly2,r););
            debug_print("A: %.2lf, A1: %.2lf, A2:  %.2lf, A1/A = %.2lf, A2/A = %.2lf, r_prev = %.2lf, r_act = %.2lf\n", A_poly, A1 , A2,  A1/A_poly, A2/A_poly, r_prev, r_act);
            debug_msg("División optima terminada\n");
            */
            
        }
        
    }
    exit(0);
    return EXIT_FAILURE;
}

int optimice2_middle_edge(int *t_original, int v_be, int *triangles, int *adj){
    
    int aux, origen,i;
    int t_incident[128];
    int t = *t_original;
    
    origen = -1;
    i = 0;
    t_incident[0] = t; i++;
    
    while (1)
    {
        //advance once triangle
        aux = t;
        t = get_adjacent_triangle_share_endpoint(t, origen, v_be, triangles, adj);
        origen = aux;
        
        //if threre is not more triangle to see
        if( t < 0) break;
        //save triangle
        t_incident[i] = t;
        i++;
    }
    //print_poly(t_incident, i);
    if(i%2 == 0){ //if the triangles surrounding the BET are even 
        i = i/2;
        *t_original = t_incident[i];
        //Choose the edge in common with the two middle triangles
        return search_prev_vertex_to_split(t_incident[i+1], v_be, t_incident[i], triangles, adj);
    }else
    {   
        //if the triangles surrounding the BET are even 
        i = (i+1)/2;
        *t_original = t_incident[i];
        //Choose any edge of the triangle in the middle; prov is choose due to this always exists
        return search_next_vertex_to_split(t_incident[i], v_be, t_incident[i-1], triangles, adj);
    }   
}