/* Funciones para manejar poligonos. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "consts.h"
#include "triang.h"
#include "polygon.h"

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define debug_block(fmt) do { if (DEBUG_TEST){ fmt }} while (0)
#define debug_print(fmt, ...) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)
#define debug_msg(fmt) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__,  __LINE__, __func__); } while (0)



//Verifica si es válido como triangulo semilla del poligono
int is_BarrierEdge(int i, int *adj, int *adj_copy, int *root_id){
    //debug_print("i = %d, root_id[i] = %d, root_id[t0_adj] = %d, root_id[t1_adj] = %d, root_id[t2_adj] = %d\n", i, root_id[i], root_id[t0_adj], root_id[t1_adj], root_id[t2_adj]  );
    //if( root_id[i] != root_id[t0_adj] && root_id[i] != root_id[t1_adj] && root_id[i] != root_id[t2_adj] ){
    int j;
    int num_fe = count_FrontierEdges(i, adj);
    debug_print("FE %d, Triangulo %d\n", num_fe, i);
    if(num_fe == 3){
        return 1;
    }else if (num_fe == 1)
    {
        for (j = 0; j < 3; j++)
        {
            if (adj_copy[3*i + j] != TRIANG_BORDER){
                if (root_id[i] !=  root_id[adj_copy[3*i + j]])
                    return 1;
            }else
            {
                return 1;
            }
            
                
        }   
    }else if(num_fe == 2)
    {
        for (j = 0; j < 3; j++)
        {
            int ind = 3*i + j;
            int ind2 = 3*i + (j + 1)%3;
            if (adj_copy[ind] != TRIANG_BORDER && adj_copy[ind2] != TRIANG_BORDER){
                if (root_id[i] !=  root_id[adj_copy[ind]] && root_id[i] !=  root_id[adj_copy[ind2]] )
                    return 1;
            }else if(adj_copy[ind] == TRIANG_BORDER){
                if (root_id[i] !=  root_id[adj_copy[ind2]] )
                    return 1;
            }else if(adj_copy[ind2] == TRIANG_BORDER){
                if (root_id[i] !=  root_id[adj_copy[ind]] )
                    return 1;
            }
        }   
    }
    return 0;
}

void remove_BarrierEdge(int *poly, int length_poly, int num_BE, int *triangles, int *adj, double *r, int tnumber, int *mesh, int *i_mesh, int *pos_poly, int *id_pos_poly, int *visited){

    debug_msg("Removiendo barrier edge de "); debug_block(print_poly(poly, length_poly););

    int *poly1 = (int *)malloc(length_poly*sizeof(int));
	int *poly2 = (int *)malloc(length_poly*sizeof(int));
	int length_poly1;
	int length_poly2;
    int num_BE_poly1;
    int num_BE_poly2;
    int v_be, v_other;
    int t1, t2;


    v_be = get_vertex_BarrierEdge(poly, length_poly);
    t1 = search_triangle_by_vertex_with_FrontierEdge(v_be, triangles, adj, tnumber);
    v_other = optimice_partition_polygon(&t1, v_be, poly, length_poly, poly1, &length_poly1, poly2, &length_poly2, num_BE, triangles, adj, r, tnumber);
    t2 = get_adjacent_triangle(t1, v_other, v_be, triangles, adj);

    debug_print("Eliminando arista de %d - %d de los triangulos  %d y %d", v_be, v_other, t1, t2);
    adj[3*t1 + get_edge(t1, v_be, v_other, triangles)] = NO_ADJ;
    adj[3*t2 + get_edge(t2, v_be, v_other, triangles)] = NO_ADJ;

    //posible bug, si se parte un BE puede omitir un vertice
    //solución, verificar si el triangulo que particiona es válido, sino se cambia.


    length_poly1 = generate_polygon(poly1, triangles, adj, r, visited, t1);
    length_poly2 = generate_polygon(poly2, triangles, adj, r, visited, t2);

    
    num_BE_poly1 = count_BarrierEdges(poly1, length_poly1);
    num_BE_poly2 = count_BarrierEdges(poly2, length_poly2);

    debug_print("num_BE_poly1 %d, num_BE_poly2 %d\n", num_BE_poly1, num_BE_poly2);




    if(num_BE_poly1 > 0 && num_BE_poly2 == 0){						
        debug_msg("Guardando poly2 y enviando recursivamente poly1\n");
        save_to_mesh(mesh, poly2, &(*i_mesh), length_poly2, pos_poly, &(*id_pos_poly));
        remove_BarrierEdge(poly1, length_poly1, num_BE_poly1, triangles, adj, r, tnumber, mesh, &(*i_mesh), pos_poly, &(*id_pos_poly), visited);
    }else if(num_BE_poly2 > 0 && num_BE_poly1 == 0){
        debug_msg("Guardando poly1 y enviando recursivamente poly2\n");
        save_to_mesh(mesh, poly1, &(*i_mesh), length_poly1, pos_poly, &(*id_pos_poly));
        remove_BarrierEdge(poly2, length_poly2, num_BE_poly2, triangles, adj, r, tnumber, mesh, &(*i_mesh), pos_poly, &(*id_pos_poly), visited);
    }else if(num_BE_poly1 > 0 && num_BE_poly2 > 0){
        debug_msg("Enviando recursivamente poly1 y poly2\n");
        remove_BarrierEdge(poly1, length_poly1, num_BE_poly1, triangles, adj, r, tnumber, mesh, &(*i_mesh), pos_poly, &(*id_pos_poly), visited);
        remove_BarrierEdge(poly2, length_poly2, num_BE_poly2, triangles, adj, r, tnumber, mesh, &(*i_mesh), pos_poly, &(*id_pos_poly), visited);
    }else{
        debug_msg("Guardando poly1 y poly2\n");
        save_to_mesh(mesh, poly1, &(*i_mesh), length_poly1, pos_poly, &(*id_pos_poly));
        save_to_mesh(mesh, poly2, &(*i_mesh), length_poly2, pos_poly, &(*id_pos_poly));
    }
    free(poly1);
    free(poly2);
}


/* Dado un poly con barrier edges
Optimiza la división de este y devuelve poly1 y poly2*/
int optimice_partition_polygon(int *t_original, int v_be, int *poly, int length_poly, int *poly1, int *length_poly1, int *poly2, int *length_poly2, int num_BE, int *triangles, int *adj, double *r, int tnumber){
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


    A1 =get_signed_area_poly(poly1, *length_poly1,r);
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
        }
        else{
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


/* Divide un poly dado un vertice e1-e2
    resultados poly1 y poly*/

void split_poly(int *original_poly, int length_poly, int *poly1, int *length_poly1, int *poly2, int *length_poly2, int e1, int e2){
    int pos1, pos2,i;
    pos1= -1;
    pos2= -1;
    for(i =0; i< length_poly; i++)
        if(original_poly[i] == e1 || original_poly[i] == e2){
            pos1 = i;
            break;
        }
    for(i =pos1 + 1; i< length_poly; i++)
        if(original_poly[i] == e1  || original_poly[i] == e2){
            pos2 = i;
            break;
        }
    debug_print("Divide pos1: %d, pos2: %d \n", pos1, pos2);
    if(pos1 == -1 || pos2 == -1){
        fprintf(stderr, "%s:%d:%s(): No se encontro pos1 o pos2\n",__FILE__,  __LINE__, __func__);
        exit(0);
    }
    
    *length_poly1 = abs(pos1-pos2) +1;
    *length_poly2 = length_poly - *length_poly1 +2;

    for (i = 0; i < *length_poly1 ; i++)
        poly1[i] = original_poly[(pos1 + i) %length_poly];

    for (i = 0; i < *length_poly2 ; i++)
        poly2[i] = original_poly[(pos2 + i) %length_poly];
}

int copy_poly(int *in, int *out, int len){
    int i;
    for (i = 0; i < len; i++)
        out[i] = in[i];
    return len;
}

void print_poly(int *poly, int length_poly){
    int i;
    fprintf(stderr,"(%d)", length_poly);
    for (i = 0; i < length_poly; i++)
        fprintf(stderr," %d", poly[i]); 
    fprintf(stderr,"\n");
}

double get_signed_area_poly(int *poly, int length_poly, double *r){
    double area = 0.0;
    double x1,y1,x2,y2;
    int i,j;
    for (i = 0; i < length_poly; i++)
    {
        j = (i+1) %length_poly;
        x1=r[2*poly[i] + 0];
        y1=r[2*poly[i] + 1];
        x2=r[2*poly[j] + 0];
        y2=r[2*poly[j] + 1];
        area += (x1 + x2)*(y2 - y1);
    }
    return area/2;
}


void save_to_mesh(int *mesh, int *poly, int *i_mesh, int length_poly, int *pos_poly, int *id_pos_poly){
    
    int i;
    for (i = 0; i < length_poly; i++){        
        mesh[*i_mesh + i] = poly[i];
    }
    *i_mesh += length_poly;
    pos_poly[*id_pos_poly] = *i_mesh;
	*id_pos_poly = *id_pos_poly + 1;
}


int count_BarrierEdges(int *poly, int length_poly){
    int count = 0;
    int x, y,i;
    for (i = 0; i < length_poly; i++)
    {
        x = i % length_poly;
        y = (i+2) % length_poly;
        if (poly[x] == poly[y])
            count++;
    }
    return count;
}

int get_vertex_BarrierEdge(int *poly, int length_poly){
    int x, y,i;
    for (i = 0; i < length_poly; i++)
    {
        x = i % length_poly;
        y = (i+2) % length_poly;
        if (poly[x] == poly[y])
            return poly[(i+1) %length_poly];
    }
    fprintf(stderr,"num_BE %d\n", count_BarrierEdges(poly, length_poly));
    fprintf(stderr,"%s:%d:%s(): No se encontro vertice BarrierEdge\n",__FILE__,  __LINE__, __func__);
    exit(0);
    return -1;
}


int generate_polygon(int * poly, int * triangles, int * adj, double *r, int * visited, int i) {
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
        poly[ind_poly] = triangles[3 * i + 0];
        ind_poly++;
        poly[ind_poly] = triangles[3 * i + 1];
        ind_poly++;
        poly[ind_poly] = triangles[3 * i + 2];
        ind_poly++;

        visited[i] = TRUE;
        return ind_poly;
    } else if(num_FrontierEdges == 2) {
        debug_print("T %d Tiene 2 Frontier edge, es oreja, se usa como semilla para generar el poly\n", i);
        /*si tiene dos FE se agregan y se empieza el ciclo*/
        for(j = 0; j<3; j++){
            ind0 = 3*i + j;
            ind1 = 3*i + (j+1)%3;
            ind2 = 3*i + (j+2)%3;
            if(adj[ind0] == NO_ADJ && adj[ind1] == NO_ADJ){
                poly[ind_poly] = triangles[ind1];
                ind_poly++;
                poly[ind_poly] = triangles[ind2];
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
                poly[ind_poly] = triangles[3*i + (j+1)%3];
                ind_poly++;
                initial_point = triangles[3*i + (j+1)%3];

                end_point = triangles[3*i + (j+2)%3];  
            }
        }
    }else {
        fprintf(stderr, "** ERROR ** num_FrontierEdges = %d", num_FrontierEdges);
    }
    
    
    /*se marca como visitado */
    visited[i] = TRUE;
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

        /*se marca el triangulo visto como visitado y se suma al area del poligono */
        
        visited[k] = TRUE;
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
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                } else if (continuous == 0) {
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];
                }
            } else if (t2 == NO_ADJ && t0 == NO_ADJ) {
                /*si endpoint es continua a t2  y t2-t0 son fe*/
                if (continuous == 0) {
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                } else if (continuous == 2) {
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                }
            } else if (t1 == NO_ADJ && t2 == NO_ADJ) {
                /*si endpoint es continua a t1 y t1-t2 son fe*/
                if (continuous == 2) {
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                } else if (continuous == 1) {
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 0];
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
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                } else if (continuous == 2) {
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                }
            } else if (t2 == NO_ADJ) {
                /*si endpoint es continua a t2  y t2 es fe*/
                if (continuous == 0) {
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                } else if (continuous == 1) {
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                }
            } else if (t1 == NO_ADJ) {
                /*si endpoint es continua a t1  y t1 es fe*/
                if (continuous == 2) {
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                } else if (continuous == 0) {
                    poly[ind_poly] = triangles[3 * k + 0];
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

/*

int generate_polygon(int * poly, int * triangles, int * adj, double *r, int * visited, int i) {
 
    while (initial_point != end_point || triangugulo_initial != k) {

        visited[k] = TRUE;
        t0 = adj[3 * k + 0];
        t1 = adj[3 * k + 1];
        t2 = adj[3 * k + 2];

        num_FrontierEdges = count_FrontierEdges(k, adj);
        debug_print("FE %d | origen %d t %d | Triangles %d %d %d | ADJ  %d %d %d\n", num_FrontierEdges, origen, k, triangles[3*k + 0], triangles[3*k + 1], triangles[3*k + 2], adj[3*k + 0], adj[3*k + 1], adj[3*k + 2]);
        if(origen == -2)
            exit(0);
        if (num_FrontierEdges == 2 && continuous != -1) {
             ///////////////////si tiene 2 frontier edge se agregan a poly //////////////////////////////////// 
            for(j = 0; j<3; j++){
                ind0 = 3*k + j;
                ind1 = 3*k + (j+1)%3;
                ind2 = 3*k + (j+2)%3;
                if(adj[ind0] == NO_ADJ && adj[ind1] == NO_ADJ){
                    if(continuous == (j+1)%3){
                        poly[ind_poly] = triangles[ind1];
                        ind_poly++;
                        poly[ind_poly] = triangles[ind2];
                        ind_poly++;

                        end_point = triangles[ind0];  
                    }else if(continuous == j){
                        poly[ind_poly] = triangles[ind0];
                        ind_poly++;
                        poly[ind_poly] = triangles[ind2];
                        ind_poly++;

                        end_point = triangles[ind1];   
                    }
                } else {
                fprintf(stderr, "** ERROR ** Adding 2 fronter edges\n");
                fprintf(stderr, "** ERROR ** k %d t %d %d %d ini %d end %d \n", k, t0, t1, t2, initial_point, end_point);
                }
            }

            aux = k;
            k = get_adjacent_triangle_share_endpoint(k, -1, end_point, triangles, adj); 
            // se le permite volver al triangulo anterior 
            continuous = is_continuous(k, end_point, triangles);
            origen = aux;

        } else if (num_FrontierEdges == 1 && continuous != -1) {
            ///////////////////si solo se tiene 1 frontier edge //////////////////////////////////// 
            //continuos es el vertice al que es continiuo dentro de un FE, el siguiente a este
            for(j = 0; j<3; j++){
                ind0 = 3*k + j;
                ind1 = 3*k + (j+1)%3;
                ind2 = 3*k + (j+2)%3;
                if(adj[ind0] == NO_ADJ){
                    if(continuous == (j+1)%3){
                        poly[ind_poly] = triangles[ind1];
                        ind_poly++;

                        end_point = triangles[ind2];  
                    }else if(continuous == (j+2)%3){
                        poly[ind_poly] = triangles[ind2];
                        ind_poly++;

                        end_point = triangles[ind1];  
                    }
                } else {
                fprintf(stderr, "** ERROR ** Adding 1 fronter edges\n");
                fprintf(stderr, "** ERROR ** k %d t %d %d %d ini %d end %d conti %d\n", k, t0, t1, t2, initial_point, end_point, continuous);
                }
            }
            //si es continuo y tiene 1 fe no puede volver, ind si se guarda  o no
            aux = k;
            k = get_adjacent_triangle_share_endpoint(k, origen, end_point, triangles, adj); // cambia el indice 
            continuous = is_continuous(k, end_point, triangles);
            origen = aux;
        } else {
            //si no es continuo no puede regresar de donde venía
            aux = k;
            k = get_adjacent_triangle_share_endpoint(k, origen, end_point, triangles, adj); // cambia el indice 
            continuous = is_continuous(k, end_point, triangles);
            origen = aux;
        }

    }
    
    return ind_poly;
}

*/