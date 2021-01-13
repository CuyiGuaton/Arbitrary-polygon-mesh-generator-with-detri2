#include <iostream>
#include <stdio.h>   
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include "detri2.h"
#include "polymesh.h"
#include <vector> 
#include <chrono>
#include <iomanip>
#include <cstdlib>



#include "io.h"
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


int main(int argc, char* argv[]){


    int nparam = 3;
    //char* params[] = {const_cast<char*> ("./detri2"), const_cast<char*> ("-z"), const_cast<char*> ("test.node")};
	//char* params[] = {const_cast<char*> ("./detri2"), const_cast<char*> ("-z"), const_cast<char*> ("506randompoints.node")};
	//char* params[] = {const_cast<char*> ("./detri2"), const_cast<char*> ("-z"), const_cast<char*> ("506equilateral.node")};
    int print_triangles = 0;
    char* ppath;
    //char* ppath = const_cast<char*> ("test");
    //TMesh *Tr = new TMesh(nparam, params);    
	auto tb_delaunay = std::chrono::high_resolution_clock::now();
	TMesh *Tr = new TMesh(argc, argv);    
	
	auto te_delaunay = std::chrono::high_resolution_clock::now();
    //Tr->print();
    
	int tnumber, pnumber, i,j;
	double *r;
	int *triangles;
	int *adj;
    int *visited;
    int *max;
	std::vector<int> border_point;

    tnumber = Tr->tnumber;
    pnumber = Tr->pnumber;

    max = (int *)malloc(tnumber*sizeof(int));
	visited = (int *)malloc(tnumber*sizeof(int));
    r = (double *)malloc(2*tnumber*sizeof(double));
    adj =(int *)malloc(3*tnumber*sizeof(int));
    triangles =   (int *)malloc(3*tnumber*sizeof(int));


	/* Se crean los poligonos */
	int *poly = (int *)malloc(tnumber*sizeof(int));
	int length_poly = 0;
	int *pos_poly = (int *)malloc(tnumber*sizeof(int));
	int id_pos_poly = 0;

	int *mesh;
	mesh = (int *)malloc(3*tnumber*sizeof(int));
	int i_mesh = 0;	
	
	int *chose_seed_triangle;
	chose_seed_triangle = (int *)malloc(tnumber*sizeof(int));
	int id_chose_seed_triangle = 0;
	int num_BE = 0;
	int est_total_be = 0;
	int est_min_triangles_be = 2147483641;
	int est_max_triangles_be = 0;
	int est_poly_with_be = 0;
	double est_ratio_be = 0;
	int tcost_be = 0;


	/* Llamada a detr2 */

    int idx =0;
    //copiar arreglo de vertices
    //std::cout<<"pnumber "<<pnumber<<std::endl;
    for (i = 0; i < Tr->trimesh->ct_in_vrts; i++) {
        if (!Tr->trimesh->io_keep_unused) { // no -IJ
            if (Tr->trimesh->in_vrts[i].typ == UNUSEDVERTEX) continue;
        }
        r[2*i + 0]= Tr->trimesh->in_vrts[i].crd[0];
        r[2*i + 1]= Tr->trimesh->in_vrts[i].crd[1];
        //std::cout<<idx<<" ("<<r[2*i + 0]<<", "<<r[2*i + 1]<<") "<<std::endl;
        Tr->trimesh->in_vrts[i].idx = idx;
        idx++;
    }
    idx = 0;
    for (int i = 0; i < Tr->trimesh->tr_tris->used_items; i++) {
        detri2::Triang* tri = (detri2::Triang *) Tr->trimesh->tr_tris->get(i);
        if (tri->is_deleted()) continue;
        if (tri->is_hulltri()) {
            tri->idx = -1;
        } else {
            tri->idx = idx;
            idx++;
        }
    }

    //std::cout<<"tnumber: "<<Tr->trimesh->tr_tris->objects - Tr->trimesh->ct_hullsize<<std::endl;
    idx = 0;
    for (int i = 0; i < Tr->trimesh->tr_tris->used_items; i++)
    {
        
        detri2::Triang* tri = (detri2::Triang *) Tr->trimesh->tr_tris->get(i);
        if (tri->is_deleted() || tri->is_hulltri()) continue;
        triangles[3*idx+0] = tri->vrt[0]->idx;
        triangles[3*idx+1] = tri->vrt[1]->idx;
        triangles[3*idx+2] = tri->vrt[2]->idx;
        adj[3*idx+ 0] = tri->nei[0].tri->idx;
        adj[3*idx+ 1] = tri->nei[1].tri->idx;
        adj[3*idx+ 2] = tri->nei[2].tri->idx;
        //std::cout<<idx<<" | "<<triangles[3*idx+0]<<" "<<triangles[3*idx+1]<<" "<<triangles[3*idx+2]<<" | ";
        //std::cout<<adj[3*idx+ 0]<<" "<<adj[3*idx+ 1]<<" "<<adj[3*idx+ 2]<<" | "<<std::endl;
        idx++;
    }
	delete Tr;

	
	auto t1 = std::chrono::high_resolution_clock::now();
	//std::cout<<tnumber<<" "<<pnumber<<" "<<std::endl;

	auto tb_label =std::chrono::high_resolution_clock::now();
	/* Etapa 1: Encontrar aristas máximas. */
	debug_msg("Etapa 1: Encontrar aristas máximas. \n");
	for(i = 0; i < tnumber; i++)
	{
		max[i] = max_edge_index(i, r, triangles);
	}
	//auto te_label =std::chrono::high_resolution_clock::now();

	//auto tb_desconect =std::chrono::high_resolution_clock::now();
	/* Etapa 2: Desconectar arcos asociados a aristas nomáx-nomáx. */
	debug_msg("Etapa 2: Desconectar arcos asociados a aristas nomáx-nomáx. \n");

	int num_terminal_edges =0;
	int num_terminal_border_edges=0;
	int num_frontier_edges=0;
	int num_frontier_border_edges=0;
	int num_interior_edges=0;


	for(i = 0; i < tnumber; i++){
		visited[i] = FALSE;
	}

	for(i = 0; i < tnumber; i++)
	{
		for(j = 0; j < 3; j++)
		{
			//Datos estadisticos
			//Borderedge and frontieredge
			if(adj[3*i +j] < 0){
				if (j == max[i])
					num_terminal_border_edges++;
				else
					num_frontier_border_edges++;
			}else
			{
				//If has frontieredge
				if(is_nomax_nomax(i, adj[3*i + j], triangles, max))
					num_frontier_edges++;
				//If has terminal_edge
				else if(is_max_max(i, adj[3*i + j], triangles, max))
					num_terminal_edges++;
				else //if is interioredge
					num_interior_edges++;
				//if(is_max_nomax(i, adj[3*i + j], triangles, max))
			}
			

			//Marcación real
			if(adj[3*i +j] < 0 || is_nomax_nomax(i, adj[3*i + j], triangles, max))
			{
				adj[3*i + j] = NO_ADJ;
			}


			if(adj[3*i +j] >= 0 && is_max_max(i, adj[3*i + j], triangles, max) && visited[adj[3*i + j]] == FALSE){
				visited[i] = TRUE;
				std::cout<<i<<" ";
			}

			//Puntos en el borde
			if(adj[3*i +j] < 0){
				border_point.push_back(triangles[3*i + (j+1 %3)]);
				border_point.push_back(triangles[3*i + (j+2 %3)]);
			}
		}
	}
	free(max);
	auto te_label =std::chrono::high_resolution_clock::now();
	 
	//debug_block(print_poly(root_id, tnumber););

	//medir tiempo
    

	debug_msg("Etapa 5: Generar poligonos\n");
	for(i = 0; i < tnumber; i++)
	{
		/*busca fronter edge en un triangulo, hacer función está wea*/

		/* si tiene 2-3 froint edge y no ha sido visitado*/
		//AGREGAR LA CONDICION DE QUE SI ROOT_ID = -1 ENTONCES NO GENERA POLY, MARCAR_ID[ALGO] == -1 DESPUES DE GENERAR C/POLY!!!
		int num_FrontierEdges = count_FrontierEdges(i, adj);
		//if(!visited[i] && num_FrontierEdges > 0){
		if(visited[i] == TRUE){

			chose_seed_triangle[id_chose_seed_triangle] = i;
			id_chose_seed_triangle++;

			length_poly = generate_polygon(poly, triangles, adj, r, visited, i);
			num_BE = count_BarrierEdges(poly, length_poly);
			
			//save_to_mesh(mesh, poly, &i_mesh, length_poly, pos_poly, &id_pos_poly);	
			
			if(num_BE>0){
				//printf("%d %d\n", num_BE, length_poly);
				est_total_be += num_BE;
				est_min_triangles_be = (num_BE < est_min_triangles_be) ? num_BE : est_min_triangles_be;
				est_max_triangles_be = (num_BE > est_max_triangles_be) ? num_BE : est_max_triangles_be;
				est_poly_with_be++;
				est_ratio_be += (float)num_BE/length_poly;

			}

			debug_msg("Poly: "); debug_block(print_poly(poly, length_poly););
			if( num_BE > 0){
				debug_print("Se dectecto %d BE\n", num_BE);
				auto tb_be = std::chrono::high_resolution_clock::now();
				remove_BarrierEdge(poly, length_poly, num_BE, triangles, adj, r, tnumber, mesh, &i_mesh, pos_poly, &id_pos_poly, visited);
				auto te_be = std::chrono::high_resolution_clock::now();
				tcost_be += std::chrono::duration_cast<std::chrono::milliseconds>( te_be - tb_be ).count();
			}else{
				debug_msg("Guardando poly\n");
				save_to_mesh(mesh, poly, &i_mesh, length_poly, pos_poly, &id_pos_poly);	
			}
		}
	}

	auto t2 = std::chrono::high_resolution_clock::now();

	
	

	
	write_geomview(r,triangles, pnumber, tnumber,i_mesh, mesh, id_pos_poly, pos_poly, print_triangles, ppath, chose_seed_triangle, id_chose_seed_triangle, border_point);

	int edges = mesh[0];
	int est_max_edges = edges;
	int est_min_edges = edges;
	int est_total_edges= edges;
	for(i = 1; i< id_pos_poly; i ++){
		edges = (pos_poly[i] - pos_poly[i-1]);
		est_max_edges = edges > est_max_edges ? edges : est_max_edges;
		est_min_edges = edges < est_min_edges ? edges : est_min_edges;
		est_total_edges += edges;
		
	}
	std::cout << std::setprecision(3) << std::fixed;
    std::cout <<"pnumber tnumber num_reg poly_be total_be min_poly_be max_poly_be ratio_be_per_poly total_edges max_edges min_edges edges_by_poly  num_terminal_edges num_terminal_border_edges num_frontier_edges num_frontier_border_edges num_interior_edges tdelaunay tlabel talgorithm ttraveandtopt ttravelalone topt"<<std::endl;
	std::cout<<pnumber<<" "<<tnumber<<" "<<id_pos_poly;
	std::cout<<" "<<est_poly_with_be<<" "<<est_total_be<<" "<<est_min_triangles_be<<" "<<est_max_triangles_be;
	std::cout<<" "<<(est_poly_with_be > 0 ? est_ratio_be/est_poly_with_be : 0.0);
	std::cout<<" "<<est_total_edges<<" "<<est_max_edges<<" "<<est_min_edges<<" "<<(float)est_total_edges/id_pos_poly;
	std::cout<<" "<<num_terminal_edges/2<<" "<<num_terminal_border_edges<<" "<<num_frontier_edges/2<<" "<<num_frontier_border_edges<<" "<<num_interior_edges/2;
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(te_delaunay - tb_delaunay).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(te_label - tb_label).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1 ).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1 ).count() - tcost_be;
	std::cout<<" "<<tcost_be<<std::endl;
	//std::cout<<3*pnumber - 3 - (num_terminal_border_edges + num_frontier_border_edges) <<" = "<<num_terminal_edges/2 + num_terminal_border_edges + num_frontier_edges/2 + num_frontier_border_edges + num_interior_edges/2<<" "<<(3*pnumber - 3 - (num_terminal_border_edges + num_frontier_border_edges) == num_terminal_edges/2 + num_terminal_border_edges + num_frontier_edges/2 + num_frontier_border_edges + num_interior_edges/2)<<std::endl;
	free(r);
	free(triangles);
	free(adj);
	free(visited );
	free(chose_seed_triangle);
	free(mesh);
	free(poly);
	free(pos_poly);
    
	return EXIT_SUCCESS;
}
    

