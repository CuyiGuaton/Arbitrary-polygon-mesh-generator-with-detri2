#include <iostream>
#include <stdio.h>   
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <vector> 
#include <chrono>
#include <iomanip>
#include <cstdlib>


#include "delaunay.h"
#include "io.h"
#include "consts.h"
#include "triang.h"
#include "polygon.h"
#include "mesh.h"
#include "BET_elimitation.h"

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 1
#endif

#define debug_block(fmt) do { if (DEBUG_TEST){ fmt }} while (0)
#define debug_print(fmt, ...) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)
#define debug_msg(fmt) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__,  __LINE__, __func__); } while (0)

//make &&  python3 rp.py 100 0 0 && ./DelaunayPolyGenerator autodata.node && geomview output/.off
//make &&  python3 rp.py 506 0 0 && ./DelaunayPolyGenerator -p -q unicorn.poly && geomview output/.off


int main(int argc, char* argv[]){

	char* ppath;
	int print_triangles = 0;

	// int nparam = 3;
    //char* params[] = {const_cast<char*> ("./detri2"), const_cast<char*> ("-z"), const_cast<char*> ("test.node")};
	//char* params[] = {const_cast<char*> ("./detri2"), const_cast<char*> ("-z"), const_cast<char*> ("506randompoints.node")};
	//char* params[] = {const_cast<char*> ("./detri2"), const_cast<char*> ("-z"), const_cast<char*> ("506equilateral.node")};
    //char* ppath = const_cast<char*> ("test");

	//number elements
	int tnumber, pnumber, i,j;
	

	//arrays inicilization
	double *r;
	int *triangles;
	int *adj;
    int *seed;
    int *max;
	int *mesh;
	int *trivertex;
	
	auto tb_delaunay = std::chrono::high_resolution_clock::now();
	generate_delaunay_from_random_points(argc, argv, pnumber,tnumber);
	auto te_delaunay = std::chrono::high_resolution_clock::now();

	r = (double *)malloc(2*pnumber*sizeof(double)); // cambiar por pnumber
    triangles = (int *)malloc(3*tnumber*sizeof(int));
	adj = (int *)malloc(3*tnumber*sizeof(int));
	seed = (int *)malloc(tnumber*sizeof(int));
    max = (int *)malloc(tnumber*sizeof(int));
	mesh = (int *)malloc(3*tnumber*sizeof(int));
	trivertex = (int *)malloc(pnumber*sizeof(int));

	copy_delaunay_arrays(tnumber, r, triangles, adj);
	
	//Asociate each vertex to an adjacent  triangle
	for(i = 0; i < pnumber; i++){
		for (j = 0; j < tnumber; j++)
		{
			if(i == triangles[3*j + 0] ||  i == triangles[3*j + 1] || i == triangles[3*j + 2]){
				trivertex[i] = j;
				break;
			}
		}
		//std::cout<<"trivertex["<<i<<"] "<<trivertex[i]<<std::endl;
	}

	//stats
	int i_mesh = 0;	
	int num_BE = 0;
	int est_total_be = 0;
	int est_min_triangles_be = 2147483641;
	int est_max_triangles_be = 0;
	int est_poly_with_be = 0;
	double est_ratio_be = 0;
	int tcost_be = 0;
	

	//initialize array
	for(i = 0; i < tnumber; i++){
		seed[i] = FALSE;
	}
	
	auto t1 = std::chrono::high_resolution_clock::now();


	auto tb_label =std::chrono::high_resolution_clock::now();
	/* Etapa 1: Encontrar aristas máximas. */
	debug_msg("Etapa 1: Encontrar aristas máximas. \n");
	for(i = 0; i < tnumber; i++)
	{
		max[i] = max_edge_index(i, r, triangles);
	}	
	
	//Marcar triangulos semilla
	for(i = 0; i < tnumber; i++){
		for(j = 0; j < 3; j++){
			if(adj[3*i +j] != -1 && is_max_max(i, adj[3*i + j], triangles, max) == TRUE ){
				if(adj[3*i + j] < i){
					seed[i] = TRUE;
					break;
				}
			}
			if(adj[3*i +j] == -1  && max[i] == (j +1)%3){
				seed[i] = TRUE;
				break;
			}
		}
	}

	/* Etapa 2: Desconectar arcos asociados a aristas nomáx-nomáx. */
	debug_msg("Etapa 2: Desconectar arcos asociados a aristas nomáx-nomáx. \n");
	
	int num_terminal_edges =0;
	int num_terminal_border_edges=0;
	int num_frontier_edges=0;
	int num_frontier_border_edges=0;
	int num_interior_edges=0;
		
	for(i = 0; i < tnumber; i++)
	{
		for(j = 0; j < 3; j++)
		{
			//Datos estadisticos, hay estadisticas que están mal, verificar
			//Borderedge and frontieredge, siento que esta wea está mala
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
			
		}
	}
	auto te_label =std::chrono::high_resolution_clock::now();
	
	free(max);
	
	
    int length_poly;

	
	debug_msg("Etapa 5: Generar poligonos\n");
	int poly[1000];	
	auto tb_travel = std::chrono::high_resolution_clock::now();
	for(i = 0; i < tnumber; i++)
	{
		if(seed[i] == TRUE ){			

			length_poly = generate_polygon(i, poly, triangles, adj, r);
			num_BE = count_BarrierEdges(poly, length_poly);
			
		
			i_mesh = save_to_mesh(mesh, poly, i_mesh, length_poly);	
					
			if(num_BE>0){
				//i_mesh = save_to_mesh(mesh, poly, i_mesh, length_poly);
				//printf("%d %d\n", num_BE, length_poly);
				est_total_be += num_BE;
				est_min_triangles_be = (num_BE < est_min_triangles_be) ? num_BE : est_min_triangles_be;
				est_max_triangles_be = (num_BE > est_max_triangles_be) ? num_BE : est_max_triangles_be;
				est_poly_with_be++;
				est_ratio_be += (float)num_BE/length_poly;
			}
			
		
			debug_msg("Poly: "); debug_block(print_poly(poly, length_poly););
			if( num_BE > 0){
				//printf("Se dectecto %d BE\n", num_BE);
				auto tb_be = std::chrono::high_resolution_clock::now();
				i_mesh = Remove_BE2(1,poly, length_poly, num_BE, triangles, adj, r, tnumber, mesh, i_mesh, trivertex);
				auto te_be = std::chrono::high_resolution_clock::now();
				tcost_be += std::chrono::duration_cast<std::chrono::milliseconds>( te_be - tb_be ).count();
				
			}else{
				debug_msg("Guardando poly\n");
				//i_mesh = save_to_mesh(mesh, poly, i_mesh, length_poly);	
				
			}
		
		}
	}
	auto te_travel = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();

	int num_region = count_regions(mesh,i_mesh);

	write_geomview(r, triangles, pnumber, tnumber, i_mesh, mesh, seed, num_region, print_triangles);

/*
	int edges = pos_poly[0];
	int est_max_edges = edges;
	int est_min_edges = edges;
	int est_total_edges= edges;
	for(i = 1; i< id_pos_poly; i ++){
		edges = (pos_poly[i] - pos_poly[i-1]);
		est_max_edges = edges > est_max_edges ? edges : est_max_edges;
		est_min_edges = edges < est_min_edges ? edges : est_min_edges;
		est_total_edges += edges;
	}*/

	int edges, est_max_edges, est_min_edges,est_total_edges;
	edges = mesh[0];
	est_max_edges = edges;
	est_min_edges = edges;
	est_total_edges= edges;
	i = 0;
	while(i < i_mesh){
        edges = mesh[i];
		est_max_edges = edges > est_max_edges ? edges : est_max_edges;
		est_min_edges = edges < est_min_edges ? edges : est_min_edges;
		est_total_edges += edges;
        i++;
        for(j=0; j < edges;j++){            
			i++;
        }
    }

	std::cout << std::setprecision(3) << std::fixed;
    std::cout <<"pnumber tnumber num_reg poly_with_be total_be min_poly_be max_poly_be ratio_be_per_poly total_edges max_edges min_edges edges_by_poly num_terminal_edges num_terminal_border_edges num_frontier_edges num_frontier_border_edges num_interior_edges tdelaunay tlabel talgorithm ttraveandtopt ttravelalone topt"<<std::endl;
	std::cout<<pnumber<<" "<<tnumber<<" "<<num_region;
	std::cout<<" "<<est_poly_with_be<<" "<<est_total_be<<" "<<est_min_triangles_be<<" "<<est_max_triangles_be;
	std::cout<<" "<<(est_poly_with_be > 0 ? est_ratio_be/est_poly_with_be : 0.0);
	std::cout<<" "<<est_total_edges<<" "<<est_max_edges<<" "<<est_min_edges<<" "<<(float)est_total_edges/num_region;
	std::cout<<" "<<num_terminal_edges/2<<" "<<num_terminal_border_edges<<" "<<num_frontier_edges/2<<" "<<num_frontier_border_edges<<" "<<num_interior_edges/2;
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(te_delaunay - tb_delaunay).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(te_label - tb_label).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1 ).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(te_travel - tb_travel).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(te_travel - tb_travel).count() - tcost_be;
	std::cout<<" "<<tcost_be<<std::endl;
	std::cout<<3*pnumber - 3 - (num_terminal_border_edges + num_frontier_border_edges) <<" = "<<num_terminal_edges/2 + num_terminal_border_edges + num_frontier_edges/2 + num_frontier_border_edges + num_interior_edges/2<<" "<<(3*pnumber - 3 - (num_terminal_border_edges + num_frontier_border_edges) == num_terminal_edges/2 + num_terminal_border_edges + num_frontier_edges/2 + num_frontier_border_edges + num_interior_edges/2)<<std::endl;
	free(r);
	free(triangles);
	free(adj);
	free(seed);
	free(mesh);    
	return EXIT_SUCCESS;
}
    

