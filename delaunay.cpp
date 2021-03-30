#include <iostream>
#include <stdio.h>   
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include "detri2.h"

detri2::Triangulation *trimesh;
//Genenerates a Delaunay triangulation with detri2 from a random point set (file.node)
//Input: arguments Detri2 x pnumber xtnumber
//output: pnumber (point number Delunay triangulation), tnumber (number of triangles of triangulation)
void generate_delaunay_from_random_points(int argc, char* argv[], int &pnumber, int &tnumber){
    trimesh = new detri2::Triangulation();
    trimesh->parse_commands(argc, argv);
    trimesh->read_nodes();
    trimesh->incremental_delaunay();

    //dd
    tnumber = trimesh->tr_tris->objects - trimesh->ct_hullsize;
    pnumber = trimesh->ct_in_vrts + (trimesh->tr_steiners != NULL ? trimesh->tr_steiners->objects : 0);              
    if (!trimesh->io_keep_unused) { // no -IJ
        pnumber -= trimesh->ct_unused_vrts;
    }

}

//Fullfil arrays with delaunay triangulation data
//Input: r (array of points), triangles(array of triangles), adj (array of neigh)
void copy_delaunay_arrays(int tnumber, double *r, int* triangles, int* adj){
    int i, idx;
    // Llamada a detr2 
    //copiar arreglo de vertices
    //std::cout<<"pnumber "<<pnumber<<std::endl;
    idx = 0;
    for (i = 0; i < trimesh->ct_in_vrts; i++) {
        if (!trimesh->io_keep_unused) { // no -IJ
            if (trimesh->in_vrts[i].typ == UNUSEDVERTEX) continue;
        }
        r[2*i + 0]= trimesh->in_vrts[i].crd[0];
        r[2*i + 1]= trimesh->in_vrts[i].crd[1];
        //std::cout<<idx<<" ("<<r[2*i + 0]<<", "<<r[2*i + 1]<<") "<<std::endl;
        trimesh->in_vrts[i].idx = idx;
        idx++;
    }
    idx = 0;
    for (int i = 0; i < trimesh->tr_tris->used_items; i++) {
        detri2::Triang* tri = (detri2::Triang *) trimesh->tr_tris->get(i);
        if (tri->is_deleted()) continue;
        if (tri->is_hulltri()) {
            tri->idx = -1;
        } else {
            tri->idx = idx;
            idx++;
        }
    }

    //std::cout<<"tnumber: "<<trimesh->tr_tris->objects - trimesh->ct_hullsize<<std::endl;
    idx = 0;
    for (int i = 0; i < trimesh->tr_tris->used_items; i++)
    {
        
        detri2::Triang* tri = (detri2::Triang *) trimesh->tr_tris->get(i);
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
    
 
	//for(i = 0; i < tnumber; i++){
	//
	//	for(int j = 0; j<3;j++){
	//		std::cout<<triangles[3*i+j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}

    delete trimesh;   
}