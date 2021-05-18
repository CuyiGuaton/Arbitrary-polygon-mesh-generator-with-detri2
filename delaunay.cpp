#include <iostream>
#include <stdio.h>   
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include "detri2.h"
#include <assert.h> 


detri2::Triangulation *trimesh;

void generate_mesh(int argc, char* argv[], detri2::Triangulation *Tr)
{
  if (argc < 2) {
    printf("Usage: detri2 [-options] filename[.node, .ele]\n");
    exit(0);    
  }

   //Tr = new detri2::Triangulation();

  // Read options.
  if (!Tr->parse_commands(argc, argv)) {
    // No input or wrong parameters.
    printf("Usage: detri2 [-options] filename[.node, .poly, .ele, .edge]\n");
    delete Tr;
    exit(0);
  }

  // Read inputs.
  if (!Tr->read_mesh()) {
    printf("Failed to read input from file %s[.poly, .node, .ele, .edge]\n",
           Tr->io_infilename);
    delete Tr;
    exit(0);
  }

  // Generate (constrained) (weighted) Delaunay detri2::Triangulation.
  if (Tr->tr_tris == NULL) {
    if (Tr->incremental_delaunay()) {
      if (Tr->tr_segs != NULL) {
        Tr->recover_segments();
        Tr->set_subdomains();
      }
    } else {
      printf("Failed to create Delaunay (regular) triangulation.\n");
      delete Tr;
      exit(0);
    }
  } else {
    Tr->reconstruct_mesh(0);
  }

  // Mesh refinement and adaptation.
  if (Tr->op_quality) { // -q
    if (Tr->io_omtfilename[0] != '\0') {
      // A background mesh is supplied.
      Tr->OMT_domain = new detri2::Triangulation();
      int myargc = 2;
      char *myargv[2];
      myargv[0] = argv[0];
      myargv[1] = Tr->io_omtfilename;
      Tr->OMT_domain->parse_commands(myargc, myargv);
      Tr->OMT_domain->read_mesh();
      // Check if this background mesh contains metrics
      if (Tr->OMT_domain->io_with_metric) {
        Tr->OMT_domain->reconstruct_mesh(0);
        Tr->set_vertex_metrics();
      } else {
        printf("Warning: Background mesh contains no metric, ignored.\n");
        delete Tr->OMT_domain;
        Tr->OMT_domain = NULL;
      }
    } else {
      assert(Tr->OMT_domain == NULL);
    }
  
    Tr->op_metric = METRIC_Euclidean_no_weight;
    
    if ((Tr->io_with_metric) ||
        (Tr->OMT_domain != NULL) ||     // -m filename
        (Tr->op_target_length > 0.)) {  // -ML=#
      Tr->coarsen_mesh();
    }
  
    Tr->delaunay_refinement();
    
    if (Tr->op_use_smoothing) { // -qS1 (Laplacian) -qS2 ()CVT
      for (int i = 0; i < Tr->op_smooth_iter; i++) { // -qI#
        Tr->smooth_vertices();
      }
    }
  } 

  // Mesh export (to files).
  if (Tr->tr_tris != NULL) {
    if ((Tr->ct_exteriors > 0) && !Tr->op_convex) { // no -c
      Tr->remove_exteriors();
    }
    Tr->save_triangulation();
    if (Tr->io_outedges) {
      Tr->save_edges();
    }
    if (Tr->io_out_ucd) { // -Iu
      Tr->save_to_ucd(0, 0);
    }
    if (Tr->io_out_voronoi) { // -Iv
      Tr->save_voronoi(Tr->io_out_ucd);
    }
  }
  

  Tr->mesh_statistics();

  //delete Tr;
  //return Tr;
}


//Genenerates a Delaunay triangulation with detri2 from a random point set (file.node)
//Input: arguments Detri2 x pnumber xtnumber
//output: pnumber (point number Delunay triangulation), tnumber (number of triangles of triangulation)
void generate_delaunay_from_random_points(int argc, char* argv[], int &pnumber, int &tnumber){
    trimesh = new detri2::Triangulation();
    //trimesh->parse_commands(argc, argv);
    //trimesh->read_nodes();
    //trimesh->incremental_delaunay();
    generate_mesh(argc, argv, trimesh);
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
        //Se le asigna a cada vértice un indice
        trimesh->in_vrts[i].idx = idx;
        idx++;
    }
    //Se le asigna a cada triangulo un indice, si es del hull es -1
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
    //Se guardan los triangulos y adj de detri2
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

/*
detri2::Triangulation *generate_mesh(int argc, char* argv[], detri2::Triangulation *Tr)
{
  if (argc < 2) {
    printf("Usage: detri2 [-options] filename[.node, .ele]\n");
    exit(0);    
  }

   Tr = new detri2::Triangulation();

  // Read options.
  if (!Tr->parse_commands(argc, argv)) {
    // No input or wrong parameters.
    printf("Usage: detri2 [-options] filename[.node, .poly, .ele, .edge]\n");
    delete Tr;
    exit(0);
  }

  // Read inputs.
  if (!Tr->read_mesh()) {
    printf("Failed to read input from file %s[.poly, .node, .ele, .edge]\n",
           Tr->io_infilename);
    delete Tr;
    exit(0);
  }

  // Generate (constrained) (weighted) Delaunay detri2::Triangulation.
  if (Tr->tr_tris == NULL) {
    if (Tr->incremental_delaunay()) {
      if (Tr->tr_segs != NULL) {
        Tr->recover_segments();
        Tr->set_subdomains();
      }
    } else {
      printf("Failed to create Delaunay (regular) triangulation.\n");
      delete Tr;
      exit(0);
    }
  } else {
    Tr->reconstruct_mesh(0);
  }

  // Mesh refinement and adaptation.
  if (Tr->op_quality) { // -q
    if (Tr->io_omtfilename[0] != '\0') {
      // A background mesh is supplied.
      Tr->OMT_domain = new detri2::Triangulation();
      int myargc = 2;
      char *myargv[2];
      myargv[0] = argv[0];
      myargv[1] = Tr->io_omtfilename;
      Tr->OMT_domain->parse_commands(myargc, myargv);
      Tr->OMT_domain->read_mesh();
      // Check if this background mesh contains metrics
      if (Tr->OMT_domain->io_with_metric) {
        Tr->OMT_domain->reconstruct_mesh(0);
        Tr->set_vertex_metrics();
      } else {
        printf("Warning: Background mesh contains no metric, ignored.\n");
        delete Tr->OMT_domain;
        Tr->OMT_domain = NULL;
      }
    } else {
      assert(Tr->OMT_domain == NULL);
    }
  
    Tr->op_metric = METRIC_Euclidean_no_weight;
    
    if ((Tr->io_with_metric) ||
        (Tr->OMT_domain != NULL) ||     // -m filename
        (Tr->op_target_length > 0.)) {  // -ML=#
      Tr->coarsen_mesh();
    }
  
    Tr->delaunay_refinement();
    
    if (Tr->op_use_smoothing) { // -qS1 (Laplacian) -qS2 ()CVT
      for (int i = 0; i < Tr->op_smooth_iter; i++) { // -qI#
        Tr->smooth_vertices();
      }
    }
  }

  // Mesh export (to files).
  if (Tr->tr_tris != NULL) {
    if ((Tr->ct_exteriors > 0) && !Tr->op_convex) { // no -c
      Tr->remove_exteriors();
    }
    Tr->save_triangulation();
    if (Tr->io_outedges) {
      Tr->save_edges();
    }
    if (Tr->io_out_ucd) { // -Iu
      Tr->save_to_ucd(0, 0);
    }
    if (Tr->io_out_voronoi) { // -Iv
      Tr->save_voronoi(Tr->io_out_ucd);
    }
  }

  Tr->mesh_statistics();

  //delete Tr;
  return Tr;
}*/


