int count_regions(int *mesh, int i_mesh){
    int i = 0;
	int num_region = 0;
    while(i < i_mesh){
		num_region++;
		i = mesh[i] + i + 1;
    }
    return num_region;
}

int save_to_mesh(int *mesh, int *poly, int i_mesh, int length_poly){    

    mesh[i_mesh] = length_poly;
    
    i_mesh++;
    for(int i = 0; i <length_poly; i++){
        mesh[i_mesh + i] = poly[i];
    }
    
    return i_mesh + length_poly;
}