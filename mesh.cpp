int count_regions(int *mesh, int i_mesh){
    int i = 0;
	int num_region = 0;
    while(i < i_mesh){
		num_region++;
		i = mesh[i] + i + 1;
    }
    return num_region;
}