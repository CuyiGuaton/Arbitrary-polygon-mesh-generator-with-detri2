/* Prototipos de funciones para manejar entrada y salida de datos. */

void write_geomview(double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, int print_triangles);
void write_VEM(double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, int print_triangles);
void write_VEM_triangles(double *r, int *triangles, int *adj, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, std::list <int> &seed_bet);

int look_triangles(int i, std::set<int> &s, int * triangles, int * adj, double *r);