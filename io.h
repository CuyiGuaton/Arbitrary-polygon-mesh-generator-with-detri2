/* Prototipos de funciones para manejar entrada y salida de datos. */

void read_qdelaunay_data(char *ppath, double **r, int **p, int **ady, int *pnumber, int *tnumuber, void *align_settings);
void read_triangles(char *cpath, int **id, int *inumber);
void write_geomview(double *r, int *p,  int pnumber,
									int tnumber, int ind_poly, 
									int *polygon, int id_start_poly,int *start_polygon, int print_triangles, char *ppath);
void read_fromfiles_data(char *ppath, double **r, int **p, int **adj, int *pnumber, int *tnumber, void **align_settings);