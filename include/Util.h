#include <stdio.h>
#include "Point.h"
#include "Vector.h"

#ifdef LINUXVERSION
#include <linux/kernel.h>
#include <sys/sysinfo.h>
#endif

#ifdef MACVERSION
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#include <sys/param.h>

#ifndef Util_h
#define Util_h

void sort_lt(int n, int list[], int f[]);
void sort_lt(int n, int list[], double f[]);
void sort_gt(int n, int list[], int f[]);
void sort_gt(int n, int list[], double f[]);
const double identity[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

int invert3X3(double m[3][3]);
double circumcenter(Point p1, Point p2, Point p3, Point p4, Point &cc);
double closest_to_edge(Point p1, Point p2, Point p, Point &ep);
double closest_to_facet(Point p1, Point p2, Point p3, Point p, Point &ep);
double distance( Point p1, Point p2);
int line_line_intersect(Point pa, Point pb, Point p1, Point p2, Point &p);
int line_facet_intersect(Point pa, Point pb, Point p1, Point p2, Point p3, Point &p,
                         double eps);
int facet_facet_intersect(Point p1, Point p2, Point p3, Point q1, Point q2, Point q3, Point &pa, Point &pb, double tol);
int tri_quad_intersect(Point t1, Point t2, Point t3, Point q1, Point q2, Point q3, Point q4, Point &pa, Point &pb, double tol);
double tetrahedral_aspect_ratio(Point p1, Point p2, Point p3, Point p4);
double tetrahedral_aspect_ratio(Point p1, Point p2, Point p3, Point p4, Point &ci, Point &cc);
double tetrahedral_volume(Point p1, Point p2, Point p3, Point p4);
double triangle_area(Point p1, Point p2, Point p3);
double triangle_aspect_ratio(Point p1, Point p2, Point p3, Point &ci, Point &cc);
double quadrilateral_area(Point p1, Point p2, Point p3, Point p4);
double quadrilateral_aspect_ratio(Point p1, Point p2, Point p3, Point p4);
double pyramid_volume(Point p1, Point p2, Point p3, Point p4, Point p5);
double pyramid_aspect_ratio(Point p1, Point p2, Point p3, Point p4, Point p5);
double prism_volume(Point p1, Point p2, Point p3, Point p4, Point p5, Point p6);
double prism_aspect_ratio(Point p1, Point p2, Point p3, Point p4, Point p5, Point p6);
double hexahedral_volume(Point p1, Point p2, Point p3, Point p4,
                         Point p5, Point p6, Point p7, Point p8);
double hexahedral_aspect_ratio(Point p1, Point p2, Point p3, Point p4,
                               Point p5, Point p6, Point p7, Point p8);
double condition_number(Point p0, Point p1, Point p2, Point p3, double w[3][3]);
double tetrahedral_condition_number(Point p0, Point p1, Point p2, Point p3, int wgt);
Vector tetrahedral_gradient(Point p1, Point p2, Point p3, Point p4,
                            double f1, double f2, double f3, double f4);
Vector pyramid_gradient(Point p1, Point p2, Point p3, Point p4, Point p5,
                        double f1, double f2, double f3, double f4, double f5);
Vector prism_gradient(Point p1, Point p2, Point p3,
                      Point p4, Point p5, Point p6,
                      double f1, double f2, double f3,
                      double f4, double f5, double f6);
Vector hexahedral_gradient(Point p1, Point p2, Point p3, Point p4,
                           Point p5, Point p6, Point p7, Point p8,
                           double f1, double f2, double f3, double f4,
                           double f5, double f6, double f7, double f8);
double solid_angle(Vector v1, Vector v2, Vector v3);
double solid_angle(double x0, double y0, double z0, double x1, double y1, double z1,
                   double x2, double y2, double z2, double x3, double y3, double z3);
void histogram(int nh, double hist[]);
double pythag(double a, double b);
int svdcmp(double **a, int m, int n, double w[], double **v);
void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
int compute_Riemannian_metric(Point p1, Point p2, Point p3, Point p4, double RT[3][3]);

template <class T>
void my_mem(int &my_dim, int &n_objs, T **obj, int my_size, int my_chunk)
{
  int new_dim;
  
  if (my_size > 0)
  {
    if (my_size < my_dim)
      new_dim = my_size;
    else
    {
      new_dim = my_dim;
      while (new_dim < my_size)
        new_dim += my_chunk;
    }
    
#ifdef LINUXVERSION
    struct sysinfo si;
    //call sysinfo
    sysinfo (&si);
    if (si.totalram < new_dim*sizeof(T))
    {
      fprintf(stderr,"\nYou are requesting more memory (%ld bytes) than the system has available (%ld bytes).", new_dim*sizeof(T), si.totalram);
      fprintf(stderr,"\nConsider making a coarser mesh, increasing spawn level to ease load balancing, adding more procs, or turning off the layers flag in P_HUGG_geometry.params and increasing g_space.");
      fflush(stderr);
    }
#endif

#ifdef MACVERSION
    int mib[2], usermem;
    size_t len;
    mib[0] = CTL_HW;
    mib[1] = HW_USERMEM;
    len = sizeof(u_int64_t);
    sysctl(mib, 2, &usermem, &len, NULL, 0);

    if (usermem < new_dim*sizeof(T))
    {
      fprintf(stderr,"\nYou are requesting more memory (%ld bytes) than the system has available (%ld bytes).", (long)new_dim*sizeof(T), (long)usermem);
      fprintf(stderr,"\nConsider making a coarser mesh, increasing spawn level to ease load balancing, adding more procs, or turning off the layers flag in P_HUGG_geometry.params and increasing g_space.");
      fflush(stderr);
    }
#endif

    if (*obj != 0)
      *obj=(T*)realloc((void*)*obj,new_dim*sizeof(T));
    else
      *obj=(T*)malloc(new_dim*sizeof(T));
    if (*obj == 0)
    {
      fprintf(stderr,"\nmy_mem(): Unable to allocate memory for object!");
      fprintf(stderr,"\n    Current dimension -> %d",my_dim);
      fprintf(stderr,"\n  Requested dimension -> %d",my_size);
#ifndef _Xcode_DEBUG
      //MPI_Abort(MPI_COMM_WORLD,0);
#endif
      exit(0);
    }
    my_dim = new_dim;
  } else
  {
    if (*obj != 0)
      free(*obj);
    *obj = 0;
    n_objs = my_dim = 0;
  }
}

#endif
