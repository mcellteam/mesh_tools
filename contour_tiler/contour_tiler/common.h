#ifndef T_COMMON
#define T_COMMON
#include "ct/ct.h" 

#ifdef i860
    #include <nx.h>
#endif

#define PLANE_VERT_DIM 2  /* no normals */ 
#define SPACE_VERT_DIM 3  /* no normals */ 
#define MAX_GROUP_NUM 400
#define MAX_LINES_NUM 4000
#define MAX_TRI_NUM 20000

typedef struct _Linklist_
{
   int index,ary_ind;
   double *ratio_ary;
   struct _Linklist_ *next;
}  Linklist;



typedef double VertType[SPACE_VERT_DIM]; /* x,y */
typedef double Normal[SPACE_VERT_DIM]; /* x,y */

typedef struct
{
   int    numpts;                /* the number of vertices                 */
   VertType *vertices;             /* the  vertices                        */
   Normal *normals;              /* the normals at the vertices            */
}  PolygonStruct;

typedef struct 
{
    double a[2],b[2],c[2]; /* two legal equations ax+by+c=0 */
    short group; /* the best tiling group */
    short ind;   /* the best tiling index */
    float dist;
    char  used;
    char  and_or; /* 0: or, 1: and */
    unsigned char  nec_polarity, nec_group;
} TileType;
    
typedef struct 
{
    short group0,group1;
    short ind0,ind1;
    char tile_dir;
} BoundaryType;


typedef struct {
    float	thres;
    int		xcube_vox,ycube_vox,zcube_vox;
    float	detz;
    int		type, beg,end,method;
    int		x1,x2,y1,y2;
    float	epsilon;
} Interpo_struct;

typedef struct Link_struct{
    short	box[6];
    struct	Link_struct *next;
} Link_struct;

/* separate pointer *prefix, *suffix, *iso_surface_name from Interpo_struct
   into another structure. So it is easier to do message passing of 
   Interpo_struct */

typedef struct {
    char	*prefix,*suffix;
    char	*output_file_name;
} Name_struct;

#endif
