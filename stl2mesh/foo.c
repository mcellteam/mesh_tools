

/* determine if a 2-D point lies within the bounds of a 2-D clipping window */
int clip_point_2D(umin,umax,vmin,vmax,u1,v1)
  double umin,umax,vmin,vmax,u1,v1;
{

  if (u1>=umin) {
    if (u1<=umax) {
      if (v1>=vmin) {
        if (v1<=vmax) {
          return(1);
        }
      }
    }
  }
  return(0);
}


/* Find which pair of partitions a point lies between. */
/*   Return index of lower partition */
/*   or return -index if point is coincident with a partition. */ 
int find_range(u,u_range,n_u_range)
double u;
double *u_range;
int n_u_range;
{
  int n_lower,n_upper,n_mid;

  n_lower=0;
  n_upper=n_u_range-1;
  while (n_lower!=(n_upper-1)) {
    n_mid=(n_lower+n_upper)/2;
    if (u>u_range[n_mid]) {
      n_lower=n_mid;
    }
    else if (u<u_range[n_mid]) {
      n_upper=n_mid;
    }
    else {
      return(-n_mid);
    }
  }
  if (u==u_range[n_lower]) {
    return(-n_lower);
  }
  if (u==u_range[n_lower+1]) {
    return(-(n_lower+1));
  }

  return(n_lower);
}

