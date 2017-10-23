
#ifndef __EIGEN__
#define __EIGEN__

#ifdef __cplusplus
extern "C" {
#endif

  long eigenCardan(int nbre, double sxx, double syy, double szz, double sxy, double sxz, double syz,
    double * vp, double * ux, double * uy, double * uz);

  void eigen_values_sdp3x3(double sxx, double sxy, double sxz, double syy, double syz, double szz, double *lambda);

  long eigenCardan2(int nbre, double sxx, double syy, double szz, double sxy, double sxz, double syz,
    double * vp, double * ux, double * uy, double * uz);


  long EqDeg2(double g_xx, double g_xy, double g_yy,
    long indice_vp, double * vp_x, double * vp_y, double * vp1, double * vp2);

#ifdef __cplusplus
}
#endif


#endif