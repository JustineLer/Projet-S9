
#include <math.h>
#include <malloc.h>

#include <util.h>
#include <eigen.h>


#define INV3  .333333333333333f
#define ROOT3 1.732050807568877f

static long eqDeg3_Cardan_Real(double a, double b, double c, double d, double * root_real)
{
  double p, q, delta, theta, norm;

  if (a == 0.0) return 0;
  p = c / a - b*b / (3.0*a*a);
  if (p > 0.0) return 0;
  q = d / a - b*c / (3.0*a*a) + 2.0*b*b*b / (27.0*a*a*a);
  delta = q*q + 4.0*p*p*p / 27.0;
  if (delta > 0.0) delta = 0.0;
  theta = atan2(sqrt(-delta), -q);
  norm = 2.0*sqrt(-p / 3.0);
  root_real[0] = norm*cos(theta / 3.0) + b / 3.0;
  root_real[1] = norm*cos((theta + PI_DOUBLE) / 3.0) + b / 3.0;
  root_real[2] = norm*cos((theta + 2.0*PI_DOUBLE) / 3.0) + b / 3.0;
  return 1;
}


long eigenCardan(int nbre, double sxx, double syy, double szz, double sxy, double sxz, double syz,
  double * vp, double * ux, double * uy, double * uz)
{
  long n, no;
  double norm, det[3], eig[3], lambda, x, y, z, a = -1.0, b, c, d;

  b = sxx + syy + szz;
  c = sxy*sxy + sxz*sxz + syz*syz - sxx*syy - sxx*szz - syy*szz;
  d = sxx*syy*szz - sxx*syz*syz - syy*sxz*sxz - szz*sxy*sxy + 2.0*sxy*sxz*syz;
  if (!eqDeg3_Cardan_Real(a, b, c, d, eig)) return 0;
  if (eig[0] < eig[1]) { lambda = eig[1]; eig[1] = eig[0]; eig[0] = lambda; }
  if (eig[0] < eig[2]) { lambda = eig[2]; eig[2] = eig[0]; eig[0] = lambda; }
  if (eig[1] < eig[2]) { lambda = eig[2]; eig[2] = eig[1]; eig[1] = lambda; }
  if (vp != NULL)
    for (n = 0; n<3; n++) vp[n] = eig[n];
  if (ux != NULL && uy != NULL && uz != NULL)
    for (n = 0; n<nbre; n++)
    {
      lambda = eig[n];
      det[0] = (sxx - lambda)*(syy - lambda) - sxy*sxy;
      det[1] = (sxx - lambda)*(szz - lambda) - sxz*sxz;
      det[2] = (syy - lambda)*(szz - lambda) - syz*syz;
      if (fabs(det[0]) > fabs(det[1]))
      {
        if (fabs(det[0]) > fabs(det[2])) no = 0;
        else if (fabs(det[1]) > fabs(det[2])) no = 1;
        else no = 2;
      }
      else
      {
        if (fabs(det[1]) > fabs(det[2])) no = 1; else no = 2;
      }
      if (det[no] == 0.0)
        x = y = z = 0.0;
      else if (no == 0)
      {
        x = ((lambda - syy)*sxz + sxy*syz) / det[no];
        y = ((lambda - sxx)*syz + sxy*sxz) / det[no];
        z = 1.0;
      }
      else if (no == 1)
      {
        x = ((lambda - szz)*sxy + sxz*syz) / det[no];
        y = 1.0;
        z = ((lambda - sxx)*syz + sxz*sxy) / det[no];
      }
      else if (no == 2)
      {
        x = 1.0;
        y = ((lambda - szz)*sxy + syz*sxz) / det[no];
        z = ((lambda - syy)*sxz + syz*sxy) / det[no];
      }
      norm = sqrt(x*x + y*y + z*z);
      if (norm != 0.0)
      {
        x /= norm;
        y /= norm;
        z /= norm;
      }
      if (y > 0.0)
      {
        x = -x;
        y = -y;
        z = -z;
      }
      ux[n] = x;
      uy[n] = y;
      uz[n] = z;
    }
  return 1;
}



void eigen_values_sdp3x3(double sxx, double sxy, double sxz, double syy, double syz, double szz, double *lambda)
{
  double c0, c1, c2, c2_div3, a_div3, mb_div2, q, magnitude, angle, inv3 = INV3, root3 = ROOT3, cs, sn;

  c0 = sxx*syy*szz + 2.0*sxy*sxz*syz - sxx*syz*syz - syy*sxz*sxz - szz*sxy*sxy;
  c1 = sxx*syy - sxy*sxy + sxx*szz - sxz*sxz + syy*szz - syz*syz;
  c2 = sxx + syy + szz;
  c2_div3 = c2 * inv3;
  a_div3 = (c1 - c2*c2_div3)*inv3;
  if (a_div3 > 0.0) a_div3 = 0.0;
  mb_div2 = .5*(c0 + c2_div3*(2.0*c2_div3*c2_div3 - c1));
  q = mb_div2 * mb_div2 + a_div3*a_div3*a_div3;
  if (q > 0.0) q = 0.0;
  magnitude = sqrt(-a_div3);
  angle = atan2(sqrt(-q), mb_div2)*inv3;
  cs = cos(angle);
  sn = sin(angle);
  lambda[0] = c2_div3 + 2.0*magnitude*cs;
  lambda[1] = c2_div3 - magnitude*(cs + root3*sn);
  lambda[2] = c2_div3 - magnitude*(cs - root3*sn);
}


long eigenCardan2(int nbre, double sxx, double syy, double szz, double sxy, double sxz, double syz,
  double * vp, double * ux, double * uy, double * uz)
{
  long n, no;
  double norm, det[3], eig[3], lambda, x, y, z;

  eigen_values_sdp3x3(sxx, sxy, sxz, syy, syz, szz, eig);
  if (eig[0] < eig[1]) { lambda = eig[1]; eig[1] = eig[0]; eig[0] = lambda; }
  if (eig[0] < eig[2]) { lambda = eig[2]; eig[2] = eig[0]; eig[0] = lambda; }
  if (eig[1] < eig[2]) { lambda = eig[2]; eig[2] = eig[1]; eig[1] = lambda; }
  if (vp != NULL)
    for (n = 0; n<3; n++) vp[n] = eig[n];
  if (ux != NULL && uy != NULL && uz != NULL)
    for (n = 0; n<nbre; n++)
    {
      lambda = eig[n];
      det[0] = (sxx - lambda)*(syy - lambda) - sxy*sxy;
      det[1] = (sxx - lambda)*(szz - lambda) - sxz*sxz;
      det[2] = (syy - lambda)*(szz - lambda) - syz*syz;
      if (fabs(det[0]) > fabs(det[1]))
      {
        if (fabs(det[0]) > fabs(det[2])) no = 0;
        else if (fabs(det[1]) > fabs(det[2])) no = 1;
        else no = 2;
      }
      else
      {
        if (fabs(det[1]) > fabs(det[2])) no = 1; else no = 2;
      }
      if (det[no] == 0.0)
        x = y = z = 0.0;
      else if (no == 0)
      {
        x = ((lambda - syy)*sxz + sxy*syz) / det[no];
        y = ((lambda - sxx)*syz + sxy*sxz) / det[no];
        z = 1.0;
      }
      else if (no == 1)
      {
        x = ((lambda - szz)*sxy + sxz*syz) / det[no];
        y = 1.0;
        z = ((lambda - sxx)*syz + sxz*sxy) / det[no];
      }
      else if (no == 2)
      {
        x = 1.0;
        y = ((lambda - szz)*sxy + syz*sxz) / det[no];
        z = ((lambda - syy)*sxz + syz*sxy) / det[no];
      }
      norm = sqrt(x*x + y*y + z*z);
      if (norm != 0.0)
      {
        x /= norm;
        y /= norm;
        z /= norm;
      }
      if (y > 0.0)
      {
        x = -x;
        y = -y;
        z = -z;
      }
      ux[n] = x;
      uy[n] = y;
      uz[n] = z;
    }
  return 1;
}


long EqDeg2(double g_xx, double g_xy, double g_yy,
            long indice_vp, double * vp_x, double * vp_y, double * vp1, double * vp2)
{
  double delta, sq, r, lambda1, lambda2, det, vx = 0.0, vy = 0.0, epsilon = 1e-10;

  delta = (g_xx - g_yy)*(g_xx - g_yy) + 4.0*g_xy*g_xy;
  if (delta < -epsilon)
    return 0;
  sq = sqrt(fabs(delta));
  r = g_xx + g_yy;
  lambda1 = (r + sq) / 2.0;
  lambda2 = (r - sq) / 2.0;
  if (vp1 != NULL)
    *vp1 = lambda1;
  if (vp2 != NULL)
    *vp2 = lambda2;
  det = sqrt(g_xy*g_xy + (g_xx - lambda1)*(g_xx - lambda1));
  if (det > 0.0)
  {
    vx = -g_xy / det;
    vy = (g_xx - lambda1) / det;
  }
  else
  {
    det = sqrt(g_xy*g_xy + (g_yy - lambda1)*(g_yy - lambda1));
    if (det > 0.0)
    {
      vx = (g_yy - lambda1) / det;
      vy = -g_xy / det;
    }
    return 0;
  }
  if (indice_vp == 1)
  {
    if (vy > 0.0)
    {
      vx = -vx;
      vy = -vy;
    }
    if (vp_x != NULL)
      *vp_x = vx;
    if (vp_y != NULL)
      *vp_y = vy;
  }
  else if (indice_vp == 2)
  {
    if (vx > 0.0)
    {
      vx = -vx;
      vy = -vy;
    }
    if (vp_x != NULL)
      *vp_x = -vy;
    if (vp_y != NULL)
      *vp_y = vx;
  }
  return 1;
}
