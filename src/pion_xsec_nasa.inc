#define MILLIBARNS_PER_GEV (1.e-27/GEV)

#define PROTON_E1   (0.3)  /* GeV */
#define PROTON_E2   (2.0)  /* GeV */
#define PROTON_E3  (50.0)  /* GeV */
#define UNUSED    DBL_MAX

typedef struct
{
   double coeff[20];
}
Pion_Xsec_Type;

static Pion_Xsec_Type Pion_Xsec_Table [] =
{
  {{UNUSED,
       6.78e-10,  -2.86, 1.82e-8, -1.92,   22.3,
          0.226,  -0.33,   -1.75, -32.1, 0.0938,
          -23.7, 0.0313,   2.5e6,  1.38,   0.25,
          -39.4,   2.88,   0.025,  0.75}},
  {{UNUSED,
       1.3e-10,  -2.86, 4.27e-9,  -2.4,   22.3,
         -1.87,   1.28,   -1.25, -33.2, 0.0938,
         -23.6, 0.0313,   2.5e6,  0.25,  60322,
          1.07,  -67.5,  UNUSED,  UNUSED}}
};

static int npf (double *a, double e_pion, double e_proton, /*{{{*/
                double *f1, double *f2)
{
   double t;

   /* e_pion, e_proton in GeV */
   *f2 = a[1] * pow (e_pion, a[2]) + a[3] * pow (e_proton, a[4]);

   t = (a[5]  + a[6] / sqrt(e_proton) + a[7] * pow (e_proton, a[8])
        + a[9] * pow (e_pion, a[10]) + a[11] * pow (e_pion, a[12]));

   *f1 = exp(t);

   return 0;
}

/*}}}*/

static int npxsec1 (double e_pion, double e_proton, double *xsec) /*{{{*/
{
   double *a = Pion_Xsec_Table[0].coeff;
   double f1, f2, t1;

   /* e_pion, e_proton in GeV,  xsec in mb/GeV */
   if (-1 == npf (a, e_pion, e_proton, &f1, &f2))
     return -1;

   t1 = exp (a[16] * sqrt(e_pion)
             + a[17] * pow(e_pion, a[18]) * pow(e_proton, a[19]));
   *xsec = (a[13]*f1/f2 + a[14]*t1) * pow (e_pion, a[15]);

   return 0;
}

/*}}}*/

static int npxsec2 (double e_pion, double e_proton, double *xsec) /*{{{*/
{
   double *b = Pion_Xsec_Table[1].coeff;
   double f1, f2, t1;

   /* e_pion, e_proton in GeV,  xsec in mb/GeV */
   if (-1 == npf (b, e_pion, e_proton, &f1, &f2))
     return -1;

   t1 = b[15] * pow (e_pion, b[16]) * exp (b[17] * sqrt(e_pion));
   *xsec = b[13] * pow (e_pion, b[14]) * f1 / f2 + t1;

   return 0;
}

/*}}}*/

static int npxsec (double e_pion, double e_proton, double *xsec) /*{{{*/
{
   int (*sigma)(double, double, double *);
   int ret;

   /*  Input:  e_pion, e_proton in erg
    * Return:  cross-section in cm^2
    */

   /* convert erg -> GeV */
   *xsec = 0.0;
   e_pion /= GEV;
   e_proton /= GEV;

   if (e_proton < PROTON_E1 || PROTON_E3 < e_proton)
     return 0;

   /* FIXME? */
   if (e_pion > 2*e_proton)
     return 0;

   if (e_proton < PROTON_E2)
     sigma = npxsec1;
   else sigma = npxsec2;

   ret = (*sigma)(e_pion, e_proton, xsec);

   /* convert mb/GeV -> cm^2/erg */
   *xsec *= MILLIBARNS_PER_GEV;

   return ret;
}

/*}}}*/
