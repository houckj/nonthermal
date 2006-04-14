require ("gsl");

% Synchrotron functions from Westfold (1959) ApJ, 130, 241

define F(x)
{
   variable s;
   try
     {
        s = synchrotron_1(x);
     }
   catch AnyError:
     {
        s = 0.0;
     }
   return s;
}

define Fp(x)
{
   variable s;
   try
     {
        s = synchrotron_2(x);
     }
   catch AnyError:
     {
        s = 0.0;
     }
   return s;
}

define F1(x)
{
   return 0.5 * (F(x) - Fp(x));
}

define F2(x)
{
   return 0.5 * (F(x) + Fp(x));
}

define Gp(x)
{
   variable k;

   try
     {
        k = bessel_Knu (1.0/3, x);
     }
   catch AnyError:
     {
        k = 0.0;
     }

   return x^(1.0/3) * k;
}

private define sG(x, gam)
{
   return ((gam+7.0/3) * Gp(x) -
     2.0 * x^(0.5*(gam-1.0)) * (F(x) - Fp(x)))/(gam+1.0);
}

define G(x, gam)
{
   if (typeof(x) == Array_Type)
     return array_map (Double_Type, &sG, x, gam);

   return sG(x,gam);
}

define westfold (nu, gam, nu_c2, nu_c1)
{
   return (G(nu/nu_c2,gam) - G(nu/nu_c1,gam)) / G(1.e-20,gam);
}

#iffalse
define main ()
{
   variable x = 10.0^[-3:3:0.01];
   variable gam = 2;

   xlog;
   plot (x, westfold(x, gam, max(x)*1.e2, min(x)*1.e-2));
}

main();
#endif
