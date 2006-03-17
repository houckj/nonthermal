require ("nonthermal");
require ("gsl");

private variable Special_Log= -50.0;
private variable Special_Val = 10.0^Special_Log;

private variable Electron_Rest_Energy = Const_m_e * Const_c^2;
private variable CBR_Temp = 2.725;
private variable Omega_CBR = (Const_k * CBR_Temp) / Electron_Rest_Energy;

% magic number 32.42 to get results similar to earlier table
private variable Omega0 = 32.42 * Omega_CBR;

% Maximum range of interest
private variable Log_X_Range = [1.0, 9.5];    % gamma
private variable Log_Y_Range = log10(10.0^[2.0,15.0] * Const_eV / Electron_Rest_Energy);  % efinal
private variable Log_X_Right = 3.0;

private variable X_Epsilon = 1.e-5 * (Log_X_Range[1] - Log_X_Range[0]);
private variable Sigma0 = 1.e-27;

private define gamma_min (omega)
{
   return 0.5 * (omega + sqrt (omega * (omega + 1.0/Omega0)));
}

private define log_to_lglg (xlg, ylg)
{
   variable xlg_left = log10(gamma_min (10.0^ylg));
   variable xx = log (xlg - xlg_left);
   return xx;
}

private define lglg_to_log (xx, ylg)
{
   variable xlg_left = log10(gamma_min (10.0^ylg));
   variable xlg = exp (xx) + xlg_left;
   return xlg;
}

private define safe_log10 (x)
{
   variable lg;

   if (typeof(x) != Array_Type)
     {
        if (x > Special_Val)
          {
             lg = log10(x);
             return lg;
          }
        return Special_Log;
     }

   lg = @x;
   variable i = where (x > Special_Val);
   lg[i] = log10(x[i]);

   % avoid comparing floats for equality
   variable j = [0:length(x)-1];
   j[i] = -1;
   lg[where (j>=0)] = Special_Log;

   return lg;
}

private define fcn (xlg, ylg)
{
   variable gam, ef;
   gam = 10.0^xlg;
   ef = 10.0^ylg;

   variable f;
   if (length(gam) == 1)
     {
        f = _invc_photon_integral(gam, ef, CBR_Temp, 0);
     }
   else
     {
        f = array_map (Double_Type, &_invc_photon_integral, gam, ef, CBR_Temp, 0);
     }

   return safe_log10 (f / Sigma0);
}

private define _fcn_lglg (xx, ylg)
{
   variable xlg = lglg_to_log (xx, ylg);
   variable f = fcn (xlg,ylg);
   return f;
}

private define fcn_lglg (xx, ylg)
{
   return array_map (Double_Type, &_fcn_lglg, xx, ylg);
}

define make_ygrid (ny)
{
   variable ymn, ymx, ygrid;
   ymn = Log_Y_Range[0];
   ymx = Log_Y_Range[1];
#iftrue
   ygrid = ymn + (ymx - ymn) * [0:ny-1]/(ny-1.0);
#else
   message ("Using non-uniform ygrid ");
   variable _fy, dty, _ty;
   _fy = [0:ny-1]/(ny - 1.0);
   dty = exp(-3*_fy^2);
   _ty = cumsum(dty) / sum(dty);
   ygrid = ymn + (ymx - ymn) * _ty;
#endif

   return ygrid;
}

private define boundary_corners (num)
{
   variable xleft, yleft;
   yleft = make_ygrid (num);
   % shift last point slightly away from the boundary
   yleft[-1] = 0.5*(yleft[-2] + yleft[-1]);
   xleft = log10(gamma_min (10.0^yleft));

   variable xx, xlg, ylg;
   xlg = xleft + X_Epsilon;
   ylg = yleft;
   xx = log_to_lglg (xlg, ylg);

   return (xx, ylg);
}

private define compute_table (file)
{
   variable nx = 1024;
   variable ny = 1024;

   variable x, y;
   (x,y) = boundary_corners (ny);

   variable l = [0:length(x)-1];

   variable xleft, yy;
   xleft = x[l];
   yy = y[l];

   variable xmn, xmx, xgrid;
   xmn = min(xleft);
   xmx = Log_X_Right;
#iftrue
   xgrid = xmn + (xmx - xmn) * [0:nx-1]/(nx-1.0);
#else
   %message ("Using non-uniform xgrid ");
   variable _fx, dtx, _tx;
   _fx = [0:nx-1]/(nx - 1.0);
   dtx = exp(-3*_fx^2);
   _tx = cumsum(dtx) / sum(dtx);
   xgrid = xmn + (xmx - xmn) * _tx;
#endif

   variable ygrid = make_ygrid (ny);

   variable xl, xr, i, j, xlg, ylg, xx;
   variable f = Double_Type[ny,nx];

   _for (0, ny-1, 1)
     {
        j = ();
        () = fprintf (stderr, "%d\r", j);
        xl = log10(gamma_min (10.0^ygrid[j]));
        xl = log_to_lglg (xl + X_Epsilon, ygrid[j]);
        xr = Log_X_Right;
        i = where (xl <= xgrid and xgrid < xr);
        if (length(i) == 0)
          continue;
        f[j,i] = fcn_lglg (xgrid[i], ygrid[j]);
     }
   () = fprintf (stderr, "\n");

   variable keys = struct
     {
        gammin, gammax, efnmin, efnmax, xepsilon, omega0, sigma0
     };
   keys.gammin = Log_X_Range[0];
   keys.gammax = Log_X_Range[1];
   keys.efnmin = Log_Y_Range[0];
   keys.efnmax = Log_Y_Range[1];
   keys.sigma0 = Sigma0;
   keys.omega0 = Omega0;
   keys.xepsilon = X_Epsilon;

   variable tf = struct {f};
   variable tx = struct {xgrid};
   variable ty = struct {ygrid};
   tx.xgrid = xgrid;
   ty.ygrid = ygrid;
   tf.f = f;

   variable fp = fits_open_file (file, "c");
   fits_write_binary_table (fp, "TABLE", tf, keys, NULL);
   fits_write_binary_table (fp, "XGRID", tx);
   fits_write_binary_table (fp, "YGRID", ty);
   fits_close_file(fp);
}

define _invc_make_table (file)
{
   compute_table(file);
}

provide ("invc_make_table");
