#define MOMENTUM_THRESHOLD 0.78  /* GeV/c */

static int npxsec (double e_pion, double e_proton, double *xsec) /*{{{*/
{
   double p_proton, e_proton_tot;

   (void) e_pion;

   /*  Input:  kinetic energies e_pion, e_proton in erg
    * Return:  cross-section in cm^2
    */

   *xsec = 0.0;

   e_proton_tot = e_proton + PROTON_REST_ENERGY;

   /* lab frame momentum */
   p_proton = (sqrt (e_proton_tot*e_proton_tot
                     - PROTON_REST_ENERGY*PROTON_REST_ENERGY)
               /GSL_CONST_CGSM_SPEED_OF_LIGHT);

   /* convert to units of GeV/c */
   p_proton /= GEV /GSL_CONST_CGSM_SPEED_OF_LIGHT;

   if (p_proton < MOMENTUM_THRESHOLD)
     return 0;

   if (p_proton < 0.96)
     {
        double m_pi2, m_chi, m_chi2, s, eta, eta2, eta4;
        m_pi2 = SQR_PION_MASS_ENERGY;
        m_chi = 2*PROTON_REST_ENERGY;
        m_chi2 = m_chi * m_chi;
        s = m_chi * (e_proton + m_chi);
        eta = 0.5 * (sqrt (((s - m_pi2 - m_chi2) - 4 * m_pi2 * m_chi2)/s)
                     / PI_ZERO_REST_ENERGY);
        eta2 = eta * eta;
        eta4 = eta2 * eta2;

        *xsec = eta2 * (0.032 + eta4 * (0.040 + eta2 * 0.047));
     }
   else if (p_proton < 1.27)
     {
        *xsec = 32.6 * pow((p_proton - 0.8), 3.21);
     }
   else if (p_proton < 8.0)
     {
        *xsec = 5.40 * pow((p_proton - 0.8), 0.81);
     }
   else /* p_proton >= 8 GeV/c */
     {
        *xsec = 32.0 * log(p_proton) + 48.5 / sqrt(p_proton) - 59.5;
     }
   
   /* convert millibarns => cm^2 */
   *xsec *= 1.e-27; 

   return 0;
}

/*}}}*/
