require ("nonthermal");

variable e_elec, e_phot, sigma_ee, sigma_ep;

e_phot = 10.0^[-3.0:9.0:0.1];   % photon energy, MeV
e_elec = 10.0^[4:9:0.5];        % electron kinetic energy, MeV

xlabel ("E [MeV]");
ylabel (latex2pg("\\sigma_{ep}, \\sigma_{ee}"));

variable x = ones(length(e_phot));

foreach (e_elec)
{
   variable ee = ();   
   
   sigma_ee = _ee_haug ((ee/0.511)*x, e_phot/0.511);
   sigma_ep = _ep_heitler ((ee/0.511)*x, e_phot/0.511);
   
   xlog;ylog;
   plot (e_phot, sigma_ep);
   oplot (e_phot, sigma_ee);
   plot_pause;
}
