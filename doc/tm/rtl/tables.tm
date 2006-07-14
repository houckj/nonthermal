\function{_sync_angular_integral}
\synopsis{Synchrotron angular integral}
\usage{_sync_angular_integral (x, interp)}
\description
  The \ifun{_sync_angular_integral} function computes the angular
  integral needed for the synchrotron spectrum model.
  The first argument, \var{x}, is defined as
#v+
     x = f / f_c
#v-
  where \var{f_c} is the synchrotron critical frequency.

  If the second argument is non-zero, the integral will be
  evaluated by spline interpolation on a pre-computed table.
  Otherwise, the value will be computed by numerical
  integration.  
\seealso{}
\done

\function{_invc_photon_integral}
\synopsis{Inverse Compton integral over incident photon energies}
\usage{_invc_photon_integral (gamma, omega, T, interp)}
\description
  The \ifun{_invc_photon_integral} function computes the 
  integral over incident photon energies needed for the 
  inverse Compton spectrum model.  The first two
  parameters are
#v+
     gamma = 1/sqrt (1 - beta^2)
     omega = E / (m c^2)
#v-
  where gamma is the Lorentz factor of the electron,
  E is the energy of the incident photon and m is
  the electron mass.  The third argument is the
  radiation field temperature [K].  

  If the fourth is non-zero, the integral will be evaluated by
  spline interpolation on a pre-computed table. Otherwise, the
  value will be computed by numerical integration.    
\seealso{}
\done

\function{_ee_haug1}
\synopsis{Electron-electron bremsstrahlung lab frame differential cross-section}
\usage{sigma = _ee_haug1 (T, omega)}
\description
  The \ifun{_ee_haug1} function computes the
  lab-frame differential cross-section for electron-electron
  bremsstrahlung by applying a numerical Lorentz transformation
  to the center of momentum frame cross-section.  The 
  parameters are
#v+
       T = incident electron kinetic energy in units of the
           electron rest energy.
   omega = scattered photon energy in units of the electron
           rest energy
#v-
\seealso{_ee_haug1_lab}
\done

\function{_ee_haug1_lab}
\synopsis{Electron-electron bremsstrahlung lab frame differential cross-section}
\usage{sigma = _ee_haug1_lab (T, omega)}
\description
  The \ifun{_ee_haug1_lab} function computes the
  lab-frame differential cross-section for electron-electron
  bremsstrahlung.  The parameters are
#v+
       T = incident electron kinetic energy in units of the
           electron rest energy.
   omega = scattered photon energy in units of the electron
           rest energy
#v-
\seealso{}
\done

\function{_ee_heitler1}
\synopsis{Bethe-Heitler cross-section}
\usage{sigma = _ee_heitler1 (T, omega)}
\description
  The \ifun{_ee_heitler1} function computes the Bethe-Heitler
  cross-section for electron-ion bremsstrahlung.  The
  parameters are  
#v+
       T = incident electron kinetic energy in units of the
           electron rest energy.
   omega = scattered photon energy in units of the electron
           rest energy
#v-
\seealso{}
\done

