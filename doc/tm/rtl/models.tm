\function{sync}
\synopsis{Synchrotron spectrum model}
\usage{sync(Int_Type i, String_Type pdf_name, Double_Type pdf_params[])}
\description
  The \ifun{sync} function computes the synchrotron spectrum
  for a specified nonthermal particle distribution function.
  The synchrotron spectrum parameters are
#v+
    B_tot = total magnetic field strength [microgauss].
     norm = normalization = A_e V / (4 \pi d^2)
#v-
where
#v+
      A_e = nonthermal electron density at 1 GeV [cm^-3 GeV^-1]
        V = emitting volume [cm^3]
        d = distance [cm]
#v-
  The first argument of the \ifun{sync} function is an integer
  index identifying a particular instance of this fit
  function.  The second and third arguments specify the
  particle distribution function (PDF) and its parameters.

  Normally, the second and third arguments are provided as
  return values by a special fit-function that acts as an
  interface to the desired PDF. For example,
#v+
      fit_fun ("sync(1, default(1))")
#v-
  causes the synchrotron function calculation to use the PDF
  named \ifun{default}, which has parameters
#v+
         index = momentum power-law index
     curvature = amount of curvature above 1 GeV
      E_cutoff = cut-off momentum [TeV].
#v-
  When this fit-function is evaluated, \ifun{default(1)} returns
  second and third parameters needed by \ifun{sync}.

  By default, the value of the angular integral is determined via
  spline interpolation in a pre-computed table. To compute the
  angular integral value via direct integration, set the
  intrinsic variable \ivar{Syn_Interpolate} to a non-zero value.
\seealso{_sync_angular_integral}
\done

\function{invc}
\synopsis{Inverse Compton spectrum model}
\usage{invc(Int_Type i, String_Type pdf_name, Double_Type pdf_params[])}
\description
  The \ifun{invc} function computes the inverse Compton
  spectrum for a specified nonthermal particle distribution
  function.  The inverse Compton spectrum parameters are
#v+
  T_photon = radiation field temperature [K]
      norm = normalization = A_e V / (4 \pi d^2)
#v-
where
#v+
      A_e = nonthermal electron density at 1 GeV [cm^-3 GeV^-1]
        V = emitting volume [cm^3]
        d = distance [cm]
#v-

  The first argument of the \ifun{invc} function is an integer
  index identifying a particular instance of this fit
  function.  The second and third arguments specify the
  particle distribution function (PDF) and its parameters.

  Normally, the second and third arguments are provided as
  return values by a special fit-function that acts as an
  interface to the desired PDF. For example,
#v+
      fit_fun ("invc(1, default(1))")
#v-
  causes the inverse Compton function calculation to use the
  PDF named \ifun{default}, which has parameters
#v+
         index = momentum power-law index
     curvature = amount of curvature above 1 GeV
      E_cutoff = cut-off momentum [TeV].
#v-
  When this fit-function is evaluated, \ifun{default(1)} returns
  second and third parameters needed by \ifun{invc}.

  By default, the value of the integral over incident photon
  energies is determined via spline interpolation in a
  pre-computed table and the value of the \var{T_photon}
  parameter is not used. To compute this integral via direct
  integration, set the intrinsic variable \ivar{IC_Interpolate}
  to a non-zero value.
\seealso{_invc_photon_integral}
\done

\function{ntbrem}
\synopsis{Nonthermal bremsstrahlung spectrum model}
\usage{ntbrem(Int_Type i, String_Type pdf_name, Double_Type pdf_params[])}
\description
  The \ifun{ntbrem} function computes the nonthermal
  bremsstrahlung spectrum for a specified nonthermal particle
  distribution function, assuming a stationary target.
  The primary nonthermal bremsstrahlung spectrum parameter 
  is the normalization,
#v+
      norm = n_t A_e V / (4 \pi d^2)
#v-
  where
#v+
      n_t = density of target nuclei [cm^-3]
      A_e = nonthermal electron density at 1 GeV [cm^-3 GeV^-1]
        V = emitting volume [cm^3]
        d = distance [cm]
#v-

  The first argument of the \ifun{ntbrem} function is an
  integer index identifying a particular instance of this fit
  function.  The second and third arguments specify the
  particle distribution function (PDF) and its parameters.

  Normally, the second and third arguments are provided as
  return values by a special fit-function that acts as an
  interface to the desired PDF.  For example,
#v+
      fit_fun ("ntbrem(1, default(1))")
#v-
  causes the nonthermal bremsstrahlung function calculation
  to use the PDF named \ifun{default}, which has parameters
#v+
         index = momentum power-law index
     curvature = amount of curvature above 1 GeV
      E_cutoff = cut-off momentum [TeV].
#v-
  When this fit-function is evaluated, \ifun{default(1)} returns
  second and third parameters needed by \ifun{ntbrem}.

  By default, the value of the electron-electron differential
  cross-section is determined via spline interpolation in a
  pre-computed table. To compute these cross-sections without
  interpolation, set the intrinsic variable
  \ivar{Ntb_Interpolate} to a non-zero value.

  The relative contributions of electron-electron and
  electron-proton bremsstrahlung are controlled by the
  \ifun{ntb_set_process_weights} function.

\seealso{ntb_set_process_weights,_ee_haug1,_ee_haug1_lab}
\done

\function{pizero}
\synopsis{Neutral-pion decay spectrum model}
\usage{pizero(Int_Type i, String_Type pdf_name, Double_Type pdf_params[])}
\description
  The \ifun{pizero} function computes the neutral-pion decay
  gamma-ray spectrum for a specified nonthermal particle
  distribution function.  The primary pion decay spectrum
  parameter is the normalization,
#v+
      norm = n_p A_p V / (4 \pi d^2)
#v-
  where
#v+
      n_p = density of target protons [cm^-3]
      A_p = nonthermal proton density at 1 GeV [cm^-3 GeV^-1]
        V = emitting volume [cm^3]
        d = distance [cm]
#v-
  The first argument of the \ifun{pizero} function is an
  integer index identifying a particular instance of this fit
  function.   The second and third arguments specify the particle
  distribution function (PDF) and its parameters.

  Normally, the second and third arguments are provided as
  return values by a special fit-function that acts as an
  interface to the desired PDF. For example,
#v+
      fit_fun ("pizero(1, default(1))")
#v-
  causes the neutral-pion decay function calculation
  to use the PDF named \ifun{default}, which has parameters
#v+
         index = momentum power-law index
     curvature = amount of curvature above 1 GeV
      E_cutoff = cut-off momentum [TeV].
#v-
  When this fit-function is evaluated, \ifun{default(1)} returns
  second and third parameters needed by \ifun{pizero}.

\seealso{}
\done

