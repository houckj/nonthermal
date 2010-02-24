#c __FILE__: ../../lib/nonthermal.sl
#c __LINE__: 309
\function{add_pdf}
\synopsis{Add a particle distribution function}
\usage{add_pdf (String_Type libname, String_Type pdfname [, param_names[], param_values[], freeze[], min[], max[]])}
\description
  The \var{add_pdf} function can be used to load a particle
  distribution function (PDF) from a shared library.  An example
  implementation of such an object is provided in the
  \var{examples} directory of the source code distribution.

  The first two arguments give the name of the shared library
  and the name of the PDF.  If the PDF has any fit parameters,
  the remaining array arguments should provide the parameter names,
  default values, default freeze/thaw state and default allowed
  minimum and maximum values.

\seealso{}
\done
#c __LINE__: 403
\function{make_sync_table}
\synopsis{Generate a synchrotron lookup table}
\usage{make_sync_table ([String_Type file])}
\description
  The \var{make_sync_table} function can be used to generate
  a new lookup table for the angular integration that arises
  in the computation of the \var{sync} model.

  The file name may be provided as a parameter.  If no
  file name is provided, the default file name is used.
  The output file is a FITS bintable.

  To customize details of the table computation, modify
  the script \var{lib/sync_make_table.sl} which is
  installed in \var{${prefix}/share/isis/nonthermal}.

\seealso{_sync_angular_integral}
\done
#c __LINE__: 449
\function{make_invc_table}
\synopsis{Generate an inverse Compton lookup table}
\usage{make_invc_table ([String_Type file])}
\description
  The \var{make_invc_table} function can be used to generate
  a new lookup table for the integral over incident photons that
  arises in the computation of the \var{invc} model.

  The file name may be provided as a parameter.  If no
  file name is provided, the default file name is used.
  The output file is a FITS bintable.

  To customize details of the table computation, modify
  the script \var{lib/invc_make_table.sl} which is
  installed in \var{${prefix}/share/isis/nonthermal}.

  A blackbody radiation field is used by default.  To use
  a different radiation field, it is necessary to modify
  the source code [replacing src/bbody.c and making any
  other necessary changes]

\seealso{_invc_photon_integral}
\done
#c __LINE__: 498
\function{make_ntbrem_table}
\synopsis{Generate an ee Bremsstrahlung lookup table}
\usage{make_ntbrem_table ([String_Type file])}
\description
  The \var{make_ntbrem_table} function can be used to generate
  a new lookup table for the electron-electron bremsstrahlung
  differential cross-sections that arise in the computation
  of the \var{ntbrem} model.

  The file name may be provided as a parameter.  If no
  file name is provided, the default file name is used.
  The output file is a FITS bintable.

  To customize details of the table computation, modify
  the script \var{lib/ntbrem_make_table.sl} which is
  installed in \var{${prefix}/share/isis/nonthermal}.

\seealso{_ee_haug1, _ee_haug1_lab}
\done
#c __LINE__: 562
\function{ntb_set_process_weights}
\synopsis{Set the relative weights of ee and ep bremsstrahlung}
\usage{ntb_set_process_weights (wt_ee, wt_ep)}
\description
  The \var{ntb_set_process_weights} function can be used to
  set the weights used for electron-electron and electron-proton
  bremsstrahlung.  The total contribution from both processes is
#v+
    S(E) = wt_ee * S_ee(E) + wt_ep * S_ep(E)
#v-
  The fluxes are scaled relative to the ambient number density
  of nuclei, \var{n}. The default weights are appropriate for a
  fully ionized gas with cosmic abundances so that most
  electrons come from hydrogen and helium.  In that case,
  the number fractions of hydrogen and helium are roughly
#v+
     X_H  = 1.0/1.1 = 0.9091,
     X_HE = 0.1/1.1 = 0.0909.
#v-
  The electron-electron bremsstrahlung contribution,
  proportional to the ambient electron density has weight
#v+
     wt_ee = X_H + 2*X_HE = 1.09091
#v-
  and the total electron-proton bremsstrahlung contribution,
  proportional to \var{n*Z^2}, has weight
#v+
     wt_ep = X_H + 4*X_HE = 1.27273
#v-

\seealso{}
\done
#c __LINE__: 602
\function{particle_info_struct}
\synopsis{Get a template structure for specifying the PDF}
\usage{Struct_Type = particle_info_struct ()}
\description
  The \var{particle_info_struct} returns a template
  structure with the following fields:
#v+
   __Name__       __Definition__
   pdf_name        name of the particle distribution function
   params          parameter values
   particle_type   electron=0, proton=1
   n_GeV           density at 1 GeV [cm^-3 GeV^-1]
   kT              temperature of thermal particles
   n_th            density of thermal particles
#v-

\seealso{}
\done
#c __LINE__: 670
\function{find_momentum_min}
\synopsis{Find the momentum at which the thermal and nonthermal PDFs intersect}
\usage{pc_min = find_momentum_min (Struct_Type p)}
\description
  The \var{find_momentum_min} function takes a structure
  defining the thermal and nonthermal components of the
  particle distribution function and computes the momentum
  at which those two distributions contribute equally.

  If the thermal PDF is everywhere below the nonthermal PDF,
  the momentum at the thermal peak is returned.

\seealso{particle_info_struct}
\done
#c __LINE__: 692
\function{nontherm_density}
\synopsis{Compute the density of nonthermal particles}
\usage{n = nontherm_density (Struct_Type p)}
\description
  The \var{nontherm_density} function takes a structure
  defining the thermal and nonthermal components of the
  particle distribution function (PDF) and computes the integral
  over the nonthermal PDF.  The lower limit of this
  integral is the momentum at which the thermal and nonthermal
  PDFs intersect. If the thermal PDF is everywhere below the
  nonthermal PDF, the momentum at the thermal peak is used.

\seealso{find_momentum_min}
\done
#c __LINE__: 715
\function{nontherm_energy_density}
\synopsis{Compute the energy density of nonthermal particles}
\usage{n = nontherm_energy_density (Struct_Type p)}
\description
  The \var{nontherm_energy_density} function takes a structure
  defining the thermal and nonthermal components of the
  particle distribution function (PDF) and computes the
  energy density in nonthermal particles by integrating over
  the nonthermal PDF.  The lower limit of this
  integral is the momentum at which the thermal and nonthermal
  PDFs intersect. If the thermal PDF is everywhere below the
  nonthermal PDF, the momentum at the thermal peak is used.

\seealso{find_momentum_min}
\done
#c __LINE__: 739
\function{force_charge_conservation}
\synopsis{Enforce charge conservation by adjusting the proton norm}
\usage{force_charge_conservation (electrons, protons [, method])}
\description
  The \var{force_charge_conservation} function take structures
  defining the electron and proton particle distribution
  functions and adjusts the normalization of the nonthermal
  proton distribution (protons.n_GeV) to enforce some definition
  of charge conservation.  The supported definitions are
#v+
   method=0 equal density of electrons and protons at some chosen
            injection kinetic energy
   method=1 equal integrated density of nonthermal electrons and protons
#v-
  The default is \var{method=0}.

\seealso{particle_info_struct}
\done
