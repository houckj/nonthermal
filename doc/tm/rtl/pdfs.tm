#i linuxdoc.tm

 \function{default}
 \synopsis{default PDF}
 \usage{default(id)}
 \description
 This distribution function provides a power-law in momentum
 with curvature above some specific particle momentum and
 with an exponential cut-off.

 The distribution function has the form
#v+
  N(p) = A (p/p0)^f(p,a) exp((p0-p)/cutoff)
#v-
 where
#v+
   f(p,a) = -index + curvature * log10(p/p0),  p>p0
   f(p,a) = -index,                            p<=p0
#v-
 The value of \tt{p0} is determined by the value
 of \ivar{Min_Curvature_Pc}=1 GeV by default.

 The fit parameters are
#v+
  index       power-law index
  curvature   curvature parameter
  cutoff      cut-off energy [TeV]
#v-

#% \example
#% \notes
 \seealso{add_pdf}
 \done

 \function{cbreak}
 \synopsis{cooling-break PDF}
 \usage{cbreak(id)}
 \description
 This distribution function is the same as the default
 distribution function except that it also has a "cooling break"
 above a particular momentum.  The cooling break is
 implemented by adding an additional multiplicative factor
 of the form
#v+
              (break/pc), for pc > break
       g(p) =
              1,          for pc < break
#v-
 where \exmp{p} is the particle momentum and \exmp{break}
 is the 'break energy'.

 The fit parameters are then
#v+
  index       power-law index
  curvature   curvature parameter
  cutoff      cut-off energy [TeV]
  break       cooling-break energy [TeV]
#v-
#% \example
#% \notes
 \seealso{add_pdf}
 \done

 \function{ke_cutoff}
 \synopsis{PDF with cutoff in particle kinetic energy}
 \usage{ke_cutoff(id)}
 \description
 This distribution function is the same as the default
 distribution function except that the exponential cutoff
 depends on kinetic energy, \exmp{T}, instead of momentum,
 \exmp{p}.

#% \example
#% \notes
 \seealso{add_pdf}
 \done

 \function{etot}
 \synopsis{Energy-dependent PDF}
 \usage{etot(id)}
 \description
 This distribution function is a power-law in total-energy,
#v+
   N(p) = A (E/E0)^(-index)
#v-
 where \exmp{E0}=1 GeV.

 The power-law index is the only fit-parameter.

#% \example
#% \notes
 \seealso{add_pdf}
 \done

 \function{dermer}
 \synopsis{PDF from Dermer (1986)}
 \usage{dermer(id)}
 \description
 This distribution function is a power-law in total-energy,
 with momentum-dependent curvature above the momentum
 specified by the intrinsic variable \ivar{Min_Curvature_Pc},
 which is 1 GeV by default.

 The fit parameters are then
#v+
  index       power-law index
  curvature   curvature parameter
#v-
#% \example
#% \notes
 \seealso{add_pdf}
 \done

 \function{mori}
 \synopsis{PDF from Mori (1997)}
 \usage{mori(id)}
 \description
 This distribution function was taken from Mori (1997)
 and is designed to represent the Galactic cosmic-ray
 proton distribution.  It has no free parameters.
#% \example
#% \notes
 \seealso{add_pdf}
\done

 \function{boltz}
 \synopsis{Non-relativistic Maxwell-Boltzmann Distribution}
 \usage{boltz(id)}
 \description
 This distribution function provides a non-relativistic
 Maxwell-Boltzmann thermal particle distribution.  The only
 fit parameter is the temperature, \ivar{kT}, in keV.
#% \example
#% \notes
 \seealso{add_pdf}
\done

 \function{rboltz}
 \synopsis{Relativistic Maxwell-Boltzmann Distribution}
 \usage{rboltz(id)}
 \description

 This distribution function provides a relativistic
 Maxwell-Boltzmann thermal particle distribution. The only fit
 parameter is the temperature, \ivar{kT}, in keV.

#% \example
#% \notes
 \seealso{add_pdf, full1}
\done

 \function{full1}
 \synopsis{RMB + Cutoff Power-Law PDF}
 \usage{full1(id)}
 \description
 This particle distribution function consists of a relativistic
 Maxwell-Boltzmann thermal distribution (\ifun{rboltz}) plus a
 nonthermal tail (\ifun{default}), extending above some
 specified momentum offset above the thermal peak. The fit
 parameters are those of \ifun{rboltz} and \ifun{default} plus
 the additional parameter, \ivar{p_offset}, giving the momentum
 offset in units of the particle momentum at the thermal peak.

 For example, setting \exmp{p_offset=2}, means that the density
 of nonthermal particles exceeds the thermal particle density
 for momenta \exmp{p > (1+p_offset)*p_peak}, where
#v+
     p_peak = gamma * m * v
     v = sqrt (2kT/m)
     beta = v/c
     gamma = 1/sqrt (1-beta^2)
#v-

 When using this distribution function, note that the overall
 spectral normalization is scaled to the thermal particle density
 rather than the nonthermal particle density at 1 GeV.  For
 example, when used in conjunction with the \ifun{ntbrem} model,
 the spectral norm should be interpreted as
#v+
      norm = n_t n0 V / (r \pi d^2)
#v-
 where \exmp{n0} is the number density of thermal electrons.

#% \example
#% \notes
 \seealso{add_pdf, rboltz, default}
\done
