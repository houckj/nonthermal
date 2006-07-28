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
   f(p,a) = -index + curvature * log10(p/p0)
#v-
 and where the value of \tt{p0} is determined by the value
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
