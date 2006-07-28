#% -*- mode: tm; mode: fold -*-

#%{{{Macros

#i linuxdoc.tm
#d it#1 <it>$1</it>

#d isis \bf{isis}
#d slang \bf{S-lang}
#d exmp#1 \tt{$1}
#d var#1 \tt{$1}

#d ivar#1 \tt{$1}
#d ifun#1 \tt{$1}
#d cvar#1 \tt{$1}
#d cfun#1 \tt{$1}
#d svar#1 \tt{$1}
#d sfun#1 \tt{$1}
#d icon#1 \tt{$1}

#d chapter#1 <chapt>$1<p>
#d preface <preface>
#d tag#1 <tag>$1</tag>

#d function#1 \sect{<bf>$1</bf>\label{$1}}<descrip>
#d variable#1 \sect{<bf>$1</bf>\label{$1}}<descrip>
#d function_sect#1 \sect{$1}
#d begin_constant_sect#1 \sect{$1}<itemize>
#d constant#1 <item><tt>$1</tt>
#d end_constant_sect </itemize>

#d synopsis#1 <tag> Synopsis </tag> $1
#d keywords#1 <tag> Keywords </tag> $1
#d usage#1 <tag> Usage </tag> <tt>$1</tt>
#d description <tag> Description </tag>
#d example <tag> Example </tag>
#d notes <tag> Notes </tag>
#d seealso#1 <tag> See Also </tag> <tt>\linuxdoc_list_to_ref{$1}</tt>
#d done </descrip><p>
#d -1 <tt>-1</tt>
#d 0 <tt>0</tt>
#d 1 <tt>1</tt>
#d 2 <tt>2</tt>
#d 3 <tt>3</tt>
#d 4 <tt>4</tt>
#d 5 <tt>5</tt>
#d 6 <tt>6</tt>
#d 7 <tt>7</tt>
#d 8 <tt>8</tt>
#d 9 <tt>9</tt>
#d NULL <tt>NULL</tt>
#d documentstyle book

#d iflatex#2 <#if output=latex2e>$1</#if><#unless output=latex2e>$2</#unless>
#d ifhtml#2 <#if output=html>$1</#if><#unless output=html>$2</#unless>

#%}}}

#d module#1 \tt{$1}
#d file#1 \tt{$1}
#d slang-documentation \
 \url{http://www.s-lang.org/doc/html/slang.html}{S-Lang documentation}

#d isis-documentation \
 \url{http://space.mit.edu/cxc/isis/manual.html}{ISIS documentation}

#d models_paper_url http://arxiv.org/abs/astro-ph/0607574
#d models_paper \ifhtml{\url{\models_paper_url}{Houck & Allen (2006)}}{Houck & Allen (2006) [ApJS, in press, (astro-ph/0607574)]}

\linuxdoc
\begin{\documentstyle}

\title Nonthermal Module Reference
\author John C. Houck, \tt{houck@space.mit.edu}
\date \__today__

\toc

\chapter{Introduction to the Nonthermal Module}

 This module provides spectral models for sources that emit
 photons via the synchrotron, inverse Compton, nonthermal
 bremsstrahlung and neutral-pion decay processes. The physical
 assumptions involved and various computational details are
 discussed in \models_paper.

 This manual describes the user interface of the module and
 presents a small number of illustrative examples.

\chapter{Examples}

 In the following examples, we assume that the nonthermal
 module has been installed in the default location so that
 \isis will find it automatically.

\sect{X-ray Synchrotron Emission}

 In this section, we describe how to fit X-ray data using the
 synchrotron model.  For a more general discussion of how to use
 \isis in data analysis, see the \isis-documentation.

 First, it is necessary to import the nonthermal module:
#v+
   require ("nonthermal");
#v-
 to make the spectral models and related functions available to
 \isis.

 Next, load the spectral data, background spectrum and
 instrument responses:
#v+
   variable pha = load_data ("pha.fits");
   assign_arf (load_arf ("arf.fits"), pha);
   assign_rmf (load_rmf ("rmf.fits"), pha);
   () = define_bgd (pha, "back.fits");
#v-
 At this point, it may be desirable to ignore selected spectral
 intervals or to rebin the data to improve the count-statistics
 in certain regions.

 Next, define the spectral model. For this example, we will
 assume that a simple model consisting of a single synchrotron
 emission component modified by line of sight absorption. To
 compute the synchrotron emission spectrum, it is also
 necessary to specify the underlying particle momentum
 distribution function.  For this example, we will use the
 \exmp{default} momentm distribution function described in
 \models_paper. To model the effects of line of sight
 absorption, we will use the \exmp{phabs} function from XSPEC.
 The spectral model is therefore:
#v+
   fit_fun ("phabs(1) * sync(1, default(1))");
#v-
 Before fitting the data, is usually a good idea to adjust the
 parameter values so that the starting model is reasonably
 close to the data.

 Once we're satisfied with the initial parameter values, we
 can perform a fit and then overplot the model and data:
#v+
   () = fit_counts;
   rplot_counts (pha);
#v-
 At this point, it is usually a good idea to investigate the
 neighborhood of the fitted parameters to make sure that the
 best possible fit has been achieved.  For example, we might
 use \exmp{conf} to compute single parameter confidence limits
 for fit-parameters of interest.

\sect{Inverse-Compton Gamma-ray Emission}

 The process of fitting gamma-ray data is essentially identical
 to that used to fit X-ray data in the previous example. The
 primary difference is that TeV gamma-ray spectra must be
 converted into a form that \isis recognizes.

 Putting HEGRA observations of Cas A into this format, we obtain
#v+
# flux in photons/s/cm^2/bin
; object   Cas A
; bintype  counts
; xunit    TeV
; exposure 1.0
#  E_lo            E_hi            flux            error  
   4.997493e-01    7.913994e-01    2.839012e-12    1.783329e-12
   7.913994e-01    1.256344e+00    3.729176e-13    2.399393e-13
   1.256344e+00    1.992807e+00    1.503555e-13    7.308948e-14
   1.992807e+00    3.168774e+00    7.704922e-14    3.311765e-14
   3.168774e+00    4.997493e+00    4.837415e-14    1.612472e-14
   4.997493e+00    7.913994e+00    1.115603e-14    1.247306e-14
   7.913994e+00    1.258409e+01    1.821113e-14    1.224027e-14
   1.258409e+01    2.001003e+01    1.804125e-15    3.373713e-15
#v-
 Details of this ASCII format are described in the
 \isis-documentation.  Lines beginning with a \exmp{#} symbol
 are ignored and may be used to insert comments.  The first two
 columns define the lower and upper edges of histogram bins.
 The bins must be in monotonic increasing order and may not
 overlap.  The next two columns contain the flux in
 \exmp{photons/s/cm^2} and the associated uncertainty, assumed
 to be symmetric. Lines beginning with a semicolon (;) define
 various keywords recognized by \isis. The \exmp{xunit} keyword
 indicates that the spectral bins are in TeV units.  The
 \exmp{exposure} keyword specifies a nominal exposure time of
 one second.  
 
 Note that the \exmp{bintype} keyword labels the data as
 ``counts'' in spite of the fact that the spectral values are
 given in flux units.  This subterfuge allows us to
 simultaneously fit X-ray data (in counts) with gamma-ray data
 (in flux units). But to make it work, we have to explicitly
 turn off one of the normal data-input validation checks. By
 default, when \isis reads counts spectra, it requires that the
 uncertainty in the number of counts in each bin be larger than
 one. If the input data do not satisfy this criterion, \isis
 replaces the input uncertainties with acceptable values.  To
 keep \isis from modifying our input uncertainties in this way,
 we set the intrinsic variable \exmp{Minimum_Stat_Err} to a
 small positive value (smaller than any of the input
 uncertainties).
 
 Having converted the Cas A HEGRA spectrum into the above
 format, we can load the data into isis:
#v+
  Minimum_Stat_Err = 1.e-30;
  variable id = load_data ("casa_hegra_data.txt");
#v-

 Once the data have been loaded, we can define and fit a model
 exactly as before.  For example, to fit a model which includes
 both inverse Compton and neutral-pion decay components, we can
 use
#v+
   fit_fun ("invc(1, default(1)) + pizero(1, default(2))");
#v-
 Note that two particle distribution functions appear in this
 spectral model. One, \exmp{default(1)}, refers to the
 nonthermal electron momentum distribution and the other,
 \exmp{default(2)}, refers to the nonthermal proton
 distribution.

\chapter{Spectral Models Reference}
#i rtl/models.tm

\chapter{Function Reference}
#i nonthermalfuns.tm
#i rtl/tables.tm

\end{\documentstyle}
