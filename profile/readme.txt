
PROFILE is a utility for manipulating and/or displaying the
coordinates and other properties of airfoils.  It processes
one or more profiles at a time in various ways (one way per
run).   The input and output profiles may be in one of four
formats described below (basically as separate surfaces, as
a wrap-around surface, or in three-column form).

If plotting of the airfoil geometry is  requested,  all  of
the profiles  in  the input  dataset will be plotted on the
same pair of axes unless the  "THREED"  option is selected.
The plots  can  be  of  the  original input data,  or  data
derived by PROFILE,  or both.   Curvature distributions may
also be plotted as may the optional pressure distributions.

Plotting of the airfoil geometries, curvature distributions
or pressure distributions is handled by a separate program,
QPLOT, which should accompany PROFILE.    Users' guides are
available for PROFILE and QPLOT.

The  four  main choices for manipulating input profiles are
to REDISTRIBUTE the abscissas using conventional splines or
parametric splines depending on whether leading or trailing
edges are rounded;  to MODIFY or perturb the  ordinates  in
some way according to shape functions added  interactively;
to REFINE the ordinates by manipulating the curvature (act-
ually 2nd derivative) distribution while seeking or retain-
ing some maximum thickness value; and to OPTIMIZE one surf-
ace of one airfoil automatically, using a predetermined set
of shape functions with some of their parameters  variable,
and a target curvature distribution to  be  matched  in the
least squares sense.

Lesser options permit the user to RECTIFY profiles which do
not have the point common to the two surfaces in the  usual
place, and to NORMALIZE or DENORMALIZE profiles. TRANSFORM-
ing between upper/lower surface representation and  camber/
thickness representation (either way) is provided for, with
decambering as an option.  Applying twist is available from
the ROTATE option.

Two options involve a "secondary" profile  (prompted for at
a lower level): an option to COMBINE profiles (one of which
may be boundary layer displacement thickness);  and  an op-
tion to LOFT linearly between two profiles.

A  "nose-job"  option permits ROUNDing or SHARPENing of the
leading edge region  -  operations which have been made  as
reversible as possible.   More generally, the SMOOTH option
smooths the entire airfoil (or just one surface) by fitting
combinations of "Wagner" functions, which are also employed
by the OPTIMIZE and MODIFY options,  or  by implicit and/or
explicit smoothing (possibly weighted nonuniformly).

Tabulation of coordinates along with derivatives and curva-
ture is provided, as well as saving the manipulated profile
as a new output file.

Saving of y" distributions is also an option,  for possible
editing and reuse in REFINE mode.

Spreadsheet-compatible output of all likely tabular quanti-
ties is also provided for the simpler operating modes. This
requires the two surfaces to have common abscissas.

NOTES:
PROFILE  has evolved considerably since its inception as  a
basic redistribute-and/or-plot program. Provision for arbi-
trary perturbations to the input geometry,  with tabulation
of the resultant coordinates and derivatives, required some
reorganization, but the structure should now serve for most
likely purposes. Some implementation considerations follow.

       *  The case of a 2-element airfoil forced the decision to plot
all input profiles on the same page  (although the "THREED"
option has since been introduced).    Normalization of such
combinations proved an awkward option to provide.  The user
should  use  the normalization option carefully.  Normaliz-
ation of 3-D wings is not available.

       *  The multiple-profile case also prevented plotting  of  more
than one frame type (such as curvature distributions in ad-
dition  to  the  airfoils)  -  hence the saving of separate
files for later plotting.

       *  Large-scale plots are feasible  (with optional  windowing),
but exact scaling cannot be guaranteed because of the vari-
ability  of  the  output devices  available.  (Plots of the
same data on the same device can vary slightly from plot to
plot.)

       *  Derivatives and curvature values are normally estimated  by
finite differences for consistency with REFINE and OPTIMIZE
modes.    It is well known that these approximations can be
poor in the presence of very small X increments and limited
precision in the Ys.

       *  An option to plot the full wrap-around curvature  distribu-
tion using parametric spline derivatives has been  provided
for a proper look at the leading edge region.  But the .ypp
file of 2nd derivatives is suppressed in this case to avoid
inappropriate use with the REFINE mode.

        *  For simplicity, each of the MODIFY,  REFINE,  and  OPTIMIZE
options assumes that the coordinates have been normalized.

 METHOD:

The basic steps are as follows:

       *  Prompt for mode of operation and the input profile file name.

      *  Set defaults for user inputs and use an input control file
to override some of them if necessary.

      *  Scan all of the input profiles, for scaling and normalizing
purposes. Use EOF to handle the unknown number of profiles.

      *  Rewind the file and process the (now known number of) profiles
one at a time, according to the selected mode.

       *  Write the following output files as requested:
       
   Revised profile coordinates in one of 4 formats

   Original and/or revised airfoil geometry for plotting
   (a QPLOT file)

   Tabulated coordinates with derivatives and curvatures

   A more complete, spreadsheet-compatible file

   Second derivatives for possible reuse by REFINE mode

   Original and revised curvature data, including target
   curvature data for OPTIMIZE mode (another QPLOT file)

   Cps estimated for original and revised airfoil (QPLOT
   format)


