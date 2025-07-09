




! The following changes to the source code were made by PDAS before
! inclusion in the library of "Public Domain Computer Programs for
! the Aeronautical Engineer"

! 2021-08-15 replaced CALL SYSTEM with CALL EXECUTE_COMMAND_LINE   RLC
! 2021-10-23 main program changed  REAL DUMMY to REAL DUMMY(1)
! 2021-10-23 Subroutine REDISTRIB changed  REAL DUMMY to REAL DUMMY(1)
! 2021-10-23 Subroutine REFINE changed REAL TARGET1,TARGET2,TARGET3 to
!                              REAL TARGET1(1),TARGET2(1),TARGET3(1)
! 2021-10-23 Subroutine WAGSMOOTH changed REAL YN to REAL YN(1)
! 2021-10-23 Added INTENT to Subroutine WAGSMOOTH
! 2021-10-23 Back in REDISTRIB, earlier change of DUMMY to DUMMY(1)
!            now requires the calls to CSFIT to use DUMMY(1)
! 2021-10-23 In PSEVAL, change TBEST to TBEST(1)
! 2021-10-23 In PSEVAL, in call to CSDVAL (2 places), change 
!            PSCOEFS(1,1,1) to (1:,1,1), (1,2,1) to (1:,2,1), and
!            (1,3,2) to (1:,3,2)
! 2021-10-24 In PSEVAL, change XLOC,XD1,XD2,YJ,YD1,YD2 to be dimension (1)
! 2021-10-24 In GETCPS, change CALL DECOMP(N,MAXN,AIC,IP) to
!            CALL DECOMP(N,MAXN,AIC,COND,IP,workarray) 
!            Also, declare COND and workarray
C+--------------------------------------------------------------------
C
      PROGRAM PROFILE
C
C  PURPOSE:
C           PROFILE is a utility for manipulating and/or displaying the
C           coordinates and other properties of airfoils.  It processes
C           one or more profiles at a time in various ways (one way per
C           run).   The input and output profiles may be in one of four
C           formats described below (basically as separate surfaces, as
C           a wrap-around surface, or in three-column form).
C
C           If plotting of the airfoil geometry is  requested,  all  of
C           the profiles  in  the input  dataset will be plotted on the
C           same pair of axes unless the  "THREED"  option is selected.
C           The plots  can  be  of  the  original input data,  or  data
C           derived by PROFILE,  or both.   Curvature distributions may
C           also be plotted as may the optional pressure distributions.
C
C           Plotting of the airfoil geometries, curvature distributions
C           or pressure distributions is handled by a separate program,
C           QPLOT, which should accompany PROFILE.    Users' guides are
C           available for PROFILE and QPLOT.
C
C           The  four  main choices for manipulating input profiles are
C           to REDISTRIBUTE the abscissas using conventional splines or
C           parametric splines depending on whether leading or trailing
C           edges are rounded;  to MODIFY or perturb the  ordinates  in
C           some way according to shape functions added  interactively;
C           to REFINE the ordinates by manipulating the curvature (act-
C           ually 2nd derivative) distribution while seeking or retain-
C           ing some maximum thickness value; and to OPTIMIZE one surf-
C           ace of one airfoil automatically, using a predetermined set
C           of shape functions with some of their parameters  variable,
C           and a target curvature distribution to  be  matched  in the
C           least squares sense.
C
C           Lesser options permit the user to RECTIFY profiles which do
C           not have the point common to the two surfaces in the  usual
C           place, and to NORMALIZE or DENORMALIZE profiles. TRANSFORM-
C           ing between upper/lower surface representation and  camber/
C           thickness representation (either way) is provided for, with
C           decambering as an option.  Applying twist is available from
C           the ROTATE option.
C
C           Two options involve a "secondary" profile  (prompted for at
C           a lower level): an option to COMBINE profiles (one of which
C           may be boundary layer displacement thickness);  and  an op-
C           tion to LOFT linearly between two profiles.
C
C           A  "nose-job"  option permits ROUNDing or SHARPENing of the
C           leading edge region  -  operations which have been made  as
C           reversible as possible.   More generally, the SMOOTH option
C           smooths the entire airfoil (or just one surface) by fitting
C           combinations of "Wagner" functions, which are also employed
C           by the OPTIMIZE and MODIFY options,  or  by implicit and/or
C           explicit smoothing (possibly weighted nonuniformly).
C
C           Tabulation of coordinates along with derivatives and curva-
C           ture is provided, as well as saving the manipulated profile
C           as a new output file.
C
C           Saving of y" distributions is also an option,  for possible
C           editing and reuse in REFINE mode.
C
C           Spreadsheet-compatible output of all likely tabular quanti-
C           ties is also provided for the simpler operating modes. This
C           requires the two surfaces to have common abscissas.
C
C  NOTES:
C           PROFILE  has evolved considerably since its inception as  a
C           basic redistribute-and/or-plot program. Provision for arbi-
C           trary perturbations to the input geometry,  with tabulation
C           of the resultant coordinates and derivatives, required some
C           reorganization, but the structure should now serve for most
C           likely purposes. Some implementation considerations follow.
C
C        *  The case of a 2-element airfoil forced the decision to plot
C           all input profiles on the same page  (although the "THREED"
C           option has since been introduced).    Normalization of such
C           combinations proved an awkward option to provide.  The user
C           should  use  the normalization option carefully.  Normaliz-
C           ation of 3-D wings is not available.
C
C        *  The multiple-profile case also prevented plotting  of  more
C           than one frame type (such as curvature distributions in ad-
C           dition  to  the  airfoils)  -  hence the saving of separate
C           files for later plotting.
C
C        *  Large-scale plots are feasible  (with optional  windowing),
C           but exact scaling cannot be guaranteed because of the vari-
C           ability  of  the  output devices  available.  (Plots of the
C           same data on the same device can vary slightly from plot to
C           plot.)
C
C        *  Derivatives and curvature values are normally estimated  by
C           finite differences for consistency with REFINE and OPTIMIZE
C           modes.    It is well known that these approximations can be
C           poor in the presence of very small X increments and limited
C           precision in the Ys.
C
C        *  An option to plot the full wrap-around curvature  distribu-
C           tion using parametric spline derivatives has been  provided
C           for a proper look at the leading edge region.  But the .ypp
C           file of 2nd derivatives is suppressed in this case to avoid
C           inappropriate use with the REFINE mode.
C
C        *  For simplicity, each of the MODIFY,  REFINE,  and  OPTIMIZE
C           options assumes that the coordinates have been normalized.
C
C  METHOD:
C
C           The basic steps are as follows:
C
C        *  Prompt for mode of operation and the input profile file name.
C
C        *  Set defaults for user inputs and use an input control file
C           to override some of them if necessary.
C
C        *  Scan all of the input profiles, for scaling and normalizing
C           purposes. Use EOF to handle the unknown number of profiles.
C
C        *  Rewind the file and process the (now known number of) profiles
C           one at a time, according to the selected mode.
C
C        *  Write the following output files as requested:
C
C              Revised profile coordinates in one of 4 formats
C
C              Original and/or revised airfoil geometry for plotting
C              (a QPLOT file)
C
C              Tabulated coordinates with derivatives and curvatures
C
C              A more complete, spreadsheet-compatible file
C
C              Second derivatives for possible reuse by REFINE mode
C
C              Original and revised curvature data, including target
C              curvature data for OPTIMIZE mode (another QPLOT file)
C
C              Cps estimated for original and revised airfoil (QPLOT
C              format)
C
C
C  MODES OF OPERATION:
C
C           MODE 0:  "Display" mode - no modifications involved.   Gen-
C                    erate requested output files,  which could include
C                    saving the coordinates in a different format. MODE
C                    <=3 is required for spreadsheet-compatible output.
C
C           MODE 1:  Rearrange or rectify the geometry data so that the
C                    common leading-edge point is indeed the  one  with
C                    minimum abscissa and shift ordinates by  an  input
C                    input amount if required.   Only  the revised pro-
C                    file may be tabulated/plotted in this case. A ver-
C                    tical shift option is also provided.
C
C           MODE 2:  Normalize profile(s) according to the total  range
C                    of x or by input chord & leading edge coordinates.
C                    A negative chord value will denormalize.  The same
C                    input values are used for each element of a multi-
C                    element airfoil.
C
C           MODE 3:  Redistribute the abscissas  and derive correspond-
C                    ing ordinates.   Conventional or parametric spline
C                    techniques are used depending on whether the lead-
C                    ing edge is sharp or rounded.  Distributions along
C                    the arc (in T rather than X) are an option.  Menu:
C
C                    -1 = Read new Xs (or Ts) from a file in standard
C                         PROFILE format (though y coordinates may be
C                         omitted if desired).
C                     0 = Distribute points uniformly.
C                     1 = Distribute points sinusoidally bunched near
C                         the leading edge.
C                     2 = Distribute  points  sinusoidally, near both
C                         the leading and trailing edges.
C                     3 = Sinusoidal bunching around an internal pt.
C                     4 = Vinokur distribution (first, last increments
C                         increments specirfied.
C
C                    A prompt  will  also  be  issued for the number of
C                    points to be generated on each surface.
C
C           MODE 4:  Perturb geometry data  according  to user-selected
C                    shape functions (interactive).
C
C           MODE 5:  Refine the airfoil,  typically modifying  (or  re-
C                    taining) its thickness while retaining (or modify-
C                    ing) its curvature distribution(s).  Numerous user
C                    inputs are prompted for in this case,  since there
C                    are several likely ways provided for  manipulating
C                    the curvature distributions. Defaults are provided
C                    where possible.  The 4 main choices:
C
C                    (1) Leave a surface unchanged  (while presumably
C                        refining the other);
C                    (2) Change the thickness with minimal changes to
C                        the existing curvature  (none of the follow-
C                        ing constraints on the  original y" values);
C                    (3) Impose a constant-2nd-derivative  constraint
C                        in some single region (to remove a bump or a
C                        spike in the curvature distribution, or per-
C                        haps  to  modify regions of flatness by con-
C                        straining second derivatives - and hence the
C                        curvature - away from zero);
C                    (4) Constrain the curvature via an input file of
C                        2nd derivative values (possibly derived from
C                        an earlier run  of  PROFILE,  or prepared by
C                        hand).  The table does not have to cover the
C                        whole abscissa range; linearly interpolating
C                        table look-ups are used.
C
C                    Brief descriptions of the inputs prompted  for  in
C                    "refine" mode follow:
C
C             *  Desired % thickness:  <CR>  retains present thickness.
C
C             *  Width param. for y:  Affects  the  nonuniform  scaling
C                                     applied  to  the ordinates  (both
C                    surfaces).   The default is 2.0.  Larger (3.0-4.0)
C                    tends to retain leading/trailing  edge  shape more
C                    while 1.0 would constrain fore and aft less.
C
C             *  Input y" table:      <CR>  means there is  none,  else
C                                     the file name is  entered.   This
C                    file should be in the standard  "PROFILE"  format.
C                    It can  cover  any  range  of  abscissas.  (Linear
C                    interpolation  is  used.)   It  may  be an  edited
C                    version  of  the  file  from  a  previous  run  of
C                    PROFILE, or it may be much cruder.  The 2nd deriv-
C                    ative values  entered  act  as constraints on  the
C                    curvature since curvature and y" are related if y'
C                    is not large.
C
C             *  Constant y" value:   <CR>  means no such constraint  -
C                                     retain existing curvature  values
C                    as much as possible.   Otherwise,  a  value  of y"
C                    entered will be sought in the  abscissa range that
C                    is prompted for next.
C
C             *  Corresp. x range:     Enter low and high  x  values on
C                                      the same line.   Allow  for  the
C                    fact that strict inequalities are  used  when  the
C                    program tests for being within this range.   E.g.:
C                    Enter  .39 .61  or  .39,.61  if you intend for the
C                    constraint to apply in [0.4,0.6].
C
C             *  Width param. for y":  Default is  3.0.   Affects  non-
C                                      uniform weighting  of  the equa-
C                    tions representing 2nd derivative  constraints  in
C                    the overdetermined system being solved.  Since the
C                    actual values of the  2nd derivatives being sought
C                    also act in  a weighting sense,  effects  of  this
C                    variable are not easy to predict. Values of 2.0 or
C                    1.0 should tend to let y" change more.
C
C             *  x for peak y" weight: The absolute values  of  Y"  are
C                                      so much bigger than those  of  Y
C                    that they all need to be scaled down in the system
C                    being solved.   If you are trying  to  flatten the
C                    curvature plot in some region,  pick the center of
C                    the region for this input. Otherwise, use the mid-
C                    chord value.
C
C             *  y" weights, x/c=0,1:  Default is 0.004. See next item.
C
C             *  peak y" weight:       Default is 0.04.  These  provide
C                                      for the fact that  the  absolute
C                    values of y" are  typically  smaller in  the  mid-
C                    section than near the  leading/trailing  edges, so
C                    they should be weighted more,  especially  in view
C                    of the fact that any  y"  constraints applied  are
C                    typically in the mid-section.  See above.
C
C           MODE 6:  Optimize one surface of one profile  using  a pre-
C                    determined set of shape functions, some parameters
C                    of which are automatically varied so as to achieve
C                    a curvature distribution matching some target cur-
C                    vatures in the least squares sense.
C
C           MODE 7:  Transform  representation  of  profile(s)  between
C                    upper/lower surface  and  camber/thickness (either
C                    way - the user is prompted for the direction).  An
C                    option to decamber a section is also provided.
C
C           MODE 8:  Rotate a profile about some point to apply twist.
C
C           MODE 9:  Combine primary profile with a  secondary  profile
C                    (read from a separate file).  This was prompted by
C                    a need to add or remove a boundary layer displace-
C                    ment thickness  (positive by definition)  but  has
C                    been arranged to handle true  addition/subtraction
C                    of distinct profiles as well.
C
C           MODE 10: Loft between primary and secondary profiles.
C
C           MODE 11: Nose-job option: Round or sharpen the leading edge
C                    region via splines.
C
C           MODE 12: Smooth either surface or both surfaces using least
C                    squares techniques (linear combination of n Wagner
C                    functions plus a "ramp" for thick trailing edges),
C                    or by implicit and/or explicit methods involving a
C                    (possibly nonuniform) weighting of y".
C
C
C  GEOMETRY INPUT:
C
C           Standard PROFILE format is shown below. The lower surface
C           is optional, but a zero must be read for NL if no lower
C           surface is included unless this is the last airfoil in the
C           file (meaning EOF can be used to indicate NL=0).  In this
C           case, a symmetrical airfoil is assumed.
C
C               TITLE                   <CHARACTER*80>
C               NU   Upper surface      <Integer # pts., first token>
C               X         Y             <Reals, first two tokens>
C               X         Y                 :
C               :         :                 :     (May be X/C, Y/C;
C               :         :                 :      Xs are increasing)
C               :         :                 :
C               NL   Lower surface      <Integer, first token> <may be 0>
C               X         Y             <Reals, first two tokens>
C               X         Y                 :
C               X         Y  ! Trailing comments are permitted
C               :         :                 :
C               :         :                 :     (Xs are increasing)
C               :         :                 :
C               ! X         Y           <Point suppressed; NL must be adjusted>
C               :         :                 :
C               :         :                 :
C
C    NOTE:  For standard format, if both surfaces are present, PROFILE
C           expects them to have the same leading edge point.  The trailing
C           edge points may differ.
C
C           The next two formats are wrap-around clockwise and wrap-around
C           counterclockwise, where the coordinates begin at the trailing
C           edge, wrap around the leading edge, and end at the trailing edge.
C           The clockwise case begins with the lower surface, and the counter-
C           clockwise case begins with the upper surface.  The format shown
C           below is essentially the same for both cases.  NPTS is the total
C           number of points on the airfoil.
C
C               TITLE                   <CHARACTER*80>
C               NPTS                    <Integer, first token>
C               X         Y             <Reals, first two tokens>
C               X         Y                 :
C               :         :                 :     (May be X/C, Y/C;
C               :         :                 :      Xs are decreasing
C               :         :                 :      until the leading
C               :         :                 :      edge, then increasing)
C
C    NOTE:  Wrap-around formats do NOT have duplicate leading edge points.
C
C           The fourth format is three-column format.  The airfoil has
C           the same abscissas for both surfaces in the 1st column and
C           ordinates for the upper and lower surfaces in the 2nd and 3rd
C           columns respectively.  Abscissas are increasing as with standard
C           format.  Here NPTS is the number of points on either surface.
C
C               TITLE                           <CHARACTER*80>
C               NPTS                            <Integer, first token>
C               X         YU        YL          <Reals, first 3 tokens>
C               X         YU        YL              :
C               :         :         :               :   (May be X/C, Y/C;
C               :         :         :               :    Xs are increasing)
C               :         :         :               :
C               :         :         :               :
C
C
C  CONTROL INPUTS:
C
C           A file containing keyword inputs and values may be used to
C           override default options.  In general, the keywords refer to
C           the airfoil plot file and other output options, and apply
C           to all modes.   Prompts are issued for inputs needed by a
C           particular mode.
C
C  KEYWORD GUIDELINES AND DEFINITIONS:
C
C           Keyword/value pairs may appear with more than one pair on a
C           line.  However, the multivalued keywords PLTLINE, CPSLINE,
C           CRVLINE, and NOFILE must not appear with other keywords on
C           the same line.
C
C           The default value in each case appears in square brackets.
C
C  KEYWORD  VALUES and synonyms     DESCRIPTION
C  -------  -------------------     -----------
C
C  FORMAT   [SAME]                  One of four formats for output profile.
C           PROFILE or STANDARD     May be in standard PROFILE format (ab-
C           CLOCKWISE or WRAPAROUND scissas increasing), clockwise wrap-around
C           COUNTERCLOCKWISE        format, counterclockwise wrap-around for-
C           THREE-COLUMN or         mat, or 3-column format.  SAME means the
C           THREE_COLUMN or         same format as the input profile.  NOTE:
C           THREECOLUMN  or         To allow easily for several synonyms for
C           TABLE                   for the THREE-COLUMN value, only the first
C                                   5 characters of the value are checked.
C
C  PLTLINE  [DEFAULT]               Controls line types of curves on profile
C           LINE                    plots.  One value may be included for
C           DASH                    each curve on the plot.  The default is
C           DOT                     symbols connected by a solid line, with
C           CHAINDASH               a different symbol type for each succes-
C           CHAINDOT                sive curve.  The first curve typically
C           THICK                   represents the original profile; the
C           SYMBOLS                 second curve represents the revised one.
C                                   Overriding the default might be desirable
C                                   when plotting multi-element airfoils or
C                                   when lines without symbols are required.
C                                   At most 20 curves are provided for.  Note:
C                                   All the line types in QPLOT are available.
C                                   SYMBOLS refers to symbols with no line
C                                   connecting them.
C
C  CPSLINE  [see PLTLINE above]     Controls line types on Cps plots in the
C                                   same manner as PLTLINE above.  
C                                   One value
C                                   per curve may be included, chosen from
C                                   the same list of values as those shown
C                                   for PLTLINE.
C
C  CRVLINE  [see PLTLINE above]     Controls line types on curvature plots
C                                   in the same way as PLTLINE and CPSLINE.
C
C  CURVATURE or [NONPARAMETRIC] or  CURVATURE and DERIVATIVES are 
C                                   synonymous
C  DERIVATIVES  [FINITE_DIFFERENCE] controls for the type of calculations
C               SPLINE     or       used for derivatives and hence 
C                                   curvature.
C               PARAMETRIC or       The default is separate-surface 
C                                   treatment
C               WRAPAROUND          using finite differences, as needed for
C                                   consistency with PROFILE's REFINE and
C                                   OPTIMIZE options.  The two surfaces 
C                                   appear
C                                   as separate frames in the curvature plot.
C                                   Otherwise, the full wrap-around 
C                                    curvature
C                                   distribution is calculated using a para-
C                                   metric spline and plotted on a single 
C                                    frame.
C
C                                   The default normally suffices except 
C                                   if the
C                                   region of interest is very near a 
C                                   rounded
C                                   leading edge.  Note that not all of the
C                                   possibilities are provided for, such as
C                                   parametric finite differences.
C
C  MINCURVATURE   [-5.]             Cutoff values for plotted curvatures.
C  MAXCURVATURE   [+5.]             Practice shows that +/-5. give useful
C                                   plot scaling by ignoring the high curv-
C                                   ature values near the leading edge.  On
C                                   the other hand, it may well be desired
C                                   to focus on the leading edge region.  
C                                   Set
C                                   both to 999. to obtain the full range.
C                                   See the CURVATURE/DERIVATIVES control.
C
C  NOFILE   [NONE]                  Used to suppress any combination of the
C           DAT                     seven output files generated by PROFILE.
C           PLT                     The values correspond to the extensions
C           TAB                     of the file names.  See elsewhere for a
C           CRV                     complete description of file contents.
C           YPP                     NONE serves only to assist leaving the
C           CPS                     NOFILE control word in the input file
C           SPREAD                  even if all outputs are desired.
C
C  PLOT     [BOTH]                  Controls plotting of original OR revised
C           ORIGINAL                profile.  The default is to plot both
C           REVISED                 original and revised (if one exists).
C
C  PRECISION   [FULL]               Controls number of digits in output
C           ENGINEERING             airfoil coordinates.  FULL gives F11.8
C                                   if possible, or E15.8 if any X >=10.
C                                   ENGINEERING gives the traditional F10.6
C                                   common to many flow solvers.
C
C  THREED   [FALSE] or [NO]         For plotting of multiple stations from
C           TRUE or YES             a 3-D wing. The default is the 2-D case.
C
C  XAXIS    [6.4]                   Length of x-axis in inches.  The default
C                                   is normally appropriate for an 8.5 x 11
C                                   page in portrait mode.
C
C           The following four keywords apply to windowing.  Any or none
C           of them may be used.
C
C  XMIN     [minima and             Minimum abscissa for desired window
C  XMAX     maxima of               Maximum abscissa for desired window
C  YMIN     the input               Minimum ordinate for desired window
C  YMAX     coordinates]            Maximum ordinate for desired window
C
C
C  SAMPLE CONTROL FILE:
C
C           A sample input file follows.  Note that keywords and values
C           may be separated with blanks, commas, colons, equals signs,
C           pr tabs. Remember, keywords with more than one value should
C           appear on separate lines.  Any keyword or text value may be
C           truncated to unambiguous leading characters.   Blank  lines
C           and trailing ! comments are ignored.
C
C
C           FORMAT = STANDARD   PRECISION = FULL
C           PLOT BOTH  THREED:NO
C           PLTLINE = SOLID, SOLID
C           CPSLINE = DOT, SYMBOLS
C           CRVLINE = SOLID, DASH, CHAINDOT
C           XAXIS = 20.
C           XMIN = 0.  XMAX 0.1
C           MAXCURVATURE = 999.   ! Both 999. means plot the full
C           MINCURVATURE = 999.   ! curvature range
C           DERIVATIVES = PARAMETRIC
C           NOFILE: YPP SPREAD
C
C
C  OUTPUT FILES:
C
C           The following seven output files are generated by  PROFILE.
C           Any  of  the  files  may  be  suppressed  using the keyword
C           NOFILE.  The user is prompted for an identifier, which will
C           become  the first part of each file name.   File extensions
C           are fixed.
C
C  <identifier>.DAT   Contains airfoil coordinates  that  have been re-
C                     vised in some way. May be in one of four formats;
C                     default is the same format as input coordinates.
C
C  <identifier>.PLT   Contains airfoil geometry coordinates  and  other
C                     information necessary for later QPLOTing.  May be
C                     a plot of original profile,  revised profile,  or
C                     both.  (Default is both, superimposed.)
C
C  <identifier>.TAB   Contains tabulated coordinates,  first and second
C                     derivatives and curvatures for the  original  and
C                     revised profile  (if  a  revised profile exists).
C                     Other diagnostics may be written here,  including
C                     a record of selected  shape  functions  from  the
C                     MODIFY option.
C
C  <identifier>.CRV   Contains curvatures of the  original  and revised
C                     profiles (if a revised profile exists)  for later
C                     QPLOTing.  Also to be used as the basis of target
C                     curvature data when OPTIMIZE option is used.  The
C                     MAXCURVATURE and MINCURVATURE keywords  determine
C                     the plot axis range, but the full surfaces appear
C                     in the file (except for the first and last pts.).
C                     See the description of these keywords for how  to
C                     obtain a full wrap-around curvature distribution.
C
C  <identifier>.YPP   Contains second derivatives  in  standard PROFILE
C                     format,  for  possible  reuse  in  "refine" mode.
C                     When a profile has been  revised,  the  file will
C                     contain  only  the  second  derivatives  of   the
C                     revised profile.   Otherwise,  the file will con-
C                     tain only the second derivatives of the  original
C                     profile.
C
C  <identifier>.CPS   Contains QPLOTable estimates of Cps for user-sup-
C                     plied alpha and free stream Mach number  (revised
C                     and/or original).
C
C  <ident>.SPREAD     Contains  spreadsheet-compatible  (tab-delimited)
C                     tabular data.   Available only if MODE <= 3,  and
C                     only if both surfaces have common abscissas.
C
C
C  LOGICAL UNIT NUMBERS:
C
C    LUNCPS      For output Cp estimates in QPLOT format.
C    LUNCRT      For prompts and diagnostics.
C    LUNCRV      For output curvature data in QPLOT format.
C    LUNDAT      For input file of geometry data.
C    LUNINP      For keyword control inputs (disk file).
C    LUNKBD      For responses to prompts.
C    LUNPLT      For airfoil geometry in QPLOT format.
C    LUNREV      For output file of revised geometry data.
C    LUNSPR      For spreadsheet-compatible output.
C    LUNTAB      For tabulations and diagnostics.
C    LUNTBL      For optional table of y" data ("refine" mode; PROFILE format).
C                Alternatively, target curvature distribution ("optimize" mode,
C                QPLOT format).
C    LUNXGR      For file for reading abscissas ("redistribute" mode),
C                or for file containing bump function info. ("optimize" mode,
C                keyword format), or for secondary profile ("combine" mode).
C    LUNYPP      For output y" data in PROFILE format.
C
C
C  EXTERNAL REFERENCES:
C
C    AFGEOM      Calculates various airfoil geometric properties
C    BOUNDS      Finds the max/min values in array(s)
C    COMBINE     Combines primary and secondary profiles
C    COPY,RVERSE Utilities for transferring data
C    CURV2D      Curvature via parametric derivatives
C    FD12K       Finite difference approx. to derivatives/curvature
C    GETCL       Computes aerodynamic coefficients
C    GETCPS      Cheap estimate of Cp distributions
C    LOFT        Lofts linearly between two profiles
C    LSTFIL      Copies a formatted file to screen or disk
C    MAXMIN      Finds maximum and minimum ordinates and abscissas
C                over upper and lower surfaces of all profiles
C    MODIFY      Interactive subroutine to modify  the  geometry  by
C                applying "bump" functions
C    NOSEJOB     Modifies the leading edge region in various ways
C    NRMLIZ      Normalizes coordinates
C    NRMSET      Prompts for the details needed by NRMLIZ
C    ROTATE      Rotates coordinates with option to renormalize
C    OPENER      File opening utility
C    OPTIMIZE    Optimizes one surface using bumps/target curvatures
C    PROTECT     Checks for duplicate points and monotonicity
C    PRREAD      Reads one airfoil profile
C    PRTAB       Tabulates coordinates, derivatives,  and curvatures
C                (which it calculates from the given derivatives)
C    PRWRIT      Saves the profile geometry data on disk in the format
C                described in PRREAD
C    PSFIT       Parametric spline fit and evaluation, used here for
C    PSTVAL      wrap-around curvature calculations.
C    QPLDAT      Entry points QPLFRM and QPLCRV permit saving of
C                QPLOTable data - used here for 2nd derivatives and Cps
C    RDKEYS      Reads keyword input control file
C    READER      Prompting utility
C    RECTIFY     Reorganizes the input geometry so that  the  common
C                leading edge point has in fact the minimum abscissa.
C                Also shifts ordinates by a user-specified amount.
C    REDISTRIB   Redistributes data points using splines
C    REFINE      Refines thickness and/or curvature distribution(s)
C    SCAN2       Scans a string for the next token
C    SELECT      Menu selection by item number or by name
C    TRANSFORM   Transforms representation of one profile, or decambers it.
C    XDERIVS     Forms X derivatives from parametric derivatives
C    XFORM       Called by TRANSFORM, but also here for spreadsheet info.
C
C  ERROR HANDLING:
C    This is mostly confined to reading the input files.  See subroutines
C    PRREAD and RDKEYS.
C
C  ENVIRONMENT:
C    DEC VMS; SGI IRIX; Fortran 90
C
C  HISTORY:
C
C  Dec. '82  DAS/RAK  Original design.
C  12/09/82    LJC    Original coding (plot/redistribute data points).
C  04/29/83    LJC    Added option to "rectify" input data points.
C  July '83    LJC    Added a 3-D capability and an option for reading
C                     a file of new abscissas in "redistribute"  mode.
C  09/27/83    LJC    Interactive MODIFY routine incorporated.
C  Oct. '83    DAS    Integrated alternative version of MODIFY  as the
C                     REFINE option; provided for saving curvature and
C                     y" values; removed REDISTRIB function from main.
C  11/09/83    LJC    Added de-normalizing to normalizing option;  in-
C                     cluded these as a separate MODE; reordered MODEs
C                     in order of complexity.
C  01/20/84    DAS    Incorporated OPTIMIZE mode.
C  04/12/84    DAS    Took advantage of QPLOT's new legend capability.
C  July '84    LJC    Added reading and writing of wraparound formats.
C  07/24/84    LJC    Added calculation of thickness for all modes.
C  Aug. '84    LJC    Changed from namelist to keyword inputs.
C  09/14/84    LJC    Changed legend entry to be read from dataset and
C                     added prompt for title.  Formerly, the title was
C                     read from dataset and the legend was hard-coded.
C  Oct. '84    LJC    Arranged for all plotting to be done  outside of
C                     PROFILE  using  QPLOT.  Formerly,  much  of  the
C                     program  was  devoted  to  plotting  the airfoil
C                     geometry with DISSPLA.
C  Dec. '84    DAS    Incorporated cheap estimates of Cp distributions
C                     using algorithm supplied by Ilan Kroo of the RAC
C                     Branch at NASA Ames.  Added Cl, Cm calculation.
C  01/24/85    LJC    Allowed for original and revised plots with MODE
C                     = 3.  Also modified to take advantage of QPLOT's
C                     new equal axis scaling option.
C  02/12/84    DAS    Added TRANSFORM option.
C  02/13/85    LJC    Added shifting of ordinates to RECTIFY mode.
C  02/19/85    DAS    Fixed REDISTRIB to handle rounded trailing edges
C                     properly; took out "BLUNT" input parameter.
C  02/28/85    LJC    Added 3-column format to PRREAD and PRWRIT.
C  03/22/85    LJC    Allowed for absent  input control file  (meaning
C                     an empty file is not needed). Also added two new
C                     keywords to RDKEYS for controlling line types on
C                     Cps and curvature plots.
C  06/17/85    DAS    Provided generalized distributions in REDISTRIB,
C                     via GETDIS/DSTRIB in place of XGRID.
C  09/05/85    DAS    Fixed minor error in QPLOTable airfoil data file.
C  09/30/85    DAS    Added COMBINE option.
C  10/09/85    DAS    Suppressed plot windowing control values from
C                     QPLOTable file if they are not in use.
C  10/21/85    DAS    Introduced LSTFIL to echo the input control file
C                     to the .TAB file for future reference.
C  11/04/85    DAS    Mixup with END CURVE/END FRAME for OPTIMIZE mode
C                     target curvatures crept in (how? used to be OK).
C  12/30/85    DAS    Added ROUND option (tentatively).
C  04/24/86    DAS    Cp at leading edge was wrong.
C  08/11/86    DAS    Added SMOOTH option (fitting of Wagner functions
C                     in linear least squares sense).
C  09/24/86    DAS    AFGEOM in place of CALCTHICK in main program;
C                     menu options by name now, as well as by number.
C  10/20/86    RAK    AFGEOM returns CAMBER and THICK as fractions; must
C                     convert to % chord at calling level.
C  02/11/87    DAS    Functionality of BOUNDS changed somewhat.
C  03/11/87    DAS    Made use of PROTECT utility (prompted by duplicate
C                     point case that dies ungracefully otherwise).
C  04/24/87    DAS    Menu items clarified and made more unique (because
C                     of SELECT's option for choosing by name or number).
C                     Greater precision in output coordinates now by default.
C                     Also: traps bad IDENT now (more than one token is
C                     presumed to be a goof).
C  04/27/87    DAS    Introduced PRECISION argument in PRWRIT.
C  08/31/87    DAS    Trailing edge angle added as AFGEOM output.
C  09/23/87    DAS    If MODE=0 but PRECISION is not the default, assume
C                     that a revised airfoil dataset is required.
C  04/29/88    DAS    Added spreadsheet-compatible output file.
C  12/01/88    DAS    Turned on spreadsheet file for MODE <= 3 now.
C  11/29/89  DAS/RGL  Plots for multi-section 3D cases with blank subtitle from
C                     the initial prompt now use the case titles as subtitles.
C                     'Original' suppressed from headers if 'Revised' does not
C                     apply.
C  02/14/90    DAS    Installed OPENER in 9 places.  File names are now in
C                     lower case to indulge the Unix community.
C  03/13/90     "     Raised MXPTS from 201 to 300.
C  06/28/90     "     Added the LOFT option.
C  06/29/90     "     Added the ROTATE option.
C  04/22/91     "     Indicated "current axes" for thickness/camber printout.
C  10/22/91     "     Replaced the ROUND option with the "nose-job" option
C                     (round or sharpen).  Providing for reversing the operation
C                     can mean letting through a non-rectified airfoil.
C  10/24/91     "     Provided for full wrap-around curvature plot (parametric
C                     derivatives).  Adjusted the tabulations accordingly.
C                     Added the control scheme description above (adapted from
C                     RDKEYS) and the geometry data description (adapted from
C                     PRREAD).
C  01/07/95     "     Original chord and leading edge are now tabulated;
C                     "normalize" mode *.crv plot file had mismatched X min/max.
C  06/13/95     "     Output of the *.ypp file is now as advertised above.
C  11/21/95     "     TRANSFORM mode now has an option to decamber a section.
C  12/19/96     "     Added SINF and four "COS" shape functions.
C  12/20/96     "     Extended SMOOTH option to allow y"-based implicit and/or
C                     explicit smoothing as an alternative to Wagner fitting.
C  10/20/99     "     Fortran 90 upgrade, mainly to eliminate use of '0' for
C                     carriage control.
C
C  AUTHORS:
C    Leslie Collins, David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C  ACKNOWLEDGMENTS:
C    Robert Kennelly (ex-Sterling, now NASA Ames) provided the inspiration and
C    many key ideas.
C
C    The Aerodynamics Division at NASA Ames funded this software under contract
C    to Sterling Software (previously known as Informatics, Inc.).
C
C-------------------------------------------------------------------------------


C     Declarations:
C     -------------

      IMPLICIT NONE

C  *  Constants:

      REAL, PARAMETER ::
     >   HALF = 0.5E+0, ONE = 1.E+0, UNDEF = 999.E+0, ZERO = 0.E+0

      CHARACTER, PARAMETER ::
     >   BLANK * 1      = ' ',
     >   CURVATURE * 10 = ' curvature',
     >   ENDCRV * 9     = 'END CURVE',
     >   ENDFRM * 9     = 'END FRAME',
     >   LO * 6         = ' lower',
     >   NEW * 3        = 'NEW',
     >   OLD * 3        = 'OLD',
     >   ORIGINAL * 9   = 'Original ',
     >   REVISED * 9    = 'Revised  ',
     >   UP * 6         = ' upper'

      INTEGER, PARAMETER ::
     >   LUNCPS=9, LUNCRT=6, LUNCRV=7,  LUNDAT=2, LUNINP=1,  LUNKBD=5,
     >   LUNPLT=8, LUNREV=3, LUNSPR=13, LUNTAB=4, LUNTBL=10, LUNXGR=11,
     >   LUNYPP=12, MAXPTS=500, MXCRVS=20, MXMENU=12, MXWRAP=2*MAXPTS+1

C     MAXPTS = Maximum number of points allowed for on one surface.
C     MXCRVS = Maximum number of curves allowed for on one plot.
C     MXMENU = Number of options on the main menu.
C     MXWRAP = 2*MAXPTS+1 (for estimation of force coefficients).

C  *  Variables:

      INTEGER
     >   CAMBCASE, FIRST, FORMAT, I, IER, INFORM, IPR, J, LAST, LUNIN,
     >   MARK, MODE, NINC, NL, NLORIG, NOSEMODE, NPR, NPTS, NU, NUORIG,
     >   NWRAP, NXTARG, PRECISION

      REAL
     >   ALFRAD, ALPHA, AREA, ARROW, CAMBER, CD, CENTRD, CHORD, CNORML,
     >   COEFS (MAXPTS, 3), CL, CM, CPB, CPT, CPBOT, CPTOP, CPSTAR,
     >   CPL (MAXPTS), CPU (MAXPTS), CRVMAX, CRVMIN, CRVNEW (MXWRAP),
     >   CRVORIG (MXWRAP), DUMMY(1), FSMACH, GAMMA, GOGM1,              ! modified RLC
     >   GP1GM1, MOMENT, PWRAP (MXWRAP), SCALE, TARGCRV (MAXPTS),
     >   TARGX (MAXPTS), TEANGLE, THICK, TOGM1,
     >   X (MAXPTS*2), XAXIS, XC (MAXPTS), XCAM,
     >   XEVAL (MAXPTS), XL (MAXPTS), XLE, XLEFT, XLORIG (MAXPTS),
     >   XMAXSV, XMAX, XMIN, XMINSV, XP (MXWRAP), XPP (MXWRAP), XRIGHT,
     >   XTH, XU (MAXPTS), XUORIG (MAXPTS), XWRAP (MXWRAP),
     >   Y (MAXPTS*2), YEVAL (MAXPTS), YKL (MAXPTS), YKLORG (MAXPTS),
     >   YKU (MAXPTS), YKUORG (MAXPTS), YL (MAXPTS), YLE,
     >   YLORIG (MAXPTS), YMAX, YMAXSV, YMIN, YMINSV, YP (MXWRAP),
     >   YPL (MAXPTS), YPP (MXWRAP), YPPL (MAXPTS), YPPU (MAXPTS),
     >   YPU (MAXPTS), YU (MAXPTS), YUORIG (MAXPTS), YWRAP (MXWRAP)

      LOGICAL
     >   CPSFIL, CRVFIL, DATFIL, DEFAULT, DISTINCT, NOMODS, NOSHOW,
     >   PLTARG, PLTFIL, PLTORIG, PLTREV, QUIT, SAVREV, SPREADFIL,
     >   TABFIL, THREED, UPPER, WRAPCRV, YPPFIL

      CHARACTER
     >   CHOICE * 16, CPSLINE (MXCRVS) * 9, CRVLINE (MXCRVS) * 9,
     >   DATAFILE * 40, ENDSTR * 9, FSTATUS * 9, HTAB * 1, IDENT * 15,
     >   INPUTFILE * 32, LEG * 8, LEGEND * 80, LEGND (MXCRVS) * 80,
     >   MENU (0:MXMENU) * 60, MNTITL * 80, PLTLINE (MXCRVS) * 9,
     >   SBTITL * 80, SBTIT2 * 80, SUBTIT (2) * 24, SURFACE (2) * 13,
     >   XLABEL * 3, YLABEL * 3

C  *  Procedures:

      EXTERNAL
     >   AFGEOM, BOUNDS, COMBINE, COPY, CURV2D, FD12K, GETCL, GETCPS,
     >   LOFT, LSTFIL, MAXMIN, MODIFY, NOSEJOB, NRMLIZ, NRMSET, OPENER,
     >   OPTIMIZE, PROTECT, PRREAD, PRTAB, PRWRIT, PSFIT, PSTVAL,
     >   QPLDAT, RDKEYS, READR, READS, READY, RECTIFY, REDISTRIB,
     >   REFINE, RVERSE, ROTATE, SCAN2, SELECT, SMOOTH, TRANSFORM,
     >   XDERIVS, XFORM

C  *  Storage:

      DATA
     >   SURFACE /'upper surface', 'lower surface'/

      DATA
     >   MENU /
     >   '  0 = DISPLAY (plot/tabulate only or alter format/precision)',
     >   '  1 = RECTIFY leading edge definition; allows vertical shift',
     >   '  2 = NORMALIZE or denormalize coordinates',
     >   '  3 = REDISTRIBUTE the abscissas',
     >   '  4 = MODIFY either surface or both (apply shape functions)',
     >   '  5 = REFINE thickness/curvature',
     >   '  6 = OPTIMIZE one surface ("bumps" + target curvatures)',
     >   '  7 = TRANSFORM YU/YL to/from camber/thickness, or decamber',
     >   '  8 = ROTATE coordinates, with option to renormalize',
     >   '  9 = COMBINE option (add or subtract profiles)',
     >   ' 10 = LOFT linearly between primary and secondary profiles',
     >   ' 11 = NOSE-JOB option: round or sharpen the leading edge',
     >   ' 12 = SMOOTH YU and/or YL: Wagner fn. fits or [im|ex]plicit' /


C     Execution:
C     ----------

      HTAB = CHAR (9)            ! Used to be a PARAMETER (non-standard).

C  *  Select mode of operation:

      WRITE (LUNCRT, 1000)
     >   BLANK,
     >   'Welcome to PROFILE from the Aerodynamics Division, NASA Ames.'

      MODE = 0
      CHOICE = 'DISPLAY'
      NOSHOW = .FALSE.

      CALL SELECT ('Select operating mode.', MXMENU+1, MENU, NOSHOW,
     >             LUNCRT, LUNKBD, MODE, CHOICE, QUIT)
      IF (QUIT) GO TO 810

C-----------------------------------------------------------------------
C         Handle input files (control file; airfoil data file).
C-----------------------------------------------------------------------

      DATAFILE = 'naca0012.dat'
      CALL OPENER (LUNCRT,
     >   'Input airfoil file? <CR>=naca0012.dat: ',
     >   LUNKBD, DATAFILE, LUNDAT, OLD)

      INPUTFILE = 'profile.inp'
      FSTATUS = 'IfPresent'
      LUNIN = LUNINP
      CALL OPENER (LUNCRT,
     >   'Input control file? <CR>=profile.inp or none: ',
     >   LUNKBD, INPUTFILE, LUNIN, FSTATUS)
      IF (FSTATUS == 'MISSING') THEN
         LUNIN = -LUNIN
         INPUTFILE = BLANK
         WRITE (LUNCRT, 1000)
     >      'No control file found - proceeding ...'
      END IF

C  *  Read input control variables (if any):

      CALL RDKEYS (LUNCRT, LUNIN, MXCRVS, UNDEF, FORMAT, PRECISION,
     >   PLTLINE, CPSLINE, CRVLINE, CRVMAX, CRVMIN, WRAPCRV,
     >   PLTORIG, PLTREV, THREED, XAXIS, XMIN, XMAX, YMIN, YMAX,
     >   DATFIL, PLTFIL, CRVFIL, TABFIL, YPPFIL, CPSFIL, SPREADFIL)

C-----------------------------------------------------------------------
C             Set up or adjust the internal control variables.
C-----------------------------------------------------------------------

      TABFIL  = TABFIL .OR.  MODE >= 4
      PLTORIG = PLTFIL .AND. PLTORIG .AND. (MODE == 0 .OR. MODE >= 3)
      PLTREV  = PLTFIL .AND. PLTREV .AND. MODE > 0
      PLTFIL  = PLTFIL .AND. (PLTORIG .OR. PLTREV)
      PLTARG  = CRVFIL .AND. MODE == 6
      CRVFIL  = CRVFIL .AND. MODE /= 7
      CPSFIL  = CPSFIL .AND. MODE /= 7
      SPREADFIL = SPREADFIL .AND. MODE <= 3
      SAVREV  = DATFIL .AND.
     >          (MODE > 0 .OR. FORMAT > 0 .OR. PRECISION /= 1)
      NOMODS  = MODE == 0 .OR. MODE == 7

      IF (WRAPCRV .AND. YPPFIL) THEN
         WRITE (LUNCRT, 1080)
         YPPFIL = .FALSE.
      END IF

      CAMBCASE = 1  ! Transform from Ys to camber/thickness

C  *  CAMBCASE applies to MODE = 7, but having it undefined for other
C     modes is inconvenient for deciding whether to invoke AFGEOM.


C-----------------------------------------------------------------------
C                        Set up the output files.
C-----------------------------------------------------------------------

C  *  Prompt for the identifier to be used for all files.  Figure that
C     more than one token, or something like 'naca0012.dat' is a goof.

  150 CONTINUE
      IDENT = 'profile'
      CALL READS (LUNCRT,
     >   'Identifier for output files? <CR>="profile": ',
     >   LUNKBD, IDENT, DEFAULT, QUIT)
      IF (QUIT) GO TO 810

      FIRST = 1
      LAST = LEN (IDENT)
      CALL SCAN2 (IDENT, ' .;/', FIRST, LAST, MARK)
      IF (LAST /= MARK) GO TO 150

C     Valid identifier may still cause a conflict with the input file.
C     (Unix systems in mind, here.)

      IF (SAVREV) THEN
         IF (IDENT (FIRST:LAST) == DATAFILE (FIRST:LAST)) THEN
            WRITE (LUNCRT, 1000)
     >         'In/out file name conflict.  Choose another identifier.'
            GO TO 150
         END IF
      END IF

      IF (SAVREV)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.dat', LUNREV, NEW)
      IF (PLTFIL)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.plt', LUNPLT, NEW)
      IF (TABFIL)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.tab', LUNTAB, NEW)
      IF (CRVFIL)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.crv', LUNCRV, NEW)
      IF (YPPFIL)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.ypp', LUNYPP, NEW)
      IF (CPSFIL)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.cps', LUNCPS, NEW)
      IF (SPREADFIL)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.spread', LUNSPR,
     >                NEW // ':153') ! Raise max. record length

C  *  Prompt for plot title and subtitle:

      MNTITL = 'PROFILE'
      CALL READS (LUNCRT,
     >   'Plot title line?  <CR> uses "PROFILE" for the title:',
     >   LUNKBD, MNTITL, DEFAULT, QUIT)
      IF (QUIT) GO TO 810

      SBTITL = BLANK
      CALL READS (LUNCRT,
     >   'Plot subtitle line?  <CR> means none:',
     >   LUNKBD, SBTITL, DEFAULT, QUIT)
      IF (QUIT) GO TO 810

C-----------------------------------------------------------------------
C                            Scan all datasets.
C-----------------------------------------------------------------------

C  *  Scan all profiles for overall data range (for axis scaling and
C     possible normalization purposes).  Count number of profiles.

      CALL MAXMIN (LUNDAT, MAXPTS, XMINSV, XMAXSV, YMINSV, YMAXSV,
     >             YLE, NPR, X, Y, XU, XL, YU, YL, LEGND (1), IER)

      IF (IER == 4) GO TO 900
      IF (IER == 3) GO TO 910
      IF (IER == 2) GO TO 920
      IF (NPR == 0) GO TO 930

C  *  Make sure both surfaces have the same leading edge point, unless
C     this is the rectify option:

      IF ((MODE /= 1) .AND.
     >    (XU (1) /= XL (1) .OR. YU (1) /= YL (1))) GO TO 850

      REWIND LUNDAT

      IF (TABFIL) THEN

         WRITE (LUNTAB, 1010)
     >      ' Program PROFILE (Applied Aerodynamics Branch, NASA Ames',
     >      ' Research Center)'
         WRITE (LUNTAB, 1020) ' Operating mode: ', MODE, ' - ', CHOICE
         WRITE (LUNTAB, 1010) ' Control file  : ', INPUTFILE
         WRITE (LUNTAB, 1010) ' Geometry file : ', DATAFILE
         WRITE (LUNTAB, 1010) ' Output file id: ', IDENT (FIRST:LAST)
         WRITE (LUNTAB, 1025) ' Profiles found: ', NPR

         IF (INPUTFILE /= BLANK) THEN

C  *        Echo control keywords to tabulation file.
C           LEGEND is as good as anything for use as a buffer.

            WRITE (LUNTAB, 1000) BLANK, 'Control file contents:', BLANK

            REWIND LUNINP
            CALL LSTFIL (LUNINP, LUNTAB, LEGEND)

         ELSE
            WRITE (LUNTAB, 1020)
     >         ' Control file was absent - all defaults taken.'
         END IF

      END IF

      IF (CPSFIL) THEN

C  *     Don't attempt to estimate Cps for more than one profile:

         CPSFIL = NPR == 1 .AND. .NOT. THREED

C  *     Save user from editing input control file too often:

         IF (CPSFIL) THEN
            CPSFIL = .FALSE.
            CALL READY (LUNCRT,
     >         'Do you really want Cp estimates? (Y/N; <CR>=No) ',
     >         LUNKBD, CPSFIL, DEFAULT, QUIT)
         END IF

         IF (.NOT. CPSFIL) THEN
            CLOSE (UNIT=LUNCPS, STATUS='DELETE')
         ELSE
            ALPHA = 0.
            CALL READR (LUNCRT,
     >         'Enter Alpha, else <CR> gives 0.: ',
     >         LUNKBD, ALPHA, DEFAULT, QUIT)
            IF (QUIT) GO TO 810

            FSMACH = 0.
            CALL READR (LUNCRT,
     >         'Enter free stream Mach number, else <CR> gives 0.: ',
     >         LUNKBD, FSMACH, DEFAULT, QUIT)
            IF (QUIT) GO TO 810

         END IF
      END IF

C-----------------------------------------------------------------------
C                     Begin main loop over all profiles:
C-----------------------------------------------------------------------

      DO IPR = 1, NPR

C  *     Read one profile.  Ignore error processing this time through.

         CALL PRREAD (LUNDAT, LEGND (IPR), MAXPTS, NU, NL, X, Y,
     >                XU, XL, YU, YL, INFORM, IER)

C  *     Normally no use in proceeding if abscissas are not monotonic:

         CALL PROTECT (NU, XU, YU, ARROW, DISTINCT)
         J = 1
         IF (ARROW /= ONE) THEN
            IF (MODE == 11) THEN     ! May need to reverse a round/sharpen
               WRITE (LUNCRT, 1070)
            ELSE IF (MODE /= 1) THEN
               GO TO 860
            END IF
         END IF

         CALL PROTECT (NL, XL, YL, ARROW, DISTINCT)
         J = 2
         IF (ARROW /= ONE) THEN
            IF (MODE == 11) THEN
               WRITE (LUNCRT, 1070)
            ELSE IF (MODE /= 1) THEN
               GO TO 860
            END IF
         END IF

C  *     Save original data before any modifications:

         CALL COPY (NU, XU, XUORIG)
         CALL COPY (NL, XL, XLORIG)
         CALL COPY (NU, YU, YUORIG)
         CALL COPY (NL, YL, YLORIG)
         NUORIG = NU
         NLORIG = NL

C-----------------------------------------------------------------------
C                           Major options.
C-----------------------------------------------------------------------

         IF (MODE == 1) THEN

C  *        Rearrange upper and lower surface points so that the
C           minimum abscissa is the leading edge point. Shift ordinates
C           by input value if required:

            CALL RECTIFY (NU, XU, YU, NL, XL, YL, X, Y, MAXPTS,
     >                    LUNCRT, LUNKBD)

         ELSE IF (MODE == 2) THEN

            IF (IPR == 1) THEN

C  *           Prompt for chord, etc., to normalize/denormalize by:

               CALL NRMSET (LUNCRT, LUNKBD, XMINSV, XMAXSV,
     >                      CNORML, XLE, YLE)
            END IF

            CALL NRMLIZ (NU, XU, XU, XLE, CNORML)
            CALL NRMLIZ (NL, XL, XL, XLE, CNORML)
            CALL NRMLIZ (NU, YU, YU, YLE, CNORML)
            CALL NRMLIZ (NL, YL, YL, YLE, CNORML)

         ELSE IF (MODE == 3) THEN

C  *        Redistribute abscissas using splines:

            CALL REDISTRIB (NU, XU, YU, NL, XL, YL, X, Y,
     >                      MAXPTS, XEVAL, YEVAL, COEFS,
     >                      LUNCRT, LUNKBD, LUNXGR)

         ELSE IF (MODE == 4) THEN

C  *        Modify profile by applying shape functions interactively:

            CALL MODIFY (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD,
     >                   LUNTAB)

         ELSE IF (MODE == 5) THEN

C  *        Refine thickness and/or curvature interactively:

            CALL REFINE (NU, XU, YU, NL, XL, YL, MNTITL,
     >                   LUNCRT, LUNKBD, LUNTAB, LUNTBL)

         ELSE IF (MODE == 6) THEN

C  *        Optimize one surface using bumps + target curvatures:

            CALL OPTIMIZE (NU, XU, YU, NL, XL, YL, MNTITL,
     >                     LUNCRT, LUNKBD, LUNTAB, LUNXGR, LUNTBL,
     >                     UPPER, NXTARG, TARGX, TARGCRV)

         ELSE IF (MODE == 7) THEN

C  *        Transform representation (upper/lower <-> camber/thickness)
C           or zero out the camber.  CAMBCASE is an output.

            CALL TRANSFORM (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD,
     >                      CAMBCASE)

         ELSE IF (MODE == 8) THEN

C  *        Rotate coordinates, with option to renormalize:

            CALL ROTATE (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD)

         ELSE IF (MODE == 9) THEN

C  *        Combine present profile with another to be prompted for.
C           Use YPU, YPPU as scratch for XU2ND, YU2ND, etc.

            CALL COMBINE (NU, XU, YU, NL, XL, YL, MAXPTS, X, Y,
     >                    YPU, YPPU, YPL, YPPL,
     >                    LUNCRT, LUNKBD, LUNXGR)

         ELSE IF (MODE == 10) THEN

C  *        Loft linearly between present profile and another to be prompted
C           for.  Use YPU, YPPU as scratch for XU2ND, YU2ND, etc.

            CALL LOFT (NU, XU, YU, NL, XL, YL, MAXPTS, X, Y,
     >                 YPU, YPPU, YPL, YPPL,
     >                 LUNCRT, LUNKBD, LUNXGR)

         ELSE IF (MODE == 11) THEN

C  *        Round or sharpen the leading edge region:

            CALL NOSEJOB (NU, XU, YU, NL, XL, YL, MAXPTS, X, Y, COEFS,
     >                    LUNCRT, LUNKBD, LUNTAB, NOSEMODE)

         ELSE IF (MODE == 12) THEN

C  *        Smooth YU and/or YL by fitting Wagner shape functions or
C           via implicit and/or explicit y"-based smoothing:

            CALL SMOOTH (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD, LUNTAB)

         END IF


C-----------------------------------------------------------------------
C                         Generate the output files.
C-----------------------------------------------------------------------


C  . . . . . . . . . Revised profile coordinates (*.DAT) . . . . . . . .

         IF (SAVREV) THEN

            IF (MODE == 7) THEN
               IF (CAMBCASE == 1) THEN
                  FORMAT = 6
               END IF
            END IF

            IF (FORMAT == 0) FORMAT = INFORM

            IF (MODE == 0) THEN
               CALL PRWRIT (LUNREV, MAXPTS, LEGND(IPR), NU, NL,
     >                      XU, XL, YU, YL, FORMAT, PRECISION)
            ELSE
               CALL PRWRIT (LUNREV, MAXPTS, MNTITL, NU, NL,
     >                      XU, XL, YU, YL, FORMAT, PRECISION)
            END IF

         END IF


C  . . . . . . . . . .  Plottable coordinates (*.PLT) . . . . . . . . . .

C  *     Multi-dataset plot frames should be self-descriptive:

         IF (NPR > 1 .AND. SBTITL == BLANK) THEN
            SBTIT2 = LEGND (IPR)
         ELSE
            SBTIT2 = SBTITL
         END IF

         IF (ABS (XU (1)  - ZERO) > 1.E-3  .OR.
     >       ABS (XU (NU) - ONE)  > 1.E-3) THEN
            XLABEL = ' x '
            YLABEL = ' y '
         ELSE
            XLABEL = 'x/c'
            YLABEL = 'y/c'
         END IF

         IF (PLTORIG) THEN

C  *        Save original airfoil coordinates in QPLOTable form. First store
C           coordinates in wrap-around format, avoiding duplicate leading edge
C           point:

            CALL RVERSE (NLORIG, XLORIG, X)
            CALL RVERSE (NLORIG, YLORIG, Y)
            CALL COPY   (NUORIG, XUORIG, X (NLORIG))
            CALL COPY   (NUORIG, YUORIG, Y (NLORIG))
            NPTS = NLORIG + NUORIG - 1

            IF (IPR == 1 .OR. THREED)
     >         CALL QPLFRM  (LUNPLT, MNTITL, SBTIT2, XLABEL, YLABEL)

            IF (.NOT. PLTREV  .AND.  IPR == NPR  .OR.  THREED) THEN
               ENDSTR = ENDFRM
            ELSE
               ENDSTR = ENDCRV
            END IF

            CALL QPLCRV (LUNPLT, NPTS, X, Y, ENDSTR, UNDEF,
     >                   'SCALE', XMIN, XMAX, YMIN, YMAX, XAXIS, UNDEF,
     >                   LEGND (IPR), PLTLINE (IPR))

         END IF

         IF (PLTREV) THEN

C  *        Now for the revised airfoil:

            CALL RVERSE (NL, XL, X)
            CALL RVERSE (NL, YL, Y)
            CALL COPY   (NU, XU, X (NL))
            CALL COPY   (NU, YU, Y (NL))
            NPTS = NL + NU - 1

            IF (.NOT. PLTORIG)
     >         CALL QPLFRM (LUNPLT, MNTITL, SBTIT2, XLABEL, YLABEL)

            CALL QPLCRV (LUNPLT, NPTS, X, Y, ENDFRM, UNDEF, 'SCALE',
     >                   XMIN, XMAX, YMIN, YMAX, XAXIS, UNDEF,
     >                   'Revised profile', PLTLINE (2))
         END IF


C  . . . . . . . . . . .  Basic tabulations (*.TAB)  . . . . . . . . . .
C                        and 2nd derivatives (*.YPP)
C                        Curvature is also required.

C  *     Don't try to calculate derivatives for a non-rectified airfoil:

         IF (MODE /= 1) THEN

            IF (.NOT. WRAPCRV) THEN

C  *           Derivatives and curvatures of original profile, using Y vs. X:

               CALL FD12K (NUORIG, XUORIG, YUORIG, YPU, YPPU, YKUORG)
               CALL FD12K (NLORIG, XLORIG, YLORIG, YPL, YPPL, YKLORG)

            ELSE

C  *           Derivatives and curvature via X vs. T, Y vs. T parametric form:

               CALL RVERSE (NLORIG, XLORIG, X)
               CALL RVERSE (NLORIG, YLORIG, Y)
               CALL COPY   (NUORIG, XUORIG, X (NLORIG))
               CALL COPY   (NUORIG, YUORIG, Y (NLORIG))
               NPTS = NLORIG + NUORIG - 1

C              Fit a conventional parametric spline (coefs. stored internally):

               CALL PSFIT (NPTS, X, Y, 'C', .FALSE., IER)
               IF (IER /= 0) GO TO 940

C              -NPTS tells PSTVAL to evaluate at the data points,
C              so there's no need to input the Ts.
C              Reuse X/YWRAP for the (unneeded) X/Y outputs.

               CALL PSTVAL (-NPTS, DUMMY, XWRAP, YWRAP, XP, YP,
     >                      XPP, YPP, NPTS, X, Y)

C              Derive the curvature from the derivatives.

               CALL CURV2D (NPTS, XP, XPP, YP, YPP, CRVORIG)

C  *           Set up the curvature and X derivatives in separate-surface
C              form for tabulation.  First, the lower surface:

               SCALE = XLORIG (NLORIG) - XLORIG (1)  ! For dx/dt ~ 0 test.

               CALL XDERIVS (NLORIG, XP, XPP, YP, YPP, YPL, YPPL, SCALE)
               CALL RVERSE  (NLORIG, YPL, YPL)
               CALL RVERSE  (NLORIG, YPPL, YPPL)
               CALL RVERSE  (NLORIG, CRVORIG, YKLORG)

C              Now the upper surface:

               I = NLORIG
               CALL XDERIVS (NUORIG, XP (I), XPP (I), YP (I), YPP (I),
     >                       YPU, YPPU, SCALE)
               CALL COPY (NUORIG, CRVORIG (I), YKUORG)

            END IF

            IF (TABFIL) THEN

C  *           Tabulate original profile with derivatives and curvatures:

               IF (NOMODS) THEN        ! Suppress 'original'
                  SUBTIT (1) = SURFACE (1)
                  SUBTIT (1) (1:1) = 'U'   ! Not 'u'
                  SUBTIT (2) = SURFACE (2)
                  SUBTIT (2) (1:1) = 'L'
               ELSE
                  SUBTIT (1) = ORIGINAL // SURFACE (1)
                  SUBTIT (2) = ORIGINAL // SURFACE (2)
               END IF

               IF (MODE == 7) THEN
                  IF (CAMBCASE == 2) THEN
                     SUBTIT (1) = 'Mean line'
                     SUBTIT (2) = 'Semi-thickness'
                  END IF
               END IF

               CALL PRTAB (LUNTAB, LEGND (IPR), SUBTIT (1), WRAPCRV,
     >                     NUORIG, XUORIG, YUORIG, YPU, YPPU, YKUORG)
               CALL PRTAB (LUNTAB, LEGND (IPR), SUBTIT (2), WRAPCRV,
     >                     NLORIG, XLORIG, YLORIG, YPL, YPPL, YKLORG)
            END IF

C  *        Calculate geometric properties, unless original profile was
C           camber/thickness:

            IF (.NOT. (MODE == 7 .AND. CAMBCASE == 2)) THEN

C  *           AFGEOM needs upper/lower surface in (*,2) form:

               CALL COPY (NUORIG, XUORIG, X (1))
               CALL COPY (NLORIG, XLORIG, X (1+MAXPTS))
               CALL COPY (NUORIG, YUORIG, Y (1))
               CALL COPY (NLORIG, YLORIG, Y (1+MAXPTS))

C              XEVAL(*) is handy for AFGEOM's YPOWER(*) work-space.

               CALL AFGEOM (NUORIG, NLORIG, MAXPTS, X, Y, AREA, CENTRD,
     >                      MOMENT, THICK, XTH, CAMBER, XCAM, TEANGLE,
     >                      XEVAL, YEVAL, COEFS (1,1), COEFS (1,2),
     >                      COEFS (1,3), IER)

               IF (NOMODS) THEN
                  LEG = BLANK
               ELSE
                  LEG = ORIGINAL
               END IF

               DO J = LUNTAB, LUNCRT, LUNCRT - LUNTAB

C  *              Don't clutter screen with original - just revised, if any:

                  IF ((J == LUNCRT .AND. NOMODS) .OR.
     >                (J == LUNTAB .AND. TABFIL)) THEN

                     IF (IER /= 0) WRITE (J, 1055) ORIGINAL, IER

                     WRITE (J, 1050) LEG, THICK * 100., XTH,
     >                               CAMBER * 100., XCAM, AREA, MOMENT
                  END IF

                  IF (J == LUNTAB .AND. TABFIL) THEN
                     CHORD = MAX (XUORIG (NUORIG), XLORIG (NLORIG)) -
     >                       XUORIG (1)
                     WRITE (J, 1060) XUORIG (1), YUORIG (1), CHORD,
     >                               TEANGLE
                  END IF

               END DO

            END IF

         END IF

         IF (YPPFIL .AND. NOMODS) THEN

C  *        Save 2nd derivatives in standard PROFILE format for possible
C           reuse by REFINE mode.  End values are not needed.

            CALL PRWRIT (LUNYPP, MAXPTS, LEGND (IPR), NUORIG-2,
     >                   NLORIG-2, XUORIG (2), XLORIG (2), YPPU (2),
     >                   YPPL (2), 5, 3)
         END IF

C  *     Repeat for the revised airfoil:
C
         IF (MODE > 0) THEN

            IF (.NOT. WRAPCRV) THEN

C  *           Derivatives and curvatures of revised profile, using Y vs. X:

               CALL FD12K (NU, XU, YU, YPU, YPPU, YKU)
               CALL FD12K (NL, XL, YL, YPL, YPPL, YKL)

            ELSE

C  *           Derivatives and curvature via X vs. T, Y vs. T parametric form:

               CALL RVERSE (NL, XL, X)
               CALL RVERSE (NL, YL, Y)
               CALL COPY   (NU, XU, X (NL))
               CALL COPY   (NU, YU, Y (NL))
               NPTS = NL + NU - 1

               CALL PSFIT (NPTS, X, Y, 'C', .FALSE., IER)
               IF (IER /= 0) GO TO 940

               CALL PSTVAL (-NPTS, DUMMY, XWRAP, YWRAP, XP, YP,
     >                      XPP, YPP, NPTS, X, Y)

               CALL CURV2D (NPTS, XP, XPP, YP, YPP, CRVNEW)

               SCALE = XL (NL) - XL (1)

               CALL XDERIVS (NL, XP, XPP, YP, YPP, YPL, YPPL, SCALE)
               CALL RVERSE  (NL, YPL, YPL)
               CALL RVERSE  (NL, YPPL, YPPL)
               CALL RVERSE  (NL, CRVNEW, YKL)

C              Now the upper surface:

               I = NL
               CALL XDERIVS (NU, XP (I), XPP (I), YP (I), YPP (I),
     >                       YPU, YPPU, SCALE)
               CALL COPY (NU, CRVNEW (I), YKU)

            END IF

            IF (TABFIL) THEN
               SUBTIT (1) = REVISED // SURFACE (1)
               SUBTIT (2) = REVISED // SURFACE (2)

               IF (MODE == 7) THEN
                  IF (CAMBCASE == 1) THEN
                     SUBTIT (1) = 'Mean line'
                     SUBTIT (2) = 'Semi-thickness'
                  END IF
               END IF

               CALL PRTAB (LUNTAB, MNTITL, SUBTIT (1), WRAPCRV,
     >                     NU, XU, YU, YPU, YPPU, YKU)
               CALL PRTAB (LUNTAB, MNTITL, SUBTIT (2), WRAPCRV,
     >                     NL, XL, YL, YPL, YPPL, YKL)
            END IF

C  *        Calculate geometric properties, unless revised profile is
C           camber/thickness:

            IF (.NOT. (MODE == 7 .AND. CAMBCASE == 1)) THEN

               CALL COPY (NU, XU, X (1))
               CALL COPY (NL, XL, X (1+MAXPTS))
               CALL COPY (NU, YU, Y (1))
               CALL COPY (NL, YL, Y (1+MAXPTS))

               CALL AFGEOM (NU, NL, MAXPTS, X, Y, AREA, CENTRD,
     >                      MOMENT, THICK, XTH, CAMBER, XCAM, TEANGLE,
     >                      XEVAL, YEVAL, COEFS (1,1), COEFS (1,2),
     >                      COEFS (1,3), IER)

               DO J = LUNTAB, LUNCRT, LUNCRT - LUNTAB

                  IF (.NOT. (J == LUNTAB .AND. .NOT. TABFIL)) THEN
                     IF (IER /= 0) WRITE (J, 1055) REVISED, IER
                     WRITE (J, 1050) REVISED, THICK * 100., XTH,
     >                               CAMBER * 100., XCAM, AREA, MOMENT
                  END IF

                  IF (J == LUNTAB .AND. TABFIL) THEN
                     CHORD = MAX (XU (NU), XL (NL)) - XU (1)
                     WRITE (J, 1060) XU (1), YU (1), CHORD, TEANGLE
                  END IF

               END DO

            END IF

            IF (YPPFIL) THEN

C  *           Save 2nd derivatives in standard PROFILE format for
C              possible reuse.

               CALL PRWRIT (LUNYPP, MAXPTS, MNTITL, NU-2, NL-2,
     >                      XU (2), XL (2), YPPU (2), YPPL (2), 5, 3)
            END IF
         END IF


C  . . . . . . . . . . . Curvature data file (*.CRV) . . . . . . . . . .

         IF (CRVFIL) THEN

C  *        Save original, revised and target curvature data (if any).
C           If this is the RECTIFY or NORMALIZE option, only revised
C           curvature data are saved.
C
            IF (MODE == 0) THEN
               LEG = BLANK
               ENDSTR = ENDFRM
            ELSE
               LEG = ORIGINAL
               ENDSTR = ENDCRV
            END IF

C  *        Normally, the two surfaces appear on separate frames.

            IF (MODE > 0) THEN   ! Fix the common normalize case
               XLEFT  = XMINSV
               XRIGHT = XMAXSV
               CHORD  = MAX (XU (NU), XL (NL)) - XU (1)

               IF (CHORD == ONE) THEN
                  XLEFT  = ZERO
                  XRIGHT = ONE
               END IF
            END IF

            IF (.NOT. WRAPCRV) THEN

C  *           First, the upper surface:

               CALL QPLFRM (LUNCRV, MNTITL, SBTIT2, XLABEL,
     >                      SURFACE (1) // CURVATURE)

               IF (MODE /= 1 .AND. MODE /= 2)
     >            CALL QPLCRV (LUNCRV, NUORIG-2, XUORIG (2), YKUORG (2),
     >                         ENDSTR, UNDEF, BLANK, XMINSV, XMAXSV,
     >                         CRVMIN, CRVMAX, UNDEF, UNDEF, LEG,
     >                         CRVLINE (1))

               IF (PLTARG .AND. UPPER)
     >            CALL QPLCRV (LUNCRV, NXTARG, TARGX, TARGCRV, ENDCRV,
     >                         UNDEF, BLANK, XMINSV, XMAXSV, CRVMIN,
     >                         CRVMAX, UNDEF, UNDEF, 'Target',
     >                         CRVLINE (3))

               IF (MODE > 0)
     >            CALL QPLCRV (LUNCRV, NU-2, XU (2), YKU (2), ENDFRM,
     >                         UNDEF, BLANK, XLEFT, XRIGHT, CRVMIN,
     >                         CRVMAX, UNDEF, UNDEF, REVISED,
     >                         CRVLINE (2))

C  *           Now for the lower surface:

               CALL QPLFRM (LUNCRV, MNTITL, SBTIT2, XLABEL,
     >                      SURFACE (2) // CURVATURE)

               IF (MODE /= 1 .AND. MODE /= 2)
     >            CALL QPLCRV (LUNCRV, NLORIG-2, XLORIG (2), YKLORG (2),
     >                         ENDSTR, UNDEF, BLANK, XMINSV, XMAXSV,
     >                         CRVMIN, CRVMAX, UNDEF, UNDEF, LEG,
     >                         CRVLINE (1))

               IF (PLTARG .AND. .NOT. UPPER)
     >            CALL QPLCRV (LUNCRV, NXTARG, TARGX, TARGCRV, ENDCRV,
     >                      UNDEF, BLANK, XMINSV, XMAXSV, CRVMIN,
     >                      CRVMAX, UNDEF, UNDEF, 'Target', CRVLINE (3))

               IF (MODE > 0)
     >            CALL QPLCRV (LUNCRV, NL-2, XL (2), YKL (2), ENDFRM,
     >                         UNDEF, BLANK, XLEFT, XRIGHT, CRVMIN,
     >                         CRVMAX, UNDEF, UNDEF, REVISED,
     >                         CRVLINE (2))

            ELSE

C  *           Wrap-around form of curvature requested:

               CALL QPLFRM (LUNCRV, MNTITL, SBTIT2, XLABEL,
     >                      'Curvature')

               IF (MODE /= 1 .AND. MODE /= 2) THEN

                  CALL RVERSE (NLORIG, XLORIG, X)
                  CALL COPY   (NUORIG, XUORIG, X (NLORIG))
                  NPTS = NLORIG + NUORIG - 1

                  WRITE (LUNCRV, 1090)
                  CALL QPLCRV (LUNCRV, NPTS, X, CRVORIG, ENDSTR, UNDEF,
     >                         BLANK, XMINSV, XMAXSV, CRVMIN, CRVMAX,
     >                         UNDEF, UNDEF, LEG, CRVLINE (1))
               END IF

               IF (MODE > 0) THEN     ! As above, for the revised profile.

                  CALL RVERSE (NL, XL, X)
                  CALL COPY   (NU, XU, X (NL))
                  NPTS = NL + NU - 1

                  WRITE (LUNCRV, 1090)
                  CALL QPLCRV (LUNCRV, NPTS, X, CRVNEW, ENDSTR, UNDEF,
     >                         BLANK, XLEFT, XRIGHT, CRVMIN, CRVMAX,
     >                         UNDEF, UNDEF, REVISED, CRVLINE (2))
                END IF
             END IF
          END IF


C  . . . . . . . .  Cp estimates in plottable form (*.CPS) . . . . . . .

         IF (CPSFIL) THEN

C  *        Generate approximate Cp distributions for original and/or
C           revised profile.
C           If this is the RECTIFY or NORMALIZE option, only the revised
C           profile is handled.

C  *        Q:  Why not modularize this to unclutter the main program?
C           A:  Too many upper/lower/original/revised arrays involved.

C  *        First, check that abscissas are same for both surfaces, as
C           required by current GETCPS:

            IF (MODE == 0) THEN

C  *           Set "revised" abscissas so just one form of check will do:

               NU = NUORIG
               NL = NLORIG
               CALL COPY (NL, XLORIG, XL)
               CALL COPY (NU, XUORIG, XU)
            END IF

            DO I = 1, MIN (NU, NL)
               IF (XU (I) /= XL (I)) CPSFIL = .FALSE.
            END DO

C  *        Moment calculations use X/C=0.25, so...

            IF (XL (1) /= ZERO) CPSFIL = .FALSE.
            IF (XL (NL) /= ONE) CPSFIL = .FALSE.

            IF (CPSFIL) THEN

               ALFRAD = ALPHA * ATAN (ONE) / 45.E+0
               IF (FSMACH == ZERO) THEN
                  CPSTAR = -999.E+0
               ELSE
                  GAMMA  = 1.4E+0
                  GOGM1  = GAMMA / (GAMMA - ONE)
                  TOGM1  = 2.E+0 / (GAMMA - ONE)
                  GP1GM1 = (GAMMA + ONE) / (GAMMA - ONE)
                  CPSTAR = (2.E+0 / (GAMMA * FSMACH**2)) *
     >               (((TOGM1 + FSMACH**2) / GP1GM1)**GOGM1 - ONE)
                  CPSTAR = MAX (-999.E+0, CPSTAR)
               END IF

               CALL QPLFRM (LUNCPS, MNTITL, SBTIT2, XLABEL, 'Cp')

               ENDSTR = ENDCRV

               IF (MODE /= 1 .AND. MODE /= 2) THEN

C  *              For the original airfoil...

                  CALL GETCPS (NUORIG, XUORIG, YUORIG, YLORIG,
     >                         ALPHA, FSMACH, XC, CPU, CPL)

C  *              XC(J), CPU(J) here correspond to center of Jth panel.
C                 Integrations require wrap-around form  (lower t.e. to
C                 upper t.e.).  This means inserting points at the true
C                 leading, trailing edges - tedious.   Just average the
C                 upper and lower surface values nearest to the leading
C                 and trailing edges.  (Would extrapolation be better?)

                  XWRAP (1) = ONE
                  YWRAP (1) = YLORIG (NLORIG)
                  PWRAP (1) = (CPL (NLORIG-1) + CPU (NLORIG-1)) * HALF

                  I = NLORIG
                  DO J = 2, NLORIG
                     XWRAP (J) = XC (I-1)
                     YWRAP (J) = (YLORIG (I-1) + YLORIG (I)) * HALF
                     PWRAP (J) = CPL (I-1)
                     I = I - 1
                  END DO

                  J = NLORIG+1
                  XWRAP (J) = ZERO
                  YWRAP (J) = YLORIG (1)
                  PWRAP (J) = (CPL (1) + CPU (1)) * HALF

                  NWRAP = NLORIG + NUORIG + 1
                  I = 2
                  DO J = NLORIG+2, NWRAP-1
                     XWRAP (J) = XC (I-1)
                     YWRAP (J) = (YUORIG (I-1) + YUORIG (I)) * HALF
                     PWRAP (J) = CPU (I-1)
                     I = I + 1
                  END DO

                  XWRAP (NWRAP) = ONE
                  YWRAP (NWRAP) = YUORIG (NUORIG)
                  PWRAP (NWRAP) = PWRAP (1)

C  *              Calculate aerodynamic coefficients.  Ignore CD though.

                  CALL GETCL (NWRAP, XWRAP, YWRAP, PWRAP, ALFRAD,
     >                        CL, CD, CM)

                  IF (NOMODS) THEN
                     LEG = BLANK
                  ELSE
                     LEG = ORIGINAL
                  END IF

                  WRITE (LEGEND, 1120) LEG, CL, CM

C  *              Note handling of Cp data range: most negative Cp
C                 determines "top" of the Cp axis, not "bottom".
C                 Also, rounded values are needed for nice scaling.
C                 Use 0.5 as a sensible increment for most Cp distrbns:

                  CPTOP = PWRAP (1)
                  CPBOT = CPTOP
                  CALL BOUNDS (NWRAP, 1, MXWRAP, PWRAP, CPTOP, CPBOT)

                  NINC  = (-CPTOP + 0.499) / HALF
                  CPTOP = -NINC * HALF
                  NINC  = (CPBOT + 0.499) / HALF
                  CPBOT =  NINC * HALF

                  CALL QPLCRV (LUNCPS, NWRAP, XWRAP, PWRAP,
     >                         ENDSTR, UNDEF, BLANK, UNDEF, UNDEF,
     >                         CPBOT, CPTOP, UNDEF, UNDEF,
     >                         LEGEND, CPSLINE (1))

                  WRITE (LUNCRT, 1150) LEG, CL, CM
                  WRITE (LUNTAB, 1100) LEG, ALPHA, FSMACH, CL, CM
                  WRITE (LUNTAB, 1110) (XC (I), CPU (I), CPL (I),
     >                                  I = 1, NUORIG-1)

               END IF

               IF (MODE > 0) THEN

C  *              For the revised airfoil...

                  CALL GETCPS (NU, XU, YU, YL,
     >                         ALPHA, FSMACH, XC, CPU, CPL)

                  XWRAP (1) = ONE
                  YWRAP (1) = YL (NL)
                  PWRAP (1) = (CPL (NL-1) + CPU (NL-1)) * HALF

                  I = NL
                  DO J = 2, NL
                     XWRAP (J) = XC (I-1)
                     YWRAP (J) = (YL (I-1) + YL (I)) * HALF
                     PWRAP (J) = CPL (I-1)
                     I = I - 1
                  END DO

                  J = NL+1
                  XWRAP (J) = ZERO
                  YWRAP (J) = YL (1)
                  PWRAP (J) = (CPL (1) + CPU (1)) * HALF

                  NWRAP = NL + NU + 1
                  I = 2
                  DO J = NL+2, NWRAP-1
                     XWRAP (J) = XC (I-1)
                     YWRAP (J) = (YU (I-1) + YU (I)) * HALF
                     PWRAP (J) = CPU (I-1)
                     I = I + 1
                  END DO

                  XWRAP (NWRAP) = ONE
                  YWRAP (NWRAP) = YU (NU)
                  PWRAP (NWRAP) = PWRAP (1)

                  CALL GETCL (NWRAP, XWRAP, YWRAP, PWRAP, ALFRAD,
     >                        CL, CD, CM)

                  WRITE (LEGEND, 1120) REVISED, CL, CM

                  CPT = PWRAP (1)
                  CPB = CPT
                  CALL BOUNDS (NWRAP, 1, MXWRAP, PWRAP, CPT, CPB)

                  NINC  = (-CPT + 0.499) / HALF
                  CPT   = MIN (CPTOP, -NINC * HALF)
                  NINC  = (CPB + 0.499) / HALF
                  CPB   = MAX (CPBOT,  NINC * HALF)

                  CALL QPLCRV (LUNCPS, NWRAP, XWRAP, PWRAP,
     >                         ENDSTR, UNDEF, BLANK, UNDEF, UNDEF,
     >                         CPB, CPT, UNDEF, UNDEF, LEGEND,
     >                         CPSLINE (2))

                  LEG = REVISED        ! Because LEG is 1 character shorter
                  WRITE (LUNCRT, 1150) LEG, CL, CM
                  WRITE (LUNTAB, 1100) REVISED, ALPHA, FSMACH, CL, CM
                  WRITE (LUNTAB, 1110) (XC (I), CPU (I), CPL (I),
     >                                  I = 1, NU - 1)
               END IF

C  *           Finally, show Cp* as a curve with a legend.
C              Use legend for Mach and alpha for now - maybe in caption later.

               XC (1)  = -0.1E+0
               XC (2)  =  HALF
               CPU (1) =  CPSTAR
               CPU (2) =  CPSTAR

               IF (CPSTAR /= -999.E+0) THEN
                  WRITE (LEGEND, 1130) CPSTAR, FSMACH, ALPHA
               ELSE
                  WRITE (LEGEND, 1140) FSMACH, ALPHA
               END IF

               CALL QPLCRV (LUNCPS, 2, XC, CPU,
     >                      ENDSTR, UNDEF, BLANK, ZERO, ONE,
     >                      UNDEF, UNDEF, UNDEF, UNDEF,
     >                      LEGEND, 'DASH')

            ELSE  ! An empty Cp file may slip through unless:

               CLOSE (UNIT=LUNCPS, STATUS='DELETE')
               WRITE (LUNCRT, 1000) BLANK,
     >            'Cannot estimate Cps - upper/lower abscissas differ.',
     >            'Use REDISTRIBUTE mode first if you want cheap Cps.'
            END IF

         END IF


C  . . . . . . . Spreadsheet-compatible output (*.SPREAD) . . . . . . .

         IF (SPREADFIL) THEN

C  *        Can't do it if abscissas aren't common to both surfaces:

            DO I = 1, MIN (NU, NL)
               IF (XU (I) /= XL (I)) SPREADFIL = .FALSE.
            END DO

            IF (SPREADFIL) THEN

C  *           We have all outputs except mean-line and thickness.
C              Use X and Y as scratch for XFORM to transform in place:

               CALL COPY (NU, YU, X)
               CALL COPY (NU, YL, Y)
               CALL XFORM (.TRUE., NU, X, Y)

               WRITE (LUNSPR, 1010) MNTITL
               IF (SBTIT2 /= BLANK) WRITE (LUNSPR, 1010) SBTIT2

               WRITE (LUNSPR, 1200) XLABEL,
     >            HTAB // YLABEL // UP, HTAB // YLABEL // LO,
     >            HTAB // 'camber',   HTAB // 'semithickness',
     >            HTAB // 'y''' // UP,  HTAB // 'y''' // LO,
     >            HTAB // 'y"' // UP,   HTAB // 'y"' // LO,
     >            HTAB // 'curvature' // UP, HTAB // 'curvature' // LO

               WRITE (LUNSPR, 1210) (XU (I),
     >            HTAB, YU (I), HTAB, YL (I), HTAB, X (I), HTAB, Y (I),
     >            HTAB, YPU (I), HTAB, YPL (I), HTAB, YPPU (I),
     >            HTAB, YPPL (I), HTAB, YKUORG (I), HTAB, YKLORG (I),
     >            I = 1, NU)

            ELSE

C              An empty spreadsheet-compatible file should be removed:

               CLOSE (UNIT=LUNSPR, STATUS='DELETE')
               WRITE (LUNCRT, 1000) BLANK,
     >           'Cannot give spreadsheet file: upper/lower Xs differ.',
     >           'Use REDISTRIBUTE mode first.'
            END IF
         END IF

      END DO  ! Next profile

C-----------------------------------------------------------------------
C                     End of main loop over all profiles.
C-----------------------------------------------------------------------


C  *  Finally, remind the user of the files generated.

      WRITE (LUNCRT, 1000)
      IF (SAVREV) WRITE (LUNCRT, 1010)
     >   '   Modified airfoil:      ', IDENT (FIRST:LAST), '.dat'
      IF (PLTFIL) WRITE (LUNCRT, 1010)
     >   '   Airfoil plot file:     ', IDENT (FIRST:LAST), '.plt'
      IF (TABFIL) WRITE (LUNCRT, 1010)
     >   '   Tabulated results:     ', IDENT (FIRST:LAST), '.tab'
      IF (CRVFIL) WRITE (LUNCRT, 1010)
     >   '   Curvature QPLOT file:  ', IDENT (FIRST:LAST), '.crv'
      IF (YPPFIL) WRITE (LUNCRT, 1010)
     >   '   2nd derivatives file:  ', IDENT (FIRST:LAST), '.ypp'
      IF (CPSFIL) WRITE (LUNCRT, 1010)
     >   '   Cp distributions file: ', IDENT (FIRST:LAST), '.cps'
      IF (SPREADFIL) WRITE (LUNCRT, 1010)
     >   '   Spreadsheet file:      ', IDENT (FIRST:LAST), '.spread'
      GO TO 999


C-----------------------------------------------------------------------
C                           Error handling.
C-----------------------------------------------------------------------

  810 WRITE (LUNCRT, 1000) 'Stopping as requested.'
      GO TO 999
  850 WRITE (LUNCRT, 1000)
     >   'Error. PROFILE does not handle differing leading edge points.'
      GO TO 999
  860 WRITE (LUNCRT, 1010) ' Error in ', SURFACE (J),
     >   ': abscissas must be monotonic except for RECTIFY mode.'
      GO TO 999
  900 WRITE (LUNCRT, 1045) NPR + 1
      GO TO 999
  910 WRITE (LUNCRT, 1030) NPR
      GO TO 999
  920 WRITE (LUNCRT, 1040) NPR
      GO TO 999
  930 WRITE (LUNCRT, 1000) 'Empty input airfoil data file - abort.'
      GO TO 999
  940 WRITE (LUNCRT, 1020)
     >   ' Bad return from PSFIT during curvature calculation.  IER: ',
     >   IER

  999 WRITE (LUNCRT, 1000)

C *** STOP ' ' ! Avoid system differences


C-----------------------------------------------------------------------
C                              Formats.
C-----------------------------------------------------------------------

 1000 FORMAT (1X, A)
 1010 FORMAT (A, A, A)
 1020 FORMAT (/, A, I2, A3, A)
 1025 FORMAT (A, I2)
 1030 FORMAT (/, ' Too many or too few data pts. found on one surface.',
     >        I3, ' profile(s) processed.')
 1040 FORMAT (/, ' Abnormal EOF or other read error.',
     >        I3, ' profile(s) processed.')
 1045 FORMAT (/, ' Missing coordinate in profile ', I2, '.')
 1050 FORMAT (/, 1X, A,
     >        /'   Maximum thickness/chord (current axes) . .', F9.4,
     >        '% at X/C =', F8.5,
     >        /'   Maximum camber/chord . . . . . . . . . . .', F9.4,
     >        '% at X/C =', F8.5,
     >        /'   Unnormalized area  . . . . . . . . . . . .', E13.6,
C****>        /'   Ordinate of centroid . . . . . . . . . . .', E13.6,
     >        /'   Moment of inertia about X axis . . . . . .', E13.6)
 1055 FORMAT (/, ' *** Geometric properties skipped for ', A,
     >        'airfoil.', /, ' IER from AFGEOM: ', I3)
 1060 FORMAT (/'   L.E.:', 2E14.6, '  Chord:', E14.6,
     >        /'   Mean-line angle at T.E. (deg.) . . . . . .', F13.2)
 1070 FORMAT (/, ' WARNING: Common leading edge point is not foremost.',
     >        /, ' Proceeding, but some results such as derivatives',
     >        ' may be affected.', /)
 1080 FORMAT (/, ' WARNING: The second derivatives file (*.ypp) is',
     >        ' being suppressed because', /,
     >        ' wrap-around (parametric) curvature was specified.', /,
     >        ' REFINE and OPTIMIZE modes require the Y vs. X form of',
     >        ' Y".', /)
 1090 FORMAT ('! The curvature distribution is in wrap-around form, ',
     >        'lower t.e. to upper t.e.', /,
     >        '! Clockwise turns along this arc correspond to ',
     >        'negative curvature.')
 1100 FORMAT (/, '1Cp distribution estimates: ', A9, 'airfoil', //,
     >        ' Alpha:', F6.2, '   Free stream Mach:', F5.2, //,
     >        ' Cl:', F7.4, '   Cm:', F8.4, //,
     >        '      X           Cpu        Cpl')
 1110 FORMAT (1X, F10.5, F11.4, F11.4)
 1120 FORMAT (A, ':  Cl =', F7.4, ',  Cm =', F8.4)
 1130 FORMAT ('Cp* =', F8.2,',  Mach =', F4.2, ',  Alpha =', F5.2)
 1140 FORMAT ('Cp* = Inf., Mach =', F4.2, ',  Alpha =', F5.2)
 1150 FORMAT (1X, A, '  Cl:', F7.4, '   Cm:', F8.4)
 1200 FORMAT (A3, 10A)
 1210 FORMAT (1P, (E13.6, 10(A1, E13.6)))

      END PROGRAM PROFILE
C+----------------------------------------------------------------------
C
      SUBROUTINE ACTIVATE (GATHER, NVARS, ACTIVE, ALLVARS,
     >                     ACTVARS, VSCALES )
C
C PURPOSE:  ACTIVATE either gathers an active set of variables  (as a
C           contiguous subset) from a complete set,  or it "scatters"
C           a given active set back into the complete set.
C
C METHOD:   Both operations are combined here at the expense  of  one
C           more argument as a switch, to keep the module count down.
C           An array of logicals switches each variable on or off; no
C           need is seen to count the active ones here.
C
C ARGUMENTS:
C    ARG       DIM TYPE I/O/S DESCRIPTION
C   GATHER      -    L    I   .TRUE.  means  "gather" or "pack", else
C                             .FALSE. means "scatter" or "unpack".
C   NVARS       -    I    I   No. of variables in the complete set.
C   ACTIVE    NVARS  L    I   ACTIVE(I)=.TRUE. means the Ith variable
C                             is active, else it is inactive.
C   ALLVARS   NVARS  R   I/O  Full set of variables; input if GATHER,
C                             and unchanged,  else updated on return.
C   ACTVARS   NVARS  R   I/O  Active subset, packed; output if GATHER
C                             else input (and unchanged).
C   VSCALES   NVARS  R    I   Multiplicative scale factors applied to
C                             packed output variables if "gather" and
C                             divided out if "scatter".   This  array
C                             should parallel the full set, ALLVARS.
C HISTORY:
C
C   01/17/84   DAS    Initial design and code.
C   01/25/84   DAS    Introduced scaling/unscaling.
C
C AUTHOR: David Saunders, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER   NVARS
      REAL      ACTVARS (NVARS), ALLVARS (NVARS), VSCALES (NVARS)
      LOGICAL   ACTIVE (NVARS), GATHER

C ... Local variables:

      INTEGER   I, J

C ... Execution:

      J = 1

      DO I = 1, NVARS

         IF (ACTIVE (I)) THEN

            IF (GATHER) THEN
               ACTVARS (J) = ALLVARS (I) * VSCALES (I)
            ELSE
               ALLVARS (I) = ACTVARS (J) / VSCALES (I)
            END IF

            J = J + 1
         END IF

      END DO

      END SUBROUTINE ACTIVATE
C+----------------------------------------------------------------------
C
      SUBROUTINE ADDBUMPS (NPTS, X, Y, MXPARM, NBUMPS, BNAMES, PARAMS,
     >                     YNEW )
C
C PURPOSE:  ADDBUMPS adds the effect of a list of "bump" functions to
C           each of the given ordinates.  It may update the ordinates
C           in-place if desired.   It is intended for perturbing air-
C           foils, one surface at a time.
C
C METHOD:   Abscissas are assumed to be in the range [0,1], for effi-
C           ciency reasons.  The original ordinates are copied to the
C           new-ordinate array, then the effect of each bump is added
C           to this output array.  Even simple scaling of Y is imple-
C           mented this way, so it is not essential for a SCALE to be
C           the FIRST bump in the list unless the caller is hoping to
C           pass the same array for YNEW as for Y.
C
C ARGUMENTS:
C    ARG    DIM    TYPE I/O/S DESCRIPTION
C   NPTS     -       I    I   Number of ordinates to be perturbed.
C   X       NPTS     R    I   Abscissas (assumed in [0,1]) and the
C   Y       NPTS     R    I   original ordinates for given surface
C   MXPARM   -       I    I   Max. # parameters  for any one bump.
C   NBUMPS   -       I    I   No. perturbing functions to be used.
C   BNAMES  NBUMPS C*(*)  I   Names of the bumps, as recognized by
C                             subroutine BEVAL (which looks at the
C                             first four characters only).
C   PARAMS  MXPARM,  R    I   PARAMS(*,J)  are parameters defining
C           NBUMPS            the Jth bump.
C   YNEW    NPTS     R    O   The updated ordinates - may  be  the
C                             same array as Y if desired,  but be-
C                             ware if SCALE is among the bump set.
C
C PROCEDURES:
C
C    BEVAL   Evaluates indicated bump function at given normalized
C            abscissas, and adds to current ordinates.
C
C HISTORY:
C
C   01/14/84  DAS  Initial design and code.
C   04/11/84   "   Now uses subroutine BEVAL instead of function BUMP.
C   04/21/97   "   The test for 'SCALE' was failing.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   MXPARM, NBUMPS, NPTS
      REAL
     >   PARAMS (MXPARM, NBUMPS), X (NPTS), Y (NPTS), YNEW (NPTS)
      CHARACTER
     >   BNAMES (NBUMPS) * (*)

C     Local constants:

      LOGICAL, PARAMETER ::
     >   ADD = .TRUE.

      INTEGER
     >   I

C     Procedures:

      EXTERNAL
     >   BEVAL

C     Execution:

      YNEW = Y

      DO I = 1, NBUMPS

         IF (BNAMES (I) (1:4) == 'SCAL') THEN

C           This is a scaling of the ordinates (independent of X(I),
C           but still arranged to be additive, not multiplicative):

            CALL BEVAL (BNAMES (I), MXPARM, PARAMS (1, I), ADD, NPTS, Y,
     >                  YNEW)
         ELSE

            CALL BEVAL (BNAMES (I), MXPARM, PARAMS (1, I), ADD, NPTS, X,
     >                  YNEW)
         END IF

      END DO

      END SUBROUTINE ADDBUMPS
C+----------------------------------------------------------------------
C
      SUBROUTINE AFGEOM (NU, NL, NMAX, X, Y, AREA, CENTRD, MOMENT,
     >                   THICK, XTH, CAMBER, XCAM, TEANGLE, YPOWER,
     >                   YINTRP, B, C, D, IER)
C
C     Description:
C
C           AFGEOM (AirFoil GEOMetry) computes several geometrical quantities
C        based on integrals over an airfoil section, including the moment of
C        inertia about the x-axis, and magnitude and location of maximum
C        thickness and camber.  Each integration is performed by summing the
C        analytically-determined integrals of the cubics defined on each
C        subinterval by the spline representation of the appropriate function.
C
C           This version also computes the angle of the mean line at the
C        trailing edge.
C
C     Arguments:
C
C        Name    Dimension   Type   I/O/S  Description
C        NU                   I     I      Number of upper surface points.
C        NL                   I     I      Number of lower surface points.
C                                          4 <= N <= NMAX for N=NU and NL.
C        NMAX                 I     I      Row dimension of X and Y arrays in
C                                          calling program.
C        X       NMAX,2       R     I      Airfoil abscissae (increasing; not
C                                          necessarily same for both surfaces).
C        Y       NMAX,2       R     I      Airfoil surface ordinates.
C                                          (*,1) = upper surface; (*,2) = lower.
C        AREA                 R       O    Cross-sectional area.
C        CENTRD               R       O    Centroid of y ordinates.
C        MOMENT               R       O    Moment of inertia about x-axis.
C        THICK                R       O    Maximum thickness/chord ratio.
C        XTH                  R       O    X coordinate of maximum thickness.
C        CAMBER               R       O    Maximum camber relative to x-axis,
C                                          as fraction of chord.
C        XCAM                 R       O    X coordinate of maximum camber.
C        TEANGLE              R       O    Trailing edge angle of mean line
C                                          with X-axis (in degrees; negative
C                                          for positive camber).
C        YPOWER  NMAX         R         S  For powers of the surface ordinates.
C        YINTRP  NMAX         R         S  For thickness/camber estimate.
C        B,C,D   NMAX         R         S  For cubic spline coeffs.
C        IER                  I       O    < 0 means NX < 4 and CSFIT will fail;
C                                            1 means non-increasing Xs detected.
C     Notes:
C
C        (1)  There are several two dimensional arrays used.  The second
C             index (*, 1) or (*, 2) means upper or lower surface, resp.
C
C        (2)  Special handling of symmetric airfoils has been omitted for
C             simplicity, while special handling of different sets of
C             abscissae (upper/lower) is included for generality. Parametric
C             splines would provide a better fit at leading edge, but are
C             harder to integrate.
C
C     External references:
C
C        CSEVAL   Needed to evaluate spline for generalized thickness calc.
C        CSFIT    Calculates conventional cubic spline coefficients.
C        CSQUAD   Integrates a conventional cubic spline efficiently.
C        PROTECT  Checks for non-increasing Xs.
C
C
C     Author:  Robert Kennelly, NASA-Ames/Sterling Software
C
C     History:
C
C        15 Dec. 1981   RAK   Initial design and coding.
C        31 Mar. 1982   RAK   Installed in FLO6QNM, with cosmetic changes.
C        18 Apr. 1986   RAK   More cosmetics, including IMPLICIT NONE.
C        22 July 1986   RAK   Only need moment of inertia of area.
C        16 Sep. 1986   DAS   Generalized for different upper/lower Xs;
C                             eliminated symmetric airfoil efficiency;
C                             CSQUAD in place of more obvious quartics. 
C        20 Oct. 1986   RAK   Minor mods. for consistency with FLO6QNM.
C                             Repaired flakey INTERP loop.  Return THICK
C                             and CAMBER as fraction of chord, not %.
C        31 Aug. 1987   DAS   Introduced trailing edge angle calculation
C                             - this is the obvious place to do it for
C                             the PROFILE application of AFGEOM.
C        22 Apr. 1991 DAS/RAK Max. camber and location of it and of max.
C                             thickness were not allowing for nonzero
C                             leading-edge coordinates.
C        21 Oct. 1999   DAS   ROTATE option can cause non-monotonic Xs.
C                             Return IER = +1 in this case.
C
C-----------------------------------------------------------------------


C     Declarations:
C     -------------

      IMPLICIT NONE

C     Constants:

      REAL, PARAMETER ::
     >   ZERO   = 0.0E+0,
     >   ONE    = 1.0E+0,
     >   TWO    = 2.0E+0,
     >   HALF   = 1.0E+0 / TWO,
     >   THIRD  = 1.0E+0 / 3.0E+0,
     >   FOURTH = 1.0E+0 / 4.0E+0,
     >   RADDEG = 57.29578E+0

C     Arguments:

      INTEGER
     >   IER, NL, NMAX, NU
      REAL
     >   AREA, B (NMAX), C (NMAX), CAMBER, CENTRD, D (NMAX), MOMENT,
     >   TEANGLE, THICK, X (NMAX, 2), XCAM, XTH, Y (NMAX, 2),
     >   YPOWER (NMAX), YINTRP (NMAX)

C     Local variables:

      INTEGER
     >   I, J, K, L, N (2), NX
      REAL
     >   ARROW, RCHORD, RESULT (3, 2), YLOWER, YUPPER
      LOGICAL
     >   INTERP, DISTINCT

C     Procedures:

      EXTERNAL
     >   CSEVAL, CSFIT, CSQUAD, PROTECT

C     Execution:
C     ----------

C     Initialize outputs in case of an early return:

      AREA    = ZERO
      CENTRD  = ZERO
      MOMENT  = ZERO
      THICK   = ZERO
      XTH     = ZERO
      CAMBER  = ZERO
      XCAM    = ZERO
      TEANGLE = ZERO

      N (1)   = NU
      N (2)   = NL

      DO J = 1, 2

         CALL PROTECT (N (J), X (1, J), Y (1, J), ARROW, DISTINCT)

         IF (ARROW /= ONE) THEN
            IER = 1
            GO TO 90
         END IF

      END DO

      RCHORD  = ONE / (MAX (X (N (1), 1), X (N (2), 2)) - X (1, 1))

C     For successive powers of the ordinates ...

      DO I = 1, 3

C        ... and for upper and lower surfaces, ...

         DO J = 1, 2

C           ... form integrand, ...

            NX = N (J)
            DO K = 1, NX
               YPOWER (K) = Y (K, J) ** I
            END DO

C           ... fit conventional cubic spline, ...

            CALL CSFIT (NX, X (1, J), YPOWER, 0, ZERO, 0, ZERO,
     >                  B, C, D, IER)
            IF (IER /= 0) GO TO 90

C           ... then integrate and save the results.  CSQUAD to overwrites
C           C for work-space reasons, but this means that thickness and
C           and camber must be calculated first to avoid refitting spline.

            IF (I == 1 .AND. J == 2) THEN

C              Need common sets of abscissae to estimate thickness easily.
C              Note that a more precise value could be computed using the
C              spline for both surfaces and a 1-dim. optimization routine.

               INTERP = (N (1) /= NX)
               IF (.NOT. INTERP) THEN
                  DO K = 1, NX
                     IF (X (K, 1) /= X (K, 2)) THEN

C                       The ordinates don't match - we'll have to interpolate.

                        INTERP = .TRUE.
                        EXIT

                     END IF
                  END DO
               END IF

               IF (INTERP) THEN

C                 Evaluate lower surface spline at upper surface Xs.

                  CALL CSEVAL (NX, X (1, 2), YPOWER, N (1), X (1, 1),
     >                         B, C, D, YINTRP)
               ELSE

C                 Copy lower surface to simplify next steps.

                  DO K = 1, NX
                     YINTRP (K) = YPOWER (K)
                  END DO
               END IF

C              The number of abscissae is now NU = N (1) either way.

               THICK  = ZERO
               CAMBER = ZERO
               XTH    = ZERO
               XCAM   = ZERO

               DO K = 1, N (1)
                  YUPPER = Y (K, 1)
                  YLOWER = YINTRP (K)
                  IF ((YUPPER - YLOWER) > THICK) THEN
                     THICK = YUPPER - YLOWER
                     XTH = X (K, 1)
                  END IF
                  IF (ABS (YUPPER + YLOWER) > ABS (CAMBER)) THEN
                     CAMBER = YUPPER + YLOWER
                     XCAM = X (K, 1)
                  END IF
               END DO

               THICK  = THICK  * RCHORD
               CAMBER = (CAMBER * HALF - Y (1, 1)) * RCHORD
               XTH  = ( XTH - X (1, 1)) * RCHORD
               XCAM = (XCAM - X (1, 1)) * RCHORD

C              This is also a convenient place to compute the angle that
C              the mean-line at the trailing edge makes with the X-axis.

               TEANGLE = ATAN2 (HALF * ((Y (NU, 1) + YINTRP (NU)) -
     >            (Y (NU - 1, 1) + YINTRP (NU - 1))),
     >            (X (NU, 1) - X (NU - 1, 1))) * RADDEG
            END IF

C           Integrate spline for Y ** I for surface J.  CSQUAD returns all
C           the subintegrals but they're not used here.

            CALL CSQUAD (NX, X (1, J), YPOWER, ZERO, C, C)
            RESULT (I, J) = C (NX)
         END DO
      END DO

      AREA   = (RESULT (1, 1) - RESULT (1, 2))
      CENTRD = (RESULT (2, 1) - RESULT (2, 2)) * HALF
      MOMENT = (RESULT (3, 1) - RESULT (3, 2)) * THIRD

C     Termination:
C     ------------

   90 RETURN

      END SUBROUTINE AFGEOM
C+----------------------------------------------------------------------
C
      SUBROUTINE BEVAL (BNAME, NP, P, ADD, NX, X, FX)
C
C PURPOSE:  BEVAL evaluates the indicated "bump" (shape function) at the
C           given abscissas, which are assumed to be normalized.  As the
C           bump names suggest,  these  functions  are commonly used for
C           perturbing airfoils.   They originated as the  "Hicks-Henne"
C           shape functions.  Some of them are also handy for generating
C           distributions of nonuniform weighting factors.
C
C           An option is provided to add the evaluations to the existing
C           values in the output array rather than just return values of
C           the current function.  This can save work-space in a calling
C           routine that is accumulating the effect of several bumps.
C
C METHOD:   The bump is selected by name rather than code number.   This
C           poses a problem of variable-length names,  resolved here  by
C           dealing with only the first four characters  -  an arbitrary
C           choice.   UPPER case is assumed  -  it seemed unnecessary to
C           handle lower case as well. Similarly, there is no attempt to
C           handle  unnormalized abscissas.
C
C           It was considered too inefficient, at this level, to attempt
C           handling of a given bump function's parameters by name. Thus
C           the ordering of the elements of P(*) IS important.
C
C           The option to ADD rather than just evaluate  is  implemented
C           by ALWAYS adding, but zeroing out first if necessary.   This
C           relieves the calling program from doing the zeroing.
C
C ARGUMENTS:
C    ARG    DIM  TYPE  I/O/S DESCRIPTION
C   BNAME    -   C*(*)   I   Name of desired  bump function.   Only
C                            the first 4 characters are  looked at,
C                            and they must be UPPER case.   Look at
C                            the code for valid names.
C    NP      -     I     I   Number of parameters (other than X) in
C                            the selected bump expression.
C    P      NP     R     I   The given values of the parameters, in
C                            a definite order.
C    ADD     -     L     I   ADD = .TRUE. means the evaluations for
C                            each X(I) are added into FX(I);
C                            ADD = .FALSE. means  FX(I)  is  zeroed
C                            out before this addition.
C    NX      -     I     I   The number of abscissas where the sel-
C                            ected function is to be evaluated.
C    X      NX     R     I   Abscissas in the range [0,1]. The case
C                            of simple scaling is an exception: the
C                            inputs and outputs involve ordinates.
C    FX     NX     R    I/O  The desired bump function values (pos-
C                            sibly added to input values; see ADD).
C
C NOTES:
C   *   The ordering of the shape function parameters is such that
C       the LAST one (P (NP)) is usually a multiplicative factor.
C
C HISTORY:
C   09/15/83   DAS   Initial design as FUNCTION BUMP.
C   01/14/84   DAS   Made simple scaling additive rather than multipli-
C                    cative: S*Y = Y + (S-1)*Y (simplifies usage).
C   02/22/84   DAS   Eliminated SQRT and SIN bump (redundant).
C   04/09/84   DAS   Adapted as SUBROUTINE BEVAL to get the loop inside
C                    rather than outside. Switched to selection by name
C                    rather than by type code; provided for  adding  as
C                    well as just evaluating; Wagner functions in-line.
C   07/17/86   DAS   Added simple "RAMP" function option.   Made Wagner
C                    functions the first choice.
C   08/11/86   DAS   Added "FLAP" and "SLAT" options (simple shearing).
C   07/26/89   RAK   ALOG changed to generic LOG.
C   11/06/93   RAK   "DROOP" function needed (1 - X) factor.
C   06/19/96   DAS   "EXP" function now has value 1. at specified X, as
C                    recommended by James Reuther;  added the symmetric
C                    forms of the modified sine function (SIN1 & SIN2),
C                    requiring use of leading 4 characters, not 3.
C   12/18/96   DAS   Added SINF, COSL, COSR, LCOS, & RCOS functions.
C   06/03/97   DAS   RADDEG name was misleading - changed it to DEGRAD.
C   10/20/99    "    Added SIN3 (fore/aft symmetry via SIN1 on [0,1])
C                    and   SIN4 (  "   "   "   "   "   "   " each half).
C
C AUTHORS: Leslie Collins, Robert Kennelly, David Saunders (Sterling);
C          Ray Hicks (NASA Ames Research Center)
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NP, NX
      REAL
     >   FX (NX), P (NP), X (NX)
      LOGICAL
     >   ADD
      CHARACTER
     >   BNAME * (*)

C     Local constants:

      REAL, PARAMETER ::
     >   DEGRAD =  0.017453292519943E+0,
     >   PI     =  3.141592653589793E+0,
     >   PIBY2  =  1.570796326794897E+0,
     >   PT5LOG = -0.693147180559945E+0,
     >   RPI    =  0.318309886183790E+0,
     >   HALF   =  5.E-1,
     >   ONE    =  1.E+0,
     >   TWO    =  2.E+0,
     >   ZERO   =  0.E+0

C     Local variables:

      INTEGER
     >   I
      REAL
     >   AEXP, BEXP, CENTER, CENTER2, N, ONEMC, RN, POWER, POWERL,
     >   POWERR, SINXI, TANGNT, THETA, XI
      CHARACTER
     >   KEY * 4

C     Statement functions:

      REAL
     >   EBUMP, SBUMP, PWR, WIDTH, XNRM

      EBUMP (WIDTH, PWR, XNRM) =
     >   XNRM ** PWR * (ONE - XNRM) * EXP (-WIDTH * XNRM)

      SBUMP (WIDTH, PWR, XNRM) =
     >   (MAX (SIN (PI * XNRM ** PWR), ZERO)) ** WIDTH

C     Execution:

C     Check for just evaluating, not adding:

      IF (.NOT. ADD) THEN
         DO I = 1, NX
            FX (I) = ZERO
         END DO
      END IF

C     Avoid comparison of different-length strings:

      KEY = BNAME (1:4)

      IF (KEY == 'WAGN') THEN   ! Wagner functions.

C        Reference:  Ramamoorthy, P. and Padmavathi, K.  "Airfoil Design
C        by Optimization" in J. Aircraft, Vol. 14 (Feb. 1977), 219-221.

         N = P (1)
         RN = ONE / N

         IF (N == ONE) THEN
            DO I = 1, NX
               THETA = TWO * ASIN (SQRT (X (I)))
               FX (I) = FX (I) + ((THETA + SIN (THETA)) * RPI -
     >            (SIN (HALF * THETA)) ** 2) * P (2)
            END DO
         ELSE
            DO I = 1, NX
               THETA = TWO * ASIN (SQRT (X (I)))
               FX (I) = FX (I) + ((SIN (N * THETA) * RN  +
     >            SIN ((N-ONE) * THETA)) * RPI) * P (2)
            END DO
         END IF

      ELSE IF (KEY == 'SINE') THEN  ! Modified "SINE" (unsymmetric):
                                      ! P1 = "center", P2 = "width"
         POWER = PT5LOG / LOG (P (1))
         DO I = 1, NX
            FX (I) = P (3) * SBUMP (P (2), POWER, X (I)) + FX (I)
         END DO

      ELSE IF (KEY == 'SINF') THEN  ! Flipped "SINE" (unsymmetric):
                                      ! P1 = "center", P2 = "width"
         POWER = PT5LOG / LOG (ONE - P (1))
         DO I = 1, NX
            FX (I) = P (3) * SBUMP (P (2), POWER, ONE - X (I)) + FX (I)
         END DO

      ELSE IF (KEY == 'SIN1') THEN  ! Modified "SINE" (symmetric):
                                      ! P1 = "center", P2 = "width"
         CENTER = P (1)
         IF (CENTER <= HALF) THEN
            POWER = PT5LOG / LOG (CENTER)
            DO I = 1, NX
               FX (I) = P (3) * SBUMP (P (2), POWER, X (I)) + FX (I)
            END DO
         ELSE
            POWER = PT5LOG / LOG (ONE - CENTER)
            DO I = 1, NX
               FX (I) = P (3) * SBUMP (P (2), POWER, ONE - X (I)) +
     >                  FX (I)
            END DO
         END IF

      ELSE IF (KEY == 'SIN2') THEN  ! Modified "SINE" (symmetric):
                                      ! P1 = "center", P2 = "width"
         CENTER = P (1)
         IF (CENTER <= HALF) THEN
            POWER = PT5LOG / LOG (ONE - CENTER)
            DO I = 1, NX
               FX (I) = P (3) * SBUMP (P (2), POWER, ONE - X (I)) +
     >                  FX (I)
            END DO
         ELSE
            POWER = PT5LOG / LOG (CENTER)
            DO I = 1, NX
               FX (I) = P (3) * SBUMP (P (2), POWER, X (I)) + FX (I)
            END DO
         END IF

      ELSE IF (KEY == 'SIN3') THEN ! Fore/aft symmetry via SIN1 on [0, 1]

C        Assume [XA, XB] = [0, 1] and CENTER in (0, 0.5]:

         CENTER = P (1)
         POWER  = PT5LOG / LOG (CENTER)

         DO I = 1, NX
            XI = X (I)
            IF (XI > HALF) XI = ONE - XI
            FX (I) = P (3) * SBUMP (P (2), POWER, XI) + FX (I)
         END DO

      ELSE IF (KEY == 'SIN4') THEN ! Fore/aft symmetry via SIN1 on each half

C        Assume [XA, XB] = [0, 0.5] and CENTER in (0, 1).  Left half:

         CENTER  = P (1)
         POWERL  = PT5LOG / LOG (CENTER)
         CENTER2 = ONE - CENTER
         POWERR  = PT5LOG / LOG (CENTER2)

         IF (CENTER <= HALF) THEN
            DO I = 1, NX
               XI = MAX (ZERO, MIN (X (I) * TWO, ONE))
               FX (I) = P (3) * SBUMP (P (2), POWERL, XI) + FX (I)
            END DO
         ELSE
            DO I = 1, NX
               XI = MAX (ZERO, MIN (X (I) * TWO, ONE))
               FX (I) = P (3) * SBUMP (P (2), POWERR, ONE - XI) + FX (I)
            END DO
         END IF

C        Now the [0.5, 0.1] half with center <-- 1 - center

         IF (CENTER2 <= HALF) THEN
            DO I = 1, NX
               XI = MAX (ZERO, MIN ((X (I) - HALF) * TWO, ONE))
               FX (I) = P (3) * SBUMP (P (2), POWERR, XI) + FX (I)
            END DO
         ELSE
            DO I = 1, NX
               XI = MAX (ZERO, MIN ((X (I) - HALF) * TWO, ONE))
               FX (I) = P (3) * SBUMP (P (2), POWERL, ONE - XI) + FX (I)
            END DO
         END IF

      ELSE IF (KEY == 'COSL') THEN  ! 1/4 cosine (or sine), peak at left

C        SIN is used instead of COS because in the alternative SFEVAL form,
C        where COSL was first installed, PI is adjusted so that SIN (PI) > 0.

         DO I = 1, NX
            SINXI = MAX (SIN (PIBY2 * (X (I) + ONE)), ZERO)
            FX (I) = P (2) * SINXI ** P (1) + FX (I)
         END DO

      ELSE IF (KEY == 'COSR') THEN  ! 1/4 cosine (or sine), peak at right

C        SIN is used instead of COS for consistency with COSL.

         DO I = 1, NX
            SINXI = MAX (SIN (PIBY2 * X (I)), ZERO)
            FX (I) = P (2) * SINXI ** P (1) + FX (I)
         END DO

      ELSE IF (KEY == 'LCOS') THEN  ! Inverted 1/4 (co)sine, peak at left

         DO I = 1, NX
            SINXI = MAX (SIN (PIBY2 * X (I)), ZERO)
            FX (I) = P (2) * (ONE - SINXI ** P (1)) + FX (I)
         END DO

      ELSE IF (KEY == 'RCOS') THEN  ! Inverted 1/4 (co)sine, peak at right

         DO I = 1, NX
            SINXI = MAX (SIN (PIBY2 * (X (I) + ONE)), ZERO)
            FX (I) = P (2) * (ONE - SINXI ** P (1)) + FX (I)
         END DO

      ELSE IF (KEY == 'EXPO') THEN  ! "EXPONENTIAL" (peak height = 1.):
                                      ! P1 = "center", P2 = "width"
         ONEMC = ONE - P (1)  
         AEXP  = P (1) * (ONE + ONEMC * P (2)) / ONEMC
         BEXP  = P (3) / EBUMP (P (2), AEXP, P (1))

         DO I = 1, NX
            FX (I) = BEXP * EBUMP (P (2), AEXP, X (I))  +  FX (I)
         END DO

      ELSE IF (KEY == 'DROO') THEN  ! "DROOP":  P1 = "width"

         DO I = 1, NX
            FX (I) = ((ONE - X (I)) * EXP (-P (1) * X (I))) * P (2)  +
     >               FX (I)
         END DO

      ELSE IF (KEY == 'LEAD') THEN  ! "LEADING"-edge:  P1 = "power"

         DO I = 1, NX
            FX (I) = ((ONE - X (I)) ** P (1)) * P (2)  +  FX (I)
         END DO

      ELSE IF (KEY == 'TRAI') THEN  ! "TRAILING"-edge:  P1 = "power"

         DO I = 1, NX
            FX (I) = (X (I) ** P (1)) * P (2)  +  FX (I)
         END DO

      ELSE IF (KEY == 'FLAP') THEN

C        "FLAP"-like function (shearing transformation only):
C        P (1) is the hinge point (fraction of chord from leading edge);
C        P (2) is the flap angle in DEGREEs (positive is "down").

         TANGNT = TAN (DEGRAD * P (2))
         DO I = 1, NX
            IF (X (I) > P (1))
     >         FX (I) = FX (I) - (X (I) - P (1)) * TANGNT
         END DO

      ELSE IF (KEY == 'SLAT') THEN

C        "SLAT"-like function (shearing transformation only):
C        P (1) is the hinge point (fraction of chord from leading edge);
C        P (2) is the flap angle in DEGREEs (positive is "down").

         TANGNT = TAN (DEGRAD * P (2))
         DO I = 1, NX
            IF (X (I) < P (1))
     >         FX (I) = FX (I) - (P (1) - X (I)) * TANGNT
         END DO

      ELSE IF (KEY == 'RAMP') THEN  ! "RAMP":  Y = P(1) * X

C        Commonly used in conjunction with Wagner functions.

         DO I = 1, NX
            FX (I) = P (1) * X (I)  +  FX (I)
         END DO

      ELSE IF (KEY == 'SCAL') THEN

C        Simple scaling - X(I) are probably ordinates.  Note that the
C        effect is arranged to be  additive,  not multiplicative,  so
C        that scaling can be treated just as for the  other functions
C        when it is included in a group for perturbing purposes.

         DO I = 1, NX
            FX (I) = X (I) * (P (1) - ONE)  +  FX (I)
         END DO

      ELSE
         WRITE (*, '(/, A)')  ' BEVAL: Illegal bump name.'
         STOP
      END IF

      END SUBROUTINE BEVAL
C+----------------------------------------------------------------------
C
      SUBROUTINE BEVAL2 (NAME, MODE, NP, P, ADD, NX, X, FX, LUNERR, IER)
C
C ONE-LINER:  Shape function utility, second collection
C
C PURPOSE:
C
C        BEVAL2 supplements the earlier BEVAL's collection of "bump"
C     (shape) functions.  The functions here are sufficiently more
C     elaborate than the standard ones to warrant separating them.
C     In fact they may require some initial calls to establish some
C     of their parameters - hence the extra MODE argument.  Otherwise
C     the usage is as for BEVAL, q.v.
C
C ARGUMENTS:
C     ARG    DIM  TYPE  I/O/S DESCRIPTION
C     NAME    -   C*(*)   I   Name of desired shape function.   Only
C                             the first 6 characters are looked at,
C                             and they must be UPPER case.  Look at
C                             the code for valid names.
C     MODE    -     I     I   MODE = 0 means initialize the function:
C                                    some element(s) of P are returned;
C                             MODE = 1 means evaluate the function at
C                                    the given abscissa(s).
C     NP      -     I     I   Number of parameters (other than X) in
C                             the shape function.
C     P      NP     R    I/O  The parameters, in a definite order - see
C                             MODE, and the NOTES below.
C     ADD     -     L     I   ADD = .TRUE. means the evaluations for
C                                   each X (I) are added into FX (I);
C                             ADD = .FALSE. means  FX (I) is zeroed
C                                   before this addition.
C     NX      -     I     I   The number of target abscissas.  NX >= 1.
C     X      NX     R     I   Abscissas in the range [0, 1].
C     FX     NX     R    I/O  The desired bump function values (possibly
C                             added to input values - see ADD).
C     LUNERR  -     I     I   LUNERR < 0 suppresses iteration printout (MODE=0);
C                             |LUNERR| is used for error messages.
C     IER     -     I     O   IER = 0 means successful initialization; ! RLC 20211025
C                             IER < 0 means a failure in the zero-finding
C                                   or in the quadrature - probably the
C                                   no. of fn. evals. limit was exceeded.
C NOTES:
C     By convention, the ordering of the shape function parameters is
C     such that the LAST one, P (NP), is a multiplicative factor.
C     See the code for the meaning of the others.
C
C PROCEDURES:
C     QUANC8RC  Reverse-communication adaptive quadrature
C     ZERORC       "      "      "    1-D zero-finder
C    
C HISTORY:
C   02/04/94  D.A.Saunders  Initial adaptation of BEVAL for some new
C                           functions proposed by James.
C
C AUTHORS: James Reuther, David Saunders, NASA Ames, Mt. View, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   MODE, NP, NX, LUNERR, IER
      REAL
     >   FX (NX), P (NP), X (NX)
      LOGICAL
     >   ADD
      CHARACTER
     >   NAME * (*)

C     Local constants:

      INTEGER, PARAMETER ::
     >   MAXFUNZ = 40

      REAL, PARAMETER ::
     >   HALF    =  5.E-1,
     >   ONE     =  1.E+0,
     >   ZERO    =  0.E+0,
     >   LOGPT5  = -0.693147180559945E+0,
     >   TWOPI   =  6.283185307179586E+0,
     >   ABSERR  = 1.E-6,   ! Quadrature's absolute error tolerance
     >   RELERR  = 1.E-5,   ! Tolerance relative to the integral size
     >   TOL     = ZERO,    ! Zero-finder's tolerance
     >   P3A     = -100.,   ! Search interval for P (3)
     >   P3B     = +100.

      CHARACTER, PARAMETER ::
     >   SUBNAME * 6 = 'BEVAL2'

C     Local variables:

      INTEGER
     >   I, I2, ISTATQ, ISTATZ, NFUNQ, NFUNZ
      REAL
     >   ERRESTQ, FLAGQ, FUNQ, FUNZ, HOLDZ (13), P1, P3, P4,
     >   POWER, SINE1, XQ
      LOGICAL
     >   TYPE1

C     Procedures:

      EXTERNAL
     >   QUANC8RC, ZERORC

C     Execution:

      IER = 0

C     Check for just evaluating, not adding:

      IF (.NOT. ADD) THEN
         IF (MODE /= 0) THEN
            FX = ZERO
         END IF
      END IF

      IF (NAME (1 : 4) == 'SINE') THEN

C        Two variations of area-conserving modified sine function.
C        Each produces a family which is symmetric about P1 = 0.5.
C
C        P1 is the zero-crossing value of X in [0, 1];
C        P2 is the power of the sine function (must be a whole number > 0);
C        P3 is the slope of the linear weight function giving zero area;
C        P4 is the overall multiplier.
C
C        'SINE1' uses P1 <= 0.5 functions (or reflections if P1 > 0.5);
C        'SINE2' uses P1 >= 0.5 functions (or reflections if P1 < 0.5).

         P1 = P (1)
         I2 = P (2)     ! (-x) ** y is undefined unless y is an integer

         IF (I2 < 1) THEN
            WRITE (ABS (LUNERR), '(/, A)') ' BEVAL2: P2 < 1.'
            GO TO 900
         END IF

         TYPE1 = NAME (5 : 5) == '1'         ! Else it is '2'
         TYPE1 = TYPE1 .AND. P1 <= HALF  .OR.
     >     .NOT. TYPE1 .AND. P1 >  HALF      ! Expression type now

         IF (TYPE1) THEN
            POWER = LOGPT5 / LOG (P1)
         ELSE
            POWER = LOGPT5 / LOG (ONE - P1)
         END IF

         IF (MODE == 0) THEN

C           Determine the area-conserving weighting factor (only).

            IF (P1 == HALF) THEN

C              Special case:  P3 = 0 is the exact solution.

               P (3) = ZERO
               GO TO 999      ! Avoids excessive indenting.

            END IF

C           We have a zero-finding iteration wrapped around an
C           adaptive quadrature iteration.  Use reverse-communication
C           to avoid passing data via common blocks.

            ISTATZ = 2        ! Initializes the zero-finder
            NFUNZ  = MAXFUNZ

   20       CONTINUE

               CALL ZERORC (P3A, P3B, P3, FUNZ, TOL, NFUNZ, SUBNAME,
     >                      LUNERR, HOLDZ, ISTATZ)

               IF (ISTATZ < 0) THEN  ! Probable application error
                  WRITE (ABS (LUNERR), 1000) ISTATZ
                  IER = ISTATZ
                  GO TO 999

               ELSE IF (ISTATZ > 0) THEN

C                 Calculate the integral for this P3 (iteratively):

                  ISTATQ = 0  ! Initialize the quadrature
                  XQ = ZERO   ! Left end of interval

   30             CONTINUE

C                    Evaluate the shape function (one of two types).
C                    Arrange for higher powers to take the sign of sin ** 1.

                     IF (TYPE1) THEN
                        SINE1 = SIN (TWOPI * XQ ** POWER)
                        FUNQ = (P3 * (XQ - P1) + ONE) *
     >                         SIGN (SINE1 ** I2, SINE1)
                     ELSE
                        SINE1 = SIN (TWOPI * (ONE-XQ) ** POWER)
                        FUNQ = -(P3 * (XQ - P1) + ONE) *
     >                         SIGN (SINE1 ** I2, SINE1)
                     END IF

                     CALL QUANC8RC (XQ, FUNQ, ZERO, ONE, ABSERR, RELERR,
     >                              FUNZ, ERRESTQ, NFUNQ, FLAGQ, ISTATQ)

                     IF (ISTATQ > 0)
     >            GO TO 30

                  IF (ISTATQ /= 0) THEN
                     WRITE (ABS (LUNERR), 1001) FLAGQ, NFUNQ, FUNZ
                     IER = ISTATQ
                     GO TO 999
                  END IF

C                 Else we have evaluated the integral successfully.

                  GO TO 20  ! Continue the search for a zero integral

C              ELSE         ! ISTATZ = 0 - a zero integral has been found
               END IF

            P (3) = P3


         ELSE    ! MODE = 1 (evaluate the shape function with its multiplier)

            P3 = P (3)
            P4 = P (4)
            IF (TYPE1) THEN
               DO I = 1, NX
                  SINE1 = SIN (TWOPI * X (I) ** POWER)
                  FX (I) = FX (I) + P4 * (P3 * (X (I) - P1) + ONE) *
     >                     SIGN (SINE1 ** I2, SINE1)
               END DO
            ELSE
               DO I = 1, NX
                  SINE1 = SIN (TWOPI * (ONE - X (I)) ** POWER)
                  FX (I) = FX (I) - P4 * (P3 * (X (I) - P1) + ONE) *
     >                     SIGN (SINE1 ** I2, SINE1)
               END DO
            END IF

         END IF

      ELSE
         WRITE (ABS (LUNERR), '(/, A)') ' BEVAL2: Bad function name.'
         GO TO 900
      END IF

      GO TO 999

  900 STOP

  999 RETURN

C     Formats:

 1000 FORMAT (/, ' *** BEVAL2 trouble (ZERORC):  status code = ', I2)
 1001 FORMAT (/, ' *** BEVAL2 trouble (QUANC8RC):', /, 5X,
     >        'FLAG = ', F7.3, '   NFUN = ', I5,  '   Integral = ', 1P,
     >        E12.3)

      END SUBROUTINE BEVAL2
C+----------------------------------------------------------------------
C
      SUBROUTINE BLENDFX (NF, NX, MAXNX, FX, WEIGHT, FBLEND)
C
C  PURPOSE:  BLENDFX combines NF discrete functions of X, presumed to be
C            defined at the same abscissas, using the linear combination
C            defined by WEIGHT(1:NF).  The original intended application
C            is to forming composite airfoils from basis airfoils.
C
C  METHOD:   Simple-minded - assume the weights are suitably normalized.
C            If these are design variables,  any normalization  is  best
C            done by the application program anyway.    If both surfaces
C            of a group of profiles are being blended, then each profile
C            would need to be in wrap-around form (one call to BLENDFX),
C            or the surfaces would have to be stored separately (meaning
C            two calls to BLENDFX).
C
C            No attempt is made to check for mismatched abscissas - they
C            are not involved here - just the ordinates.
C
C  ARGUMENTS:
C   ARG      DIM    TYPE I/O/S DESCRIPTION
C    NF       -       I    I   Number of discrete functions to blend
C    NX       -       I    I   Number of data points (same for all fns.)
C  MAXNX      -       I    I   Row dimension of FX(*,*) in calling prog.
C    FX    MAXNX,NF   R    I   The different functions stored by columns
C  WEIGHT    NF       R    I   WEIGHT(J) is applied to FX(I,J) for all I
C  FBLEND    NX       R    O   Blended function desired
C
C  DEVELOPMENT HISTORY:
C  DATE    INITIALS    DESCRIPTION
C  10/24/85   DAS      Initial implementation - after a comment by
C                      Ray Hicks that linear combinations of airfoils
C                      have been used effectively in airfoil design by
C                      optimization.
C
C  AUTHORS:  David Saunders, Informatics, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   NF, NX, MAXNX

      REAL
     >   FX (MAXNX, NF), WEIGHT (NF), FBLEND (NX)

C ... Local variables:

      INTEGER
     >   I, J

C ... Execution:

C ... Considerations:

C     (1)  NX is likely to be much greater than NF.
C     (2)  Why bother to zero out the sums?

      DO I = 1, NX
         FBLEND (I) = WEIGHT(1) * FX (I, 1)
      END DO

      DO J = 2, NF
         DO I = 1, NX
            FBLEND (I) = WEIGHT(J) * FX (I, J) + FBLEND (I)
         END DO
      END DO

      END SUBROUTINE BLENDFX
C+----------------------------------------------------------------------
C
      SUBROUTINE CALCTHICK (NL, NU, XL, XU, YL, YU, THICKNESS,
     >                      XATMAX, B, C, D, YEVAL)
C
C  PURPOSE: CALCTHICK  determines the maximum thickness of an airfoil
C           and the corresponding abscissa. The thickness is returned
C           as a percentage of the chord (which may or may not be 1).
C
C  METHOD:  If the  upper and lower surfaces are not  defined at  the
C           same  set  of abscissas,  a spline  is  fit to  the lower
C           surface and evaluated using the upper surface  abscissas.
C           Otherwise the thickness is determined using the  original
C           coordinates.  No attempt is made to estimate the location
C           BETWEEN data points where the thickness appears greatest.
C
C  ARGUMENTS: These are obvious, except for B, C, D, and YEVAL, which
C             are  work-space for spline coefficents for fitting just
C             one surface,  plus  evaluation of the spline at the ab-
C             scissas of the other surface.
C
C  HISTORY:
C     Oct. '83  DAS  Original design and coding
C     01/25/84  LJC  Added spline fitting for unlike upper and lower
C                    surface abscissas
C     03/16/84  DAS  Arranged for it not to matter which surface is
C                    upper and which is lower
C     07/20/84  LJC  Replaced SPLINE and SEVAL with CSFIT and CSEVAL
C     01/29/90  DAS  Removed tabs, underscores, and END DOs
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C     Arguments:

      INTEGER   NL, NU
      REAL      XATMAX, THICKNESS,
     >          B(NL), C(NL), D(NL), YEVAL(NU),
     >          XL(NL), XU(NU), YL(NL), YU(NU)

C     Local variables:

      INTEGER   I, IER
      REAL      DY
      LOGICAL   FITSPLINE

C     Procedures:

      EXTERNAL  CSEVAL, CSFIT

C     Execution:

      THICKNESS = 0.E+0
      FITSPLINE = NU/=NL
    
      IF (.NOT. FITSPLINE) THEN
         DO I = 1, NU
            IF (XU(I) /= XL(I)) FITSPLINE = .TRUE.
         END DO
      END IF
       
      IF (FITSPLINE) THEN

C        Spline lower surface and evaluate at upper surface abscissas:

         CALL CSFIT (NL, XL, YL, 2, 0., 2, 0., B, C, D, IER)
         IF ( IER/=0 ) GO TO 900

         CALL CSEVAL (NL, XL, YL, NU, XU, B, C, D, YEVAL)
      END IF
   
      DO I = 1, NU
         IF (FITSPLINE) THEN
            DY = YU(I) - YEVAL(I)
         ELSE
            DY = YU(I) - YL(I)
         END IF
         IF (THICKNESS < ABS (DY)) THEN
            THICKNESS = ABS (DY)
            XATMAX = XU(I)
         END IF
      END DO

      THICKNESS = THICKNESS * 100.E+0 / (XU(NU) - XU(1))

      RETURN

  900 WRITE (*, '(/, A)') ' CALCTHICK: Error in fitting spline.'

      END SUBROUTINE CALCTHICK
C+----------------------------------------------------------------------
C
      SUBROUTINE CFDISTRIB (NOPTVARS, OPTVARS, SUMSQS)
C
C ACRONYM: Cf. (compare) distributions
C          --            -------
C PURPOSE: CFDISTRIB  computes  the sum-of-squares type of objective
C          function associated with two 1-dimensional distributions.
C          The argument list is that expected by the QNMDIF optimiz-
C          ing algorithm.   The target and current distributions are
C          NOT assumed to be defined at the same abscissas.  Provis-
C          is made for weighting the elements of the sum of squares.
C
C          This routine is application-dependent,  but the structure
C          may be appropriate in other contexts.
C
C          This version of CFDISTRIB  deals with curvature distribu-
C          tions for one airfoil surface. The curvature distribution
C          is a function of "bump" functions applied to the surface.
C          The active variables/parameters for these  bump functions
C          are optimized by QNMDIF to give a perturbed airfoil surf-
C          ace with curvature as close as possible to a  target dis-
C          tribution - an approach to designing or refining  airfoil
C          profiles.
C
C METHOD:  The target distribution is assumed already available, and
C          may cover any part or all of the range of the  calculated
C          distribution, which is generated here with each call. The
C          bulk of the quantities involved are passed through COMMON
C          because of the restricted calling sequence above. An out-
C          line of the steps required for this application  follows.
C
C          *  "Scatter" (unpack) the given active variables.   (This
C             handles cases where certain shape function  parameters
C             are fixed but must be present during evaluation.)
C
C          *  Apply the shape functions to each of the ordinates de-
C             fining the original surface (NOT in-place).
C
C          *  Calculate the resulting curvature distribution.
C
C          *  Evaluate this calculated distribution at each  of  the
C             target abscissas using table-look-up/linear interpola-
C             tion, and accumulate the desired sum of squares.
C
C          *  If necessary, impose constraints, probably in the form
C             of a penalty function constraining the thickness.
C
C ARGUMENTS:
C    ARG       DIM    TYPE I/O/S DESCRIPTION
C  NOPTVARS     -       I    I   Number of variables being optimized
C                                (passed for compatibility with QNM-
C                                DIF, but not actually used directly
C                                here).
C  OPTVARS   NOPTVARS   R    I   Current values of the variables.
C  SUMSQS       -       R    O   Corresponding value of the function
C                                being minimized.
C
C PARAMETER CONSTANTS:
C
C  LUNERR      Logical unit number for error messages.
C  MXBUMPS     Maximum no. of bump functions provided for.
C  MXOPT       Maximum no. of optimization variables allowed.
C  MXPARM      Maximum no. of parameters associated with any of
C              the bump functions (including a multiplier).
C  MXSURF      Maximum no. of data points handled per surface.
C  MXTARG      Maximum no. of target curvature values handled.
C
C COMMON BLOCKS USED:
C
C   /ACTUAL/  (Current values corresponding to input opt. variables)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C   YCALC     MXSURF     R    S   Perturbed airfoil surface
C  Y1CALC     MXSURF     R    S   1st and 2nd derivatives - by-prod-
C  Y2CALC     MXSURF     R    S   ucts of the curvature calculations
C  YKCALC     MXSURF     R    S   Curvature distribution correspond-
C                                 to the current optimization vars.
C
C   /BCHARS/  (Bump names - wish they could be in /BUMPS/)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C  BNAMES    MXBUMPS   C*11   I   List of bump names
C
C   /BUMPS /  (List of bump function parameters, etc.)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C  NBUMPS       -        I    I   Number of bump functions active
C  PARAMS MXPARM,MXBUMPS R  I/O/S Complete set of parameters for the
C                                 active bump functions,   including
C                                 some possible unused ones  present
C                                 because a 2-dim. array approach is
C                                 used for storing them  (as opposed
C                                 to packing the variable numbers of
C                                 parameters associated with differ-
C                                 ent shape functions).  See ACTIVE.
C  ACTIVE MXPARM,MXBUMPS L    I   Marks the bump function parameters
C                                 as active or inactive.  The unused
C                                 ones must be marked inactive along
C                                 with the fixed (but used) ones.
C  VSCALE MXPARM,NBUMPS  R    I   Scale factors needed by optimizing
C                                 algorithm for active variables.
C
C   /ORIGNL/  (Copies of the original airfoil surface data)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C    NXY        -        I    I   Number of points defining  surface
C   XORIG      NXY       R    I   Abscissas for the surface
C   YORIG      NXY       R    I   Ordinates for the unperturbed srf.
C
C   /TARGET/  (Target distribution quantities)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C   NTARG       -        I    I   Number of target data points
C   XTARG     NTARG      R    I   Abscissas for the target data
C  YKTARG     NTARG      R    I   Target curvature values
C  WEIGHT     NTARG      R    I   Weights to be applied (multiplica-
C                                 tively) to each element of the sum
C                                 of squares.  May be all 1s.
C FILES USED:
C    LUN      I/O/S DESCRIPTION
C   LUNERR      O   Error messages
C
C EXTERNAL REFERENCES:
C  ACTIVATE   For inserting active bump variables into complete set
C  ADDBUMPS   For evaluating and applying the bump functions
C  FD12K      For evaluating the perturbed curvature distribution
C  TABLE1     For matching target and current curvature distributions
C
C HISTORY:
C   01/17/84   DAS   Initial design.
C   03/16/84   DAS   Provided for constraining thickness (penalty fun.)
C   02/07/90   DAS   Removed list-directed I/O.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER  NOPTVARS
      REAL     OPTVARS(NOPTVARS), SUMSQS

C ... Parameter constants and
C ... Common blocks:

C ... Quantities passed to objective function routine via COMMON
C     because of fixed calling sequence expected by QNMDIF.  Any
C     changes affecting these COMMONs should also be made in one
C     other routine - OPTIMIZE.
C
      INTEGER, PARAMETER ::
     >   MXBUMPS = 20,  MXPARM = 3, MXOPT = MXPARM*MXBUMPS,
     >   MXSURF  = 180, MXTARG = MXSURF, LUNERR = 6

      REAL            YCALC, Y1CALC, Y2CALC, YKCALC
      COMMON /ACTUAL/ YCALC(MXSURF), Y1CALC(MXSURF), Y2CALC(MXSURF),
     >                YKCALC(MXSURF)

      CHARACTER*11    BNAMES
      COMMON /BCHARS/ BNAMES(MXBUMPS)

      INTEGER         NBUMPS
      REAL            PARAMS, VSCALE
      LOGICAL         ACTIVE
      COMMON /BUMPS / PARAMS(MXPARM,MXBUMPS), VSCALE(MXOPT),
     >                ACTIVE(MXPARM,MXBUMPS), NBUMPS

      INTEGER         NOTHER, NXY
      REAL            XORIG, XOTHER, YORIG, YOTHER
      COMMON /ORIGNL/ XORIG(MXSURF), YORIG(MXSURF),
     >                XOTHER(MXSURF), YOTHER(MXSURF), NOTHER, NXY

      INTEGER         NTARG
      REAL            OPTWRK, PENLTY, TARGTH, XTARG, YKTARG, WEIGHT
      COMMON /TARGET/ XTARG(MXTARG), YKTARG(MXTARG), WEIGHT(MXTARG),
     >                TARGTH, PENLTY, OPTWRK(4*MXSURF), NTARG

C ... Local variables:

      INTEGER  I, IER, INDEX
      REAL     THICK, XATMAX

C ... Procedures:

      EXTERNAL ACTIVATE, ADDBUMPS, CALCTHICK, FD12K, TABLE1
      REAL     TABLE1


C ... Execution:

C ... Unpack the active (optimization) variables so that the bump
C ... functions can be evaluated. (Some parameters may be fixed.)

      CALL ACTIVATE (.FALSE., MXPARM*NBUMPS, ACTIVE, PARAMS,
     >               OPTVARS, VSCALE)

C ... Evaluate and apply the bump functions to the original surface.

      CALL ADDBUMPS (NXY, XORIG, YORIG, MXPARM,
     >               NBUMPS, BNAMES, PARAMS, YCALC)

C ... Generate the corresponding curvature distribution.

      CALL FD12K (NXY, XORIG, YCALC, Y1CALC, Y2CALC, YKCALC)

C ... Check for constraining the thickness using a penalty function:

      IF (PENLTY > 0.E+0) THEN

         CALL CALCTHICK (NXY, NOTHER, XORIG, XOTHER, YCALC, YOTHER,
     >                   THICK, XATMAX, OPTWRK, OPTWRK(MXSURF+1),
     >                   OPTWRK(2*MXSURF+1), OPTWRK(3*MXSURF+1))

C ...    (Note: CALCTHICK doesn't care which surface is which.)

      ELSE
         THICK = TARGTH
      END IF

C ... Compute the objective function (the penalty part being optional):

      SUMSQS = PENLTY * (THICK - TARGTH) ** 2
      INDEX = 1

      DO I = 1, NTARG

C ...    Interpolate the current distribution to each target abscissa.
C        Assume the target abscissas are ordered (to use updated INDEX).

         SUMSQS = (TABLE1 ( NXY, XORIG, YKCALC, INDEX, XTARG(I), IER)
     >            -  YKTARG(I)) ** 2  * WEIGHT(I)  +  SUMSQS

         IF (IER /= 0) GO TO 900

      END DO

      GO TO 999

C ... Error handling:

  900 WRITE (LUNERR, '(/, A, 2I6)')
     >   ' CFDISTRIB: Table look-up error. IER, I: ', IER, I

  999 RETURN

      END SUBROUTINE CFDISTRIB
C+----------------------------------------------------------------------
C
      SUBROUTINE COMBINE (NU, XU, YU, NL, XL, YL,
     >                    MAXPTS, X, Y, XU2ND, YU2ND, XL2ND, YL2ND,
     >                    LUNCRT, LUNKBD, LUN2ND)
C
C  PURPOSE:  COMBINE combines two "airfoils" by addition or subtraction.
C            Adding or removing boundary layer displacement thickness is
C            the original rationale, but there may be other uses.
C
C  METHOD:   COMBINE is intended to be one of PROFILE's high-level modes
C            of operation.   Therefore a "primary" airfoil is assumed to
C            be read by PROFILE in the usual way,  and a secondary "air-
C            foil" (which may actually represent a boundary layer)  will
C            be prompted for here.   The latter should be passed through
C            PROFILE's REDISTRIBUTE mode first if necessary, so that the
C            abscissas match.   No attempt is made to do the redistribu-
C            tion here - one look at REDISTRIB explains why.
C
C  ARGUMENTS:
C   ARG      DIM    TYPE I/O/S DESCRIPTION
C  NU,NL      -       I   I/O  No. of upper/lower surface pts. (primary)
C  XU,XL    NU,NL     R   I/O  Abscissas, upper/lower (primary)
C  YU,YL    NU,NL     R   I/O  Corresponding ordinates, upper/lower.
C                              Surfaces must have common leading edge.
C  MAXPTS     -       I    I   Max. no. pts. provided for on 1 surface.
C  X,Y     2*MAXPTS   R    S   Buffers for reading secondary airfoil.
C  XU2ND,YU2ND MAXPTS R    S   Secondary abscissas found.
C  XL2ND,YL2ND MAXPTS R    S   Secondary ordinates found.
C  LUNCRT     -       I    I   Logical unit for prompts.
C  LUNKBD     -       I    I   Logical unit for responses.
C  LUN2ND     -       I    I   Logical unit for secondary coordinates
C                              in one of the PROFILE formats.
C
C  PROCEDURES:
C    OPENER   File opening utility
C    PRREAD   For reading secondary profile
C    READER   Prompting utility
C
C  HISTORY:
C  09/30/85   DAS   Adapted from REDISTRIB.
C  01/29/90   DAS   Removed the END DOs; installed OPENER.
C
C  AUTHOR:  David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUN2ND, MAXPTS, NL, NU
      REAL
     >   X(MAXPTS*2), XL(NL), XU(NU), XL2ND(MAXPTS), XU2ND(MAXPTS),
     >   Y(MAXPTS*2), YL(NL), YU(NU), YL2ND(MAXPTS), YU2ND(MAXPTS)

C ... Local variables:

      INTEGER
     >   FORMAT, I, IER, NL2ND, NU2ND
      REAL
     >   SIGNL, SIGNU
      LOGICAL
     >   ADDED, BLAYER, DEFAULT, QUIT
      CHARACTER
     >   TITLE*80, XY2ND*48

C ... Procedures:

      EXTERNAL
     >   OPENER, PRREAD, READY

C ... Execution:

C ... Read secondary coordinates (in any of the PROFILE formats):

      CALL OPENER (LUNCRT,
     >             'Enter file name for secondary coordinates: ',
     >             LUNKBD, XY2ND, LUN2ND, 'OLD')

      CALL PRREAD (LUN2ND, TITLE, MAXPTS, NU2ND, NL2ND, X, Y,
     >             XU2ND, XL2ND, YU2ND, YL2ND, FORMAT, IER)
      IF (IER /= 0) GO TO 830

C ... Check for mismatched abscissas:

      IF (NU2ND /= NU) GO TO 840
      IF (NL2ND /= NL) GO TO 840

      DO I = 1, NU2ND
         IF (XU2ND(I) /= XU(I)) GO TO 840
      END DO

      DO I = 1, NL2ND
         IF (XL2ND(I) /= XL(I)) GO TO 840
      END DO

C ... Check for adding or subtracting:
C     (Hard to avoid distinguishing boundary layer operation from
C     two-true-airfoils case here.)

      BLAYER = .TRUE.
      CALL READY (LUNCRT,
     >   'Do secondary coordinates represent displacement thickness? '//
     >   '<CR>=Y: ',
     >   LUNKBD, BLAYER, DEFAULT, QUIT)

      ADDED = .TRUE.
      CALL READY (LUNCRT,
     >   'Are they to be ADDed? <CR>=Y=added; N=subtracted: ',
     >   LUNKBD, ADDED, DEFAULT, QUIT)

      IF (ADDED) THEN
         SIGNU = +1.E+0
      ELSE
         SIGNU = -1.E+0
      END IF

      IF (BLAYER) THEN
         SIGNL = -SIGNU
      ELSE
         SIGNL = +SIGNU
      END IF

      DO I = 1, NU2ND
         YU(I) = YU(I) + SIGNU * YU2ND(I)
      END DO

      DO I = 1, NL2ND
         YL(I) = YL(I) + SIGNL * YL2ND(I)
      END DO

      RETURN

C ... Error handling:

  830 WRITE (LUNCRT, 1001)
     >   ' COMBINE: Error reading secondary coordinates - quitting.'
      GO TO 990

  840 WRITE (LUNCRT, 1001)
     >   ' COMBINE: The two sets of coordinates must have the same Xs',
     >   '          for corresponding surfaces.',
     >   '          Use PROFILE''s REDISTRIBUTE option then try again.'
  990 STOP

C ... Formats:

 1001 FORMAT (/, (A))

      END SUBROUTINE COMBINE
C+----------------------------------------------------------------------
C
      SUBROUTINE DUMMY (IENTRY, N, NLDIM, X, F, G, H, L, D, NFTOTL,
     >                  NITER, NTYPE, CONTIN)
C
C
C     Description and Usage:
C
C           This is a null version of the user-supplied routine that is
C        required by QNMDIF.  A more meaningful routine can be used, if
C        necessary, to monitor the progress of the optimization.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C
C                  (see QNMDIF header...)
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     History:
C
C        22 Mar. 1983    RAK    Initial coding.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

      LOGICAL
     >   CONTIN

      INTEGER
     >   IENTRY, N, NLDIM, NFTOTL, NITER, NTYPE

      REAL
     >   X, F, G, H, L, D

      END SUBROUTINE DUMMY
C+----------------------------------------------------------------------
C
      SUBROUTINE GETBUMPS (LUNCRT, LUNKBD, LUNBMP, MXPARM, MXBUMPS,
     >                     NBUMPS, BNAMES, PARAMS, PNAMES, ACTIVE,
     >                     NACTIVE, SCALES)
C
C PURPOSE: GETBUMPS returns a user-specified "bump" set description  in
C          the form of a list of bump (shape function) names, a corres-
C          ponding list of values for all of the parameters  associated
C          with the set, and an indication of which of these parameters
C          are to be considered variable (active) or fixed  (inactive).
C          This is useful when  numerical  optimization  techniques are
C          employed to perturb airfoil profiles.
C
C METHOD:  EITHER -
C          The first N Wagner functions may be invoked:   N is prompted
C          for; scale factors are assigned empirically, of order 1000.;
C          and the initial multipliers are all set to  zero, meaning no
C          "bumps" on the first function evaluation.
C
C          OR -
C          A previously-prepared file is invoked, for keyword-style in-
C          put of any of the bump functions known to subroutine  BEVAL.
C          The format of this file is something like this:
C
C          BUMP: SINE
C          CENTER: 0.5   STATUS: INACTIVE
C          WIDTH:  3.0   STATUS: INACTIVE
C          MULTIPLIER: 0.  STATUS: ACTIVE  SCALE: 100.
C
C          BUMP = EXP
C          POWER = 15.  STATUS = FIXED
C          WIDTH = 10.  STATUS = FIXED   SCALE = 1.
C          MULT = .001  STATUS = VARIABLE  SCALE: 10.
C
C             .
C             .
C             .
C
C          Points to note:
C
C            *  Free format, with several possible delimiters.
C            *  Blank lines are optional.
C            *  Keywords need to be long enough to be unambiguous.
C            *  Any BUMP keyword must be the first for that bump,  and
C               on its own.
C            *  The ordering of subsequent lines describing  a  bump's
C               variables or parameters is unimportant,  but  if  some
C               are omitted,  they will be detected as undefined,  and
C               execution will halt.  There  is  NO ATTEMPT TO DEFAULT
C               values (yet).
C            *  Variable names for a given bump  (e.g. WIDTH)  must be
C               the FIRST keyword on a line, one per line.
C            *  Either ordering of STATUS  and SCALE  within a line is
C               acceptable.
C            *  SCALE is optional if STATUS=FIXED/INACTIVE/CONSTANT.
C            *  SCALE is defaulted if STATUS  is  ACTIVE/FREE/VARIABLE
C               and no entry is given.  (1.0 is the likely default.)
C            *  EOF is used to signal end of data.   Some special key-
C               word would probably be necessary on a mainframe.
C
C ARGUMENTS:
C    ARG    DIM    TYPE I/O/S DESCRIPTION
C   LUNCRT   -       I    I   Logical unit for screen prompts.
C   LUNKBD   -       I    I   Logical unit for keyboard entries.
C   LUNBMP   -       I    I   Logical unit for reading bump info.
C   MXPARM   -       I    I   Max. # parameters  for any one bump.
C   MXBUMPS  -       I    I   Max. # bumps allowed for.
C   NBUMPS   -       I    O   Actual # bumps returned.
C   BNAMES  NBUMPS  C*(*) O   Bump names found.
C   PARAMS  MXPARM,  R    O   PARAMS(1:?,J) are the parameters defining
C           NBUMPS            the Jth bump.
C   PNAMES  MXPARM, C*(*) O   Names of these parameters.
C           NBUMPS
C   ACTIVE  MXPARM,  L    O   ACTIVE(I,J)=.TRUE. if the Ith parameter of
C           NBUMPS            of the Jth bump is to be treated as active
C                             (variable).  Otherwise, the parameter is
C                             either to remain at the value returned here
C                             (fixed), or bump J has fewer than I parameters.
C   NACTIVE  -       I    O   Number of active parameters found.
C   SCALES  MXPARM,  R    O   Scale factors, to be applied multiplicatively 
C           NBUMPS            to the active parameters returned here (and
C                             divided out prior to any evaluation of the bumps).
C
C FILES USED: See argument list.
C
C SIGNIFICANT LOCAL VARIABLES:
C    VAR     DIM    TYPE  DESCRIPTION
C   ATTRS     3     C*6   Dictionary of names of "attributes" of a variable
C   BUMPS   NNAMES  C*11  Dictionary of bump function names
C   DICT      4     C*10  Variable dictionary of valid variable names
C   FACTORS MXWAGNER R    Empirical factors used to derive scale factors
C                         for Wagner functions N=2,3,4,.. from the value
C                         entered by the user for N=1
C   LINE      -     C*80  Buffer for one record of bump-description file
C   STATI   NSTATI  C*17  Dictionary of "STATUS" values (e.g., 'ACTIVE')
C
C EXTERNAL REFERENCES:
C   GETLINE  Reads one line; handles comments.
C   LOOKUP   Looks up given key in given dictionary - allows partial
C            matches, and checks for ambiguities.
C   OPENER   File opening utility.
C   PAIRS    Parses a string of keywords + values text into pair(s).
C   READER   Prompting utility with multiple entry points.
C
C ERROR HANDLING:
C     Execution halts with a diagnostic if more bumps or more parameters
C     are found than have been provided for by the calling program. Com-
C     prehensive handling of keyword-driven bump descriptions includes
C     checking for invalid keywords, ambiguous partial keywords, and
C     undefined parameters.
C
C HISTORY:
C   01/20/84   DAS   Initial implementation (Wagner functions only).
C   01/25/84   DAS   Introduced scaling of the variables.
C   02/23/84   DAS   Introduced keyword entry of bump description.
C   03/16/84   DAS   Scale factors for Wagner functions are now pseudo-
C                    variable.  (One can always resort to the keyword-
C                    driven scheme if greater flexibility is needed.)
C   09/24/86   DAS   Clarified use of '|' dictionary entry; minor clean-up.
C   10/24/88   DAS   Removed a scale-factor prompt that showed up undesirably
C                    in the BPLOT application (Wagner fns.; no great loss).
C   02/14/90   DAS   Installed OPENER.
C   12/03/91   DAS   Installed GETLINE; updated coding style.
C   06/19/96   DAS   Added SIN1 and SIN2 functions; EXPONENTIAL function
C                    now expects CENTER, not POWER.
C   12/19/96   DAS   Added SINF and four "COS" functions.
C   10/20/99   DAS   Added SIN3 and SIN4 functions.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNBMP, MXBUMPS, MXPARM, NACTIVE, NBUMPS

      REAL
     >   PARAMS (MXPARM, MXBUMPS), SCALES (MXPARM, MXBUMPS)

      LOGICAL
     >   ACTIVE (MXPARM, MXBUMPS)

      CHARACTER
     >   BNAMES (MXBUMPS) * (*), PNAMES (MXPARM, MXBUMPS) * (*)

C     Local constants:

      INTEGER, PARAMETER ::
     >   MXWAGNER = 20

      REAL, PARAMETER ::
     >   UNDEF = 999.E+0

      LOGICAL, PARAMETER ::
     >   ORDERED = .TRUE.

      CHARACTER, PARAMETER ::
     >   COMMENT * 1 = '!'

C     Local variables:

      INTEGER
     >   I, IENTRY, IK, IOS, IP, J, LAST, N, NEEDED, NKEYS, NP

      REAL
     >   FACTORS (MXWAGNER)

      LOGICAL
     >   DEFAULT, QUIT

      CHARACTER
     >   DICT (3) * 10, FILENAME * 48, KEYS (3) * 10, LINE * 80,
     >   MULTIP * 10, VALUES (3) * 20

C     Procedures:

      EXTERNAL
     >   GETLINE, LOOKUP, OPENER, PAIRS, READI

C     Local character arrays used as static dictionaries:

      INTEGER, PARAMETER ::
     >   NNAMES = 16, NSTATI = 6

      CHARACTER
     >   ATTRS (2) * 6, BUMPS (NNAMES) * 11, STATI (NSTATI) * 17

C ... Helpful comments for future users of LOOKUP (this being the first
C     application):
C
C     The dictionaries must be CAPITALIZED, and LEFT-JUSTIFIED, with no
C     entry exceeding the length indicated in  the  declaration.  Where
C     possible, entries should be ALPHABETIZED for efficiency.  The use
C     of dictionary lookups permits, among other things, partial match-
C     ing and trapping of ambiguities and invalid keywords in a modular
C     sort of way.  The application code runs the risk of becoming NON-
C     mnemonic - it tends to deal with integer item numbers in the dic-
C     tionary rather than with quoted strings. However, the application
C     programmer has the option of staying with strings for readability
C     if efficiency is not an issue, by comparing  DICT (ENTRY) against
C     the (full) string of interest rather than comparing ENTRY against
C     the corresponding item number of interest.    For dictionaries of
C     known length (as in the applications here), the terminating entry
C     '|' is unnecessary.  See LOOKUP for use of '|' to end searches.
C
C     Synonyms may be handled in two ways.   The way that leads to just
C     one result from  LOOKUP  (regardless of which synonym the key re-
C     presents) is illustrated by the "STATI" dictionary.  Briefly, one
C     keyword is treated as fundamental;  any synonyms are entered with
C     this fundamental key appended using any of the delimiters  recog-
C     nized by SCANNR.   This can mean the fundamental string shows  up
C     multiple times, but handling the return from LOOKUP is facilitat-
C     ed.  (The alternative involves checking more than one possibility
C     on return from LOOKUP.)

      DATA
     >   ATTRS  /'SCALE', 'STATUS'/,
     >   BUMPS  /'COSL', 'COSR', 'DROOP', 'EXPONENTIAL',
     >           'LCOS', 'LEADING', 'RCOS', 'SCALE',
     >           'SIN1', 'SIN2', 'SIN3', 'SIN4', 'SINE', 'SINF',
     >           'TRAILING', 'WAGNER'/,
     >   MULTIP /'MULTIPLIER'/,
     >   STATI  /'ACTIVE',
     >           'CONSTANT:INACTIVE',
     >           'FIXED   :INACTIVE',
     >           'FREE    :ACTIVE',
     >           'INACTIVE',
     >           'VARIABLE:ACTIVE'/

      DATA
     >   FACTORS/5*1.E+0, 5*3.E+0, 5*1.E+1, 5*3.E+1/

C     These empirical factors for the relative scaling of the Wagner
C     functions N = 1:MXWAGNER should be tuned some day.


C     Execution:

C     Initialize all parameters as inactive - some will never be used:

      NACTIVE = 0

      DO J = 1, MXBUMPS
         DO I = 1, MXPARM
            ACTIVE (I, J) = .FALSE.
            PARAMS (I, J) = UNDEF
            SCALES (I, J) = 1.E+0
         END DO
      END DO

C     Provide for using the first N Wagner functions easily:

  200 CALL READI (LUNCRT, 'Enter N to use Wagner functions 1:N, ' //
     >            'or <CR> to use an input file: ',
     >            LUNKBD, N, DEFAULT, QUIT )

      IF (.NOT. DEFAULT) THEN

         IF (N > MIN (MXBUMPS, MXWAGNER) .OR. N < 1) GO TO 200

            DO J = 1, N
               BNAMES (J)   = 'WAGNER'
               PNAMES (1, J) = 'N'
               PNAMES (2, J) = MULTIP
               PARAMS (1, J) = J
               PARAMS (2, J) = 0.E+0
               ACTIVE (2, J) = .TRUE.
               SCALES (2, J) = FACTORS (J)
            END DO

            NACTIVE = N
            NBUMPS  = N

      ELSE

C        Bump set has been previously prepared:

         FILENAME = 'bplot.inp'
         CALL OPENER (LUNCRT,
     >      'Enter file name for bumps. Default is bplot.inp: ',
     >      LUNKBD, FILENAME, LUNBMP, 'OLD')

         NBUMPS = 0
  400    CONTINUE

C           Read a line, expecting either 'BUMP', EOF, or empty line:
C
            CALL GETLINE (LUNBMP, COMMENT, LINE, LAST, IOS)

            IF (IOS  < 0)  GO TO 999     ! EOF
            IF (IOS  > 0)  GO TO 800     ! System-dependent error code
            IF (LAST == 0) GO TO 400     ! Empty line

            WRITE (LUNCRT, 1002) LINE (1 : LAST)

C           Organize the line as a keyword/value pair:

            NKEYS = 1
            CALL PAIRS (LINE (1 : LAST), NKEYS, KEYS, VALUES)

            IF (KEYS (1) /= 'BUMP') GO TO 801

C           Identify the bump function name:

            CALL LOOKUP (NNAMES, BUMPS, ORDERED, VALUES (1), IENTRY)
            IF (IENTRY < 1) GO TO 803

            NBUMPS = NBUMPS + 1
            BNAMES (NBUMPS) = BUMPS (IENTRY)

C           Set up a (short) dictionary of variable names for the
C           appropriate bump.  In order to handle many bumps the same
C           way, the order of this dictionary MUST be the order expected
C           by subroutine BEVAL.  This means the dictionary cannot also
C           be alphabetized - hence the two ways of searching.  Note that
C           the items in the CASE statement must be in the same order
C           as they are in the BUMPS dictionary (which IS alphabetized).

            GO TO (410, 420, 430, 440, 450, 460, 470, 480, 490, 500,
     >             510, 520, 530, 540, 550, 560) IENTRY

  410          CONTINUE  ! 'COSL' bump:  (1/4 cosine, peak at left) ** POWER

               NEEDED   = 2
               DICT (1) = 'POWER'
               DICT (2) = MULTIP
               GO TO 600

  420          CONTINUE  ! 'COSR' bump:  (1/4 cosine, peak at right) ** POWER

               NEEDED   = 2
               DICT (1) = 'POWER'
               DICT (2) = MULTIP
               GO TO 600

  430          CONTINUE  ! 'DROOP' bump:  (1 - x) * exp (-WIDTH*x) * MULT

               NEEDED   = 2
               DICT (1) = 'WIDTH'
               DICT (2) = MULTIP
               GO TO 600

  440          CONTINUE  
C              'EXPONENTIAL' bump:  MULT * x**P (1-x) exp (-WIDTH*x) /
C                                   CENTER**P (1-CENTER) exp (-WIDTH*CENTER)
               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  450          CONTINUE  ! 'LCOS' bump:  (1/4 cosine, peak at left) ** POWER

               NEEDED   = 2
               DICT (1) = 'POWER'
               DICT (2) = MULTIP
               GO TO 600

  460          CONTINUE  ! 'LEADING' bump:  (1-x)**POWER * MULT

               NEEDED   = 2
               DICT (1) = 'POWER'
               DICT (2) = MULTIP
               GO TO 600

  470          CONTINUE  ! 'RCOS' bump:  (1/4 cosine, peak at right) ** POWER

               NEEDED   = 2
               DICT (1) = 'POWER'
               DICT (2) = MULTIP
               GO TO 600

  480          CONTINUE  ! 'SCALE' bump:  y = FACTOR * y

               NEEDED   = 1
               DICT (1) = 'FACTOR'
               GO TO 600

  490          CONTINUE  ! 'SIN1' bump (L/L form of SINE)

               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  500          CONTINUE  ! 'SIN2' bump (R/R form of SINE)

               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  510          CONTINUE  ! 'SIN3' bump (fore/aft symmetry)

               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  520          CONTINUE  ! 'SIN4' bump (fore/aft symmetry)

               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  530          CONTINUE  ! 'SINE' bump (unsymmetric, L/R):
                         !  SIN (PI*x** (LOG (.5)/LOG (CENTER))**WIDTH * MULT

               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  540          CONTINUE  ! 'SINF' bump (flipped SINE, unsymmetric, R/L)

               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  550          CONTINUE  ! 'TRAILING' bump:  x**POWER * MULT

               NEEDED   = 2
               DICT (1) = 'POWER'
               DICT (2) = MULTIP
               GO TO 600

  560          CONTINUE  ! 'WAGNER' function:

               NEEDED   = 2
               DICT (1) = 'N'
               DICT (2) = MULTIP
               GO TO 600

  600       CONTINUE
            NP = 1

  610       CONTINUE

C              Read a line, expecting a variable name as first keyword:

               CALL GETLINE (LUNBMP, COMMENT, LINE, LAST, IOS)

               IF (IOS  < 0)  GO TO 805     ! EOF
               IF (IOS  > 0)  GO TO 800     ! System-dependent error code
               IF (LAST == 0) GO TO 610     ! Empty line

               WRITE (LUNCRT, 1002) LINE (1 : LAST)

               NKEYS = 3
               CALL PAIRS (LINE (1 : LAST), NKEYS, KEYS, VALUES)

               CALL LOOKUP (NEEDED, DICT, .NOT. ORDERED, KEYS (1),
     >                      IENTRY)
               IF (IENTRY < 1) GO TO 807

C              Get the value of the variable into PARAMS (*,*) in the
C              order implied by the dynamic dictionary:

               READ (VALUES (1), 1020, ERR=809) PARAMS (IENTRY, NBUMPS)
               PNAMES (IENTRY, NBUMPS) = KEYS (1)

C              Now process any remaining keywords from this line -
C              normally STATUS and SCALE (either order). Special
C              cases: if STATUS is INACTIVE, SCALE is OPTIONAL, and
C              if STATUS is omitted, it defaults to INACTIVE.

               IP = IENTRY
               IK = 2
  700          IF (IK <= NKEYS) THEN

C                 What "attribute" does the next keyword represent?

                  CALL LOOKUP (2, ATTRS, ORDERED, KEYS (IK), IENTRY)
                  IF (IENTRY < 1) GO TO 811

                  IF (ATTRS (IENTRY) == 'STATUS') THEN

C                    Several synonymous values are permitted:

                     CALL LOOKUP (NSTATI, STATI, ORDERED, VALUES (IK),
     >                            IENTRY)
                     IF (IENTRY < 1) GO TO 813

                     IF (STATI (IENTRY) == 'ACTIVE') THEN
                        ACTIVE (IP, NBUMPS) = .TRUE.
                     ELSE  ! Parameter is fixed - ignore any scale factor:
                        IK = 3
                     END IF

                  ELSE IF (ATTRS (IENTRY) == 'SCALE') THEN

C                    Decode the value of the scale factor:

                     READ (VALUES (IK), 1020, ERR = 815)
     >                  SCALES (IP, NBUMPS)

                  END IF

C                 Check for another keyword on this line (3 at most):

                  IK = IK + 1
                  GO TO 700
               END IF

C              Check for more variables needed:

               NP = NP + 1
               IF (NP <= NEEDED)
     >      GO TO 610

C           Still no guarantee the same variable wasn't entered twice.
C           Check for undefined variables for this bump function, and
C           update the number of active variables identified so far:

            DO NP = 1, NEEDED
               IF (PARAMS (NP, NBUMPS) == UNDEF) GO TO 817
               IF (ACTIVE (NP, NBUMPS)) NACTIVE = NACTIVE + 1
            END DO

C           Look for further bumps:

            IF (NBUMPS < MXBUMPS)
     >   GO TO 400

      END IF

      GO TO 999

C     Error Handling:

  800 WRITE (LUNCRT, 1005) IOS
      GO TO 990
  801 WRITE (LUNCRT, 1003) 'Bad keyword where "BUMP" expected: ',
     >                     KEYS (1)
      GO TO 990
  803 WRITE (LUNCRT, 1003) 'Invalid bump name: ', VALUES (1)
      GO TO 990
  805 WRITE (LUNCRT, 1002) 'Unexpected end of data.'
      GO TO 890
  807 WRITE (LUNCRT, 1003) 'Invalid variable name for this bump: ',
     >                     KEYS (1)
      GO TO 890
  809 WRITE (LUNCRT, 1003) 'Bad value for this variable: ', VALUES (1)
      GO TO 890
  811 WRITE (LUNCRT, 1003) 'Bad keyword: ', KEYS (IK)
      GO TO 890
  813 WRITE (LUNCRT, 1003) 'Bad value for STATUS keyword: ', 
     >                     VALUES (IK)
      GO TO 890
  815 WRITE (LUNCRT, 1003) 'Bad value for SCALE: ', VALUES (IK)
      GO TO 890
  817 WRITE (LUNCRT, 1002) 'Undefined variable detected:'
      GO TO 890

  890 WRITE (LUNCRT, 1004) NP, NBUMPS
      WRITE (LUNCRT, 1002) 'Reqd. variables:', (DICT (I), I = 1, NEEDED)

  990 STOP

  999 RETURN

C     Formats:

 1002 FORMAT (1X, A)
 1003 FORMAT (1X, A, A)
 1004 FORMAT (' Parameter line', I2, ' of bump number', I3, '.')
 1005 FORMAT (//' *** GETBUMPS: Read error.  IOS = ', I4)
 1020 FORMAT (BN, F20.0)

      END SUBROUTINE GETBUMPS
C+----------------------------------------------------------------------
C
      SUBROUTINE GETCL (NX, XOC, YOC, CP, ALPHA, CL, CD, CM)
C
C
C     Description and Usage:
C
C           Calculate aerodynamic force coefficients for an airfoil by
C        integrating pressure coefficients using the trapezoid method.
C        Adapted from SUBROUTINE FORCF of FLO6 by A. Jameson.
C
C           In this version, airfoil coordinates (XOC, YOC) are assumed
C        to be normalized to the chord, angle of attack is in radians,
C        and reference for the moment coefficient is the quarter chord.
C
C
C     Arguments:
C
C        Name    Dim.     Type   I/O/S   Description
C        NX       -         I      I     No. of airfoil coordinates
C        XOC      NX        R      I     Normalized coords. in wraparound
C        YOC      "         "      "     form (lower t.e. to upper t.e.)
C        CP       NX        R      I     Corresponding Cp values
C        ALPHA    -         R      I     Associated angle of attack
C        CL       -         R      O     Lift coefficient
C        CD       -         R      O     Drag coefficient
C        CM       -         R      O     Moment coefficient (1/4-chord)
C
C
C     Author:  Robert Kennelly, Informatics Inc., 12 Feb. 1982
C
C-----------------------------------------------------------------------


C     Arguments:
C     ----------

      INTEGER NX
      REAL    XOC(NX), YOC(NX), CP(NX), ALPHA, CL, CD, CM

C     Execution:
C     ----------

      XREF = .25E+0
      CL = 0.E+0
      CD = 0.E+0
      CM = 0.E+0

      DO I = 1, NX - 1
         DX = (XOC(I+1) - XOC(I))
         DY = (YOC(I+1) - YOC(I))
         XA = .5E+0 * (XOC(I+1) + XOC(I))
         YA = .5E+0 * (YOC(I+1) + YOC(I))
         CPA = .5E+0 * (CP(I+1) + CP(I))
         DCL = -CPA * DX
         DCD = CPA * DY
         DCM = DCD * YA - DCL * (XA - XREF)
         CL = CL + DCL
         CD = CD + DCD
         CM = CM + DCM
      END DO

C     Express forces in "lab frame" - rotate by ALPHA.

      CLTEMP = CL * COS (ALPHA) - CD * SIN (ALPHA)
      CD = CL * SIN (ALPHA) + CD * COS (ALPHA)
      CL = CLTEMP

      END SUBROUTINE GETCL
C+----------------------------------------------------------------------
C
      SUBROUTINE GETCPS (NPTS, X, YU, YL, ALPHA, FSMACH, XC, CPU, CPL)
C
C ACRONYM: GET (approximate) Cps for an airfoil
C
C PURPOSE: GETCPS is the subroutine version of program FOIL, implemented
C          by Ilan Kroo (NASA Ames/Stanford U.) as a cheap way of gener-
C          ating Cp estimates for a thick airfoil at an angle of attack,
C          using a vortex and source method.    This version was adapted
C          for use by program PROFILE, and may find other uses.
C
C METHOD:  The same abscissas are expected for upper and lower surfaces,
C          for simplicity here.   The coordinates are not necessarily in
C          normalized form.
C
C          The airfoil is modeled with discrete sources and vortices  on
C          the X-axis at 1/4-panel locations.   The thickness effect  is
C          computed from linearized theory with Riegel's correction.  An
C          N*N dense system (N=NPTS-1) is solved for vortex strengths at
C          each coordinate, from which Cp estimates are derived.   These
C          are then corrected for the given free stream Mach number.
C
C          The number of scratch arrays is large enough that no  attempt
C          has been made to pass work-space through the  argument  list.
C          This version is largely self-contained as a result.
C
C REFERENCES:
C          *  Ilan Kroo, Advanced Aerodynamic Concepts Branch (RAC),
C             NASA Ames Research Center, Moffett Field, CA 94035:
C             Forthcoming note on the basic algorithm (1984)
C          *  Schlichting, H. and Truckenbrodt E.,  "Aerodynamics of the
C             Airplane", 1979, McGraw-Hill: Thickness effects (following
C             Riegels).
C
C ARGUMENTS:
C    ARG    DIM     TYPE I/O/S DESCRIPTION
C
C   NPTS     -        I    I   Number of coordinates for each surface
C     X     NPTS      R    I   Abscissas common to both surfaces
C   YU,YL   NPTS      R    I   Ordinates for upper and lower surfaces
C   ALPHA    -        R    I   Desired angle of attack (degrees)
C   FSMACH   -        R    I   Desired free stream Mach number
C    XC     NPTS-1    R    O   Mid-points of panels
C  CPU,CPL  NPTS-1    R    O   Cp estimates at mid-point locations
C
C PARAMETER CONSTANTS:
C
C   MAXN      I     Maximum size of N = NPTS-1 handled by local arrays
C
C EXTERNAL REFERENCES:
C  DECOMP     LU-decomposition of square matrix
C  SOLVE      Solution of factorized linear system
C
C HISTORY:
C   12/03/84  D.A.Saunders   Adapted as subroutine for use in PROFILE.
C   02/21/85    "     "      Switched to version of DECOMP that estimates
C                            the matrix condition number
C   10/23/91    "     "      Switched back to FORSYTHE version of DECOMP
C                            because of a clash with another PROFILE option.
C
C AUTHOR: Ilan Kroo, NASA Ames/Stanford University, Calif.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C ... Local constants:

      INTEGER, PARAMETER ::
     >   MAXN = 201

      REAL, PARAMETER ::
     >   HALF = 0.5, ONE = 1.0, R14 = 0.25, R34 = 0.75, ZERO = 0.

C ... Arguments:

      INTEGER
     >   NPTS
      REAL
     >   X (NPTS), YU (NPTS), YL (NPTS), ALPHA, FSMACH,
     >   XC (NPTS-1), CPU (NPTS-1), CPL (NPTS-1)

C ... Local variables:

      INTEGER
     >   I, J, N
      INTEGER
     >   IP (MAXN)
      REAL
     >   CONST, COSA, DR, PI, SINA, SUM, VTANL, VTANU,
     >   AIC (MAXN, MAXN), DX (MAXN), DYDXM (MAXN),
     >   DYDXT (MAXN), GAMMA (MAXN), RF (MAXN), VTS (MAXN),
     >   XCTL (MAXN), XVORT (MAXN), YM (MAXN+1), YT (MAXN+1)
      REAL cond, workarray(maxn)                          ! RLC 20211024 
      LOGICAL
     >   NONZRO

C ... External references:

      EXTERNAL
     >   DECOMP, SOLVE

C ... Execution:

      N = NPTS - 1
      IF (N > MAXN) THEN
         WRITE (*, '(/, A)')
     >      ' *** GETCPS: Too many airfoil coordinates. ***'
         GO TO 99
      END IF

      PI = 4.E+0 * ATAN (ONE)
      DR = PI / 180.E+0

C ... Generate the vortex locations (1/4-panel) and normal vectors:

      DO I = 1, N
         XVORT (I) = X (I) * R34 + X (I+1) * R14
         XCTL (I)  = X (I) * R14 + X (I+1) * R34
         XC (I)    = X (I) * HALF+ X (I+1) * HALF
         DX (I) =  (X (I+1) - X (I))
         YM (I) =  (YU (I) + YL (I)) * HALF
         YT (I) =  (YU (I) - YL (I)) * HALF
      END DO

      YM (N+1) =  (YU (N+1) + YL (N+1)) * HALF
      YT (N+1) =  (YU (N+1) - YL (N+1)) * HALF

      DO I = 1, N
         DYDXT (I) =  (YT (I+1) - YT (I)) / DX (I)
         DYDXM (I) =  (YM (I+1) - YM (I)) / DX (I)
      END DO

C ... The thickness effect is computed from linearized theory
C     with Riegel's correction:

      DO I = 1, N
         SUM = ZERO
         DO J = 1, N
            IF ((XC (I) - X (J)) * (XC (I) - X (J+1)) <= ZERO)
     >         CYCLE
            SUM = SUM + DYDXT (J) *
     >         LOG ((XC (I) - X (J)) / (XC (I) - X (J+1)))
         END DO
         RF (I)  = ONE / SQRT (DYDXT (I) ** 2 + ONE)
         VTS (I) = (ONE + SUM / PI) * RF (I)
      END DO

C ... Set up the influence coefficients of the vortices.  The
C     flow tangency conditions define the right-hand-side:

      SINA = SIN (ALPHA * DR)
      COSA = COS (ALPHA * DR)
      NONZRO = .FALSE.

      DO I = 1, N
         DO J = 1, N
            AIC (I, J) = HALF / (PI * (XCTL (I) - XVORT (J)))
         END DO
         AIC (I, I) = AIC (I, I) + HALF * DYDXT (I) / DX (I)
         GAMMA (I) = SINA - DYDXM (I) * VTS(I) * COSA
         IF (GAMMA (I) /= ZERO) NONZRO = .TRUE.
      END DO

      IF (NONZRO) THEN

C ...    Solve the system for the vortex strengths (else they are zeros):

         CALL DECOMP (N, MAXN, AIC, cond, IP, workarray)   ! RLC 20211028

         IF (IP (N) == 0) THEN
            WRITE (*, '(/, A)') ' *** GETCPS: Matrix is singular ***'
            GO TO 99
         END IF

         CALL SOLVE (N, MAXN, AIC, GAMMA, IP)

      END IF

C ... Cp is given by  1 - (Vlocal/U)^2.
C     The Karman-Tsien correction for subsonic Mach numbers is applied.

      CONST = SQRT (ONE - FSMACH ** 2)
      DO I = 1, N
         VTANU = VTS (I) * COSA + HALF * RF (I) * GAMMA (I) / DX (I)
         VTANL = VTS (I) * COSA - HALF * RF (I) * GAMMA (I) / DX (I)
         CPU (I) = ONE - VTANU ** 2
         CPL (I) = ONE - VTANL ** 2
         CPU (I) = CPU (I) / (CONST + HALF * (ONE - CONST) * CPU (I))
         CPL (I) = CPL (I) / (CONST + HALF * (ONE - CONST) * CPL (I))
      END DO

   99 RETURN

      END SUBROUTINE GETCPS
C+----------------------------------------------------------------------
C
      SUBROUTINE GETDIS (LUNCRT, LUNKBD, MODE, NP, P, NPTS)
C
C  PURPOSE:  GETDIS does the prompting needed for most uses of DSTRIB,
C            and possibly other 1-D grid utilities.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S   DIM     DESCRIPTION
C    LUNCRT   I     I      -      Logical unit for screen prompts.
C    LUNKBD   I     I      -        "   "   "   "  keyboard responses.
C    MODE     I     O      -      MODE with which DSTRIB is to be used
C                                 (-1, 0, 1, 2, or 3 - see DSTRIB - or
C                                 4 for the Vinokur distribution, or
C                                 5 for the sinusoid + quadratic). If
C                                 MODE = -99, calling program can quit.
C    NP       I    I/O     -      Number of distribution params. needed.
C                                 Input with maximum room provided (at
C                                 least 1);  output with actual length
C                                 of P (*) if MODE > 0, else not set.
C    P        R     O      NP     Parameters prompted for.  See NP.
C    NPTS     I     O      -      No. of points in desired distribution.
C
C  PROCEDURES:
C    READER      Prompting utility
C
C  HISTORY:
C  05/26/85    DAS    Original design and code (DSTRIB options only).
C  06/17/85    DAS    Switched to MODE=-99, not NPTS<0, to mean "quit".
C  05/29/90    DAS    Clarified internal point prompt (X, not I).
C  05/04/92    DAS    Provided for the Vinokur distribution.
C  04/01/95    DAS    Provided for FOILGRID.
C  10/18/96    DAS    Replaced FOILGRID with FOILGRD.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, MODE, NP, NPTS
      REAL
     >   P (NP)

C     Local variables:

      CHARACTER
     >   YESNO * 1
      LOGICAL
     >   DEFAULT, QUIT

C     Procedures:

      EXTERNAL
     >   READI, READR

C     Execution:

      WRITE (LUNCRT, '(A)') ' ',
     >' Options (5 with indicated defaults is suggested for airfoils):',
     >  '   -1:  Xs or Ts are to be read from a file, not generated',
     >  '    0:  Uniform distribution',
     >  '    1:  Sinusoidal bunching towards the lower end',
     >  '    2:  Symmetric sinusoidal bunching towards both ends',
     >  '    3:  Sinusoidal bunching around an internal point',
     >  '    4:  Vinokur distribution (1st & last increment specified)',
     >  '    5:  Linear + Quadratic + Sine + Cosine combination',
     >  ' '

      MODE = 5
      CALL READI (LUNCRT,
     >   'Enter distribution choice. <CR> = 5; ^Z (or ^D) = quit: ',
     >   LUNKBD, MODE, DEFAULT, QUIT)

      IF (QUIT) THEN
         MODE = -99
         GO TO 999
      END IF

      IF (MODE == -1) GO TO 999

      NPTS = 100
      CALL READI (LUNCRT, 'How many points?  <CR> gives 100: ',
     >            LUNKBD, NPTS, DEFAULT, QUIT)

      IF (MODE == 0) THEN

C  *     Uniform distribution - no other parameters needed:

      ELSE IF (MODE >= 1 .AND. MODE <= 3) THEN

C  *     Sinusoidal-type distributions require a "WIDTH"-type exponent:

         WRITE (LUNCRT, '(A)') ' ',
     >      ' Exponents > 1.0 give bunching in further regions;',
     >      ' fractional exponents in the range (0.,1.] do not.',
     >      ' '

         P (1) = 1.E+0
         CALL READR (LUNCRT, 'Enter exponent. (Default is 1.) ',
     >               LUNKBD, P (1), DEFAULT, QUIT)

         IF (MODE == 3) THEN

            NP = 2

C           Default here is misleading - data not necessarily in [0.,1.]

            P (2) = 0.5E+0
            CALL READR (LUNCRT,
     >     'Internal pt. (X or T) about which pts. are to be bunched: ',
     >         LUNKBD, P (2), DEFAULT, QUIT)

         ELSE

            NP = 1

         END IF

      ELSE IF (MODE == 4) THEN   ! Vinokur distribution

         NP = 2
  400    P (1) = -0.2
         WRITE (LUNCRT, '(A)')
     >      ' First increment?  +ve is absolute; -ve is relative;'
         CALL READR (LUNCRT,
     >      '-r means r% of range; default = 0.2% of range: ',
     >      LUNKBD, P (1), DEFAULT, QUIT)
         IF (P (1) == 0.) GO TO 400

         P (2) = P (1) + P (1)
         CALL READR (LUNCRT,
     >      'Last increment?  Default = twice the first: ',
     >      LUNKBD, P (2), DEFAULT, QUIT)
         IF (P (2) == 0.) GO TO 400

      ELSE IF (MODE == 5) THEN   ! Linear + quadratic + sine + cosine

         NP = 4
         P (1) = 0.04
         P (2) = 0.0
         P (3) = 0.3
         P (4) = 0.66

         WRITE (LUNCRT, '(A)')
     >      ' Defaults for L, Q, S, C terms are 0.04, 0.0, 0.3, 0.66.'
         YESNO = 'Y'
         CALL READC (LUNCRT, 'Take the defaults? (Y/N; <CR>=Y): ',
     >      LUNKBD, YESNO, DEFAULT, QUIT)

         IF (YESNO /= 'Y') THEN
            CALL READR (LUNCRT, 'Weight on    LINEAR term?  [0.04]: ',
     >         LUNKBD, P (1), DEFAULT, QUIT)
            CALL READR (LUNCRT, 'Weight on QUADRATIC term?  [0.00]: ',
     >         LUNKBD, P (2), DEFAULT, QUIT)
            CALL READR (LUNCRT, 'Weight on      SINE term?  [0.30]: ',
     >         LUNKBD, P (3), DEFAULT, QUIT)
            CALL READR (LUNCRT, 'Weight on    COSINE term?  [0.66]: ',
     >         LUNKBD, P (4), DEFAULT, QUIT)
         END IF

      END IF

  999 RETURN

      END SUBROUTINE GETDIS
C+----------------------------------------------------------------------
C
      SUBROUTINE GETSHAPE (LUNCRT, LUNKBD, PROMPT, SUPRES, IBUMP, BNAME,
     >                     NPARAM, PARAMS, PNAMES, BMULT)
C
C  PURPOSE:  GETSHAPE interactively defines a Hicks-Henne-type shape
C            function or "bump", and its parameters.  A multiplier is
C            also returned if argument BMULT is input as zero.  (The
C            application may already have an intended multiplier, so
C            the prompt needs to be suppressible.)  A "quit" request
C            is signalled by NPARAM = 0 on output.
C
C  PROCEDURES:
C     READER  Prompting utility
C     SELECT  Menu utility
C
C  HISTORY:
C     12/21/96  DAS  Initial adaptation of PROFILE's MODIFY/MODSRF module,
C                    including introduction of SELECT menu utility, as
C                    needed for the implicit/explicit smoothing option.
C     12/24/96   "   Added PNAMES argument for use by MODIFY/MODSRF.
C     12/26/96   "   MODIFY's choice of "done" as the default forced adding
C                    this option to the menu.
C     12/18/96   "   Added SINF, COSL, COSR, LCOS, & RCOS functions.
C     10/20/99   "   Added SIN3 (fore/aft symmetry via SIN1 on [0,1])
C                    and   SIN4 (  "   "   "   "   "   "   " each half).
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C ------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER   LUNCRT, LUNKBD  ! I   Screen & keyboard
      CHARACTER PROMPT * (*)    ! I   Prompt to appear first on screen;
                                !     see usage described in SELECT
      LOGICAL   SUPRES          ! I   SUPRES = T means don't display the
                                !     menu (unless the user asks for it)
      INTEGER   IBUMP           ! I/O "Bump" number matching BNAME if
                                !     BNAME is non-blank (reqd. by SELECT
                                !     for defaulting); IBUMP = 0 on output
                                !     means "no more shape functions"
      CHARACTER BNAME * (*)     ! I/O "Bump" name as expected by BEVAL;
                                !     blank on input means no default;
                                !     must match IBUMP if non-blank
      INTEGER   NPARAM          ! O   # shape fn. parameters, including
                                !     a multiplier; NPARAM <= 3 in BEVAL;
                                !     NPARAM = 0 means user aborted
      REAL      PARAMS (*)      ! O   Parameter values; PARAMS (NPARAM)
                                !     is (generally) the multiplier,
                                !     copied from BMULT if BMULT |= 0.
      CHARACTER PNAMES (*) * 10 ! O   Parameter names for the chosen bump
                                !     as needed for writing to a log file
      REAL      BMULT           ! I   Multiplier, prompted for if BMULT is
                                !     0., else copied to PARAMS (NPARAM)
C-------------------------------------------------------------------------

C  *  Local constants:

      INTEGER, PARAMETER ::
     >   MXBUMP = 19,    ! Max. # bump functions available
     >   MXPARM = 3      ! Max. # parameters per function incl. multiplier

      CHARACTER, PARAMETER ::
     >   BLANK * 1 = ' '

C  *  Local variables:

      INTEGER
     >   I, IPARM (MXBUMP, 3), IPNAME, NPARMS (MXBUMP), NPROMPT

      CHARACTER
     >   LPROMPT * 43, MENU (0:MXBUMP) * 63, PARMS (6) * 10

      LOGICAL 
     >   CR, QUIT

C  *  Procedures:

      EXTERNAL
     >   READR, SELECT

C  *  Storage:

      DATA
     >   MENU /
     >' 0.  Done.                                                     ',
     >' 1.  SCALE:  Y <- Y * P1                                       ',
     >' 2.  RAMP:   Y <- X * P1                                       ',
     >' 3.  FLAP:   Y <- Y - (X - P1) * TAN (P2 deg.)  for X > P1     ',
     >' 4.  SLAT:   Y <- Y - (P1 - X) * TAN (P2 deg.)  for X < P1     ',
     >' 5.  TRAIL:  X ** P1                                           ',
     >' 6.  DROOP:  (1 - X) * EXP (-P1 * X)                           ',
     >' 7.  LEAD:   (1 - X) ** P1                                     ',
     >' 8.  EXPO:   (X^P (1-X) EXP (-P2 X)) / (P1^P (1-P) EXP (-P2 P))',
     >' 9.  SINE:   SIN ** P2  of  pi * X ** (LOG (0.5) / LOG (P1))   ',
     >'10.  SINF:   Flipped form of 9 [SINE] - left & right swapped   ',
     >'11.  SIN1:   Symmetric form of 9 [SINE] - left half            ',
     >'12.  SIN2:   Symmetric form of 9 [SINE] - right half           ',
     >'13.  SIN3:   Ensures fore/aft symmetry (SIN1 on [0, 1])        ',
     >'14.  SIN4:   Ensures fore/aft symmetry (SIN1 on each half)     ',
     >'15.  COSL:   1/4 cosine, peak at left; power P1                ',
     >'16.  COSR:   1/4 cosine, peak at left; power P1                ',
     >'17.  LCOS:   1 - 1/4 cosine, peak at left; power P1            ',
     >'18.  RCOS:   1 - 1/4 cosine, peak at left; power P1            ',
     >'19.  WAGNER: Wagner function; P1 = order of term N             '/

      DATA
     >   NPARMS
     >      /1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2/

      DATA
     >   IPARM
     >      /6, 6, 1, 1, 2, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 5,
     >       0, 0, 4, 4, 6, 6, 6, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6,
     >       0, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0/

      DATA
     >   PARMS
     >      /'CENTER', 'POWER', 'WIDTH', 'ANGLE', 'ORDER', 'MULTIPLIER'/

C     IPARM (I, J) selects PARMS (J) for shape function I.

      SAVE
     >   IPARM, MENU, PARMS

C  *  Execution:

C     Prompt for a shape function by number or name.

      CALL SELECT (PROMPT, MXBUMP + 1, MENU, SUPRES, LUNCRT, LUNKBD,
     >             IBUMP, BNAME, QUIT)

      IF (QUIT) GO TO 900
      IF (IBUMP == 0) GO TO 950

      NPARAM  = NPARMS (IBUMP)
      NPROMPT = NPARAM

      IF (BMULT /= 0.) THEN
         IF (IPARM (IBUMP, NPARAM) == 6) THEN
            PARAMS (NPARAM) = BMULT
            NPROMPT = NPROMPT - 1
         END IF
      END IF

      DO I = 1, NPROMPT

C  *     Prompt for and read shape function parameters required:

         IPNAME = IPARM (IBUMP, I)
  200    WRITE (LPROMPT, '(5A)')
     >      'Shape function ', MENU (IBUMP) (6 : 13),
     >      '  Enter ', PARMS (IPNAME), ': '

         CALL READR (LUNCRT, LPROMPT, LUNKBD, PARAMS (I), CR, QUIT)

         IF (QUIT) GO TO 900
         IF (CR) GO TO 200

         PNAMES (I) = PARMS (IPNAME)
      END DO

      GO TO 999


  900 NPARAM = 0  ! Quit
      GO TO 999

  950 NPARAM = 1  ! Distinguish "quit (start over)" from "done"

  999 RETURN

      END SUBROUTINE GETSHAPE
C+-----------------------------------------------------------------------
C
      SUBROUTINE IMPSMOOTH (NPTS, X, Y, ISRF, LUNCRT, LUNKBD, LUNOUT)
C
C  PURPOSE:  IMPSMOOTH applies implicit and/or explicit smoothing to one
C            surface of an airfoil, which need not be normalized.
C
C  METHOD:   Basically,
C
C               _
C               y = y + dx**2 eps y"  (explicit; dx**2 for stability)  or
C               _       _
C               y - eps y" = y        (implicit/stable)
C
C            where y" is the second derivative w.r.t. arc length, s,
C            and eps may be a function of s, using a choice of standard
C            shape functions.
C
C            Suggested usage:  implicit with eps ~ 0.01  and/or
C                              explicit with eps ~ 0.1 (1+ iterations)
C
C  ARGUMENTS:
C
C  VAR   DIM   I/O/S   DESCRIPTION
C  NPTS   -    I       Number of points on current surface
C  X     NPTS  I       Surface coordinates
C  Y     NPTS  I/O     (not necessarily normalized)
C  ISRF   -    I       1 means uppper surface;
C                      2 means lower surface
C  LUNCRT -    I       Logical units for screen, keyboard, and
C  LUNKBD -    I       log file
C  LUNOUT -    I
C
C  PROCEDURES:
C
C  BEVAL     Evaluates numerous shape functions on [0., 1]
C  CHORDS2D  Normalized chord lengths
C  READER    Prompting utility
C  TRDIAG    Tridiagonal solver
C
C  HISTORY:
C
C  12/20/96  DAS  Initial form of scheme outlined by James Reuther,
C                 adapted from WAGSMOOTH (formerly SMSURF).
C  12/24/96  DAS  GETSHAPE needed IBUMP and PNAMES arguments.
C  01/31/97  DAS  Steve Edwards pointed to the need for dT ~ dS**2 to
C                 stabilize the explicit method.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C-------------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   ISRF, LUNCRT, LUNKBD, LUNOUT, NPTS

       REAL
     >   X (NPTS), Y (NPTS)

C  *  Local constants:

      INTEGER, PARAMETER ::
     >   MXPTS = 256   ! Max. # points per surface handled

      REAL, PARAMETER ::
     >   ONE   = 1.E+0,
     >   ZERO  = 0.E+0

      CHARACTER, PARAMETER ::
     >   BLANK * 1 = ' '

C  *  Local variables:

      INTEGER
     >   I, IBUMP, ITER, NITER, NPARAM

      REAL
     >   A (MXPTS), B (MXPTS), C (MXPTS), EPS (MXPTS), S (MXPTS),
     >   PARAMS (3), AI, BI, CI, BMULT, EPSEXP, EPSIMP, HL, HR, TERM,
     >   XLE, XTE, YLE, YTE

      CHARACTER
     >   BNAME * 4, PNAMES (3) * 10, SURFCE (2) * 5

      LOGICAL 
     >   CR, EOF, SAME, SUPRES, UNIFORM

C  *  Procedures:

      EXTERNAL
     >   BEVAL, CHORDS2D, READI, READR, READY, TRDIAG

C  *  Storage:

      SAVE
     >   BNAME, EPSEXP, EPSIMP, IBUMP, NITER, NPARAM, PARAMS

      DATA
     >   SURFCE / 'UPPER', 'LOWER' /

C  *  Execution:

      IF (NPTS > MXPTS) THEN
         WRITE (LUNCRT, '(/, (A, I4))') ' IMPSMOOTH: ', NPTS,
     >      ' points exceeds local limit: ', MXPTS
         STOP
      END IF

      IF (ISRF == 1) THEN

         WRITE (LUNCRT, '(/, 2A)')
     >      ' Enter smoothing controls for the UPPER surface.',
     >      ' (0. = leave it alone.)'
         EPSIMP = 0.01
         CALL READR (LUNCRT,
     >      '[Peak] IMPLICIT smoothing coefficient? (>=0.; [0.01]) ',
     >      LUNKBD, EPSIMP, CR, EOF)

         EPSEXP = 0.1
         CALL READR (LUNCRT,
     >   '[Peak] EXPLICIT smoothing coefficient? (>=0.; [0.10*dS**2]) ',
     >      LUNKBD, EPSEXP, CR, EOF)

         NITER = 0
         IF (EPSEXP /= ZERO) THEN
            NITER = 1
            CALL READI (LUNCRT, '# explicit smoothing iterations? [1] ',
     >      LUNKBD, NITER, CR, EOF)
         END IF

         UNIFORM = .FALSE.
         CALL READY (LUNCRT,
     >      'Uniform smoothing along the chord? (Y/N; [N]) ',
     >      LUNKBD, UNIFORM, CR, EOF)

         IF (.NOT. UNIFORM) THEN

            SUPRES = .FALSE.  ! Don't suppress the menu
            BNAME  = BLANK    ! Suppresses defaulting
            IBUMP  = 0        ! Matches BNAME
            BMULT  = ONE      ! Non-zero suppresses prompt for a multiplier

            CALL GETSHAPE (LUNCRT, LUNKBD,
     >         'Choose a weighting function.',
     >         SUPRES, IBUMP, BNAME, NPARAM, PARAMS, PNAMES, BMULT)

         ELSE
            BNAME = 'N/A '
         END IF

      ELSE  ! ISRF = 2

         IF (EPSIMP == ZERO .AND. EPSEXP == ZERO) THEN

            WRITE (LUNCRT, '(/, A)')
     >         ' Enter smoothing controls for the LOWER surface.'
            SAME   = .FALSE.
            SUPRES = .FALSE.
            BNAME  = BLANK
            IBUMP  = 0

         ELSE

            SAME = .TRUE.
            CALL READY (LUNCRT,
     >         'Smooth LOWER surface the same way? (Y/N; [Y]) ',
     >         LUNKBD, SAME, CR, EOF)

            SUPRES = SAME
         END IF

         IF (.NOT. SAME) THEN
            EPSIMP = 0.01
            CALL READR (LUNCRT,
     >         '[Peak] IMPLICIT smoothing coefficient? (>=0.; [0.01]) ',
     >         LUNKBD, EPSIMP, CR, EOF)

            EPSEXP = 0.01
            CALL READR (LUNCRT,
     >         '[Peak] EXPLICIT smoothing coefficient? (>=0.; [0.01]) ',
     >         LUNKBD, EPSEXP, CR, EOF)

            NITER = 0
            IF (EPSEXP /= ZERO) THEN
               NITER = 1
               CALL READI (LUNCRT,
     >            '# explicit smoothing iterations? [1] ',
     >            LUNKBD, NITER, CR, EOF)
            END IF

            UNIFORM = .FALSE.
            CALL READY (LUNCRT,
     >         'Uniform smoothing along the chord? (Y/N; [N]) ',
     >         LUNKBD, UNIFORM, CR, EOF)

            IF (.NOT. UNIFORM) THEN

               BMULT = ONE
               IF (BNAME == 'N/A ') BNAME = BLANK
               IF (BNAME == BLANK)  IBUMP = 0

               CALL GETSHAPE (LUNCRT, LUNKBD,
     >            'Choose a weighting function.',
     >            SUPRES, IBUMP, BNAME, NPARAM, PARAMS, PNAMES, BMULT)

            ELSE
               BNAME = 'N/A '
            END IF

         END IF

      END IF

      IF (EOF .OR. (EPSIMP == ZERO .AND. EPSEXP == ZERO)) GO TO 900


C     Log the details:

      WRITE (LUNOUT, 1100) SURFCE (ISRF), EPSIMP, EPSEXP, NITER

      IF (.NOT. UNIFORM) THEN
         PARAMS (NPARAM) = ONE  ! May be left over as eps from a previous call
         PNAMES (NPARAM) = 'MULTIPLIER'

         WRITE (LUNOUT, 1200)
     >      BNAME, (PNAMES (I), PARAMS (I), I = 1, NPARAM)
      END IF


C     Calculate normalized chord-lengths.  TERM = total arc is not used.

      CALL CHORDS2D (NPTS, X, Y, .TRUE., TERM, S)

      XLE = X (1)
      YLE = Y (1)
      XTE = X (NPTS)
      YTE = Y (NPTS)

C     Implicit smoothing (for the lower frequencies):

      IF (EPSIMP /= ZERO) THEN

         IF (UNIFORM) THEN
            DO I = 1, NPTS
               EPS (I) = EPSIMP
            END DO
         ELSE
            PARAMS (NPARAM) = EPSIMP
            CALL BEVAL (BNAME, NPARAM, PARAMS, .FALSE., NPTS, S, EPS)
         END IF

         A (1) = ZERO
         B (1) = ONE
         C (1) = ZERO

         DO I = 2, NPTS - 1
            TERM  = -2. * EPS (I) / (S (I + 1) - S (I - 1))
            A (I) = TERM / (S (I) - S (I - 1))
            C (I) = TERM / (S (I + 1) - S (I))
            B (I) = ONE - (A (I) + C (I))
         END DO

         A (NPTS) = ZERO
         B (NPTS) = ONE
         C (NPTS) = ZERO

         CALL TRDIAG (A, B, C, Y, Y, NPTS)

      END IF


C     Explicit smoothing (for the higher frequencies):

      IF (EPSEXP /= ZERO) THEN

         IF (UNIFORM) THEN
            DO I = 1, NPTS
               EPS (I) = EPSEXP
            END DO
         ELSE
            PARAMS (NPARAM) = EPSEXP
            CALL BEVAL (BNAME, NPARAM, PARAMS, .FALSE., NPTS, S, EPS)
         END IF

         DO ITER = 1, NITER

            DO I = 1, NPTS
               A (I) = Y (I)
            END DO

            DO I = 2, NPTS - 1
               HL    = S (I) - S (I - 1)
               HR    = S (I + 1) - S (I)
               TERM  = 2. * EPS (I) * ((MIN (HL, HR) ** 2) / (HL + HR))
               AI    = TERM / HL
               CI    = TERM / HR
               BI    = ONE - (AI + CI)
               Y (I) = AI * A (I - 1) + BI * A (I) + CI * A (I + 1)
            END DO

         END DO

      END IF

      X (1)    = XLE
      X (NPTS) = XTE
      Y (1)    = YLE
      Y (NPTS) = YTE

  900 RETURN

C     Formats:

 1100 FORMAT (//, ' Details of smoothing for ', A5, ' surface:',
     >        //, ' [Peak] implicit coefficient:', F8.3,
     >        /,  ' [Peak] explicit coefficient:', F8.3,
     >        /,  ' Number of explicit iterations: ', I3)
 1200 FORMAT (    ' Nonuniform weighting function: ', A,
     >        /,  ' Weighting function parameters:',
     >        /,  (4X, A10, F10.6))

      END SUBROUTINE IMPSMOOTH
C+----------------------------------------------------------------------
C
      SUBROUTINE LOFT (NU, XU, YU, NL, XL, YL,
     >                 MAXPTS, X, Y, XU2ND, YU2ND, XL2ND, YL2ND,
     >                 LUNCRT, LUNKBD, LUN2ND)
C
C  PURPOSE:  LOFT performs the simplest type of lofting between two wing
C            sections by linear interpolation in the third dimension.
C
C  METHOD:   LOFT is intended to be one of PROFILE's high-level modes of
C            operation.   Therefore a "primary" airfoil is assumed to be
C            read by PROFILE in the usual way, and a "secondary" airfoil
C            is prompted for here.   The latter should be passed through
C            PROFILE's REDISTRIBUTE mode first if necessary, so that the
C            numbers of abscissas match the corresponding numbers on the
C            primary section's upper and lower surfaces.   No attempt is
C            made to do the redistribution here - one look at  REDISTRIB
C            explains why.
C
C  ARGUMENTS:
C   ARG      DIM    TYPE I/O/S DESCRIPTION
C  NU,NL      -       I   I/O  No. of upper/lower surface pts. (primary)
C  XU,XL    NU,NL     R   I/O  Abscissas, upper/lower (primary)
C  YU,YL    NU,NL     R   I/O  Corresponding ordinates, upper/lower.
C                              Surfaces must have common leading edge.
C  MAXPTS     -       I    I   Max. no. pts. provided for on 1 surface.
C  X,Y     2*MAXPTS   R    S   Buffers for reading secondary airfoil.
C  XU2ND,YU2ND MAXPTS R    S   Secondary abscissas found.
C  XL2ND,YL2ND MAXPTS R    S   Secondary ordinates found.
C  LUNCRT     -       I    I   Logical unit for prompts.
C  LUNKBD     -       I    I   Logical unit for responses.
C  LUN2ND     -       I    I   Logical unit for secondary coordinates
C                              in one of the PROFILE formats.
C
C  PROCEDURES:
C    LINTRP   Modularizes constant interpolation between multiple pairs
C    NRMLIZ   Normalizes/denormalizes coordinates
C    OPENER   File opening utility
C    PRREAD   For reading secondary profile
C    RDREALS  Prompts for one or more real values
C
C  HISTORY:
C  06/28/90   DAS   Adapted from COMBINE option.
C  10/21/93   DAS   Lofting of twisted airfoils meant to safeguard of
C                   checking for equal relative Xs was not appropriate.
C                   Just require the same numbers of points.
C
C  AUTHOR:  David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUN2ND, MAXPTS, NL, NU
      REAL
     >   X (MAXPTS*2), XL (NL), XU (NU), XL2ND (MAXPTS), XU2ND (MAXPTS),
     >   Y (MAXPTS*2), YL (NL), YU (NU), YL2ND (MAXPTS), YU2ND (MAXPTS)

C ... Local constants:

      REAL, PARAMETER ::
     >   ONE  = 1.E+0,
     >   ZERO = 0.E+0,
     >   TOL  = 1.E-6   ! For comparing normalized Xs

C ... Local variables:

      INTEGER
     >   FORMAT, I, IER, NL2ND, NU2ND, NVALS
      REAL
     >   CHORD, XLE, Z (3)
      CHARACTER
     >   TITLE*80, XY2ND*48

C ... Procedures:

      EXTERNAL
     >   LINTRP, NRMLIZ, OPENER, PRREAD, RDREALS

C ... Execution:

C ... Read secondary coordinates (in any of the PROFILE formats):

      CALL OPENER (LUNCRT,
     >   'Enter file name for secondary coordinates: ',
     >   LUNKBD, XY2ND, LUN2ND, 'OLD')

      CALL PRREAD (LUN2ND, TITLE, MAXPTS, NU2ND, NL2ND, X, Y,
     >   XU2ND, XL2ND, YU2ND, YL2ND, FORMAT, IER)
      IF (IER /= 0) GO TO 830

C ... Check for mismatched numbers of points:

      IF (NU2ND /= NU) GO TO 840
      IF (NL2ND /= NL) GO TO 840

C ... Programmer note:  Use local variables NU2ND, NL2ND everywhere now -
C     it may save a little code.

C     The RELATIVE distributions of abscissas should also match.
C*****NO: not for twisted sections.  Abandon this check.
C     Use X (*) for normalizing the primary Xs; Y (*) for the secondary Xs:

C*****CALL NRMLIZ (NU2ND, XU, X, XU (1), XU (NU2ND) - XU (1))
C*****CALL NRMLIZ (NU2ND, XU2ND, Y, XU2ND (1),
C****>             XU2ND (NU2ND) - XU2ND (1))

C*****DO 410, I = 1, NU2ND
C*****   IF ( ABS (X (I) - Y (I)) > TOL) GO TO 840
C*410 CONTINUE

C*****CALL NRMLIZ (NL2ND, XL, X, XL (1), XL (NL2ND) - XL (1))
C*****CALL NRMLIZ (NL2ND, XL2ND, Y, XL2ND (1),
C****>             XL2ND (NL2ND) - XL2ND (1))

C*****DO 420, I = 1, NL2ND
C*****   IF (ABS (X (I) - Y (I)) > TOL) GO TO 840
C*420 CONTINUE

C ... Allow for the possibility of starting with normalized sections,
C     which may or may not need to be denormalized:

      IF (XU (1) == ZERO  .AND.  YU (1) == ZERO  .AND.
     >    XU (NU2ND) == ONE) THEN

  510    NVALS = 3
         WRITE (LUNCRT, 1001)
     >      ' Enter primary section''s unnormalized (Xle, Yle)' //
     >      ' and chord.'
         CALL RDREALS (LUNCRT, '$(3 values; <CR> = leave normalized): ',
     >      LUNKBD, NVALS, Z)
         IF (NVALS == 1 .OR. NVALS == 2) GO TO 510

         IF (NVALS == 3) THEN
            CALL NRMLIZ (NU2ND, XU, XU, Z (1), -Z (3))
            CALL NRMLIZ (NU2ND, YU, YU, Z (2), -Z (3))
            CALL NRMLIZ (NL2ND, XL, XL, Z (1), -Z (3))
            CALL NRMLIZ (NL2ND, YL, YL, Z (2), -Z (3))
         END IF
      END IF

      IF (XU2ND (1) == ZERO  .AND.  YU2ND (1) == ZERO  .AND.
     >    XU2ND (NU2ND) == ONE) THEN

  520    NVALS = 3
         WRITE (LUNCRT, 1001)
     >      ' Enter secondary section''s unnormalized (Xle, Yle)' //
     >      ' and chord.'
         CALL RDREALS (LUNCRT, '$(3 values; <CR> = leave normalized): ',
     >      LUNKBD, NVALS, Z)
         IF (NVALS == 1 .OR. NVALS == 2) GO TO 520

         IF (NVALS == 3) THEN
            CALL NRMLIZ (NU2ND, XU2ND, XU2ND, Z (1), -Z (3))
            CALL NRMLIZ (NU2ND, YU2ND, YU2ND, Z (2), -Z (3))
            CALL NRMLIZ (NL2ND, XL2ND, XL2ND, Z (1), -Z (3))
            CALL NRMLIZ (NL2ND, YL2ND, YL2ND, Z (2), -Z (3))
         END IF
      END IF

C ... Now for the third dimension...

  600 NVALS = 3
      WRITE (LUNCRT, '(A)')
      CALL RDREALS (LUNCRT,
     >   ' Enter span stations of primary, secondary, and desired ' //
     >   'sections: ', LUNKBD, NVALS, Z)
      IF (NVALS /= 3) GO TO 600

C ... Do the spanwise linear interpolation for X, Y on each surface,
C     overwriting the primary section data with the results:

      CALL LINTRP (NU2ND, XU, XU2ND, Z (1), Z (2), Z (3), XU)
      CALL LINTRP (NU2ND, YU, YU2ND, Z (1), Z (2), Z (3), YU)
      CALL LINTRP (NL2ND, XL, XL2ND, Z (1), Z (2), Z (3), XL)
      CALL LINTRP (NL2ND, YL, YL2ND, Z (1), Z (2), Z (3), YL)

      RETURN

C ... Error handling:

  830 WRITE (LUNCRT, 1001)
     >   ' LOFT: Error reading secondary coordinates - quitting.'
      GO TO 990

  840 WRITE (LUNCRT, 1001)
     >   ' LOFT: The two sections must have the same RELATIVE point',
     >   '       distributions on corresponding surfaces.',
     >   '       PROFILE''s REDISTRIBUTE and NORMALIZE/DENORMALIZE',
     >   '       options can help.  Quitting...'

  990 STOP

C ... Formats:

 1001 FORMAT (/, (A))

      END SUBROUTINE LOFT
C+---------------------------------------------------------------------
C
      SUBROUTINE MAXMIN (LUNRD, MAXPTS, XMIN, XMAX, YMIN, YMAX,
     >                   YLE, NPR, X, Y, XU, XL, YU, YL, LEGND, IER)
C
C  PURPOSE:  MAXMIN calls PRREAD to scan all of the profiles in the
C            given file.  It counts the number of profiles read and
C            finds the maximum and minimum abscissas and  ordinates
C            over all profiles.  (There may be only one profile.)
C
C  ARGUMENTS:
C   ARG    TYPE   I/O/S    DIM     DESCRIPTION
C  LUNRD     I      I       -      Logical unit used in PRREAD
C  MAXPTS    I      I       -      Maximum number of points allowed
C                                  for on a surface
C  XMIN      R      O       -      Minimum abscissa of profile(s)
C  XMAX      R      O       -      Maximum abscissa of profile(s)
C  YMIN      R      O       -      Minimum ordinate of profile(s)
C  YMAX      R      O       -      Maximum ordinate of profile(s)
C  NPR       I      O       -      Number of profiles read
C  XU        R      S     MAXPTS   Upper surface abscissas
C  XL        R      S     MAXPTS   Lower surface abscissas (if any)
C  YU        R      S     MAXPTS   Upper surface ordinates
C  YL        R      S     MAXPTS   Lower surface ordinates (if any)
C  YLE       R      O       -      Ordinate of leading edge point
C  LEGND     C     S/O      -      Legend entry associated with profile(s).
C                                  The last one found may be used in
C                                  the initial plot setup unless the
C                                  run is a "THREED" case.
C  IER       I      O       -      Error return code.  See PRREAD.
C
C  EXTERNAL REFERENCES:
C  MODULE    DESCRIPTION
C  BOUNDS    Determines maximum and minimum values of array(s)
C  PRREAD    Reads one profile
C
C  AUTHOR: Leslie Collins, Informatics General Corporation, Palo Alto, CA
C
C  HISTORY:
C
C    02/15/83   LJC   Original design and coding
C    03/03/83   LJC   Included indefinite scanning loop to determine
C                     overall data range
C    08/03/84   LJC   Added YLE as an argument (previously determined
C                     in main program)
C    02/11/87   DAS   Functionality of BOUNDS changed somewhat; took
C                     out nonstandard DO WHILE.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   I, IER, LUNRD, MAXPTS, NPR

       REAL
     >    XMAX, XMIN, XL(MAXPTS), XU(MAXPTS), YLE, YMIN, YMAX,
     >    YL(MAXPTS), YU(MAXPTS)

      REAL:: X(*),Y(*)     ! changed by RLC 20210816
      
      CHARACTER
     >   LEGND * 80 

C  *  Local variables:

      INTEGER
     >   FORMAT, NL, NU
      
C  *  Set up indefinite loop to scan all profiles for overall data range 
C     (for axis scaling and possible normalization purposes).  Count the
C     number of profiles as we go:

C  *  Initialize saved minimums and maximums for first comparison:

      XMIN =  1.E10
      XMAX = -XMIN
      YMIN =  XMIN
      YMAX =  XMAX
      NPR  =  0

  200 CONTINUE

C  *     Read a profile:

         CALL PRREAD (LUNRD, LEGND, MAXPTS, NU, NL, X, Y, XU, XL,
     >                YU, YL, FORMAT, IER)
         IF (IER /= 0) GO TO 999

         NPR = NPR+1

C  *     Find minimum and maximum X and Y so far.
C        Note that BOUNDS expects input MAX/MINs.

         CALL BOUNDS (NU, 1, MAXPTS, XU, XMIN, XMAX)
         CALL BOUNDS (NU, 1, MAXPTS, YU, YMIN, YMAX)

         IF (NL /= 0) THEN
            CALL BOUNDS (NL, 1, MAXPTS, XL, XMIN, XMAX)
            CALL BOUNDS (NL, 1, MAXPTS, YL, YMIN, YMAX)
         END IF
         
C  *     Retrieve Y coordinate of point with minimum X (clumsy):

         DO I = 1, NU
            IF (XU(I) == XMIN) THEN
               YLE = YU(I)
               GO TO 900
            END IF
         END DO

         DO I = 1, NL
            IF (XL(I) == XMIN) THEN
               YLE = YL(I)
               GO TO 900
            END IF
         END DO

  900    CONTINUE

      write(*,*) 'DEBUG: completed maxmin, xmin,xmax,ymin,ymax',
     >   xmin,xmax,ymin,ymax 


C  *     Look for another profile:

      GO TO 200
      
  999 CONTINUE    
      write(*,*) 'DEBUG: return from minmax,ier=', ier
      RETURN

      END SUBROUTINE MAXMIN
C+----------------------------------------------------------------------
C
      SUBROUTINE MODIFY (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD, LUNOUT)
C
C  PURPOSE:  MODIFY drives the interactive application of a variety of
C            shape functions to airfoil surfaces.    It simply invokes
C            MODSRF once for each surface.   (Other versions of MODIFY
C            may need to treat the profile as a whole; hence the above
C            calling sequence.)
C            
C  ARGUMENTS:
C     VAR   DIM   I/O/S   DESCRIPTION
C     NU     -      I     Number of points on upper surface
C     XU     NU     I     Abscissas for upper surface
C     YU     NU   I/O     Ordinates for upper surface
C     NL     -      I     Number of points on lower surface
C     XL     NL     I     Abscissas for lower surface
C     YL     NL   I/O     Ordinates for lower surface
C     LUNCRT -      I     Logical unit for screen
C     LUNKBD -      I     Logical unit for keyboard
C     LUNOUT -      I     Logical unit for printed output
C
C  PROCEDURES:
C
C     MODSRF   Perturbs an airfoil surface using shape functions
C
C  HISTORY:
C     09/27/83   LJC   Initial coding
C     10/28/83   DAS   Added LUNs as arguments, for consistency
C                        
C  AUTHOR: Leslie Collins, Informatics, Palo Alto, CA
C
C-------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER  LUNCRT, LUNKBD, LUNOUT, NU, NL
      REAL     XU (NU), XL (NL), YU (NU), YL (NL)
      EXTERNAL MODSRF

      CALL MODSRF (NU, XU, YU, 1, LUNCRT, LUNKBD, LUNOUT)
      CALL MODSRF (NL, XL, YL, 2, LUNCRT, LUNKBD, LUNOUT)
    
      END SUBROUTINE MODIFY
C+-----------------------------------------------------------------------
C
      SUBROUTINE MODSRF (NPTS, X, Y, ISRF, LUNCRT, LUNKBD, LUNOUT)
C
C  PURPOSE:  MODSRF allows the user to perturb one airfoil surface by
C            selecting shape functions interactively.  The airfoil is
C            assumed to be normalized.
C
C  METHOD:   DO WHILE more shape functions desired
C               Prompt for and read no. of shape function from menu
C               DO I = 1, number of parameters for this shape function
C                  Prompt for and read function parameters
C               END DO
C               Prompt for and read function multiplier
C            END DO
C
C            Evaluate each function at all (normalized) abscissas,
C            using stored parameters, and add to ordinates.
C            
C  ARGUMENTS:
C     VAR   DIM   I/O/S   DESCRIPTION
C     NPTS   -      I     Number of points on current surface
C     X     NPTS    I     Abscissas of current surface (normalized)
C     Y     NPTS  I/O     Ordinates of current surface
C     ISRF   -      I     1 means uppper surface;
C                         2 means lower surface.
C     LUNCRT -      I     Logical unit for screen
C     LUNKBD -      I     Logical unit for keyboard
C     LUNOUT -      I     Logical unit for printed output
C
C  FILES USED:  See arguments.
C
C  PROCEDURES:
C     ADDBUMPS   Applies selected bumps to airfoil ordinates
C     GETSHAPE   Modularization of bump selection (used for smoothing too)
C     READER     Prompting utility
C
C  HISTORY:
C     08/05/83   LJC   Initial design and coding
C     09/07/83   LJC   Added echoing of inputs
C     09/27/83   LJC   Replaced statement functions with subroutine BUMP
C     10/04/83   LJC   Introduced READER routine for accepting inputs
C     01/14/84   DAS   Introduced ADDBUMPS in place of in-line code
C     02/22/84   DAS   Eliminated SQRT and SIN bumps (redundant)
C     04/13/84   LJC   Modified to handle ADDBUMPS which expects bump
C                      names now, not integer code numbers
C     12/27/85   DAS   Streamlined prompts and other I/O
C     08/12/86   DAS   Added RAMP, FLAP, and SLAT options
C     01/29/90   DAS   Removed END DOs and DO WHILE
C     11/09/93   DAS   RAK noticed the DROOP fn. was missing a (1-X) factor
C     06/19/96   DAS   The EXPONENTIAL function now has peak 1. at
C                      specified X; introduced symmetric forms of
C                      the modified sine function.
C     12/18/96   DAS   Installed SINF, COSL, COSR, LCOS, & RCOS functions;
C                      no more concatenations in I/O lists.
C     12/24/96   DAS   Introduced GETSHAPE after it was needed for doing
C                      implicit/explicit smoothing. This keeps the shape
C                      function menu in one place (apart from GETBUMPS),
C                      and allows selecting bumps by name.
C                        
C  AUTHOR: Leslie Collins, Informatics/NASA Ames, Mt. View, CA
C
C-------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   ISRF, LUNCRT, LUNKBD, LUNOUT, NPTS
       REAL
     >   X (NPTS), Y (NPTS)

C     Local constants:

      INTEGER, PARAMETER ::
     >   MXNSEL = 20,   ! Max. # function selections per surface
     >   MXPARM = 3     ! Max. # parameters per function incl. multiplier

      CHARACTER, PARAMETER ::
     >   BLANK * 1 = ' '

C     Local variables:

      INTEGER
     >   I, IBUMP, NBUMP, NPARAM
      REAL
     >   BMULT, PARAMS (MXPARM), PARM (MXPARM, MXNSEL)
      CHARACTER
     >   BNAME * 4, BNAMES (MXNSEL) * 4, PNAMES (MXPARM) * 10,
     >   SURFCE (2) * 13
      LOGICAL 
     >   CR, EOF, RETRY, SUPRES

C     Procedures:

      EXTERNAL
     >   ADDBUMPS, GETSHAPE, READY

C     Storage:

      DATA  
     >   SURFCE /'UPPER surface', 'LOWER surface'/

      SAVE
     >   SURFCE

C     Execution:

      BMULT = 0.     ! Non-zero suppresses prompt for a multiplier
      RETRY = .FALSE.

  100 IF (ISRF == 1 .OR. RETRY) THEN
         SUPRES = .FALSE.  ! Don't suppress the shape function menu
         WRITE (LUNCRT, '(A)')
      END IF

      IBUMP = 0
      BNAME = 'DONE'  ! Default is "no more shape functions"
      NBUMP = 0

      CALL GETSHAPE (LUNCRT, LUNKBD, 
     >   'First shape function to add to the ' // SURFCE (ISRF) // '?',
     >   SUPRES, IBUMP, BNAME, NPARAM, PARAMS, PNAMES, BMULT)

      IF (NPARAM == 0) GO TO 900  ! Abort
      IF (IBUMP  == 0) GO TO 999  ! Done

      WRITE (LUNOUT, '(/, 3A)' )
     >   ' Shape functions applied to the ', SURFCE (ISRF), ':'

C     Loop over packing of bumps and prompting for further shape functions:

  200 CONTINUE

         NBUMP = NBUMP + 1
         BNAMES (NBUMP) = BNAME
         DO I = 1, NPARAM
            PARM (I, NBUMP) = PARAMS (I)
         END DO

C        Echo selection to output file:

         WRITE (LUNOUT, '(/, 4X, 2A)') BNAME, ':'
         WRITE (LUNOUT, '(4X, A, F10.6)')
     >      (PNAMES (I), PARAMS (I), I = 1, NPARAM)

         IF (NBUMP < MXNSEL) THEN

C           Prompt for the next shape function:

            SUPRES = .TRUE.  ! Suppress the shape function menu
            IBUMP  = 0
            BNAME  = 'DONE'

            CALL GETSHAPE (LUNCRT, LUNKBD, 
     >         'Next shape function for the ' // SURFCE (ISRF) // '?',
     >         SUPRES, IBUMP, BNAME, NPARAM, PARAMS, PNAMES, BMULT)

            IF (NPARAM == 0) GO TO 900  ! Abort
            IF (IBUMP  >  0) GO TO 200

         END IF


C     Apply the selected bumps to this surface, in-place:

      CALL ADDBUMPS (NPTS, X, Y, MXPARM, NBUMP, BNAMES, PARM, Y)

      GO TO 999
      

  900 CONTINUE  ! Exit or start over:

      RETRY = .TRUE.
      CALL READY (LUNCRT, 
     >   'Do you want to start this surface over? (Y/N/EOF; <CR>=Yes) ',
     >   LUNKBD, RETRY, CR, EOF)
      IF (EOF)   GO TO 999
      IF (RETRY) GO TO 100


  999 WRITE (LUNOUT, '(A)')

      END SUBROUTINE MODSRF
C+------------------------------------------------------------------------------
C
      SUBROUTINE NOSEJOB (NU, XU, YU, NL, XL, YL, MAXPTS, X, Y, 
     >                    COEFS, LUNCRT, LUNKBD, LUNTAB, MODE)
C
C  PURPOSE:
C
C        NOSEJOB performs various modifications to airfoil leading edges.
C     Initially, it has the following options (others may arise):
C
C        > Round off a sharp leading edge.
C        > Sharpen a rounded leading edge.
C
C        After the leading edge is modified, the original chord and thickness
C     could in principle be retrieved here through in-line reuse of PROFILE's
C     REFINE module - maybe some day.
C
C  METHOD:
C
C        A menu is presented for the various options.  The input coordinates
C     are overwritten by the modified coordinates.  The number of points on
C     each surface is held the same.  In the modified nose region, the relative
C     point distribution in terms of arc length is also held the same.  The
C     rounding and sharpening options are thus as reversible as possible,
C     especially if the option to retrieve original chord and thickness is not
C     used.  The user is asked to specify the points I1, I2 for each surface at
C     which the modified surface should blend with the original.  
C
C                                      x
C                      I2     x
C                        *
C                     x+
C                   x +
C                 x  +
C                   x +
C                     x+
C                        *
C                      I1     x
C                                      x
C
C        Rounding is achieved by rotating the X axis to be parallel to the
C     bisector of the angle between the tangents at I1 and I2, then fitting a
C     conventional spline to the region near the nose where the "abscissas"
C     (rotated Y coordinates) are monotonic.  Evaluating this spline at points
C     spaced along the arc in the same relative way as the original points is
C     awkward!  Conservatively careful approximations are achieved by evaluating
C     the rounded surface at ~30 points spaced uniformly in the abscissa,
C     refitting those points parametrically, and working with cumulative chord
C     lengths from there.  Local storage is used for this set of gyrations.
C
C        Sharpening the nose is achieved by determining where the extrapolated
C     surfaces intersect, offering that point to the user as the default for
C     the new leading edge, then allowing an alternative point to be entered
C     interactively.  Much the same tedious steps as above are then taken to
C     preserve the original point distribution along the modified portions of
C     each surface, although local spline LCSFIT is used instead of PSFIT,
C     because it should be quite adequate for sharp noses.
C
C  ARGUMENTS:
C   ARG      DIM    TYPE I/O/S DESCRIPTION
C  NU,NL      -       I    I   Number of upper/lower surface pts.
C  XU,XL    NU,NL     R   I/O  Abscissas, upper/lower.
C  YU,YL    NU,NL     R   I/O  Corresponding ordinates, upper/lower.
C  X,Y     NU+NL-1    R    S   Storage for wraparound form of airfoil.
C  MAXPTS     -       I    I   Max. no. of pts. provided for on 1 surface.
C  COEFS   MAXPTS*3   R    S   Coefficients used by conventional spline.
C  LUNCRT     -       I    I   Logical unit for prompts.
C  LUNKBD     -       I    I   Logical unit for responses.
C  LUNTAB     -       I    I   Logical unit for printable record of the
C                              iterations performed by the REFINE option,
C                              if it is ever installed.
C  MODE       -       I    O   The menu choice, in case it is needed by
C                              the calling program.
C
C  PROCEDURES:
C    CHORD         Chord-length utility
C    COPY, RVERSE  Data transfer utilities         
C    CSEVAL        Evaluates the spline previously fit by CSFIT
C    CSFIT         Fits a conventional cubic interpolating spline
C    FDCNTR        1st and 2nd derivative by central differencing
C    INTSEC2       Used to find where extrapolated surfaces meet
C    LCSFIT        Local spline, used for the "sharpen" option
C    PSFIT         Parametric form of CSFIT
C    PSTVAL        Evaluates PSFIT's spline at specified values of T
C    ROTATE2D      Rotates point(s) (X, Y) about (P, Q)
C    READER        Prompting utility
C    TSUBJ         Gives arc-length T associated with Jth data point for PSFIT
C    XGRID         Used to generate a uniform distribution
C
C  HISTORY:
C  10/21/91  DAS  Initial implementation, starting with a copy of REDISTRIB.
C                 No use of REFINE yet.  Difficulties with the leading edge
C                 (no longer necessarily the foremost point) would require
C                 RECTIFY as well, and this conflicts with reversing the
C                 operation.  Pursue it again some day, perhaps.
C  03/29/92  DAS  Use with I1=I2=20 revealed that TORIG(NLOCAL=30) needed to
C                 be TORIG(2*NLOCAL) for the "round" option.
C  08/01/93  DAS  INTSEC2 had CALCT1, CALCT2 arguments added.
C  04/03/95  DAS  INTSEC2 had TOL argument added.
C 
C  AUTHORS: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNTAB, MAXPTS, MODE, NL, NU

      REAL
     >   COEFS (MAXPTS, 3), X (MAXPTS * 2), XL (NL), XU (NU),
     >   Y (MAXPTS * 2), YL (NL), YU (NU)

C     Local constants:

      INTEGER
     >   MXMENU, NEND, NLOCAL
      REAL
     >   TOL
      CHARACTER
     >   BLANK * 1
      LOGICAL
     >   CALCT1, CALCT2
      PARAMETER
     >  (BLANK   = ' ',
     >   CALCT1  = .TRUE., ! Switches for INTSEC2
     >   CALCT2  = .TRUE.,
     >   TOL     = 1.E-6,  ! Used by INTSEC2 relative to the curve lengths
     >   MXMENU  = 2,
     >   NEND    = 4,      ! Number of original points included at each
                           ! end of the NLOCAL points used to refit the
                           ! nose region parametrically (round option)
     >   NLOCAL  = 40)     ! Number of evaluations in local storage of
                           ! the X vs. Y spline used for rounding, in
                           ! order to translate to arc-lengths carefully
C     Local variables:

      INTEGER
     >   I, I1, I1R, I2, IER, J, J1, J2, K, K1, K2, M, NMONO, NNOSE,
     >   NPTS, NTEMP

      REAL
     >   DUMMY, S1, S2, SCALE, THETA1, THETA2, THETAM, TOFFSET,
     >   XLE, YLE,
     >   TLOCAL (NLOCAL), TORIG (2*NLOCAL), XLOCAL (2*NLOCAL),
     >   YLOCAL (2*NLOCAL)

      LOGICAL
     >   DEFAULT, QUIT

      CHARACTER
     >   MENU (0 : MXMENU) * 37

C     Procedures:

      REAL
     >   CHORD, TSUBJ

      EXTERNAL
     >   CHORD, COPY, CSEVAL, CSFIT, FDCNTR, INTSEC2, LCSFIT, PSFIT,
     >   PSTVAL, READI, READR, ROTATE2D, RVERSE, TSUBJ, XGRID

C     Storage:

      DATA
     >   MENU /
     >   ' PROFILE''s "nose-job" options are:',
     >   '   (1) Round off a sharp leading edge',
     >   '   (2) Sharpen a rounded leading edge'/


C     Execution:

      WRITE (LUNCRT, '(/, A, /, A, A, /)') MENU
  100 M = 0
      CALL READI (LUNCRT, 'What''ll it be? ', LUNKBD, M, DEFAULT, QUIT)
      IF (QUIT) GO TO 999
      IF (M < 1 .OR. M > MXMENU) GO TO 100
      MODE = M

C     Round and sharpen options both need points specified at which to blend.

      I2 = 0
  110 WRITE (LUNCRT, '(A)')
      CALL READI (LUNCRT,
     >   'Index of first UPPER surface point to retain: ',
     >   LUNKBD, I2, DEFAULT, QUIT)
      IF (QUIT) GO TO 999
      IF (I2 <= 1 .OR. I2 >= NU) GO TO 110

      I1 = 0
  120 CALL READI (LUNCRT,
     >   'Index of first LOWER surface point to retain: ',
     >   LUNKBD, I1, DEFAULT, QUIT)
      IF (QUIT) GO TO 999
      IF (I1 <= 1 .OR. I1 >= NL) GO TO 120

      IF (I1 > NLOCAL .OR. I2 > NLOCAL) GO TO 800


C                              ---------------
C                              Rounding option
C                              ---------------

      IF (M == 1) THEN

C        Determine the original relative locations along the arc:

         TORIG (1) = 0.
         J = I1
         DO I = 2, I1
            TORIG (I) = CHORD (XL, YL, J, J - 1) + TORIG (I - 1)
            J = J - 1
         END DO

         I = I1
         DO J = 2, I2
            I = I + 1
            TORIG (I) = CHORD (XU, YU, J - 1, J) + TORIG (I - 1)
         END DO
            
         NNOSE = I      ! Number of nose points from I1 to I2 inclusive

C        Determine the gradients at I1 and I2 (finite differencing):

         CALL FDCNTR (I1, XL, YL, S1, S1)
         CALL FDCNTR (I2, XU, YU, S2, S2)

         THETA1 = ATAN (S1)
         THETA2 = ATAN (S2)
         THETAM = (THETA1 + THETA2) * 0.5  ! Mean angle, for rotating X axis to
         THETAM = THETAM * 45.0 / ATAN (1.0)

C        Set up the airfoil as a single wraparound curve, eliminating the
C        indicated leading edge region:

         CALL RVERSE (NL, XL, X)
         CALL RVERSE (NL, YL, Y)

         I1R   = NL - I1 + 1
         NTEMP = NU - I2 + 1
         CALL COPY (NTEMP, XU (I2), X (I1R + 1))
         CALL COPY (NTEMP, YU (I2), Y (I1R + 1))

         NPTS = NU + NL - I1 - I2 + 2   ! No. of points left after removing
                                        ! the points to be changed

C        Rotate the airfoil so the X axis is parallel to the angle bisector.
C        The leading edge will do for the center of rotation (arbitrary).

         CALL ROTATE2D (NPTS, X, Y, -THETAM, XU (1), YU (1))

C        Locate the nose portion where the transformed Ys are monotonic:

         DO I = I1R, 2, -1              ! Lower surface
            IF (Y (I - 1) >= Y (I)) THEN
               J1 = I
               GO TO 220
            END IF
         END DO
         J1 = 1

  220    DO I = I1R + 1, NPTS - 1      ! Upper
            IF (Y (I + 1) <= Y (I)) THEN
               J2 = I
               GO TO 240
            END IF
         END DO
         J2 = NPTS

  240    NMONO = J2 - J1 + 1

C        Spline transformed X vs. transformed Y for the nose region:

         CALL CSFIT (NMONO, Y (J1), X (J1), 0, DUMMY, 0, DUMMY,
     >               COEFS (1, 1), COEFS (1, 2), COEFS (1, 3), IER)
         IF (IER /= 0) GO TO 810

C        In order to work with arc lengths and with decent resolution,
C        evaluate this spline at sufficient points for reliable parametric
C        reinterpolation.  We need to include some of the unchanged points,
C        yet we need to include points I1 and I2 exactly in order to find
C        the new arc length (for redistributing it as for the original points).
C        We can't be certain of the number of points beyond J1, J2 that are
C        monotonic in transformed Y.  This is messy!

         K1 = MAX (I1R - NEND, J1)     ! +/-4 should control ends adequately.
         K2 = MIN (I1R + 1 + NEND, J2) ! There may not always be this many.

         J = 0
         DO I = K1, I1R - 1            ! Low end
            J = J + 1
            YLOCAL (J) = Y (I)
         END DO

         K = NLOCAL + 1
         DO I = K2, I1R + 2, -1        ! High end
            K = K - 1
            YLOCAL (K) = Y (I)
         END DO

C        Now for the middle: uniform in the Y range of transformed points
C        I1R and I1R + 1.

         NTEMP = NLOCAL - J - (NLOCAL + 1 - K)
         J = J + 1    ! These are now the YLOCAL (*) indices corresponding
         K = K - 1    ! to data points I1R, I1R + 1

         CALL XGRID (NTEMP, 0, Y (I1R), Y (I1R + 1), YLOCAL (J))

C        Evaluate the conventional spline at these Ys.

         CALL CSEVAL (NMONO, Y (J1), X (J1), NLOCAL, YLOCAL,
     >                COEFS (1, 1), COEFS (1, 2), COEFS (1, 3), XLOCAL)

C        Now refit the rounded curve parametrically:

         CALL PSFIT (NLOCAL, YLOCAL, XLOCAL, 'C', .FALSE., IER)
         IF (IER /= 0) GO TO 820

C        ... and evaluate it at the same relative arc lengths as the originals:

         TOFFSET = TSUBJ (J)  ! TSUBJ extracts an arc length from PSFIT's COMMON

         SCALE = (TSUBJ (K) - TOFFSET) / TORIG (NNOSE)
         DO I = 2, NNOSE
            TORIG (I) = TORIG (I) * SCALE + TOFFSET
         END DO

         CALL PSTVAL (NNOSE, TORIG, Y (I1R), X (I1R), Y (I1R), X (I1R),
     >                Y (I1R), X (I1R), NLOCAL, YLOCAL, XLOCAL)

C        Rotate the new points back to the original coordinate system:

         CALL ROTATE2D (NNOSE, X (I1R), Y (I1R), THETAM, XU (1), YU (1))

C        Finally, overwrite the original points with the rounded ones:

         J = I1R
         DO I = I1 - 1, 1, -1
            J = J + 1
            XL (I) = X (J)
            YL (I) = Y (J)
         END DO

         DO I = 1, I2 - 1
            XU (I) = X (J)
            YU (I) = Y (J)
            J = J + 1
         END DO

         IF (XL (2) <= XL (1) .OR. XU (2) <= XU (1)) THEN
            WRITE (LUNCRT, 1005)
         END IF

C                             -----------------
C                             Sharpening option
C                             -----------------

      ELSE IF (M == 2) THEN

C        Determine where the extrapolated surfaces meet.  INTSEC2's use
C        of local splines means the first four points per curve suffice.

         J1 = 1  ! Estimate of index near point of intersection
         J2 = 1
         CALL INTSEC2 (4, XL (I1), YL (I1), TORIG (1), J1, CALCT1,
     >                 4, XU (I2), YU (I2), TORIG (5), J2, CALCT2, TOL,
     >                 XLE, YLE, -LUNCRT, IER)
         IF (IER /= 0) GO TO 830

         WRITE (LUNCRT, 1004) XLE, YLE
         CALL READR (LUNCRT, 'Different X?  <CR> = above value: ',
     >               LUNKBD, XLE, DEFAULT, QUIT)
         IF (QUIT) GO TO 999
         CALL READR (LUNCRT, 'Different Y?  <CR> = above value: ',
     >               LUNKBD, YLE, DEFAULT, QUIT)
         IF (QUIT) GO TO 999


C        Lower surface first:

C        Determine arc-length distribution of points to be replaced.

         TORIG (1) = 0.
         DO I = 2, I1
            TORIG (I) = CHORD (XL, YL, I - 1, I) + TORIG (I - 1)
         END DO

C        The new leading edge and the first 3 retained points suffice
C        for interpolation in the extrapolated region by local methods.

         XLOCAL (1) = XLE
         YLOCAL (1) = YLE
         TLOCAL (1) = 0.
         I = I1
         DO J = 2, 4
            XLOCAL (J) = XL (I)
            YLOCAL (J) = YL (I)
            TLOCAL (J) = CHORD (XLOCAL, YLOCAL, J - 1, J) + TLOCAL (J-1)
            I = I + 1
         END DO

C        For a sharp leading edge, one chord is near enough to the arc length:

         SCALE = TLOCAL (2) / TORIG (I1)
         DO I = 2, I1
            TORIG (I) = TORIG (I) * SCALE
         END DO

C        Interpolate at the equivalent arc lengths:

         CALL LCSFIT (4, TLOCAL, XLOCAL, .TRUE., 'B', I1 - 1, TORIG,
     >                XL, XL)
         CALL LCSFIT (4, TLOCAL, YLOCAL, .TRUE., 'B', I1 - 1, TORIG,
     >                YL, YL)


C        Repeat for upper surface:

C        Determine arc-length distribution of points to be replaced.

         DO I = 2, I2
            TORIG (I) = CHORD (XU, YU, I - 1, I) + TORIG (I - 1)
         END DO

         I = I2
         DO J = 2, 4
            XLOCAL (J) = XU (I)
            YLOCAL (J) = YU (I)
            TLOCAL (J) = CHORD (XLOCAL, YLOCAL, J - 1, J) + TLOCAL (J-1)
            I = I + 1
         END DO

         SCALE = TLOCAL (2) / TORIG (I2)
         DO I = 2, I2
            TORIG (I) = TORIG (I) * SCALE
         END DO

         CALL LCSFIT (4, TLOCAL, XLOCAL, .TRUE., 'B', I2 - 1, TORIG,
     >                XU, XU)
         CALL LCSFIT (4, TLOCAL, YLOCAL, .TRUE., 'B', I2 - 1, TORIG,
     >                YU, YU)

      END IF


      RETURN

C     Error handling:

  800 WRITE (LUNCRT, 1002) I1, I2, NLOCAL
      GO TO 110
  810 WRITE (LUNCRT, 1003) 'CSFIT', IER
      GO TO 999
  820 WRITE (LUNCRT, 1003) 'PSFIT', IER
      GO TO 999
  830 WRITE (LUNCRT, 1003) 'INTSEC2', IER
C*****GO TO 999

  999 WRITE (LUNCRT, 1001) ' Stopping in NOSEJOB.'
      STOP

C     Formats:

 1001 FORMAT (/, A)
 1002 FORMAT (/, ' Index ', I2, ' and/or ', I2,
     >        ' exceeds local limit of ', I2, '.  Try again.')
 1003 FORMAT (/, ' NOSEJOB:  Bad return from ', A, '.  IER: ', I3)
 1004 FORMAT (/, ' Estimate of leading edge by extrapolation:', /,
     >        ' X = ', G13.6, '   Y = ', G13.6, /)
 1005 FORMAT (/, ' WARNING:  The common leading edge point is no',
     >        ' longer foremost.', /,
     >        ' This may be necessary if you intend to reverse the',
     >        ' operation.', /,
     >        ' But some results such as derivatives are affected. ',
     >        ' Proceeding ...')

      END SUBROUTINE NOSEJOB
C+----------------------------------------------------------------------
C
      SUBROUTINE NRMLIZ (NPTS, XYIN, XYOUT, XYLE, CHORD)
C
C  PURPOSE:
C     NRMLIZ normalizes X or Y coordinates (presumably from airfoils)
C     according to the given chord and leading edge coordinate.  It will
C     DE-normalize if the given chord is negative.  In-place is safe.
C     The relevant formulas are indicated by
C
C        X(norm)   = (X - X(l.e.)) / Chord                 and
C        X(denorm) = X(l.e.) + X(norm) * Chord
C
C  ARGUMENTS:
C   ARG    TYPE  I/O/S   DIM     DESCRIPTION
C   NPTS    I      I      -      Number of coordinates (X or Y)
C   XYIN    R      I    NPTS     Data to be normalized/de-normalized
C   XYOUT   R      O    NPTS     Results. May be same array as XYIN.
C   XYLE    R      I      -      X or Y coordinate of leading edge
C   CHORD   R      I      -      CHORD>0 means normalize by this chord;
C                                CHORD<0 means de-normalize by -CHORD.
C
C  AUTHORS: Various (Informatics, 1983).
C
C   10/21/88  D. Saunders  Ensured 1.0 exactly, where appropriate.
C
C-----------------------------------------------------------------------

C     Arguments:

      REAL XYIN(NPTS), XYOUT(NPTS)

C     Constants:

      REAL, PARAMETER ::
     >   ONE = 1., ZERO = 0.

C     Execution:

C  *  Use local copy of leading edge coordinate, to protect use of
C     an XYIN(*) element passed as this argument:

      XYLEAD = XYLE

      IF (CHORD > ZERO) THEN  ! Normalize coordinates:

         RCHORD = ONE / CHORD
         DO I = 1, NPTS
            XYOUT(I) = (XYIN(I) - XYLEAD) * RCHORD
         END DO

         IF (ABS (XYOUT(NPTS) - ONE) < 1.E-6) XYOUT(NPTS) = ONE

      ELSE  ! De-normalize coordinates:

         DO I = 1, NPTS
            XYOUT(I) = XYLEAD - XYIN(I) * CHORD
         END DO

      END IF

      END SUBROUTINE NRMLIZ
C+----------------------------------------------------------------------
C
      SUBROUTINE NRMSET (LUNCRT, LUNKBD, XMINSV, XMAXSV,
     >                   CNORML, XLE, YLE)
C
C  PURPOSE: NRMSET determines values to normalize/denormalize the
C           profile(s) by - introduced to keep the prompting out of
C           the main program.
C           
C  ARGUMENTS:
C   ARG    TYPE  I/O/S   DIM  DESCRIPTION
C   LUNCRT  I    I        -   Logical unit for prompts (screen).
C   LUNKBD  I    I        -   Logical unit for responses (keyboard).
C   XMINSV, R    I        -   Minimum and maximum abscissas found over
C   XMAXSV  R    I        -   all profiles.
C   CNORML  R      O      -   Chord; either input at the terminal or
C                             determined by the data range.
C   XLE     R      O      -   Abscissa of leading edge point; either
C                             input from the terminal or determined from
C                             the data range checking.
C   YLE     R    I/O      -   Ordinate of leading edge point; input as
C                             Y corresp. to XMINSV; may be updated here.
C
C  PROCEDURES:
C  READER   Prompts for and reads integer, real, etc.
C   
C  AUTHOR: Leslie Collins, Informatics General, Palo Alto, CA
C
C  HISTORY:
C  08/07/84   LJC   Initial design and coding.
C  10/09/85   DAS   Eliminated normalization of plot scale info.
C  02/07/90   DAS   Eliminated list-directed I/O.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE
      
C  *  Arguments:

      INTEGER
     >   LUNCRT, LUNKBD
      REAL
     >   CNORML, XLE, XMAXSV, XMINSV, YLE

C  *  Local variables:

      LOGICAL
     >   DEFAULT, QUIT
      REAL
     >   CHORD

C  *  Procedures:

      EXTERNAL
     >   READR

C  *  Execution:

      CHORD = XMAXSV - XMINSV

      WRITE (LUNCRT, 1001) 'Current chord: ', CHORD,
     >   'Enter <CR> to normalize by this chord, or enter a'
      CALL READR (LUNCRT,
     >   'different chord; (a negative value de-normalizes): ',
     >   LUNKBD, CNORML, DEFAULT, QUIT)
      IF (QUIT) GO TO 999

      IF (DEFAULT) THEN

         CNORML = CHORD
         XLE = XMINSV

      ELSE ! YLE is already defined by MAXMIN called from main program.

         WRITE (LUNCRT, 1001)
         CALL READR (LUNCRT, 'Enter abscissa of leading edge point: ',
     >      LUNKBD, XLE, DEFAULT, QUIT)
         IF (QUIT) GO TO 999

         CALL READR (LUNCRT, 'Enter ordinate of leading edge point: ',
     >      LUNKBD, YLE, DEFAULT, QUIT)
         IF (QUIT) GO TO 999

      END IF

      RETURN

  999 WRITE (LUNCRT, 1001) 'Stopping as requested.'
      STOP

 1001 FORMAT (1X, A, G13.6)

      END SUBROUTINE NRMSET
C+----------------------------------------------------------------------
C
      SUBROUTINE OPTIMIZE (NU, XU, YU, NL, XL, YL, TITLE, 
     +                     LUNCRT, LUNKBD, LUNPRT, LUNBUMPS, LUNTARG,
     +                     UPPER, NXTARG, TARGX, TARGCRV)
C
C  PURPOSE:  OPTIMIZE  perturbs  one  surface of an airfoil by applying
C            given shape functions to the surface and  optimizing  some
C            of the bump parameters so as to match a given target curv-
C            ature distribution in the least squares sense.
C
C  METHOD:   No attempt is made to optimize both surfaces  in  the same
C            run - separate runs are likely to be more manageable, from
C            both the user's standpoint and the programmer's.
C
C            The chosen bump function set is read from an editable disk
C            file rather than entered interactively.   Repeated running
C            would be too busy otherwise.   Values for ALL of the shape
C            function parameters must be in this file, some as starting
C            guesses for the optimizing algorithm, the others inactive.
C
C            Provision is made for using the first N  Wagner  functions
C            more easily than by reading a text file, since this common
C            case is easily generated.
C            
C            The basic procedure follows:
C
C            *  Prompt for which surface to optimize.
C
C            *  Prompt for and read the previously-prepared bump set.
C
C            *  Prompt for and read the target curvature  distribution,
C               assumed to be in simple QPLOT format.
C
C            *  Prompt for any additional constraints (thickness?...).
C
C            *  Set up optimizer  QNMDIF  (initial function evaluation,
C               tolerances, and estimation of optimal finite difference
C               intervals, etc.).
C
C            *  Minimize the sum-of-squares-type objective function.
C
C            *  Update the airfoil surface using the optimal bumps, in-
C               place as expected by the calling program which performs
C               any plotting/tabulating/saving of results requested.
C
C  ARGUMENTS:
C  ARG    DIM  TYPE I/O/S DESCRIPTION
C  NU      -    I     I   Number of upper surface points
C  XU     NU    R     I   Upper surface abscissas
C  YU     NU    R    I/O  Lower surface ordinates
C  NL      -    I     I   Number of lower surface points
C  XL     NL    R     I   Lower surface abscissas
C  YL     NL    R    I/O  Lower surface ordinates
C  TITLE   *    C     I   Title for printout generated by this run (LUNPRT)
C  LUNCRT       I     I   Logical unit for screen
C  LUNKBD       I     I   Logical unit for keyboard
C  LUNPRT       I     I   Logical unit for printout from OPTIMIZE & QNMDIF
C  LUNBUMPS     I     I   Logical unit for previously-prepared bump data
C  LUNTARG      I     I   Logical unit for target curvature distribution
C  UPPER   -    L     O   TRUE if target curvature applies to upper surface
C  NXTARG  -    I     O   Number of target curvature points (for plotting)
C  TARGX NXTARG R     O   Target curvature distribution (needed in COMMON,
C  TARGCRV " "  R     O   but returned to calling program for plotting - an
C                         approach that does not clutter the calling program
C                         with the localized COMMON block).
C
C PARAMETER CONSTANTS:
C    PARAM   TYPE   DESCRIPTION
C  LUNERR      I    Logical unit number for error messages.
C  MXBUMPS     I    Maximum no. of bump functions provided for.
C  MXOPT       I    Maximum no. of optimization variables allowed.
C  MXPARM      I    Maximum no. of parameters associated with any of
C                   the bump functions (including a multiplier).
C  MXSURF      I    Maximum no. of data points handled per surface.
C  MXTARG      I    Maximum no. of target curvature values handled.
C
C COMMON BLOCKS USED:
C
C   /ACTUAL/  (Current values corresponding to input opt. variables)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C   YCALC     MXSURF     R    S   Perturbed airfoil surface
C  Y1CALC     MXSURF     R    S   1st and 2nd derivatives - by-prod-
C  Y2CALC     MXSURF     R    S   ucts of the curvature calculations
C  YKCALC     MXSURF     R    S   Curvature distribution correspond-
C                                 ing to the current bump variables
C
C   /BCHARS/  (Bump function names - can't go in /BUMPS/ but should)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C  BNAMES     MXBUMPS  C*11   S   Names of bumps found by GETBUMPS
C
C   /BUMPS /  (List of bump function parameters, etc.)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C  NBUMPS       -        I    S   Number of bump functions active
C  PARAMS MXPARM,MXBUMPS R    S   Complete set of parameters for the
C                                 active bump functions,   including
C                                 some possible unused ones  present
C                                 because a 2-dim. array approach is
C                                 used for storing them  (as opposed
C                                 to packing the variable numbers of
C                                 parameters associated with differ-
C                                 ent shape functions).  See ACTIVE.
C  ACTIVE MXPARM,MXBUMPS L    S   Marks the bump function parameters
C                                 as active or inactive.  The unused
C                                 ones must be marked inactive along
C                                 with the fixed (but used) ones.
C  VSCALE MXPARM,MXBUMPS R    S   Scale factors needed by the optim-
C                                 izing algorithm.
C
C   /ORIGNL/  (Copies of the original airfoil surface data)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C    NXY        -        I    S   Number of points defining  surface
C   XORIG      NXY       R    S   Abscissas for the surface
C   YORIG      NXY       R    S   Ordinates for the unperturbed srf.
C
C   /TARGET/  (Target distribution quantities)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C   NTARG       -        I    S   Number of target data points
C   XTARG     NTARG      R    S   Abscissas for the target data
C  YKTARG     NTARG      R    S   Target curvature values
C  WEIGHT     NTARG      R    S   Weights to be applied (multiplica-
C                                 tively) to each element of the sum
C                                 of squares.  May be all 1s.
C
C  SIGNIFICANT LOCAL VARIABLES:
C  VAR       DIM     TYPE   DESCRIPTION
C  D        MXOPT      R    Diagonal factor of Hessian approximation
C  G        MXOPT      R    Objective function gradient approximation
C  H        MXOPT      R    Finite difference intervals
C  L MXOPT*(MXOPT-1)/2 R    Lower triangle factor of Hessian approx.
C  OPTVARS  MXOPT      R    Active bump function variables
C  PNAMES   MXPARM*  C*6    Names of bump function parameters found
C           MXBUMPS
C
C  FILES USED:
C  LUN      I/O/S   DESCRIPTION
C  LUNBUMPS   I     Previously-prepared bump function data
C  LUNCRT     O     User prompts
C  LUNKBD     I     User responses
C  LUNPRT     O     Optimizer printout
C  LUNTARG    I     Target curvature distribution (QPLOT format)
C
C  EXTERNAL REFERENCES:
C  ACTIVATE   Packs or unpacks active variables; scales/unscales too.
C  CALCTHICK  Calculates airfoil thickness.
C  CENDIF     Estimates optimal finite difference intervals.
C  CFDISTRIB  Function needed by QNMDIF - compares distributions and
C             returns corresponding sum of squares.
C  COPY       Copies one array to another.
C  DUMMY      "USER" routine called by QNMDIF - do nothing in this case.
C  GETBUMPS   Preprocesses bump function data read from disk file.
C  LSTFIL     Echoes target curvature file to printable file.
C  OPENER     File opening utility.
C  PRTBUMPS   Prints bump function data.
C  PRWRIT     Used to echo target curvature distribution found.
C  QNMDIF     General purpose minimizer not requiring derivatives.
C  READER     Prompting utility.
C  RDQPL      Reads target curvature distribution in QPLOT format.
C
C  HISTORY:
C  01/18/84   DAS   Initial design and coding.
C  01/25/84   DAS   Added scaling, printing of bump data, etc.
C  01/27/84   DAS   Target data now returned as arguments for plotting
C  01/30/84   LJC   Provided for new CALCTHICK (handles NU/=NL)
C  03/16/84   DAS   Provided for constraining thickness (penalty fun.)
C  10/24/85   DAS   Used LSTFIL in place of PRWRIT for echoing target data
C  08/19/86   DAS   Calling sequence to CENDIF was revised
C  01/29/90   DAS   Removed END DOs, underscores, list-directed I/O, and
C                   concatenations in I/O lists (for IRIS 4D purposes);
C                   installed OPENER.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C ... Arguments:

      INTEGER   LUNBUMPS, LUNCRT, LUNKBD, LUNPRT, LUNTARG, NL, NU,
     >          NXTARG
      REAL      TARGX(*), TARGCRV(*), XL(NL), XU(NU), YL(NL), YU(NU)
      LOGICAL   UPPER
      CHARACTER TITLE*(*)

C ... Parameter constants and global variables:

C ... Quantities passed to objective function routine via COMMON
C     because of fixed calling sequence expected by QNMDIF.  Any
C     changes affecting these COMMONs should also be made in one
C     other routine - CFDISTRIB.

      INTEGER, PARAMETER ::
     >   MXBUMPS = 20,  MXPARM = 3, MXOPT = MXPARM*MXBUMPS,
     >   MXSURF  = 180, MXTARG = MXSURF, LUNERR = 6

      REAL            YCALC, Y1CALC, Y2CALC, YKCALC
      COMMON /ACTUAL/ YCALC(MXSURF), Y1CALC(MXSURF), Y2CALC(MXSURF),
     >                YKCALC(MXSURF)

      CHARACTER*11    BNAMES
      COMMON /BCHARS/ BNAMES(MXBUMPS)

      INTEGER         NBUMPS
      REAL            PARAMS, VSCALE
      LOGICAL         ACTIVE
      COMMON /BUMPS / PARAMS(MXPARM,MXBUMPS), VSCALE(MXOPT),
     >                ACTIVE(MXPARM,MXBUMPS), NBUMPS

      INTEGER         NOTHER, NXY
      REAL            XORIG, XOTHER, YORIG, YOTHER
      COMMON /ORIGNL/ XORIG(MXSURF), YORIG(MXSURF),
     >                XOTHER(MXSURF), YOTHER(MXSURF), NOTHER, NXY

      INTEGER         NTARG
      REAL            OPTWRK, PENLTY, TARGTH, XTARG, YKTARG, WEIGHT
      COMMON /TARGET/ XTARG(MXTARG), YKTARG(MXTARG), WEIGHT(MXTARG),
     >                TARGTH, PENLTY, OPTWRK(4*MXSURF), NTARG

C ... Local variables (more below):

      INTEGER   I, IER, NOPTVARS
      REAL      D(MXOPT), G(MXOPT), H(MXOPT), L( MXOPT*(MXOPT-1)/2 ),
     >          OPTVARS(MXOPT), THICKNESS, XATMAX
      LOGICAL   DEFAULT, QUIT
      CHARACTER FILENAME*50, PNAMES(MXPARM, MXBUMPS)*6,
     >          SURFACE*1, TEXT*80, YESNO*1

C ... The following locals are required by QNMDIF:

      INTEGER   NFCEN, NFTOTL, NITER, NLDIM, NTYPE
      REAL      EPSMCH, EPSOBJ, ETA, SSQ, SSQMIN, STEPMX, TOL
      LOGICAL   UNITL, LOCAL, CONV, CONTIN, RESCUE, PRINT

C ... Procedures:

      EXTERNAL  ACTIVATE, CALCTHICK, CENDIF, CFDISTRIB, COPY,
     >          DUMMY, GETBUMPS, LSTFIL, OPENER, PRTBUMPS, QNMDIF,
     >          RDQPL, READC, READR, READY

C ... Execution:


      IF (MAX (NU, NL) > MXSURF) THEN
         WRITE (LUNCRT, '(/, A)')
     >      ' Not enough work-space. Recompile OPTIMIZE, CFDISTRIB.'
      END IF

      WEIGHT = 1.E+0

C ... Prompt for the surface being perturbed, and copy it to where
C     the objective function routine can get at it:

  200 CALL READC (LUNCRT, 'Enter U(pper) or L(ower) to indicate'//
     +            ' the surface being optimized: ',
     +            LUNKBD, SURFACE, DEFAULT, QUIT)

      IF (QUIT) THEN
         WRITE (LUNCRT, '(/, A)') ' Stopping at user request.'
         STOP
      END IF

      IF (SURFACE == 'L') THEN
         UPPER = .FALSE.
         NXY = NL
         CALL COPY (NL, XL, XORIG)
         CALL COPY (NL, YL, YORIG)
         NOTHER = NU
         CALL COPY (NU, XU, XOTHER)
         CALL COPY (NU, YU, YOTHER)
      ELSE IF (SURFACE == 'U') THEN
         UPPER = .TRUE.
         NXY = NU
         CALL COPY (NU, XU, XORIG)
         CALL COPY (NU, YU, YORIG)
         NOTHER = NL
         CALL COPY (NL, XL, XOTHER)
         CALL COPY (NL, YL, YOTHER)
      ELSE
         WRITE (LUNCRT, 1001) 'Invalid response. Try again.'
         GO TO 200
      END IF

C ... Get the original thickness ratio, and a target thickness if any:

      CALL CALCTHICK (NL, NU, XL, XU, YL, YU, THICKNESS,
     >                XATMAX, OPTWRK, OPTWRK(MXSURF+1),
     >                OPTWRK(2*MXSURF+1), OPTWRK(3*MXSURF+1))

      WRITE (LUNCRT, 1001) 'Original thickness:', THICKNESS

      TARGTH = THICKNESS
      CALL READR (LUNCRT, 'Enter desired % thickness ' //
     >            '(<CR> means keep same or don''t care): ',
     >            LUNKBD, TARGTH, DEFAULT, QUIT)

C ... Is thickness going to be constrained?

      PENLTY = 0.E+0
      CALL READR (LUNCRT, 'Enter thickness penalty parameter ' //
     >            '(<CR> means no constraint): ',
     >            LUNKBD, PENLTY, DEFAULT, QUIT)

      WRITE (LUNPRT, 1002) 'Case: ', TITLE
      WRITE (LUNPRT, 1001) 'Original thickness:', THICKNESS,
     >                     'Corresponding  x/c:', XATMAX,
     >                     'Target thickness  :', TARGTH,
     >                     'Penalty parameter :', PENLTY

C ... Get the complete set of bump functions, previously prepared:

      CALL GETBUMPS (LUNCRT, LUNKBD, LUNBUMPS, MXPARM, MXBUMPS,
     >               NBUMPS, BNAMES, PARAMS, PNAMES,
     >               ACTIVE, NOPTVARS, VSCALE)

      CALL PRTBUMPS (LUNPRT, MXPARM, NBUMPS, BNAMES, PARAMS,
     >               PNAMES, ACTIVE, NOPTVARS, VSCALE)

C ... Get the target data (assumed to be in simple QPLOT format).
C     It is needed in COMMON for CFDISTRIB, but also by the calling
C     program for plotting.  Arguments were added for this purpose,
C     so as to keep the main program free of COMMON blocks.

      CALL OPENER (LUNCRT,
     >             'Enter target curvature distribution file name: ',
     >             LUNKBD, FILENAME, LUNTARG, 'OLD')

      CALL RDQPL (MXSURF, TEXT, NTARG, XTARG, YKTARG, LUNTARG)

      WRITE (LUNPRT, 1003)
     >   ' ', 'Target curvature distribution found:', ' '

      REWIND LUNTARG
      CALL LSTFIL (LUNTARG, LUNPRT, TEXT)

      NXTARG = NTARG
      CALL COPY (NXTARG,  XTARG, TARGX)
      CALL COPY (NXTARG, YKTARG, TARGCRV)

C ... Prepare for the optimizing algorithm.
C     First, extract the active variables as contiguous starting guesses:

      CALL ACTIVATE (.TRUE., MXPARM * NBUMPS, ACTIVE, PARAMS, OPTVARS,
     >               VSCALE)

      ETA = 1.0E-1
      TOL = 1.0E-3
      SSQMIN = 0.0E+0
      STEPMX = 1.0E+10
      EPSMCH = 5.0E-8
      EPSOBJ = 1.0E+2 * EPSMCH
      PRINT  = .TRUE.
      LOCAL  = .FALSE.
      CONTIN = .TRUE.
      RESCUE = .FALSE.
      UNITL  = .FALSE.
      NITER  = 100
      NTYPE  = 1
      NLDIM  = MAX (NOPTVARS * (NOPTVARS - 1) / 2, 1)

      DO I = 1, NLDIM
         L(I) = 0.0E+0
      END DO
      DO I = 1, NOPTVARS
         H(I) = -1.0E-3
      END DO

      PRINT = .FALSE.
      CALL READY (LUNCRT, 'Do you want full optimization printout? '//
     >            '(Y/N; <CR>=No): ',
     >            LUNKBD, PRINT, DEFAULT, QUIT)

C ... Compute initial value of objective function:

      CALL CFDISTRIB (NOPTVARS, OPTVARS, SSQ)

      WRITE (LUNCRT, 1001) 'Initial function value:', SSQ
      WRITE (LUNPRT, 1001)
      WRITE (LUNPRT, 1001) 'Initial function value:', SSQ, ' '

C ... Estimate good finite difference step-sizes:

      CALL CENDIF (NOPTVARS, OPTVARS, SSQ, EPSOBJ, H, G, D, NFCEN,
     >             LUNPRT, CFDISTRIB)
      NFTOTL = 1 + NFCEN

C ... Minimize the objective function:

      CALL QNMDIF (NOPTVARS, NLDIM, NFTOTL, NITER, NTYPE, LUNPRT,
     >             OPTVARS, SSQ, SSQMIN, G, H, L, D, ETA, TOL, STEPMX,
     >             EPSMCH, EPSOBJ, UNITL, LOCAL, CONV, CONTIN, RESCUE,
     >             PRINT, CFDISTRIB, DUMMY)

C ... Regenerate best result found:

      CALL CFDISTRIB (NOPTVARS, OPTVARS, SSQ)

      WRITE (LUNCRT, 1001) '  Final function value:', SSQ
      WRITE (LUNPRT, 1001)
      WRITE (LUNPRT, 1001) 'Repeat of best function evaluation:', SSQ

      CALL PRTBUMPS (LUNPRT, MXPARM, NBUMPS, BNAMES, PARAMS,
     >               PNAMES, ACTIVE, NOPTVARS, VSCALE)

C ... Update the airfoil permanently:

      IF (UPPER) THEN
         CALL COPY (NU, YCALC, YU)
      ELSE
         CALL COPY (NL, YCALC, YL)
      END IF

C ... Compute the modified thickness achieved:

      CALL CALCTHICK (NL, NU, XL, XU, YL, YU, THICKNESS,
     >                XATMAX, OPTWRK, OPTWRK(MXSURF+1),
     >                OPTWRK(2*MXSURF+1), OPTWRK(3*MXSURF+1))

      WRITE (LUNPRT, 1001)
      WRITE (LUNPRT, 1001) 'Modified % thickness:  ', THICKNESS,
     >                     'Corresponding abscissa:', XATMAX
      WRITE (LUNCRT, 1001) 'Modified % thickness:  ', THICKNESS,
     >                     ' Corresponding abscissa:', XATMAX

C ... Formats:

 1001 FORMAT (' ', A, G13.6)
 1002 FORMAT (/, '1', A, A, //)
 1003 FORMAT (' ', A, A)

      END SUBROUTINE OPTIMIZE
C+------------------------------------------------------------------------------
C
      SUBROUTINE PRREAD (LUNRD, TITLE, MAXPTS, NU, NL, X, Y, XU, XL,
     >                   YU, YL, FORMAT, IER)
C
C  ACRONYM: PRofile: READ one
C
C  PURPOSE: PRREAD reads one airfoil profile per call from a file that is
C           assumed to be open (LUNRD).  Four file formats are supported.
C
C           PRREAD returns the data as XU, YU, XL, YL, matching the standard
C           "PROFILE" format and returns a flag indicating which format the
C           values read were found in.  Standard PROFILE format is shown below.
C           The lower surface values are optional, but a zero must be read for
C           NL if no lower surface is included unless this is the last airfoil
C           in the file (meaning EOF can be used to indicate NL = 0).
C           In this case, PRREAD returns a symmetrical airfoil.
C
C               TITLE                   <Variable-length title>
C               NU   Upper surface      <Integer, first token>
C               X         Y             <Two reals>
C               X         Y                 .
C               .         .                 .     (May be X/C, Y/C;
C               .         .                 .      Xs are increasing.)
C               .         .                 .
C               NL   Lower surface      <Integer, first token> <may be 0 or EOF>
C               X         Y             <Two reals>
C               X         Y                 .
C               X         Y  ! Trailing comments are permitted
C               .         .                 .
C               .         .                 .     (Xs are increasing.)
C               .         .                 .
C               ! X         Y           <Point suppressed; NL must be adjusted>
C               .         .                 .
C
C    NOTE:  If FORMAT = 1 and both surfaces are present, PROFILE expects
C           them to have the same leading edge point.  The trailing edge
C           points may differ.  However, checking for the same leading
C           edge point is left to the calling program in case PRREAD is
C           being used to read other upper/lower surface-type distributions.
C
C           The next two formats are wraparound clockwise and wraparound
C           counterclockwise, where the coordinates begin at the trailing
C           edge, wrap around the leading edge, and end at the trailing edge.
C           The clockwise case begins with the lower surface, and the counter-
C           clockwise case begins with the upper surface.  The format shown
C           below is essentially the same for both cases. NPTS is the total
C           number of points on the airfoil.
C
C               TITLE                   <CHARACTER*80>
C               NPTS                    <Integer, first token>
C               X         Y             <Two reals>
C               X         Y                 .
C               .         .                 .     (May be X/C, Y/C;
C               .         .                 .      Xs are decreasing
C               .         .                 .      until the leading
C               .         .                 .      edge, then increasing)
C
C    NOTE:  Wraparound formats do NOT have duplicate leading edge points.
C
C           The fourth format is called 3-column format.  The airfoil is
C           represented in three columns, with the same abscissas for both
C           surfaces in the 1st column and ordinates for the upper and lower
C           surfaces in the 2nd and 3rd columns respectively.  Abscissas are
C           increasing as with standard format.  Here NPTS is the number of
C           points on either surface.  Further columns may be present in
C           this case only (as for Cp distributions - only the first three
C           columns of the data proper will be read here).
C
C               TITLE                           <CHARACTER*80>
C               NPTS                            <Integer, first token>
C               X         YU        YL     [Cp] <Reals, first 3 tokens;
C               X         YU        YL     [Cp]  further columns are ignored.>
C               .         .         .      [..]     .
C               .         .         .               .   (May be X/C, Y/C;
C               .         .         .               .    Xs are increasing.)
C               .         .         .               .
C               .         .         .               .
C
C  ARGUMENTS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   LUNRD    -      I     I     Logical unit number for file being read.
C   MAXPTS   -      I     I     Maximum number of coordinates on any one
C                               surface expected by calling program.
C   TITLE    *      C     O     Variable-length title for this profile.
C   NU,NL    -      I     O     Number of data points found for upper and
C                               lower surfaces.  NL is set to NU if NL = 0
C                               (or if EOF is encountered) on the read.
C   X,Y   MAXPTS*2  R    S/O    Buffers for reading coordinates; returned
C                               with the airfoil in clockwise wraparound form
C                               (NU + NL - 1 points) in case that is more
C                               convenient than the form in XU, YU, XL, YL.  
C   XU     MAXPTS   R     O     Upper surface abscissas found.
C   XL     MAXPTS   R     O     Lower surface abscissas.
C   YU,YL  NU,NL    R     O     Corresponding ordinates.
C   FORMAT   -      I     O     FORMAT=1 means data found in standard format
C                                     =2 means clockwise wraparound format
C                                     =3 means counterclockwise wraparound
C                                     =4 means 3-column format
C   IER      -      I     O     IER=0 means one profile found normally;
C                                  =1 means EOF encountered on the first
C                                     read -- normal unless this was the
C                                     first call to PRREAD.
C                                  =2 means an unexpected EOF or other
C                                     read error encountered -- fatal.
C                                  =3 means NU or NL  was found out of
C                                     range (> MAXPTS, or too small).
C                                  =4 means a coordinate was missing.
C
C  PROCEDURES:
C   COPY           Copies a vector
C   GETLINE        Reads a line as text; handles suppressed points & comments
C   RVERSE         Reverses a vector
C   TOKENS         Tokenizes a string
C
C  HISTORY:
C  03/19/82    DAS    Original implementation.
C  06/29/84    LJC    Added reading of wraparound formats.
C  09/14/84    LJC    Changed "legend" entry to be read from the dataset
C                     and added a prompt for the plot title at a higher
C                     level. (Formerly the title was read from the
C                     dataset and the legend (TITLE here) was hard-coded.)
C  10/18/84    DAS    Took out check for duplicate leading edge point -
C                     not wanted if other distributions are being read.
C                     (Left to the calling program where appropriate.)
C  02/27/85    LJC    Added 3-column format.
C  09/18/87    DAS    Handled as few as 2 points on a surface properly;
C                     allowed for EOF on reading NL -- treat as NL = 0.
C  02/06/90  DAS/RAK  GETLINE introduced to permit trailing comments or
C                     suppression of coordinates.
C  10/23/91    DAS    Should use LINE (1:LAST) everywhere, not LINE.
C  09/18/92    DAS    LENTOK = 20 failed on VAX double-precision numbers.
C  11/06/92    DAS    3-column format can now ignore additional columns.
C                     Too few coordinates on a line is now trapped.
C  11/16/93    DAS    Returned additional clockwise wraparound form of
C                     the airfoil in the available X, Y arrays in case
C                     that is more convenient for some applications.
C
C  AUTHOR:     David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   FORMAT, IER, LUNRD, MAXPTS, NU, NL

      CHARACTER
     >   TITLE * 80

      REAL
     >   X (MAXPTS * 2), XU (MAXPTS), XL (MAXPTS), Y (MAXPTS * 2),
     >   YU (MAXPTS), YL (MAXPTS)

C  *  Local constants:

      INTEGER, PARAMETER ::
     >   LENTOK  = 24

      CHARACTER, PARAMETER ::
     >   BLANK * 1   = ' ',
     >   COMMENT * 1 = '!',
     >   IFMT * 8    = '(BN,I24)',
     >   RFMT * 10   = '(BN,F24.0)'

C  *  Local variables:

      INTEGER
     >   I, IOS, LAST, MIDSRF, LE, NPTS, NUMBER

      CHARACTER
     >   LINE * (3 * LENTOK), LIST (3) * (LENTOK)

C  *  Procedures:

      EXTERNAL
     >   COPY, GETLINE, RVERSE, TOKENS

C  *  Execution:

      IER = 0

C  *  Look for descriptive text for a profile.  There may be no more profiles.

      CALL GETLINE (LUNRD, BLANK, TITLE, LAST, IOS)
      IF (IOS  < 0) GO TO 800                    ! Normal EOF
      IF (IOS /= 0) GO TO 900

C  *  Look for a count of the points to follow:

  100 CALL GETLINE (LUNRD, COMMENT, LINE, LAST, IOS)
      IF (IOS  /= 0) GO TO 900
      IF (LAST == 0) GO TO 100    ! Insignificant line

      NUMBER = 1
      CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)

      READ (LIST (1), IFMT, ERR=900) NPTS

      IF (NPTS <= MAXPTS*2 .AND. NPTS > 1) THEN

C  *     Check number of columns in first line of data:

  110    CALL GETLINE (LUNRD, COMMENT, LINE, LAST, IOS)
         IF (IOS  /= 0) GO TO 900
         IF (LAST == 0) GO TO 110
         NUMBER = 3
         CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)

         IF (NUMBER == 2) THEN

C  *        Only two columns of data found. Process one of the first 3 formats.
C           No need to backspace any more.  (But a DO loop is inconvenient.)

            I = 1
  200       CONTINUE
               READ (LIST (1), RFMT, ERR=900) X (I)
               READ (LIST (2), RFMT, ERR=900) Y (I)
               I = I + 1

               IF (I <= NPTS) THEN

  210             CALL GETLINE (LUNRD, COMMENT, LINE, LAST, IOS)

                  IF (IOS  /= 0) GO TO 900
                  IF (LAST == 0) GO TO 210

                  CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)

                  IF (NUMBER < 2) GO TO 910
                  GO TO 200
               END IF

C  *        Distinguish standard and wrap-around formats automatically:

            MIDSRF = MAX (1, NPTS / 4)
            IF (X (MIDSRF) < X (MIDSRF + 1)) THEN

C  *           Abscissas are increasing; "PROFILE" format assumed.
C              NPTS must be less than MAXPTS (not MAXPTS*2) now:

               IF (NPTS <= MAXPTS) THEN
                  FORMAT = 1
                  NU = NPTS

C  *              Copy X and Y into upper surface arrays:

                  CALL COPY (NU, X, XU)
                  CALL COPY (NU, Y, YU)                        

C  *              Look for a lower surface point count (may be 0 or EOF):

  300             CALL GETLINE (LUNRD, COMMENT, LINE, LAST, IOS)

                  IF (IOS < 0) THEN          ! Normal EOF
                     NL = 0
                  ELSE IF (IOS  /= 0) THEN
                     GO TO 900
                  ELSE IF (LAST == 0) THEN
                     GO TO 300
                  ELSE
                     NUMBER = 1
                     CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)
                     READ (LIST (1), IFMT, ERR=900) NL
                  END IF

                  IF (NL > 0) THEN

                     IF (NL <= MAXPTS) THEN
                        NUMBER = 2
                        DO I = 1, NL
  310                      CALL GETLINE (LUNRD, COMMENT, LINE, LAST,IOS)
                           IF (IOS  /= 0) GO TO 900
                           IF (LAST == 0) GO TO 310
                           CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)
                           IF (NUMBER < 2) GO TO 910
                           READ (LIST (1), RFMT, ERR=900) XL (I)
                           READ (LIST (2), RFMT, ERR=900) YL (I)
                        END DO
                     ELSE
C  *                    Number of lower surface points found out of range:
                        IER = 3
                     END IF

                  ELSE

C  *                 Generate symmetrical lower surface points:
  
                     DO I = 1, NU
                        XL (I) = XU (I)
                        YL (I) =-YU (I)
                     END DO
                     NL = NU

                  END IF

               ELSE

C  *              Number of upper surface points found out of range:

                  IER = 3

               END IF

            ELSE

C  *           Process wraparound airfoil:           

               DO I = 1, NPTS-1  ! Search for leading edge:

                  IF (X (I) < X (I+1)) THEN
                     LE = I
                     GO TO 510
                  END IF
               END DO

C  *           No leading edge found:

               GO TO 900
 
  510          CONTINUE

C  *           Arrange coordinates in standard PROFILE format:

               IF (Y (MIDSRF) < Y (NPTS - MIDSRF)) THEN

C  *              The lower surface is first; this is the clockwise case:

                  FORMAT = 2
                  NU = NPTS - LE + 1
                  NL = LE
                  CALL RVERSE (NL, X, XL)
                  CALL RVERSE (NL, Y, YL)
                  CALL COPY (NU, X (NL), XU)
                  CALL COPY (NU, Y (NL), YU)

               ELSE  ! Counterclockwise case:

                  FORMAT = 3
                  NL = NPTS - LE + 1
                  NU = LE
                  CALL COPY (NL, X (NU), XL)
                  CALL COPY (NL, Y (NU), YL)
                  CALL RVERSE (NU, X, XU)
                  CALL RVERSE (NU, Y, YU)

               END IF
  
            END IF

         ELSE IF (NUMBER == 3) THEN  ! Process 3-column format:

            FORMAT = 4
            I = 1
  600       CONTINUE
               READ (LIST (1), RFMT, ERR=900) XU (I)
               READ (LIST (2), RFMT, ERR=900) YU (I)
               READ (LIST (3), RFMT, ERR=900) YL (I)
               I = I + 1
               IF (I <= NPTS) THEN
  610             CALL GETLINE (LUNRD, COMMENT, LINE, LAST, IOS)
                  IF (IOS  /= 0) GO TO 900
                  IF (LAST == 0) GO TO 610
                  CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)
                  IF (NUMBER < 3) GO TO 910
                  GO TO 600
               END IF

            CALL COPY (NPTS, XU, XL)
            NU = NPTS
            NL = NPTS
    
         ELSE  ! Inappropriate number of columns found:

            GO TO 900

         END IF
      
      ELSE  ! Number of points found out of range:

         IER = 3

      END IF

C  *  Return a second copy of the coordinates in clockwise wraparound form
C     in case that is more convenient for some application other than PROFILE:

      IF (IER == 0) THEN
         CALL RVERSE (NL, XL, X)
         CALL COPY   (NU, XU, X (NL))
         CALL RVERSE (NL, YL, Y)
         CALL COPY   (NU, YU, Y (NL))
      END IF

      GO TO 999


  800 CONTINUE

C  *  Normal EOF encountered on first read -- no more profiles:

      IER = 1
      GO TO 999

  900 CONTINUE

C  *  Abnormal EOF or other read error -- fatal:

      IER = 2
      GO TO 999

  910 CONTINUE

C  *  Too few values on a line:

      IER = 4
!     GO TO 999

  999 RETURN

      END SUBROUTINE PRREAD
C+------------------------------------------------------------------------
C
      SUBROUTINE PRTAB (LUNTAB, TITLE, SUBTITLE, WRAPCRV,
     >                  N, X, Y, YP, YPP, YK)
C
C  PURPOSE: PRTAB tabulates geometrical properties of an airfoil profile,
C           one surface per call.
C
C  ARGUMENTS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   LUNTAB   -      I     I     Logical unit for file being written to
C   TITLE    *      C     I     Main title for tabulation (probably the
C                               title from the file of airfoil coordinates)
C   SUBTITLE *      C     I     Additional subtitle that can help make the
C                               tabulation more self-descriptive
C   WRAPCRV  -      L     I     .TRUE. means the derivatives and curvatures
C                               were calculated parametrically using splines;
C                               .FALSE. means they were finite difference
C                               approximations with surfaces treated separately.
C   N        -      I     I     Number of data points on the surface
C   X        N      R     I     Abscissas for the surface
C   Y        N      R     I     Corresponding ordinates
C   YP       N      R     I     1st derivatives  at given coordinates
C   YPP      N      R     I     2nd derivatives  at given coordinates
C   YK       N      R     I     Curvature values at given coordinates
C
C  HISTORY:
C  01/20/83   LJC   Initial coding
C  04/05/83   LJC   Added calculation of curvatures.
C  04/08/83   LJC   Handle each surface separately.
C  10/17/83   DAS   Removed curvature calculn. - now done by FD12K.
C  10/29/83   DAS   Introduced subtitle.
C  10/28/91   DAS   Introduced WRAPCRV.
C
C  AUTHOR:   Leslie Collins, Informatics, Palo Alto, CA.
C
C----------------------------------------------------------------------

C     Arguments:

      INTEGER       LUNTAB, N
      REAL          X (N), Y (N), YP (N), YPP (N), YK (N)
      LOGICAL       WRAPCRV
      CHARACTER*(*) TITLE, SUBTITLE

C     Local variables:

      INTEGER       I

C     Execution:

      WRITE (LUNTAB, 1001) '1', TITLE, SUBTITLE
      WRITE (LUNTAB, 1002) N
      IF (WRAPCRV) THEN
         WRITE (LUNTAB, 1003) 'parametric spline'
      ELSE
         WRITE (LUNTAB, 1003) 'nonparametric finite differencing'
      END IF
      WRITE (LUNTAB, 1004)
      WRITE (LUNTAB, 1005) (X(I), Y(I), YP(I), YPP(I), YK(I), I = 1, N)

C     Formats:

 1001 FORMAT (/, A1, A, //, 1X, A)
 1002 FORMAT (/, ' Number of points: ', I3)
 1003 FORMAT (/, ' Derivatives and curvature calculated by ',
     >        A, '.')
 1004 FORMAT (//, T9, 'X',  T25, 'Y', T39, 'Y''', T54, 'Y"',
     >        T66, 'CURVATURE' )
 1005 FORMAT (2E15.7, 3E15.6)

      END SUBROUTINE PRTAB
C+----------------------------------------------------------------------
C
      SUBROUTINE PRTBUMPS (LUNPRT, MXPARM, NBUMPS, BNAMES,
     >                     PARAMS, PNAMES, ACTIVE, NACTIVE, SCALES)
C
C PURPOSE: PRTBUMPS prints a description of a given bump set.  It is
C          intended for echoing an initial set prior to optimization
C          of the active variables, and for tabulating the optimized
C          results in unscaled form.
C
C ARGUMENTS:
C    ARG     DIM   TYPE I/O/S DESCRIPTION
C   LUNPRT    -      I    I   Logical unit for tabulations.
C   MXPARM    -      I    I   Max. # parameters  for any one bump.
C   NBUMPS    -      I    I   Number of bump functions involved.
C   BNAMES  NBUMPS  C*(*) I   Names of the bump functions.
C   PARAMS  MXPARM,  R    I   PARAMS(1:?,J) are the parameters de-
C           NBUMPS            fining the Jth bump.
C   PNAMES  MXPARM, C*(*) I   Names of function parameters
C           NBUMPS
C   ACTIVE  MXPARM,  L    I   ACTIVE(I,J)=.TRUE. if the Ith param-
C           NBUMPS            eter of the Jth bump is to be treat-
C                             ed as active (variable).  Otherwise,
C                             the parameter is either to remain at
C                             the value returned here (fixed),  or
C                             bump J has fewer than I parameters.
C   NACTIVE  -       I    I   Number of active parameters.
C   SCALES  MXPARM,  R    I   Scale factors, included here for in-
C           NBUMPS            formation only.
C
C HISTORY:
C   01/25/84   DAS    Initial design and code.
C   02/23/84   DAS    Now prints bump names instead of code numbers.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER   LUNPRT, MXPARM, NACTIVE, NBUMPS
      REAL      PARAMS(MXPARM,NBUMPS), SCALES(MXPARM,NBUMPS)
      LOGICAL   ACTIVE(MXPARM,NBUMPS)
      CHARACTER BNAMES(NBUMPS)*(*), PNAMES(MXPARM,NBUMPS)*(*)

C     Local variables:

      INTEGER   I, J
      REAL      UNDEF

      DATA      UNDEF / 999.E+0 / ! This preset value should match GETBUMPS.

C     Execution:

      WRITE (LUNPRT, '(//, (A))') ' Shape Function Description',
     >                            ' --------------------------'

      DO J = 1, NBUMPS
         WRITE (LUNPRT, '(//, A, I3, 2A, //, A, /)')
     >      ' Function', J, ':  Name: ', BNAMES(J),
     >      ' Parameter    Value       Scale      Active?'

         DO I = 1, MXPARM
            IF (PARAMS(I,J) == UNDEF) CYCLE
               WRITE (LUNPRT, '( 1X, A, F12.8, F12.5, L10 )')
     >            PNAMES(I,J), PARAMS(I,J), SCALES(I,J), ACTIVE(I,J)
         END DO
      END DO

      WRITE (LUNPRT, '(//, A, I4)')
     >   ' Total number of active parameters:', NACTIVE

      END SUBROUTINE PRTBUMPS
C+---------------------------------------------------------------------
C
      SUBROUTINE PRWRIT (LUN, MAXPTS, TITLE, NU, NL, XU, XL, YU, YL,
     >                   FORMAT, PRECISION)
C
C  ACRONYM: Program PRofile: WRITe one airfoil profile
C                   --       ----
C  PURPOSE: PRWRIT writes one airfoil profile to the indicated file, in
C           one of four formats - "standard",  wraparound (either way),
C           or "three-column" - described in PRREAD.
C
C           PRWRIT  may also be used to save camber/thickness distribu-
C           tions,  or to save second derivative information for reuse.
C           Both of these extended uses employ  "standard"  format (not
C           wrap-around), with suitable labeling.
C
C           This version allows for three types of precision,  prompted
C           by the need to retain more digits for manipulating the all-
C           important leading edge region effectively, and to deal with
C           large magnitudes such as when the coordinates are in milli-
C           meters, or when the "ordinates" are really 2nd derivatives.
C
C  METHOD:  Values are passed to PRWRIT as separate surfaces,  and  are
C           returned untouched.  Note that the wrap-around formats omit
C           one of the two leading edge points, which are assumed to be
C           the same point.
C           
C  ARGUMENTS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   LUN      -      I     I     Logical unit for file being written to.
C   MAXPTS   -      I     I     Max. no. of points on any one surface.
C   TITLE    -     C*(*)  I     Variable-length title for profile.
C   NU,NL    -      I     I     Number of data points for upper and lower
C                               surfaces (NU > 0; NL=0 is OK if FORMAT=1;
C                               NU=NL=no. of pts. in 3-column format and
C                               the camber/thickness distributions - FORMATs
C                               4 & 5).
C   XU     MAXPTS   R     I     Upper surface abscissas.
C   XL     MAXPTS   R     I     Lower surface abscissas (if any).
C   YU,YL  MAXPTS   R     I     Corresponding ordinates, y" values, or
C                               camber/thickness values, accdg. to FORMAT.
C   FORMAT   -      I     I     Requested format for output where:
C                               = 1 means standard PROFILE format (coords.)
C                               = 2 means clockwise wrap-around format
C                               = 3 means counterclockwise wrap-around
C                               = 4 means 3-column format
C                               = 5 means 2nd derivatives (standard format)
C                               = 6 means camber/thickness (standard format)
C   PRECISION -     I     I     Controls number of digits in saved values:
C                               = 1 means "full" (see PURPOSE above);
C                               = 2 means "engineering" or "flow code" (F10.6);
C                               = 3 means "low" - appropriate for y" values.
C  FILES USED:
C   UNIT     I/O/S    DESCRIPTION
C   LUN        O      File (assumed open) to contain one or more datasets
C
C  HISTORY:
C  12/23/82   LJC   Coding adapted from PRREAD.
C  10/30/83   DAS   Introduced alternative E-format so that PRWRIT
C                   can be used to save 2nd derivatives similarly.
C  07/09/84   LJC   Added writing of wraparound formats.
C  02/11/84   DAS   Handled thickness/camber - hard to avoid dupli-
C                   cation of code now.  But PRWRIT is still handy
C                   for getting this stuff out of the main program.
C  02/28/85   LJC   Added 3-column format. (FORMAT=4)
C  04/09/87   DAS   Values except X are written with F10.7 format
C                   instead of F10.6.   (May help difficulties in
C                   critical leading edge region; will not affect
C                   formatted 2-column reads, but may be a little
C                   misleading.)
C  04/23/87   DAS   Erred further in direction of more precision:
C                   use E format except on basically normalized
C                   data; go to F12.8 for normalized airfoils.
C  04/24/87   DAS   Retained F10.6 option after all for old flow
C                   codes (FORMAT=7; standard PROFILE format only.
C  04/27/87   DAS   Reluctantly introduced PRECISION argument after
C                   the above failed to handle FORMAT="flowcode" and
C                   COUNTERCLOCKWISE both.
C  11/03/87   DAS   Switched to G formats for "full" and "low" precision.
C  12/12/91   DAS   Had to "comment out" the text following NU, NL
C                   because of how RDXYZ works now.  (PROFILE uses
C                   PRWRIT to save 2nd derivatives, which are read back
C                   via RDXYZ.)
C
C  AUTHOR:   Leslie Collins, Informatics, Palo Alto, CA.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      CHARACTER
     >   TITLE * (*)
      INTEGER
     >   FORMAT, LUN, MAXPTS, NL, NU, PRECISION
      REAL
     >   XU (MAXPTS), XL (MAXPTS), YU (MAXPTS), YL (MAXPTS)

C  *  Local variables:

      INTEGER
     >   I
      CHARACTER
     >   FMT * 13

C  *  Execution:

      IF (PRECISION == 1) THEN        ! "Full" precision.
         FMT = '(2(1X,G15.7))'          ! Gw.d needs w >= d + 8 (basically).

      ELSE IF (PRECISION == 2) THEN   ! "Engineering" precision.
         FMT = '(2F10.6)     '          ! Traditional for many flow codes.

      ELSE                              ! "Low" precision.
         FMT = '(2(1X,G13.5))'          ! Appropriate for y" values (FORMAT=5).
      END IF

      IF (FORMAT == 4) FMT (2:2) = '3'

      WRITE (LUN, '(A)') TITLE

      IF (FORMAT == 1) THEN           ! Standard PROFILE format.

         WRITE (LUN, 1001) NU, 'Upper Coordinates'
         WRITE (LUN, FMT) (XU (I), YU (I), I = 1, NU)

         IF (NL > 0) THEN
            WRITE (LUN, 1001) NL, 'Lower Coordinates'
            WRITE (LUN, FMT) (XL (I), YL (I), I = 1, NL)
         END IF            

      ELSE IF (FORMAT == 2) THEN      ! Wrap-around (lower surface first).

         WRITE (LUN, 1001) NU + NL - 1, 'Coordinates Clockwise'
         WRITE (LUN, FMT) (XL (I), YL (I), I = NL, 1, -1)
         WRITE (LUN, FMT) (XU (I), YU (I), I = 2, NU)

      ELSE IF (FORMAT == 3) THEN      ! Wrap-around (upper surface first).

         WRITE (LUN, 1001) NU + NL - 1, 'Coordinates Counter-clockwise'
         WRITE (LUN, FMT) (XU (I), YU (I), I = NU, 1, -1)
         WRITE (LUN, FMT) (XL (I), YL (I), I = 2, NL)

      ELSE IF (FORMAT == 4) THEN      ! Three-column format.

         WRITE (LUN, 1001) NU, 'Coordinates per surface'
         WRITE (LUN, FMT) (XU (I), YU (I), YL (I), I = 1, NU)

      ELSE IF (FORMAT == 5) THEN      ! 2nd derivatives (both surfaces).
                                        ! Use standard format.
         WRITE (LUN, 1001) NU, 'Upper 2nd Derivatives'
         WRITE (LUN, FMT) (XU (I), YU (I), I = 1, NU)

         WRITE (LUN, 1001) NL, 'Lower 2nd Derivatives'
         WRITE (LUN, FMT) (XL (I), YL (I), I = 1, NL)

      ELSE  ! FORMAT = 6: Camber and Thickness distributions in Standard format.

         WRITE (LUN, 1001) NU, 'Camber'
         WRITE (LUN, FMT) (XU (I), YU (I), I = 1, NU)

         WRITE (LUN, 1001) NL, 'Thickness'
         WRITE (LUN, FMT) (XL (I), YL (I), I = 1, NL)

      END IF        

C  *  Formats:

 1001 FORMAT (I4, ' ! ', A)

      END SUBROUTINE PRWRIT
C+----------------------------------------------------------------------
C
      SUBROUTINE RDKEYS (LUNCRT, LUNIN, MXCRVS, UNDEF, FORMAT,
     >                   PRECISION, PLTLINE, CPSLINE, CRVLINE, MAXCRV,
     >                   MINCRV, WRAPCRV, PLTORIG, PLTREV, THREED,
     >                   XAXIS, XMIN, XMAX, YMIN, YMAX, DATFIL, PLTFIL,
     >                   CRVFIL, TABFIL, YPPFIL, CPSFIL, SPREADFIL)
C
C  PURPOSE:
C
C     RDKEYS reads PROFILE's keyword input control file.   Most of these
C     inputs are common to all modes of operation.  Other mode-dependent
C     inputs are entered interactively in the appropriate subprogram.
C
C  METHOD:
C
C     Read a line of input and break up the line into pairs of keywords and
C     values.  For each pair, look up the keyword in a dictionary and then
C     look up the value in a separate dictionary associated with that keyword.
C     Assign a value to the corresponding internal control variable via the
C     argument list.
C
C     Some keywords may have more than one corresponding value.  In this case,
C     all of the tokens following the first token are treated as values and
C     not pairs.  A full description of the keywords follows.
C
C  KEYWORD GUIDELINES AND DEFINITIONS:
C
C     Keyword/value pairs may appear with more than one pair on a line.
C     However, the multivalued keywords PLTLINE, CPSLINE, CRVLINE, and
C     NOFILE must not appear with other keywords on the same line.
C
C     The default value in each case appears in square brackets.
C
C  KEYWORD  VALUES and synonyms     DESCRIPTION
C  -------  -------------------     -----------
C
C  FORMAT   [SAME]                  One of four formats for output profile.
C           PROFILE or STANDARD     May be in standard PROFILE format (ab-
C           CLOCKWISE or WRAPAROUND scissas increasing), clockwise wraparound
C           COUNTERCLOCKWISE        format, counterclockwise wraparound for-
C           THREE-COLUMN or         mat, or 3-column format.  SAME means the
C           THREE_COLUMN or         same format as the input profile.  NOTE:
C           THREECOLUMN  or         To allow easily for several synonyms for
C           TABLE                   for the THREE-COLUMN value, only the first
C                                   5 characters of the value are checked.
C
C  PRECISION   [FULL]               Controls number of digits in output
C           ENGINEERING             airfoil coordinates.  FULL gives F11.8
C                                   if possible, or E15.8 if any X >=10.
C                                   ENGINEERING gives the traditional F10.6
C                                   common to many flow solvers.
C
C  PLTLINE  [DEFAULT]               Controls line types of curves on profile
C           LINE                    plots.  One value may be included for
C           DASH                    each curve on the plot.  The default is
C           DOT                     symbols connected by a solid line, with
C           CHAINDASH               a different symbol type for each succes-
C           CHAINDOT                sive curve.  The first curve typically
C           THICK                   represents the original profile; the
C           SYMBOLS                 second curve represents the revised one.
C                                   Overriding the default might be desirable
C                                   when plotting multi-element airfoils or
C                                   when lines without symbols are required.
C                                   At most 20 curves are provided for.  Note:
C                                   All the line types in QPLOT are available. 
C                                   SYMBOLS refers to symbols with no line
C                                   connecting them.
C
C  CPSLINE  [see PLTLINE above]     Controls line types on Cps plots in the
C                                   same manner as PLTLINE above.  One value
C                                   per curve may be included, chosen from
C                                   the same list of values as those shown
C                                   for PLTLINE.
C
C  CRVLINE  [see PLTLINE above]     Controls line types on curvature plots
C                                   in the same way as PLTLINE and CPSLINE.
C
C  CURVATURE or [NONPARAMETRIC] or  CURVATURE and DERIVATIVES are synonymous
C  DERIVATIVES  [FINITE_DIFFERENCE] controls for the type of calculations
C               SPLINE     or       used for derivatives and hence curvature.
C               PARAMETRIC or       The default is separate-surface treatment
C               WRAPAROUND          using finite differences, as needed for
C                                   consistency with PROFILE's REFINE and
C                                   OPTIMIZE options.  The two surfaces appear
C                                   as separate frames in the curvature plot.
C                                   Otherwise, the full wrap-around curvature
C                                   distribution is calculated using a para-
C                                   metric spline and plotted on a single frame.
C
C                                   The default normally suffices except if
C                                   the region of interest is very near a
C                                   rounded leading edge.  Note that not all
C                                   of the possibilites are provided for, such
C                                   as parametric finite differences.
C
C  MINCURVATURE   [-5.]             Cutoff values for plotted curvatures.
C  MAXCURVATURE   [+5.]             Practice shows that +/-5. give useful
C                                   plot scaling by ignoring the high curv-
C                                   ature values near the leading edge.  On
C                                   the other hand, it may well be desired
C                                   to focus on the leading edge region.  Set
C                                   both to 999. to obtain the full range.
C                                   See the CURVATURE/DERIVATIVES control.
C
C  NOFILE   [NONE]                  Used to suppress any combination of the
C           DAT                     seven output files generated by PROFILE.
C           PLT                     The values correspond to the extensions
C           TAB                     of the file names.  See elsewhere for a
C           CRV                     complete description of file contents.
C           YPP                     NONE serves only to assist leaving the
C           CPS                     NOFILE control word in the input file
C           SPREAD                  even if all outputs are desired.
C
C  PLOT     [BOTH]                  Controls plotting of original OR revised
C           ORIGINAL                profile.  The default is to plot both
C           REVISED                 original and revised (if one exists).
C 
C  THREED   [FALSE] or [NO]         For plotting of multiple stations from
C           TRUE or YES             a 3-D wing. The default is the 2-D case.
C 
C  XAXIS    [6.4]                   Length of x-axis in inches.  The default
C                                   is normally appropriate for an 8.5 x 11
C                                   page in portrait mode.
C
C           The following four keywords apply to windowing.  Any or none
C           of them may be used.
C 
C  XMIN     [minima and             Minimum abscissa for desired window
C  XMAX     maxima of               Maximum abscissa for desired window
C  YMIN     the input               Minimum ordinate for desired window
C  YMAX     coordinates]            Maximum ordinate for desired window
C
C
C  SAMPLE CONTROL FILE:
C
C     A sample input file follows.  Note that keywords and values
C     may be separated with blanks, commas, colons, equals signs,
C     or tabs. Remember, keywords with more than one value should
C     appear on separate lines.  Any keyword or text value may be
C     truncated to unambiguous leading characters.   Blank  lines
C     and trailing ! comments are ignored.
C
C
C           FORMAT = STANDARD   PRECISION = FULL
C           PLOT BOTH  THREED:NO
C           PLTLINE = SOLID, SOLID
C           CPSLINE = DOT, SYMBOLS
C           CRVLINE = SOLID, DASH, CHAINDOT
C           XAXIS = 20.
C           XMIN = 0.  XMAX 0.1
C           MAXCURVATURE = 999.   ! Both 999. means plot the full
C           MINCURVATURE = 999.   ! curvature range
C           DERIVATIVES = PARAMETRIC
C           NOFILE: YPP SPREAD
C 
C
C  ERROR HANDLING:
C
C     The dictionary lookup routine returns a pointer indicating which entry
C     matched the key.  When the pointer is negative, meaning no match was
C     found, RDKEYS prints an error message along with the string in question
C     and stops.
C
C  EXTERNAL REFERENCES:
C
C  GETLINE     Reads one line, and handles trailing comments.
C  LOOKUP      Performs dictionary lookups.
C  PAIRS       Breaks up alternating fields into pairs of keywords and values.
C  SCANNR      Looks for non-blank fields in a string.
C  TOKENS      Separates groups of contiguous characters into an array.
C
C  HISTORY:
C
C  08/01/84    LJC    Initial design and coding.
C  12/05/84    DAS    Added CPSFIL handling.
C  02/19/85    DAS    Took out SHARP/BLUNT control - now in REDISTRIB.
C  03/22/85    LJC    Replaced LINE keyword with PLTLINE, CPSLINE and CRVLINE,
C                     to control line types on all plots.
C  10/09/85    DAS    Introduced UNDEF argument because of revised QPLOT.
C  04/24/87    DAS    Added FORMAT=FLOWCODE after PRWRIT was changed to
C                     provide more precision than F10.6.
C  04/27/87    DAS    The above didn't handle FLOWCODE and COUNTERCLOCKWISE
C                     needed by one user - introduced PRECISION instead.
C                     Allowed FORMAT=WRAPAROUND [=CLOCKWISE - arbitrary,
C                     but it had to be one way or the other].
C  05/08/87    DAS    Bugs: VALUES () * 10 is too short for ENGINEERING;
C                     Dimension of LIST () has to be MAX (MXCRVS,NFILDC-1),
C                     except NFILDC is all that will fit on a line anyway
C                     for CPSLINE, CRVLINE, and PLTLINE.
C  04/29/88    DAS    Added spreadsheet-compatible output file control.
C  06/18/91    DAS    Made MXLIST = MXCRVS + 1 = 21 for 3D case of many
C                     sections.  (Plot line types CAN be abbreviated, so
C                     the comment of 05/08/87 is not really right.)
C  10/23/91    DAS    Introduced GETLINE to allow commented control files.
C  10/28/91    DAS    Introduced a means of specifying full wrap-around
C                     curvature distribution (DERIVATIVES/CURVATURE keywords).
C
C  AUTHOR:     Leslie Collins, Informatics General, Palo Alto, CA
C
C------------------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   FORMAT, LUNCRT, LUNIN, MXCRVS, PRECISION
      REAL
     >   MAXCRV, MINCRV, UNDEF, XAXIS, XMAX, XMIN, YMAX, YMIN
      CHARACTER
     >   CPSLINE (MXCRVS) * 9, CRVLINE (MXCRVS) * 9,
     >   PLTLINE (MXCRVS) * 9
      LOGICAL
     >   CPSFIL, CRVFIL, DATFIL, PLTFIL, PLTORIG, PLTREV, SPREADFIL,
     >   TABFIL, THREED, WRAPCRV, YPPFIL

C  *  Local constants:

      INTEGER, PARAMETER ::
     >   MXCHARS = 81, MXLIST = 21, NDERDC = 5, NFILDC = 8, NFMTDC = 8,
     >   NKEYDC  = 17, NLINDC = 9,  NPLTDC = 3, NPREDC = 2, NTHRDC = 4

      CHARACTER, PARAMETER ::
     >   BLANK * 1 = ' '

C  *  Local variables:

      INTEGER
     >   FIRST, I, IOS, J, LAST, MARK, NUMBER, NPAIRS, POINTR
      LOGICAL
     >   CR
      CHARACTER
     >   VAL*3,
     >   BUFFER * (MXCHARS), DERDIC (NDERDC) * 13, KEYDIC (NKEYDC) * 12,
     >   FILDIC (NFILDC) * 6, FMTDIC (NFMTDC) * 11, KEYS (NKEYDC) * 12,
     >   LINDIC (NLINDC) * 9, LIST (MXLIST) * 10, PLTDIC (NPLTDC) * 8,
     >   PREDIC (NPREDC) * 11, THRDIC (NTHRDC) * 8, VALUES (NKEYDC) * 12

C  *  Procedures:

      EXTERNAL
     >   GETLINE, LOOKUP, PAIRS, SCANNR, TOKENS

C  *  Dictionaries:

C     Guidelines:
C     >  All calls to LOOKUP assume that the dictionaries are alphabetical.
C     >  Use of synonyms either leads to long dictionary entries or forces
C        (safe) shortening of the target string to (say) 5 characters.  Any
C        characters beyond that in the other dictionary entries are strictly
C        for programmer readability.

C     DERIVATIVES/CURVATURE dictionary:

      DATA DERDIC
     >   /'FINITE=NONPAR', 'NONPARAMETRIC', 'PARAM=SPLINE',
     >    'SPLINE', 'WRAPA=SPLINE'/

C     Main control keyword dictionary:

      DATA KEYDIC
     >   /'CPSLINE', 'CRVLINE',      'CURVATURE',    'DERIV=CURVA',
     >    'FORMAT',  'MAXCURVATURE', 'MINCURVATURE', 'NOFILE',
     >    'PLOT',    'PLTLINE',      'PRECISION',    'THREED',
     >    'XAXIS',   'XMAX', 'XMIN', 'YMAX', 'YMIN'/

C     NOFILE dictionary:

      DATA FILDIC
     >   /'CPS', 'CRV', 'DAT', 'NONE', 'PLT', 'SPREAD', 'TAB', 'YPP'/

C     FORMAT dictionary:

      DATA FMTDIC
     >   /'CLOCK', 'COUNT', 'PROFI', 'SAME', 'STAND:PROFI',
     >    'TABLE:THREE', 'THREE', 'WRAPA:CLOCK'/

C     Line-type dictionary:

      DATA LINDIC
     >   /'CHAINDASH', 'CHAINDOT', 'DASH', 'DEFAULT', 'DOT',
     >    'LONGDASH', 'SOLID', 'SYMBOLS', 'THICK'/

C     PLOT dictionary:

      DATA PLTDIC
     >   /'BOTH', 'ORIGINAL', 'REVISED'/

C     PRECISION dictionary:

      DATA PREDIC
     >   /'ENGINEERING', 'FULL'/

C     THREED dictionary:

      DATA THRDIC
     >   /'FALSE', 'NO:FALSE', 'TRUE', 'YES:TRUE'/


C  *  Execution:

C  *  Set defaults:

      FORMAT = 0
      PRECISION = 1

      DO J = 1, MXCRVS
         CPSLINE (J) = BLANK
         CRVLINE (J) = BLANK
         PLTLINE (J) = BLANK
      END DO

      MAXCRV = +5.
      MINCRV = -5.
      WRAPCRV = .FALSE.
      PLTORIG = .TRUE.
      PLTREV = .TRUE.
      THREED = .FALSE.
      XAXIS = 6.4
      XMIN = UNDEF
      XMAX = UNDEF
      YMIN = UNDEF
      YMAX = UNDEF
      DATFIL = .TRUE.
      PLTFIL = .TRUE.
      CRVFIL = .TRUE.
      TABFIL = .TRUE.
      YPPFIL = .TRUE.
      CPSFIL = .TRUE.
      SPREADFIL = .TRUE.

C  *  No control file was found if LUNIN < 0.  Return with defaults:

      IF (LUNIN < 0) GO TO 950


    5 CONTINUE

C  *     Read one line of input and break into pairs of keywords and values:

         CALL GETLINE (LUNIN, '!', BUFFER, LAST, IOS)
         IF (IOS  <  0) GO TO 950   ! Normal EOF
         IF (IOS /=  0) GO TO 980   ! Read error
         IF (LAST == 0) GO TO 5

         NPAIRS = NKEYDC
         FIRST = 1

         CALL PAIRS (BUFFER (1 : LAST), NPAIRS, KEYS, VALUES)

         DO I = 1, NPAIRS

C  *        Look for a keyword in the dictionary.

            KEYS (I) (6:) = BLANK        ! Only first 5 are significant
            CALL LOOKUP (NKEYDC, KEYDIC, .TRUE., KEYS (I), POINTR)
            IF (POINTR <= 0) GO TO 960

C  *        Read the accompanying value(s) into an internal control
C           variable/array, or find a text value in the appropriate
C           secondary dictionary and assign accordingly:

            GO TO (10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120,
     >             130, 140, 150, 160, 170) POINTR

   10          CONTINUE  !  --CPSLINE--

C  *           Here's a special case. We need to find the beginning of the
C              first value and separate the rest of the tokens in the buffer.

               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               FIRST = MARK + 2
               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               NUMBER = MXLIST
               CALL TOKENS (BUFFER (FIRST:LAST), NUMBER, LIST)

               DO J = 1, NUMBER

C  *              Find each value in the dictionary and assign accordingly:
   
                  CALL LOOKUP (NLINDC, LINDIC, .TRUE., LIST (J), POINTR)
                  IF (POINTR <= 0) GO TO 970
                  CPSLINE (J) = LIST (J)
               END DO
               GO TO 5

   20          CONTINUE   !  --CRVLINE--

C  *           Another special case. (See CPSLINE above)

               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               FIRST = MARK + 2
               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               NUMBER = MXLIST
               CALL TOKENS (BUFFER (FIRST:LAST), NUMBER, LIST)

               DO J = 1, NUMBER

C  *              Find each value in the dictionary and assign accordingly:
   
                  CALL LOOKUP (NLINDC, LINDIC, .TRUE., LIST (J), POINTR)
                  IF (POINTR <= 0) GO TO 970
                  CRVLINE (J) = LIST (J)
               END DO
               GO TO 5

   30          CONTINUE   !  --CURVATURE--
   40          CONTINUE   !  --DERIVATIVES--

C  *           These two are synonymous controls for derivative estimates.

               VALUES (I) (6:) = BLANK
               CALL LOOKUP (NDERDC, DERDIC, .TRUE., VALUES (I), POINTR)
               IF (POINTR <= 0) GO TO 970

               WRAPCRV = VALUES (I) == 'SPLINE'
               CYCLE


   50          CONTINUE   !  --FORMAT--

C  *           Pass only the first 5 characters of the value in this case.
C              We want to recognize any of a number of synonyms:

               VALUES (I) (6:) = BLANK
               CALL LOOKUP (NFMTDC, FMTDIC, .TRUE., VALUES (I), POINTR)
               IF (POINTR <= 0) GO TO 970
               IF (VALUES (I) == 'PROFI') THEN
                  FORMAT = 1
               ELSE IF (VALUES (I) == 'CLOCK') THEN
                  FORMAT = 2
               ELSE IF (VALUES (I) == 'COUNT') THEN
                  FORMAT = 3
               ELSE IF (VALUES (I) == 'THREE') THEN
                  FORMAT = 4
               ELSE
                  FORMAT = 0
               END IF
               CYCLE

   60          CONTINUE   !  --MAXCURVATURE--
               READ (VALUES (I), 1000) MAXCRV
               CYCLE

   70          CONTINUE   !  --MINCURVATURE--
               READ (VALUES (I), 1000) MINCRV
               CYCLE

   80          CONTINUE   !  --NOFILE--

C  *           Another special case (see CPSLINE):

               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               FIRST = MARK + 2
               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               NUMBER = NFILDC
               CALL TOKENS (BUFFER (FIRST:LAST), NUMBER, LIST)

               DO J = 1, NUMBER
                  CALL LOOKUP (NFILDC, FILDIC, .TRUE., LIST (J), POINTR)
                  IF (POINTR <= 0) GO TO 970

                  VAL = LIST (J) (1:3)
                  IF (VAL == 'CPS') THEN
                     CPSFIL = .FALSE.
                  ELSE IF (VAL == 'CRV') THEN
                     CRVFIL = .FALSE.
                  ELSE IF (VAL == 'DAT') THEN
                     DATFIL = .FALSE.
                  ELSE IF (VAL == 'PLT') THEN
                     PLTFIL = .FALSE.
                  ELSE IF (VAL == 'TAB') THEN
                     TABFIL = .FALSE.
                  ELSE IF (VAL == 'YPP') THEN
                     YPPFIL = .FALSE.
                  ELSE IF (VAL == 'SPR') THEN
                     SPREADFIL = .FALSE.
                  ELSE IF (VAL == 'NON') THEN
C                    Do nothing.  'None' allows NOFILE keyword to be present
C                    even if all output files are desired.
                  END IF
               END DO
               GO TO 5

   90          CONTINUE   !  --PLOT--

               CALL LOOKUP (NPLTDC, PLTDIC, .TRUE., VALUES (I), POINTR)
               IF (POINTR <= 0) GO TO 970
               IF (VALUES (I) == 'ORIGINAL') THEN
                  PLTREV = .FALSE.
               ELSE IF (VALUES (I) == 'REVISED') THEN
                  PLTORIG = .FALSE.
               END IF
               CYCLE

  100          CONTINUE   !  --PLTLINE--

C  *           Another special case. (See CPSLINE above.)

               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               FIRST = MARK + 2
               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               NUMBER = MXLIST
               CALL TOKENS (BUFFER (FIRST:LAST), NUMBER, LIST)

               DO J = 1, NUMBER

C  *              Find each value in the dictionary and assign accordingly:
   
                  CALL LOOKUP (NLINDC, LINDIC, .TRUE., LIST (J), POINTR)
                  IF (POINTR <= 0) GO TO 970
                  PLTLINE (J) = LIST (J)
               END DO
               GO TO 5

  110          CONTINUE   !  --PRECISION--

               CALL LOOKUP (NPREDC, PREDIC, .TRUE., VALUES (I), POINTR)
               IF (POINTR <= 0) GO TO 970
               IF (VALUES (I) == 'ENGINEERING') THEN
                  PRECISION = 2
               END IF
               CYCLE

  120          CONTINUE   !  --THREED--

               CALL LOOKUP (NTHRDC, THRDIC, .TRUE., VALUES (I), POINTR)
               IF (POINTR <= 0) GO TO 970
               IF (VALUES (I) == 'TRUE') THEN
                  THREED = .TRUE.
               END IF
               CYCLE

  130          READ (VALUES (I), 1000) XAXIS
               CYCLE

  140          READ (VALUES (I), 1000) XMAX
               CYCLE

  150          READ (VALUES (I), 1000) XMIN
               CYCLE

  160          READ (VALUES (I), 1000) YMAX
               CYCLE

  170          READ (VALUES (I), 1000) YMIN

         END DO ! Next keyword/value pair

C  *     Look for another line of keywords:

      GO TO 5


  950 RETURN


C  *  Error handling:

  960 WRITE (LUNCRT, 1020) 'keyword', BUFFER
      GO TO 999

  970 WRITE (LUNCRT, 1020) 'value', BUFFER
      GO TO 999

  980 WRITE (LUNCRT, 1030) IOS

  999 STOP

C  *  Formats:

 1000 FORMAT (BN, F10.0)
 1020 FORMAT (' Abnormal termination - invalid ', A,
     >         ' in the following line:', /, 1X, A)
 1030 FORMAT (/, ' System error reading a line of keyword text.  IOS: ',
     >        I6)
      END SUBROUTINE RDKEYS
C+-----------------------------------------------------------------------

      SUBROUTINE RECTIFY (NU, XU, YU, NL, XL, YL, X, Y, MAXPTS, LUNCRT,
     >                    LUNKBD )
C
C  PURPOSE:  RECTIFY reorganizes airfoil geometry data so that the common
C            leading edge point is in fact the minimum abscissa.   It can
C            also shift the ordinates to ensure a user-specified leading-
C            edge ordinate, if desired.
C 
C  METHOD:
C   *  Merge the separate surfaces as a single wrap-around surface.
C   *  Identify the subscript corresponding to minimum abscissa.
C   *  Separate into two surfaces again, with truly monotonically increasing
C      abscissas.   
C   *  Shift ordinates if required.
C
C  ARGUMENTS:
C   ARG    DIM  I/O/S   DESCRIPTION
C   NU      -     I     Number of points on upper surface before and after
C                       rearranging
C   XU   MAXPTS  I/O    Upper surface abscissas in ascending order
C   YU   MAXPTS  I/O    Upper surface ordinates in ascending order
C   NL      -    I      Number of points on lower surface before and after
C                       rearranging
C   XL   MAXPTS  I/O    Lower surface abscissas in ascending order
C   YL   MAXPTS  I/O    Lower surface ordinates in ascending order
C   X   MAXPTS*2    S   Abscissas of both surfaces (wrap-around order)
C   Y   MAXPTS*2    S   Ordinates of both surfaces (wrap-around order)
C  MAXPTS   -    I      Maximum number of points allowed for on a surface
C  LUNCRT   -    I      Logical unit for prompts (screen)
C  LUNKBD   -    I      Logical unit for responses (keyboard)
C
C  EXTERNAL REFERENCES:
C
C  COPY      Copies one array to another
C  RVERSE    Reverses order of lower surface points for wrap-around ordering
C  READER    Prompts for/reads integer, real, etc.; handles <CR>
C
C  HISTORY:
C
C    04/29/83    LJC     Original design and coding
C    02/13/85    LJC     Added shifting of ordinates
C    02/26/85    LJC     Added protection against shifting thicknesses
C
C  AUTHOR: Leslie Collins, Informatics General Corporation, Palo Alto, CA
C
C----------------------------------------------------------------------

      IMPLICIT NONE
    
C  *  Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, MAXPTS, NL, NU

      REAL
     >   XU(MAXPTS), YU(MAXPTS), XL(MAXPTS), YL(MAXPTS),  X(MAXPTS*2),
     >   Y(MAXPTS*2)

C  *  Local variables:

      INTEGER
     >   I, NPTS

      REAL
     >   LE, SHIFT

      CHARACTER
     >   YESNO

      LOGICAL
     >   DEFAULT, QUIT
     
C  *  Execution:

C  *  Store coordinates of both surfaces in one wrap-around array:

      CALL RVERSE (NL, XL, X)
      CALL RVERSE (NL, YL, Y)
      CALL COPY   (NU, XU, X(NL))
      CALL COPY   (NU, YU, Y(NL))
      NPTS = NU + NL - 1

C  *  Search for subscript of minimum abscissa:

      DO I = 2, NPTS
         IF (X(I) < X(I-1)) NL = I
      END DO

      NU = NPTS - NL + 1

C  *  Reverse lower surface points so that abscissas are in ascending 
C  *  order:

      CALL RVERSE (NL, X, XL)
      CALL RVERSE (NL, Y, YL)

C  *  Store remaining elements of X and Y in XU and YU arrays:

      CALL COPY (NU, X(NL), XU)
      CALL COPY (NU, Y(NL), YU)

C  *  Permit vertical shift if desired:

      WRITE (LUNCRT, 1000)
     >   ' You have an option to shift the airfoil vertically.'
      WRITE (LUNCRT, 1001)
     >   ' The current leading edge ordinate is ', YU(1)
      CALL READR (LUNCRT,
     >   'Enter a new leading edge ordinate or <CR> to leave as is: ',
     >   LUNKBD, LE, DEFAULT, QUIT)

      IF (.NOT. DEFAULT) THEN

         SHIFT = YU(1) - LE
         DO I = 1, NU
            YU(I) = YU(I) - SHIFT
         END DO

C  *     Protect user from shifting a thickness distribution:

         WRITE (LUNCRT, 1000)
     >      ' Is this a camber/thickness distribution? (Y/N, <CR>=N)'
         CALL READC (LUNCRT,
     >      '(If so, the thickness will not be shifted.) ',
     >      LUNKBD, YESNO, DEFAULT, QUIT)

         IF (YESNO == 'N' .OR. DEFAULT) THEN
            DO I = 1, NL
               YL(I) = YL(I) - SHIFT
            END DO
         END IF

      END IF

 1000 FORMAT (/, A)
 1001 FORMAT (A, F10.6)

      END SUBROUTINE RECTIFY
C+------------------------------------------------------------------------------
C
      SUBROUTINE REDISTRIB (NU, XU, YU, NL, XL, YL, X, Y, MAXPTS,
     >                      XEVAL, YEVAL, COEFS, LUNCRT, LUNKBD, LUNXGR)
C
C  PURPOSE:  REDISTRIB redistributes the data points defining a given
C            airfoil profile, using one of several possible methods 
C            (uniform, sinusoidal, or Vinokur-type distributions).
C            If the leading edge is rounded, the distributions may be
C            specified as being along the arc if this is preferred over
C            the standard chord-wise distribution.
C
C  METHOD:   A menu is used to determine how the revised coordinates are
C            produced.  The new Xs (or Ts) are either generated over the
C            original data range in a choice of ways, or read from a disk
C            file in standard PROFILE format (TITLE and NU followed by
C            a leading-to-trailing-edge distribution, then either NL
C            and a second leading-to-trailing-edge distribution, or NL=0
C            or, equivalently, EOF).  Thus the disk file may contain points 
C            for one or two surfaces.  (A single distribution read will be
C            reused if necessary.)  Any columns after the first are ignored.
C
C            If the airfoil is treated as a single curve, the parametric
C            conventional spline of PSFIT is applied.  Otherwise, several
C            choices of spline are offered.  The original coordinates are
C            overwritten by the redistributed points.
C
C  ARGUMENTS:
C   ARG      DIM    TYPE I/O/S DESCRIPTION
C  NU,NL      -       I   I/O  Number of upper/lower surface pts.
C  XU,XL    NU,NL     R   I/O  Given abscissas, upper/lower.
C  YU,YL    NU,NL     R   I/O  Corresponding ordinates, upper/lower.
C                              The surfaces must have common leading edge
C                              points if they are to be treated as rounded,
C                              and common trailing edge points if they are
C                              to be treated as closed smoothly.
C  X,Y     NU+NL-1    R    S   Workspace for the wraparound form of the
C                              airfoil, if a parametric spline is used.
C  MAXPTS     -       I    I   Max. no. of pts. provided for on 1 surface.
C  XEVAL,YEVAL MAXPTS R    S   Workspace for new abscissas/ordinates.
C  COEFS   MAXPTS*3   R    S   Coefficients used by nonparametric splines;
C                              handy as workspace for the parametric spline.
C  LUNCRT     -       I    I   Logical unit for prompts.
C  LUNKBD     -       I    I   Logical unit for responses.
C  LUNXGR     -       I    I   Logical unit for disk file of input
C                              abscissas or Ts in standard PROFILE format.
C                              (Ordinates are ignored if present.)
C  PROCEDURES:
C    COPY         Copies an array
C    CSEVAL       Evaluates the spline previously fit by CSFIT
C    CSFIT        Fits a conventional cubic interpolating spline
C    DSTRIB       Generates indicated distribution
C    FOILGRD      Linear + Quadratic + Sine + Cosine distribution
C    GETDIS       Prompts for details of generating a distribution
C    GETLINE      File reading utility
C    GETXFORM     Shift/scale factors for transforming [A, B] to [P, Q]
C    HTDIS2       Vinokur's 1-D distribution method
C    LCSFIT       Local cubic spline utility with monotonic and
C                 piecewise linear options as well
C    OPENER       File opening utility
C    PSFIT        Fits an parametric conventional or monotonic cubic spline
C    PSEVAL       Evaluates the PSFIT spline at specified Xs
C    PSTVAL       Evaluates the PSFIT spline at specified Ts
C    READER       Prompting utility
C    RVERSE       Reverses the order of an array
C    TOKENS       Tokenizes a string
C    TSUBJ        Extracts parameter value T associated with data point J
C    USESCALE     Applies the transformation from GETXFORM
C
C  HISTORY:
C  10/28/83  DAS  Code extracted as a module from early PROFILE.
C  07/20/84  LJC  Replaced calls to SPLINE and SEVAL with CSFIT and CSEVAL.
C  08/02/84  LJC  Moved opening of abscissa file from main program and added
C                 prompts for NUEVAL and NLEVAL (originally in namelist).
C  10/20/84  LJC  Modified calls to new versions of PSFIT and PSEVAL.
C  02/19/85  DAS  Handled rounded trailing edges using PSFIT's CLOSED option;
C                 confined rounded leading-edge handling to this routine,
C                 where BLUNT had originally been an input.
C  06/17/85  DAS  Incorporated GETDIS/DSTRIB in place of XGRID and in-line
C                 prompts. (DSTRIB offers generalized distributions -
C                 fractional powers of cosine allowed.)
C  07/23/85  RGL  Allowed for possibility that lower surface only
C                 would be redistributed from a file of abscissas.
C  09/05/85  DAS  Adjusted [T1,T2] definition to reflect modified
C                 parametric spline routines, which use cumulative
C                 arc length now for the parametric variable.
C  08/25/88  DAS  Provided for choosing a monotonic spline in place of
C                 the conventional spline, for nonparametric cases.
C  02/07/90  DAS  Introduced OPENER and GETLINE.
C  03/07/90  DAS  Bug fix: GETDIS return of "quit" was not acted upon.
C  03/23/91  DAS  Replaced MSFIT with LCSFIT in order to provide a
C                 piecewise linear option for wedge-type sections.
C  10/16/91  DAS  PSFIT's calling sequence changed (METHOD); TOL changed
C                 to 1.E-6 instead of 1.E-5 since it is used relatively now.
C  05/04/92  DAS  Provided Vinokur's distribution as an option; provided
C                 for distributions in terms of arc length as well as X;
C                 reused upper surface points read if none found for the
C                 lower surface.
C  08/24/93  DAS  Arc-based redistribution now uses monotonic spline for
C                 X vs. T to avoid possible excursion at the leading edge.
C                 Y vs. T remains conventional.
C  04/01/95  DAS  Installed FOILGRID option.
C  10/18/96  DAS  Replaced FOILGRID with FOILGRD.
C
C  AUTHORS: David Saunders/Leslie Collins, Sterling/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNXGR, MAXPTS, NL, NU

      REAL
     >   COEFS (MAXPTS, 3), X (MAXPTS*2), XEVAL (MAXPTS), XL (NL),
     >   XU (NU), Y (MAXPTS*2), YEVAL (MAXPTS), YL (NL), YU (NU)

C     Local constants:

      REAL, PARAMETER ::
     >   PERCENT = 1.E-2,
     >   TOL     = 1.E-6,
     >   ZERO    = 0.E+0

      CHARACTER, PARAMETER ::
     >   BLANK * 1   = ' ',
     >   COMMENT * 1 = '!',
     >   IFMT * 8    = '(BN,I20)',
     >   RFMT * 10   = '(BN,F20.0)'

C     Local variables:

      INTEGER
     >   I, IER, IOS, LAST, MODE, NLEVAL, NUEVAL, NUMBER, NP, NPTS

      REAL
     >   A, B, DUMMY(1), P (4), SCALE(1), SHIFT(1), T1, T2

      LOGICAL
     >   ARC, BLUNT, CLOSED, DEFAULT, OPENED, QUIT, REUSE, YES

      CHARACTER
     >   METHOD * 1, TOKEN * 25, XFILE * 50

C     Procedures:

      REAL
     >   TSUBJ

      EXTERNAL
     >   COPY, CSEVAL, CSFIT, DSTRIB, FOILGRD, GETDIS, GETLINE,
     >   GETXFORM, HTDIS2, LCSFIT, OPENER, PSFIT, PSEVAL, PSTVAL,
     >   READC, READY, RVERSE, TOKENS, TSUBJ, USESCALE

C     Execution:

C     Determine whether airfoil is to be treated as one curve or two:

      WRITE (LUNCRT, 1001)
     >   '    PROFILE treats the airfoil as one curve if the leading ',
     >   '    edge is rounded. In this case, if the trailing edge is ',
     >   '    also rounded, the curve is closed smoothly. Otherwise, ',
     >   '    the trailing edge is sharp or open (no continuity).',
     >   BLANK

      BLUNT = .TRUE.
      CALL READY (LUNCRT,
     >   'Is the  LEADing edge rounded? (Y/N; <CR> means Y(es)): ',
     >   LUNKBD, BLUNT, DEFAULT, QUIT)
      IF (QUIT) GO TO 999

      CLOSED = .FALSE.

      IF (BLUNT) THEN
         CALL READY (LUNCRT,
     >      'Is the TRAILing edge rounded? (Y/N; <CR> means NO): ',
     >      LUNKBD, CLOSED, DEFAULT, QUIT)
         IF (QUIT) GO TO 999

      ELSE  ! Don't try to handle a sharp leading edge/rounded trailing edge

C        Non-rounded airfoils may require a choice of interpolation methods:

         WRITE (LUNCRT, 1001)
     >     '    Conventional splines are normally appropriate, but the',
     >     '    monotonic spline or piecewise linear interpolation may',
     >     '    be preferable for angular sections.',
     >     BLANK

  100    METHOD = 'C'
         CALL READC (LUNCRT,
     > 'Conventional spline, Monotonic, or Linear? (C/M/L; <CR> = C): ',
     >      LUNKBD, METHOD, DEFAULT, QUIT)
         IF (QUIT) GO TO 999
         IF (METHOD /= 'C' .AND.
     >       METHOD /= 'M' .AND. METHOD /= 'L') GO TO 100
      END IF

      ARC = .FALSE.          ! For nonparametric cases

      IF (BLUNT) THEN        ! Parametric spline permits two choices:

         METHOD = 'S'
         CALL READC (LUNCRT,
     >      'Distribute along the chord or the arc? (X/S; <CR> = S): ',
     >      LUNKBD, METHOD, DEFAULT, QUIT)
         IF (QUIT) GO TO 999
         ARC = METHOD /= 'X'

C        Set up airfoil as a single wraparound curve, eliminating the
C        duplicate leading edge point.

         IF (XU (1) /= XL (1)  .OR.  YU (1) /= YL (1)) GO TO 830
         IF (CLOSED  .AND.
     >      (XU (NU) /= XL (NL) .OR. YU (NU) /= YL (NL))) GO TO 840
       
         CALL RVERSE (NL, XL, X)
         CALL RVERSE (NL, YL, Y)
         CALL COPY (NU - 1, XU (2), X (NL + 1))
         CALL COPY (NU - 1, YU (2), Y (NL + 1))
         NPTS = NL + NU - 1

C        Fit a parametric spline, using the "monotonic" option for X vs. T
C        (and a conventional spline for Y vs. T), to avoid possible trouble
C        at the leading edge of blunt airfoils.

         CALL PSFIT (NPTS, X, Y, 'MC', CLOSED, IER)
         IF (IER /= 0) GO TO 900

      END IF

      OPENED = .FALSE.   ! For possible input file of Xs or Ts.


C     Upper surface:
C     -------------

  200 CONTINUE

      YES = .TRUE.
      CALL READY (LUNCRT,
     >   'Do you want to redistribute the UPPER surface? (<CR>=Y): ',
     >   LUNKBD, YES, DEFAULT, QUIT)
      IF (QUIT) GO TO 999

      IF (YES) THEN

         IF (ARC) THEN
            A = TSUBJ (NL)    ! Cumulative chord length limits from PSFIT
            B = TSUBJ (NPTS)  ! This is the range of T to be distributed
         ELSE
            A = XU (1)        ! This is the X range to be distributed
            B = XU (NU)
         END IF

C        Determine the details for this redistribution:

         NP = 2
         CALL GETDIS (LUNCRT, LUNKBD, MODE, NP, P, NUEVAL)

         IF (MODE < -1) THEN

C           User must have changed his mind - quit:

            GO TO 999

         ELSE IF (MODE == -1) THEN

C           Read from file with new Xs (or Ts; standard PROFILE format):

            XFILE = 'xgrid.dat'
            CALL OPENER (LUNCRT,
     >         'Enter filename for new Xs or Ts. <CR> uses xgrid.dat:',
     >         LUNKBD, XFILE, LUNXGR, 'OLD')

            OPENED = .TRUE.

            READ (LUNXGR, 1000, ERR = 870)     ! Skip title

C           Get the no. of pts.  XFILE is a handy character buffer.

            CALL GETLINE (LUNXGR, BLANK, XFILE, LAST, IOS)
            IF (IOS /= 0) GO TO 870
            NUMBER = 1
            CALL TOKENS (XFILE, NUMBER, TOKEN)
            READ (TOKEN, IFMT, ERR = 870) NUEVAL

            IF (NUEVAL > MAXPTS) GO TO 850

C           Pick off the first item on each line:

            DO 320, I = 1, NUEVAL
  310          CALL GETLINE (LUNXGR, COMMENT, XFILE, LAST, IOS)
               IF (IOS /= 0) GO TO 870
               IF (LAST == 0) GO TO 310

               CALL TOKENS (XFILE, NUMBER, TOKEN)
               READ (TOKEN, RFMT, ERR = 870) XEVAL (I)
  320       CONTINUE

            IF (.NOT. ARC) THEN   ! Require a full range of abscissas

               IF (XEVAL (1) /= XU (1) .OR.
     >             XEVAL (NUEVAL) /= XU (NU)) GO TO 880

            ELSE  ! Distribution must be l.e. to t.e. but may be relative

               CALL GETXFORM (XEVAL (1), XEVAL (NUEVAL), A, B,
     >                        SCALE(1), SHIFT(1))

               CALL USESCALE ('D', 1, NUEVAL, XEVAL, DUMMY, DUMMY,
     >                        SCALE, SHIFT, IER)

               XEVAL (1) = A      ! Avoid round-off at the end points
               XEVAL (NUEVAL) = B

            END IF

         ELSE     ! Generate the specified distribution

            IF (NUEVAL > MAXPTS) GO TO 850

            IF (MODE <= 3) THEN   ! Uniform or sinusoidal-type

               CALL DSTRIB (MODE, NP, P, NUEVAL, A, B, XEVAL)

            ELSE IF (MODE == 4) THEN   ! Vinokur distribution

C              Relative or absolute?

               IF (P (1) < ZERO) P (1) = -PERCENT * (B - A) * P (1)
               IF (P (2) < ZERO) P (2) = -PERCENT * (B - A) * P (2)

               CALL HTDIS2 (.TRUE., A, B, P (1), P (2), NUEVAL, XEVAL,
     >                      -LUNCRT, IER)
               IF (IER /= 0) GO TO 855

            ELSE IF (MODE == 5) THEN   ! Linear + Quadratic + Sine + Cosine

               CALL FOILGRD (NUEVAL, A, B, P (1), P (2), P (3), P (4),
     >                       XEVAL)
            END IF

         END IF


C        Now evaluate the new upper surface coordinates:

         IF (BLUNT) THEN

            IF (ARC) THEN     ! Evaluate X and Y at the indicated Ts

C              Make use of COEFS (*, 1) for the Ts, and COEFS (*, 2)
C              and COEFS (*, 3) for the unneeded 1st & 2nd derivatives.

               CALL COPY (NUEVAL, XEVAL, COEFS)

               CALL PSTVAL (NUEVAL, COEFS (1, 1), XEVAL, YEVAL,
     >                      COEFS (1, 2), COEFS (1, 3),
     >                      COEFS (1, 2), COEFS (1, 3), NPTS, X, Y)

            ELSE              ! Evaluate Y at the indicated Xs

               T1 = TSUBJ (NL)   ! Define the relevant curve that is
               T2 = TSUBJ (NPTS) ! monotonic in X

               CALL PSEVAL (NUEVAL, XEVAL, YEVAL, COEFS, COEFS, T1,
     >                      T2, NPTS, X, Y, TOL, IER)
               IF (IER /= 0) GO TO 910

            END IF

         ELSE

C           Separate spline for each surface (conventional, monotonic or linear)
  
            IF (METHOD == 'C') THEN

               CALL CSFIT (NU, XU, YU, 0, DUMMY(1), 0, DUMMY(1), 
     >             COEFS (1,1), COEFS (1,2), COEFS (1,3), IER)
               IF (IER /= 0) GO TO 890

               CALL CSEVAL (NU, XU, YU, NUEVAL, XEVAL, COEFS (1,1),
     >                      COEFS (1,2), COEFS (1,3), YEVAL)
            ELSE

               CALL LCSFIT (NU, XU, YU, .TRUE., METHOD, NUEVAL, XEVAL,
     >                      YEVAL, YEVAL)
            END IF

         END IF

C        Overwrite the original coordinates with the redistributed ones:

         CALL COPY (NUEVAL, XEVAL, XU)
         CALL COPY (NUEVAL, YEVAL, YU)
         NU = NUEVAL

      END IF


C     Repeat for lower surface:
C     ------------------------

  400 CONTINUE

      YES = .TRUE.
      CALL READY (LUNCRT,
     >   'Do you want to redistribute the LOWER surface? (<CR>=Y): ',
     >   LUNKBD, YES, DEFAULT, QUIT)
      IF (QUIT) GO TO 999

      IF (YES) THEN

C        Interval to be distributed:

         IF (ARC) THEN
            A = ZERO
            B = TSUBJ (NL)
         ELSE
            A = XL (1)
            B = XL (NL)
         END IF

C        Determine the details for this redistribution:

         NP = 2
         CALL GETDIS (LUNCRT, LUNKBD, MODE, NP, P, NLEVAL)

         IF (MODE < -1) THEN

C           User must have changed his mind - quit:

            GO TO 999

         ELSE IF (MODE == -1) THEN

            IF (.NOT. OPENED) THEN

C              Open file of new Xs or Ts:

               XFILE = 'xgrid.dat'
               CALL OPENER (LUNCRT,
     >            'Enter filename for new abscissas. <CR>=xgrid.dat:',
     >            LUNKBD, XFILE, LUNXGR, 'OLD')

C              Skip upper surface (unused but may be present).
C              HOWEVER: Store it in case only one surface is present.

               READ (LUNXGR, 1000, ERR = 870)  ! Skip title

               CALL GETLINE (LUNXGR, BLANK, XFILE, LAST, IOS)
               IF (IOS /= 0) GO TO 870
               NUMBER = 1
               CALL TOKENS (XFILE, NUMBER, TOKEN)
               READ (TOKEN, IFMT, ERR = 870) NUEVAL

               DO 490, I = 1, NUEVAL
  480             CALL GETLINE (LUNXGR, COMMENT, XFILE, LAST, IOS)
                  IF (IOS /= 0) GO TO 870
                  IF (LAST == 0) GO TO 480

                  CALL TOKENS (XFILE, NUMBER, TOKEN)
                  READ (TOKEN, RFMT, ERR=870) XEVAL (I)
  490          CONTINUE
            END IF

C           Read number of lower surface abscissas (if any):

            CALL GETLINE (LUNXGR, BLANK, XFILE, LAST, IOS)
            IF (IOS < 0) THEN
               REUSE = .TRUE.     ! EOF means same as NL = 0
            ELSE IF (IOS /= 0) THEN
               GO TO 870
            ELSE
               NUMBER = 1
               CALL TOKENS (XFILE, NUMBER, TOKEN)
               READ (TOKEN, IFMT, ERR = 870) NLEVAL

               REUSE = NLEVAL == 0
               IF (NLEVAL > MAXPTS) GO TO 860
            END IF

            IF (REUSE) THEN

               NLEVAL = NUEVAL

               IF (ARC) THEN
                  IF (OPENED) THEN

C                    Recover the relative Ts from the upper surface.
C                    They have been transformed, but transforming then
C                    again is OK apart from round-off.

                     CALL COPY (NLEVAL, COEFS, XEVAL)
                     WRITE (LUNCRT,1001) ' Reusing upper surface Ts ...'
                  END IF
               ELSE IF (OPENED) THEN
                  WRITE (LUNCRT, 1001) ' Reusing upper surface Xs ...'
               END IF               

            ELSE  ! Read lower surface Xs (first item on each valid line)

               DO 520, I = 1, NLEVAL
  510             CALL GETLINE (LUNXGR, COMMENT, XFILE, LAST, IOS)
                  IF (IOS /= 0) GO TO 870
                  IF (LAST == 0) GO TO 510

                  CALL TOKENS (XFILE, NUMBER, TOKEN)
                  READ (TOKEN, RFMT, ERR = 870) XEVAL (I)
  520          CONTINUE

            END IF

            IF (.NOT. ARC) THEN   ! Require a full range of abscissas

               IF (XEVAL (1) /= XL (1) .OR.
     >             XEVAL (NLEVAL) /= XL (NL)) GO TO 880

            ELSE  ! Transform relative arc distribution to desired interval

               CALL GETXFORM (XEVAL (1), XEVAL (NLEVAL), A, B,
     >                        SCALE(1), SHIFT(1))

               CALL USESCALE ('D', 1, NLEVAL, XEVAL, DUMMY, DUMMY,
     >                        SCALE, SHIFT, IER)

               XEVAL (1) = A      ! Avoid round-off
               XEVAL (NLEVAL) = B

            END IF

         ELSE     ! Generate the specified distribution

            IF (NLEVAL > MAXPTS) GO TO 860

            IF (MODE <= 3) THEN  ! Uniform or sinusoidal-type

               CALL DSTRIB (MODE, NP, P, NLEVAL, A, B, XEVAL)

            ELSE IF (MODE == 4) THEN  ! Vinokur distribution

C              Relative or absolute?

               IF (P (1) < ZERO) P (1) = -PERCENT * (B - A) * P (1)
               IF (P (2) < ZERO) P (2) = -PERCENT * (B - A) * P (2)

               CALL HTDIS2 (.TRUE., A, B, P (1), P (2), NLEVAL, XEVAL,
     >                      -LUNCRT, IER)
               IF (IER /= 0) GO TO 865

            ELSE IF (MODE == 5) THEN   ! Linear + Quadratic + Sine + Cosine

               CALL FOILGRD (NLEVAL, A, B, P (1), P (2), P (3), P (4),
     >                       XEVAL)
            END IF

         END IF


C        Now evaluate the new lower surface coordinates:

         IF (BLUNT) THEN
 
            IF (ARC) THEN     ! Evaluate X and Y at the indicated Ts

C              [0, T(l.e.)] has been treated as though it were
C              [T(l.e.), T(t.e.)], so the Ts need to be transformed.
C              (We want to come out with X/YEVAL(*) going from l.e. to t.e.)

               CALL GETXFORM (XEVAL (1), XEVAL (NLEVAL), B, A,
     >                        SCALE(1), SHIFT(1))

               CALL COPY (NLEVAL, XEVAL, COEFS) ! Use COEFS (*, 1) for Ts

               CALL USESCALE ('D', 1, NLEVAL, COEFS, DUMMY, DUMMY,
     >                        SCALE, SHIFT, IER)

               COEFS (1, 1) = B      ! Avoid round-off
               COEFS (NLEVAL, 1) = A

               CALL PSTVAL (NLEVAL, COEFS (1, 1), XEVAL, YEVAL,
     >                      COEFS (1, 2), COEFS (1, 3),
     >                      COEFS (1, 2), COEFS (1, 3), NPTS, X, Y)

            ELSE              ! Evaluate Y at the indicated Xs

               T1 = ZERO      ! Define monotonic-in-X sub-curve
               T2 = TSUBJ (NL)

               CALL PSEVAL (NLEVAL, XEVAL, YEVAL, COEFS, COEFS, T1,
     >                      T2, NPTS, X, Y, TOL, IER)
               IF (IER /= 0) GO TO 910

            END IF

         ELSE

            IF (METHOD == 'C') THEN

               CALL CSFIT (NL, XL, YL, 0, DUMMY(1), 0, DUMMY(1), 
     >            COEFS (1,1), COEFS (1,2), COEFS (1,3), IER)
               IF (IER /= 0) GO TO 890

               CALL CSEVAL (NL, XL, YL, NLEVAL, XEVAL, COEFS (1,1),
     >                      COEFS (1,2), COEFS (1,3), YEVAL)
            ELSE

               CALL LCSFIT (NL, XL, YL, .TRUE., METHOD, NLEVAL, XEVAL,
     >                      YEVAL, YEVAL)
            END IF

         END IF

C        Overwrite the original coordinates with the redistributed ones:

         NL = NLEVAL
         CALL COPY (NLEVAL, XEVAL, XL)
         CALL COPY (NLEVAL, YEVAL, YL)

      END IF

      RETURN


C     Error handling:

  830 WRITE (LUNCRT, 1001)
     >   ' The two surfaces have different leading points - quitting.'
      GO TO 999
  840 WRITE (LUNCRT, 1001)
     >   ' The two surfaces have different trailing points - quitting.'
      GO TO 999
  850 WRITE (LUNCRT, 1002) MAXPTS
      GO TO 200
  855 WRITE (LUNCRT, 1003) 'HTDIS2', IER
      GO TO 200
  860 WRITE (LUNCRT, 1002) MAXPTS
      GO TO 400
  865 WRITE (LUNCRT, 1003) 'HTDIS2', IER
      GO TO 400
  870 WRITE (LUNCRT, 1001) ' Error reading file of input abscissas.'
      GO TO 999
  880 WRITE (LUNCRT, 1001) ' Must supply full range of abscissas.'
      GO TO 999
  890 WRITE (LUNCRT, 1003) 'CSFIT', IER
      GO TO 999
  900 WRITE (LUNCRT, 1003) 'PSFIT', IER
      GO TO 999
  910 WRITE (LUNCRT, 1003) 'PSEVAL', IER
C*****GO TO 999

  999 WRITE (LUNCRT, 1001) ' Stopping in REDISTRIB.'
      STOP ' '

C     Formats:

 1000 FORMAT (A)
 1001 FORMAT (/, (A))
 1002 FORMAT (/, ' Can''t handle so many abscissas. Limit = ', I3)
 1003 FORMAT (/, ' Error in subroutine ', A, '.  IER: ', I3)

      END SUBROUTINE REDISTRIB
C+----------------------------------------------------------------------
C
      SUBROUTINE REFINE (NU, XU, YU, NL, XL, YL, TITLE, 
     +                   LUNCRT, LUNKBD, LUNTAB, LUNYPP)
C
C  PURPOSE:  REFINE  perturbs the surface(s) of the given airfoil in  a
C            way that strives for a specified thickness while retaining
C            the original curvature (actually 2nd derivative) distribu-
C            tion as much as possible, specially near the leading edge.
C            It also provides for refining  the curvature distribution,
C            offering a variety of ways to indicate target y" values.
C
C  METHOD:   Solve the overdetermined system obtained from the  simple-
C            minded scaling,
C
C                         ynew(i) = yscale(i) * yold(i)    (i = 2, n-1)
C
C            along with constraints obtained from the finite difference
C            expression for the 2nd derivative at each  interior point:
C
C                         ynew"(i) = old or desired y"(i)  (i = 2, n-1)
C
C            Provision is made for weighting (down) the y" portion.
C
C            Note:   ynew(i) = yscale(i)*yold(i) exactly (i = 1 and n)
C
C            because the leading edge point should come out EXACTLY the
C            same for both surfaces (not just close), and the same  MAY
C            be true of a sharp trailing edge.
C
C            The very sparse system involved was originally solved by a
C            dense method (HDESOL, still shown but commented out); this
C            has since been replaced with a specialized  sparse  method
C            (DTDLSQ).
C            
C            Since the desired thickness is not likely to  be  obtained
C            exactly by this linear least squares approach,  the method
C            is iterated until the thickness is arbitrarily close.
C
C            Any target y" values read from a file are assumed to be in
C            the standard PROFILE format for (x,y) pairs, although NU=0
C            or NL=0 is permitted, meaning leave that surface alone.
C
C  ARGUMENTS:
C  VAR   DIM  TYPE I/O/S DESCRIPTION
C  NU     -    I     I   Number of upper surface points
C  XU     NU   R     I   Upper surface abscissas (assumed increasing and
C  YU     NU   R    I/O  Lower surface ordinates             normalized)
C  NL     -    I     I   Number of lower surface points
C  XL     NL   R     I   Lower surface abscissas
C  YL     NL   R    I/O  Lower surface ordinates
C  TITLE  *    C     I   Tabulations title
C  LUNCRT -    I     I   Logical unit for screen
C  LUNKBD -    I     I   Logical unit for keyboard
C  LUNTAB -    I     I   Logical unit for tabulations
C  LUNYPP -    I     I   Logical unit for optional table of y" values
C
C  PARAMETER CONSTANTS:
C  PARAM   TYPE   DESCRIPTION:
C  MAXITER   I    Max. no. iterations for achieving target thickness  
C  MAXN      I    MAXSRF-2, because leading and trailing edge points
C                 are not solved for with the interior points
C  MAXSRF    I    Max. no. pts. on one surface allowed by local arrays
C
C  SIGNIFICANT LOCAL VARIABLES:
C  VAR       DIM     TYPE  DESCRIPTION
C**A       2*N,N+1     R   For setting up 2N*N overdetermined system in
C**                        dense form (N=MAXSRF-2)
C**C         2*N       R   Work-space for HDESOL
C    NOTE: A, B, C, D, are of length MAXSRF, not MAXN, because they are
C          reused in the thickness calculation for spline coefs., etc.
C  A       MAXSRF      R   For sparse solution of overdetermined system
C  B       MAXSRF      R    "    "    "    "    "    "    "    "    "
C  C       MAXSRF      R    "    "    "    "    "    "    "    "    "
C  D       MAXSRF      R    "    "    "    "    "    "    "    "    "
C  R          N        R    "    "    "    "    "    "    "    "    "
C  S          N        R    "    "    "    "    "    "    "    "    "
C  YSCALE   MAXSRF     R   Scale factors applied to each y(I), possibly
C                          functions of I
C  YLORIG,  MAXSRF     R   Needed for refining result by repeating the
C  YUORIG                  analysis with an extrapolated target thickness
C
C  FILES USED:
C  LUN       I/O/S     DESCRIPTION
C  LUNCRT      O       User prompts
C  LUNKBD      I       User responses
C  LUNTAB      O       Tabulations
C  LUNYPP      I       Optional input table of y" targets
C
C  PROCEDURES:
C  CALCTHICK   Calculates airfoil thickness
C  COPY        Copies an array
C  LINTRP      Linear interpolation (for thickness iteration)
C  LSTFIL      Echoes the table of y" values to the tabulation file
C  OPENER      File opening utility
C  READER      Prompting utility
C  REFSRF      Performs the above algorithm on one surface at a time
C  RDXYZ       Used to (x, y") pairs - one surface per call
C
C  HISTORY:
C  06/12/95   DAS   Suppressed printing details of linear system in .TAB file.
C                   Thickness iterations still go to .TAB and screen.
C  07/30/93   DAS   Global SAVE introduced.
C  12/12/91   DAS   RDXYZ's calling sequence changed.
C  07/19/90   DAS   No! PRREAD doesn't permit 0 for either surface's
C                   target y" values.  Introduced RDXYZ instead.
C  02/07/90   DAS   Introduced OPENER; removed concatenation from I/O;
C                   revived use of PRREAD to read the y" values.
C  02/10/87   DAS   XATMAXTHICKNESS=1.0 caused failure (YBUMPPARAMS(1)).
C                   Safeguarded this abnormal case. Also safeguarded use
C                   on unnormalized airfoil - shape functions fail then.
C  09/15/86   DAS   PRREAD no longer used for target y" values - it does
C                   not handle NU=0 or NL=0 properly for this context.
C  04/16/86   DAS   Removed YP, YPP, YK arrays and the call to FD12K -
C                   the finite differencing is done in-line in REFSRF.
C  03/19/86   DAS   Streamlined the diagnostic I/O and control of the
C                   prompting done in REFSRF first time in/each surface.
C  10/24/85   DAS   Used LSTFIL in place of PRWRIT for echoing y" table.
C  11/05/84   DAS   Introduced DTDLSQ in place of HDESOL to save storage
C                   and CPU time.
C  10/25/84   DAS   Optional on-line help added; prompts clarified.
C  04/11/84   DAS   Eliminated handling of non-normalized airfoil.
C  01/30/84   LJC   Incorporated new calling sequence for CALCTHICK.
C  12/21/83   DAS   Made use of PRWRIT to print the table of y" values
C                   in the tabulation output (desirable because it is
C                   easy to be working with a table that is something
C                   other than what you thought it was).
C  11/03/83   DAS   Implemented no-y"-constraint case (as opposed to
C                   leave-surface-alone case - shoot for original y").
C  10/28/83   DAS   Renamed REFINE and integrated into PROFILE as an
C                   option distinct from MODIFY.
C  10/26/83   DAS   Incorporated table look-up for y" distributions.
C  10/20/83   DAS   Provided for retaining original curvature (without
C                   knowing what it is ahead of time).
C  10/19/83   DAS   Revised thickness refinement to handle case where
C                   desired thickness = original thickness.
C  10/12/83   DAS   Put y" weighting and constraining into REFSRF so
C                   that surfaces may be treated differently.
C  10/07/83   DAS   Changed back to true y" for each equation in the
C                   overdetermined system (rather than multiplying
C                   throughout by functions of DXI, DXIM1); bounded
C                   RHS y" values in center region (to flatten the
C                   upper curvature distribution somewhat); rethought
C                   the thickness-refinement iteration.
C  09/30/83   DAS   Added refinement step for more precise thickness;
C                   nonlinear weighting of y"; curvature data plotting.
C  09/27/83   DAS   Separated out REFSRF for treating both surfaces.
C  09/26/83   DAS   Put in fancier scaling using sine "bump" (R.Kennelly).
C  09/23/83   DAS   y" before and after from finite differencing now,
C                   for greater consistency.
C  09/20/83   DAS   Initial design and coding (y" from SPLINE; upper
C                   surface only).
C
C  AUTHOR: David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

      SAVE

C ... Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNTAB, LUNYPP, NL, NU

      REAL
     >   XL (NL), XU (NU), YL (NL), YU (NU)

      CHARACTER
     >   TITLE * (*)

C ... Local constants:

      INTEGER, PARAMETER ::
     >   MAXITER   = 10,
     >   MAXSRF    = 301,
     >   MAXN      = MAXSRF - 2

      REAL, PARAMETER ::
     >   ONE       = 1.E+0,
     >   TOLERANCE = 5.E-5,
     >   ZERO      = 0.E+0

      CHARACTER, PARAMETER ::
     >   BLANK * 1 = ' '

C ... Local variables:

      INTEGER
     >   COLUMNS (2), I, IER, ITER, NPTSYPPTABLE (2)

C*****REAL      A (2*MAXN, MAXN+1), C (2*MAXN),  ! For dense proof-of-concept

      REAL
     >   A (MAXSRF), B (MAXSRF), C (MAXSRF), D (MAXSRF), R (MAXN),
     >   S (MAXN), DELTATHICKNESS, DESIREDTHICKNESS, BETTER1, BETTER2,
     >   ORIGINALTHICKNESS, TARGET1(1), TARGET2(1), TARGET3(1),
     >   XATMAXTHICKNESS,
     >   XYPPTABLE (MAXSRF, 2), YPPTABLE (MAXSRF, 2), YLORIG (MAXSRF),
     >   YPPWEIGHT (MAXSRF), YSCALE (MAXSRF), YUORIG (MAXSRF),
     >   YBUMPPARAMS (3), YBUMPWIDTH

      LOGICAL
     >   DEFAULT, DETAILS, FINIS, KEEPSAMETHICKNESS, LEAVEALONE (2),
     >   NOYPPTABLE (2), PROMPT, QUIT, YES

      CHARACTER
     >   BUFFER * 80, YPPTABLENAME * 50

C ... Procedures:

      EXTERNAL
     >   CALCTHICK, COPY, LINTRP, LSTFIL, OPENER, RDXYZ, READR, READY,
     >   REFSRF

C ... Execution:

      IF (XU (1) < ZERO .OR. XU (NU) > ONE .OR.
     >    XL (1) < ZERO .OR. XL (NL) > ONE) THEN
         WRITE (LUNCRT, 1005)
     >      ' "Refine" mode requires Xs in [0,1]. Use "Normalize."'
         GO TO 990
      END IF

      IF (MAX (NU, NL) > MAXSRF) THEN
         WRITE (LUNCRT, 1005)
     >      ' Not enough work-space. Recompile module REFINE.'
         GO TO 990
      END IF

C ... The following were moved from REFSRF because of problems with
C     SAVE and DATA statements:

      PROMPT = .TRUE.
      LEAVEALONE (1) = .FALSE.
      LEAVEALONE (2) = .FALSE.

C ... On-line help may encourage first-time users:

      YES = .FALSE.
      CALL READY (LUNCRT,
     >   'Do you want an explanation of the prompts to follow?   '//
     >   '(Y/N; <CR>=No): ',
     >   LUNKBD, YES, DEFAULT, QUIT)

      IF (YES) THEN
         WRITE (LUNCRT, 1001)
         CALL READY (LUNCRT,
     >      '                   <Hit RETURN for more.>',
     >      LUNKBD, YES, DEFAULT, QUIT)
         WRITE (LUNCRT, 1002)
      END IF

      DETAILS = .FALSE.        ! Suppress output to .TAB file
C****      CALL READY (LUNCRT,
C****     >   'Do you want full details of the refinement iterations? '//
C****     >   '(Y/N; <CR>=No): ',
C****     >   LUNKBD, DETAILS, DEFAULT, QUIT)

C ... Calculate the original thickness.  Use available workspace
C     for the spline coefficients and spline evaluation.

      CALL CALCTHICK (NL, NU, XL, XU, YL, YU, ORIGINALTHICKNESS,
     >                XATMAXTHICKNESS,
C****>                A (1, 1), A (1, 2), A (1, 3), A (1, 4))
     >                A, B, C, D)

      WRITE (LUNCRT, 1020) '0Original % thickness:', ORIGINALTHICKNESS  

      DESIREDTHICKNESS = ORIGINALTHICKNESS
      CALL READR (LUNCRT,
     >   'Enter desired % thickness, or <CR> to keep same: ',
     >   LUNKBD, DESIREDTHICKNESS, KEEPSAMETHICKNESS, QUIT)

      DELTATHICKNESS = ABS (DESIREDTHICKNESS - ORIGINALTHICKNESS)
      IF (KEEPSAMETHICKNESS .OR. DELTATHICKNESS < 1.E-02) THEN

C ...    Scales on all the ordinates will be essentially 1.0, so
C        pick a reasonable width for the sine bump and don't prompt:

         YBUMPWIDTH = ONE
      ELSE IF (DELTATHICKNESS < 0.2E+0) THEN
         YBUMPWIDTH = 1.5E+0
      ELSE
         YBUMPWIDTH = 2.0E+0
         CALL READR (LUNCRT,
     >      'Enter power of sine to use for y scaling (<CR> = 2.): ',
     >      LUNKBD, YBUMPWIDTH, DEFAULT, QUIT)
      END IF

      WRITE (LUNTAB, 1005)
      WRITE (LUNTAB, 1010) '1Case: ', TITLE, BLANK
      WRITE (LUNTAB, 1020)
     >   ' Original thickness  :', ORIGINALTHICKNESS,
     >   ' Corresponding  x/c  :', XATMAXTHICKNESS,
     >   ' Desired  thickness  :', DESIREDTHICKNESS,
     >   ' y-scaling sine power:', YBUMPWIDTH

      YPPTABLENAME = BLANK
      CALL OPENER (LUNCRT,
     >   'Enter file name for target 2nd derivatives; <CR> if none: ',
     >   LUNKBD, YPPTABLENAME, LUNYPP, 'OLD')

      IF (YPPTABLENAME == BLANK) THEN
         NOYPPTABLE (1) = .TRUE.
         NOYPPTABLE (2) = .TRUE.
      ELSE
         WRITE (LUNTAB, 1010) BLANK, BLANK,
     >      ' Table of y" values used: ', YPPTABLENAME,
     >      ' Values of y" found in table follow:', BLANK, BLANK

         CALL LSTFIL (LUNYPP, LUNTAB, BUFFER)
         REWIND LUNYPP

C ...    Skip the y" file's title:

         READ (LUNYPP, '(A)')

C ...    Look for target y" data one surface at a time:

         COLUMNS (1) = 1
         COLUMNS (2) = 2
         DO I = 1, 2
            CALL RDXYZ (2, LUNYPP, LUNCRT, COLUMNS, MAXSRF,
     >                  NPTSYPPTABLE (I), XYPPTABLE (1, I),
     >                  YPPTABLE (1, I), YPPTABLE (1, I), FINIS, IER)
            IF (IER /= 0) GO TO 920

            NOYPPTABLE (I) = NPTSYPPTABLE (I) == 0
         END DO

      END IF

      CALL COPY (NL, YL, YLORIG)
      CALL COPY (NU, YU, YUORIG)

      YBUMPPARAMS (2) = YBUMPWIDTH
      YBUMPPARAMS (1) = XATMAXTHICKNESS
      IF (XATMAXTHICKNESS == ZERO .OR.
     >    XATMAXTHICKNESS == ONE) YBUMPPARAMS (1) = 0.5E+0

C ... Begin thickness/second-derivative computation - iterated because
C     the desired thickness is not achieved exactly by solution of the
C     overdetermined systems involved.  The target thickness is varied
C     above or below the desired thickness till the computed thickness
C     is arbitrarily close to the desired one.

      ITER = 0
      TARGET3(1) = DESIREDTHICKNESS

  200 CONTINUE

         YBUMPPARAMS (3) = TARGET3(1) / ORIGINALTHICKNESS - ONE

C ...    Modify each surface separately:

         CALL REFSRF (1, NU, XU, YU, YSCALE,
     >                YBUMPPARAMS, YPPWEIGHT,
C****>                A, C,
     >                A, B, C, D, R, S,
     >                NOYPPTABLE (1), NPTSYPPTABLE (1),
     >                XYPPTABLE (1, 1), YPPTABLE (1, 1),
     >                DETAILS, PROMPT, LEAVEALONE (1),
     >                LUNCRT, LUNKBD, LUNTAB)

         CALL REFSRF (2, NL, XL, YL, YSCALE,
     >                YBUMPPARAMS, YPPWEIGHT,
C****>                A, C,
     >                A, B, C, D, R, S,
     >                NOYPPTABLE (2), NPTSYPPTABLE (2),
     >                XYPPTABLE (1, 2), YPPTABLE (1, 2),
     >                DETAILS, PROMPT, LEAVEALONE (2),
     >                LUNCRT, LUNKBD, LUNTAB)

C ...    Reestimate the actual thickness achieved:

         IF (ITER == 0) THEN
            WRITE (LUNTAB, 1010)
            WRITE (LUNCRT, 1010)
         ELSE
            BETTER1 = BETTER2
         END IF

         CALL CALCTHICK (NL, NU, XL, XU, YL, YU, BETTER2,
     >                   XATMAXTHICKNESS,
C****>                   A (1, 1), A (1, 2), A (1, 3), A (1, 4))
     >                   A, B, C, D)

         WRITE (LUNTAB, 1040) ' Itn.:', ITER,
     >      'Current % thickness:', BETTER2,
     >      'Corresponding abscissa:', XATMAXTHICKNESS
         WRITE (LUNCRT, 1040) ' Itn.:', ITER,
     >      'Current % thickness:', BETTER2,
     >      'Corresponding x/c:', XATMAXTHICKNESS

C ...    Check for refining this new thickness:

         ITER = ITER + 1
         DELTATHICKNESS = BETTER2 - DESIREDTHICKNESS
         IF (ABS (DELTATHICKNESS) > TOLERANCE) THEN

            IF (ITER < MAXITER) THEN

               IF (ITER == 1) THEN

C ...             Generate second calibration point prior to
C                 linear interpolation:

                  TARGET1 = ORIGINALTHICKNESS
                  TARGET2 = DESIREDTHICKNESS - DELTATHICKNESS
                  TARGET3 = TARGET2

               ELSE IF (ITER >= 3) THEN

C ...             Do a normal iteration, whereas passes 1, 2 were special:

                  TARGET1 = TARGET2
                  TARGET2 = TARGET3
               END IF

               IF (ITER >= 2) THEN
                  CALL LINTRP (1, TARGET1, TARGET2, BETTER1, BETTER2,
     >                         DESIREDTHICKNESS, TARGET3)
               END IF

               CALL COPY (NL, YLORIG, YL)
               CALL COPY (NU, YUORIG, YU)
               GO TO 200

            ELSE
               WRITE (LUNCRT, 1003)
            END IF

         END IF

      GO TO 999

  920 WRITE (LUNCRT, 1030)
     >   ' REFINE: Error reading y" table - aborting. ',
     >   ' IER from RDXYZ: ', IER

  990 STOP

  999 RETURN

C ... Formats:

 1001 FORMAT
     >  (/' REFINE works with  2N  equations in  N  unknowns (the ys):',
     >  //'        N of the form         y (new) = scale * y (old)',
     >   /'    and N of the form   wt * y" (new) = wt * y" (desired)',
     >  //' where  wt  represents weighting of the  second  derivative',
     >   /' equations to equilibrate the two halves of the system.',
     >  //' The first half enables thickness to be changed,  while the',
     >   /' second half enables the second derivatives  (and hence the',
     >   /' curvature distribution) to be smoothed.    The scaling and',
     >   /' the weighting use  "sine"  shape functions which  must  be',
     >   /' controlled by you the user  -  hence the series of prompts',
     >   /' to be discussed next.',//)
 1002 FORMAT
     > (//' Thickness ratio refinement:',
     >  //'    The nonlinear y scaling is intended to preserve as much',
     >   /'    as possible the curvature near the leading and trailing',
     >   /'    edges.    The sine function is centered at the point of',
     >   /'    maximum thickness.    Higher powers of the sine tend to',
     >   /'    preserve the leading/trailing-edge properties better.',
     >  //' Curvature smoothing:',
     >  //'    Typical weighting of the y" equations varies from 0.004',
     >   /'    at the leading and trailing edges to 0.04 at the center',
     >   /'    of the region of interest  (where most of the smoothing',
     >   /'    is sought).  Use a power of 3. or 4. to smooth out some',
     >   /'    NARROW region of noisy curvature,  else  a lesser power',
     >   /'    (1., 1.5, or 2.) if BROAD smoothing is sought (probably',
     >   /'    in conjunction with a table of target 2nd derivatives).',
     >   /)
 1003 FORMAT
     >  (/' REFINE: The T/C refinement iteration has not converged.',
     >   /' You are probably asking for too much.',
     >   /' Either try again, seeking less, or REFINE this result.')

 1005 FORMAT (/, A)
 1010 FORMAT (A, A)
 1020 FORMAT (A, F10.5)
 1030 FORMAT (/, 2A, I5, /)
 1040 FORMAT (A, I3, 4X, A, F10.5, 4X, A, F10.5)

      END SUBROUTINE REFINE
C+----------------------------------------------------------------------
C
      SUBROUTINE REFSRF (ISRF, NSRF, X, Y, YSCALE,
     >                   YBUMPPARAMS, YPPWEIGHT,
C****>                   A, C,
     >                   A, B, C, D, R, S,
     >                   NOYPPTABLE, NPTSYPPTABLE,
     >                   XYPPTABLE, YPPTABLE,
     >                   DETAILS, PROMPT, LEAVEALONE,
     >                   LUNCRT, LUNKBD, LUNTAB)
C
C  PURPOSE:  REFSRF was introduced so that REFINE can treat both surf-
C            aces without repeated code.  See REFINE for more details.
C
C  ARGUMENTS:
C  VAR   DIM   TYPE I/O/S DESCRIPTION
C  ISRF   -     I     I   ISRF=1 means upper surface; 2 means lower
C  NSRF   -     I     I   Number of points on this surface
C  X     NSRF   R     I   Surface abscissas
C  Y     NSRF   R    I/O  Surface ordinates (probably modified on return)
C YSCALE NSRF   R     S   For nonuniform scaling of the ordinates
C YBUMPPARAMS   3  R  I   Parameters for nonuniform y scaling bump
C YPPWEIGHT   NSRF R  S   For nonuniform weighting of y" equations
C**A   2*N,N+1  R     S   For setting up overdetermined system (N=NSRF-2)
C**C     2*N    R     S   For solution of overdetermined system
C  A      N     R     S   For sparse solution of overdetermined system
C  B     N+1    R     S    "    "    "    "    "    "    "    "    "
C  C      N     R     S    "    "    "    "    "    "    "    "    "
C  D      N     R     S    "    "    "    "    "    "    "    "    "
C  R      N     R     S    "    "    "    "    "    "    "    "    "
C  S      N     R     S    "    "    "    "    "    "    "    "    "
C  NOYPPTABLE   L     I   T means no y" table look-up for this surface
C  NPTSYPPTABLE I     I   Obvious input
C  XYPPTABLE, YPPTABLE    Abscissas and ordinates of y" table if present
C  DETAILS -    L     I   Controls detailed writes to LUNTAB
C  PROMPT  -    L    I/O  Since REFSRF is called for either surface and
C                         doesn't know about the outside iteration,
C                         control of the prompts common to both surfaces
C                         is a little awkward.  PROMPT needs to be TRUE
C                         for each surface for the first iteration, and
C                         turned off thereafter.  It also controls some
C                         of the detailed printing to LUNTAB (needed
C                         once for each surface, but not each iteration).
C  LEAVEALONE  L    I/O   Similarly, this logical must be FALSE for the
C                         first time in for each surface, but may be set
C                         TRUE for one of the surfaces at the prompt here.
C  LUNCRT -     I     I   Logical unit for screen
C  LUNKBD -     I     I   Logical unit for keyboard
C  LUNTAB -     I     I   Logical unit for tabulations
C
C  PROCEDURES:
C  BEVAL       For the bumps needed by the y scaling and y" weighting
C  DTDLSQ      Alternative to HDESOL, for special structure involved
C**HDESOL      Linear least squares of dense system by Householder
C**            transformations
C  READER      Prompting utility
C  TABLE1      1-D table look-up (linear interpolation)
C
C  HISTORY:
C  07/30/93   DAS   Introduced global SAVE after all, to be certain.
C  02/07/90   DAS   Removed concatenation from I/O; replaced ('$...')
C                   with (A,$) for IRIS 4D purposes.
C  03/19/86   DAS   Eliminated the END DO VAX dependency for PC compilers
C                   and arranged for suppression of details in .TAB file.
C                   Also avoided problems with SAVE and DATA statements
C                   by making DETAILS, PROMPT, and LEAVEALONE arguments.
C  11/05/84   DAS   Introduced DTDLSQ to avoid large 2-D array.
C                   Left dense-method solution in, but commented out.
C  10/25/84   DAS   Prompts tidied up in view of new on-line help.
C  04/11/84   DAS   Eliminated handling of non-normalized airfoil; put
C                   subroutine BEVAL in place of function BUMP.
C  11/26/83   DAS   Simplified set-up of 2nd derivative constraints  -
C                   prompted by difficulty when I=1.
C  10/28/83   DAS   Revised for integration with REFINE in PROFILE, as
C                   a permanent new capability.
C  10/26/83   DAS   First shot at doing table look-up for  target y"s.
C  10/20/83   DAS   Suppressed prompts for case where a surface is re-
C                   quired to be left unperturbed.
C  10/12/83   DAS   y" weightings and mid-section constant constraints
C                   may be different on each surface now.
C  10/07/83   DAS   Added ISRF to permit different things done to each
C                   surface; enabled bounding of y" in some x/c region
C                   on the upper surface only.
C  09/30/83   DAS   Added saving of plottable curvature data.
C  09/27/83   DAS   Modularized as REFSRF, for treating both surfaces.
C  09/26/83   DAS   Nonuniform y-scaling (modified sine bump, suggest-
C                   by Rob Kennelly).  Both surfaces treated same way.
C  09/23/83   DAS   y"  before and after from finite differencing now,
C                   for greater consistency.
C  09/20/83   DAS   Initial design/coding of special version of MODIFY
C                   (y" from spline; upper surface only).
C
C  AUTHOR: David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

      SAVE

C ... Arguments:

      INTEGER
     >   ISRF, LUNCRT, LUNKBD, LUNTAB, NPTSYPPTABLE, NSRF

C*****REAL A (2*(NSRF-2), NSRF-1), C (2*(NSRF-2)),

      REAL
     >   A (NSRF-2), B (NSRF-1), C (NSRF-2), D (NSRF-2), R (NSRF-2),
     >   S (NSRF-2), YSCALE (NSRF), X (NSRF), Y (NSRF),
     >   YPPWEIGHT (NSRF), XYPPTABLE (*), YPPTABLE (*), YBUMPPARAMS (3)

      LOGICAL
     >   DETAILS, LEAVEALONE, NOYPPTABLE, PROMPT

C ... Local variables:

      INTEGER
     >   I, IER, INDEX, IOS, J, N

      REAL
     >   INNERWEIGHT (2), OUTERWEIGHT (2), SBUMPPARAMS (3),
     >   XATMAXWEIGHT (2), XHI (2), XLO (2), YPPBUMPWIDTH (2),
     >   YPPCONST (2), CI, CIM1, CIP1, DXI, DXIM1, DXSUM, SSQMIN,
     >   WEIGHT, YPPINTERP

      LOGICAL
     >   DEFAULT, NOYPPCONSTRAINT (2), QUIT

      CHARACTER
     >   SRF (2) * 15

C ... Procedures:

      REAL
     >   TABLE1

      EXTERNAL
     >   BEVAL, DTDLSQ, READR, READY, TABLE1

C ... Storage:

      DATA
     >   SRF / ' upper surface ', ' lower surface ' /


C ... Execution:

      IF (LEAVEALONE) GO TO 999

      IF (PROMPT) THEN

         IF (NOYPPTABLE) THEN
            WRITE (LUNCRT, 1010)
     >         ' No target y" values have been read for the',
     >         SRF (ISRF), '-'
            CALL READY (LUNCRT, 'leave the' // SRF (ISRF) //
     >         'unchanged? (Y/N; <CR>=No): ',
     >         LUNKBD, LEAVEALONE, DEFAULT, QUIT)
         END IF

C ...    <Default is FALSE from calling program.>

         IF (LEAVEALONE) THEN
            WRITE (LUNTAB, 1010) ' The', SRF (ISRF), 'is unchanged.'
            PROMPT = ISRF == 1
            GO TO 999
         ELSE
            WRITE (LUNTAB, 1010) ' Controls for', SRF (ISRF), 'follow:'
            WRITE (LUNCRT, 1010) '    Smoothing of', SRF (ISRF),
     >                           'is controlled by the following:'

            IF (NOYPPTABLE) THEN
               CALL READR (LUNCRT,
     > '   Enter y" value sought over some range. <CR>=no constraint: ',
     >            LUNKBD, YPPCONST (ISRF), NOYPPCONSTRAINT (ISRF), QUIT)

               IF (.NOT. NOYPPCONSTRAINT (ISRF)) THEN
C ...             Permit simple flattening of curvature in some interval:

  200             WRITE (LUNCRT, '(A)', ADVANCE = 'NO')
     >              '    Corresp. x/c-range? (2 values):               '
                  READ  (LUNKBD, *, IOSTAT = IOS) XLO (ISRF), XHI (ISRF)
                  IF (IOS /= 0) THEN
                     WRITE (LUNCRT, '(A)')
                     GO TO 200
                  END IF
                  WRITE (LUNTAB, 1020) 'y" constraint and x range:',
     >               YPPCONST (ISRF), XLO (ISRF), XHI (ISRF)
               END IF
            END IF

            XATMAXWEIGHT (ISRF) = 0.5E+0
            CALL READR (LUNCRT,
     >         '   Center of smoothed region? <CR> gives x/c=0.5: ',
     >         LUNKBD, XATMAXWEIGHT (ISRF), DEFAULT, QUIT)

            YPPBUMPWIDTH (ISRF) = 3.E+0
            CALL READR (LUNCRT,
     >         '   Power of sine for y" weights?  <CR> gives 3.0: ',
     >         LUNKBD, YPPBUMPWIDTH (ISRF), DEFAULT, QUIT)

            WRITE (LUNTAB, 1020) ' x/c at peak y" weight    :',
     >                           XATMAXWEIGHT (ISRF)
            WRITE (LUNTAB, 1020) ' Sine power for y" weights:',
     >                           YPPBUMPWIDTH (ISRF)

            OUTERWEIGHT (ISRF) = 0.004E+0
            CALL READR (LUNCRT,
     >         '   Weight at x/c=0,1 for y"? <CR> gives .004: ',
     >         LUNKBD, OUTERWEIGHT (ISRF), DEFAULT, QUIT)
            WRITE (LUNTAB, 1020) ' y" weight at x/c=0 and 1 :',
     >                           OUTERWEIGHT (ISRF)

            INNERWEIGHT (ISRF) = 0.04E+0
            CALL READR (LUNCRT,
     >         '   Peak weight for y"?       <CR> gives .040: ',
     >         LUNKBD, INNERWEIGHT (ISRF), DEFAULT, QUIT)
            WRITE (LUNTAB, 1020) ' Peak y" weight           :',
     >                           INNERWEIGHT (ISRF)
            WRITE (LUNCRT, 1005)
         END IF
      END IF

C ... Set up scale factors for the ordinates (nonlinear thinning):

      DO I = 1, NSRF
         YSCALE (I) = 1.E+0
      END DO

      CALL BEVAL ('SINE', 3, YBUMPPARAMS, .TRUE., NSRF, X, YSCALE)

      IF (DETAILS) THEN
         WRITE (LUNTAB, 1005)
     >      ' ', ' y scale factors sought at each x:', ' '
         WRITE (LUNTAB, 1030) (X (I), YSCALE (I), I = 1, NSRF)
      END IF

C ... Generate the y" weights for this surface:

      SBUMPPARAMS (1) = XATMAXWEIGHT (ISRF)
      SBUMPPARAMS (2) = YPPBUMPWIDTH (ISRF)
      SBUMPPARAMS (3) = INNERWEIGHT (ISRF) - OUTERWEIGHT (ISRF)

      DO I = 1, NSRF
         YPPWEIGHT (I) = OUTERWEIGHT (ISRF)
      END DO

      CALL BEVAL ('SINE', 3, SBUMPPARAMS, .TRUE., NSRF, X, YPPWEIGHT)

      IF (PROMPT) THEN
         PROMPT = ISRF == 1

         IF (DETAILS) THEN
            WRITE (LUNTAB, 1010) ' Weights used for y" at each',
     >                           SRF (ISRF), 'x/c:', ' '
            WRITE (LUNTAB, 1030) (X (I), YPPWEIGHT (I), I = 1, NSRF)
         END IF
      END IF

C ... Set up overdetermined system.
C     The original solution as a dense system is left in comment form
C     as an aid to understanding the solution by a specialized least
C     squares algorithm.
C     First, zero out the system (treated as dense):
C
C     N = NSRF - 2
C
C     DO J = 1, N+1
C        DO I = 1, 2*N
C           A (I, J) = 0.E+0
C        END DO
C     END DO
C
C ... Next, the simple scaling part:
C
C     DO I = 1, N
C        A (I, I)   = 1.E+0
C        A (I, N+1) = Y (I+1)*YSCALE (I+1)
C     END DO
C
C ... Next, the 2nd derivative "constraints":
C
C     DO I = 1, N
C        DXI   = X (I+2) - X (I+1)
C        DXIM1 = X (I+1) - X (I)
C        DXSUM = X (I+2) - X (I)
C        WEIGHT= YPPWEIGHT (I+1)*2.E+0
C        CIM1  = WEIGHT / (DXIM1*DXSUM)
C        CI    =-WEIGHT / (DXIM1*DXI)
C        CIP1  = WEIGHT / (DXI*DXSUM)
C        IF (I>1) A (I+N, I-1) = CIM1
C        A (I+N, I)   = CI
C        A (I+N, I+1) = CIP1
C        A (I+N, N+1) = CIM1*Y (I) + CI*Y (I+1) + CIP1*Y (I+2)
C
C        IF (NOYPPTABLE) THEN
C           IF (.NOT. NOYPPCONSTRAINT (ISRF)) THEN
C ...          Flatten the curvature plot in given x interval:
C
C              IF (X (I+1)>XLO (ISRF) .AND. X (I+1)<XHI (ISRF))
C    >            A (I+N, N+1) = 0.5E+0*WEIGHT*YPPCONST (ISRF)
C           ELSE
C ...          Don't constrain the y" target - shoot for original.
C              CONTINUE
C           END IF
C        ELSE
C           INDEX = 1
C           YPPINTERP = TABLE1 (NPTSYPPTABLE, XYPPTABLE,
C    >                          YPPTABLE, INDEX, X (I+1), IER)
C           IF (IER==0) THEN
C              A (I+N, N+1) = 0.5E+0*WEIGHT*YPPINTERP
C           ELSE IF (IER<=4) THEN
C              WRITE (LUNCRT, *) 'Table look-up error. IER:', IER
C              WRITE (LUNCRT, *) 'NTABLE:', NPTSYPPTABLE,
C    >                           '  ISRF:', ISRF
C              WRITE (LUNCRT, *) 'Abscissa:', X (I+1)
C              STOP
C           ELSE
C ...          Abscissa was out of table range - use original y".
C              CONTINUE
C           END IF
C        END IF
C     END DO
C
C ... Adjust RHS for known leading and trailing edge values:
C
C     Y (1)    = YSCALE (1)*Y (1)
C     Y (NSRF) = YSCALE (NSRF)*Y (NSRF)
C     A (1+N, N+1) = A (1+N, N+1) - 2.E+0*YPPWEIGHT (2)*Y (1) /
C    >                        ( (X (2)-X (1))* (X (3)-X (1)))
C     A (N+N, N+1) = A (N+N, N+1) - 2.E+0*YPPWEIGHT (N)*Y (NSRF) /
C    >                                               (DXI*DXSUM)
C
C     WRITE (LUNTAB, *)
C     WRITE (LUNTAB, *)
C    >   'RHS vector, of length 2N, where N = NSRF-2 =', N, ':'
C     WRITE (LUNTAB, *)
C     WRITE (LUNTAB, ' (1X, 1P, 10E13.4)') (A (I, N+1), I = 1, N)
C     WRITE (LUNTAB, ' (1X, 1P, 10E13.4)') (A (I, N+1), I = N+1, 2*N)
C     WRITE (LUNTAB, *)
C
C ... Solve the system:
C
C     CALL HDESOL (2*N, 2*N, N+1, A, C, SSQMIN)
C
C ... Extract the solution:
C
C     DO I = 1, N
C        Y (I+1) = C (I)
C     END DO


C ... Set up overdetermined system.
C     DTDLSQ was written specially for this diagonal+tridiagonal structure.

      N = NSRF - 2

C ... First, the simple scaling part:

      DO I = 1, N
         D (I) = 1.E+0
         R (I) = Y (I+1) * YSCALE (I+1)
      END DO

C ... Next, the 2nd derivative "constraints":

      DO I = 1, N

         DXI   = X (I+2) - X (I+1)
         DXIM1 = X (I+1) - X (I)
         DXSUM = X (I+2) - X (I)
         WEIGHT= YPPWEIGHT (I+1) * 2.E+0
         CIM1  = WEIGHT / (DXIM1*DXSUM)
         CI    =-WEIGHT / (DXIM1*DXI)
         CIP1  = WEIGHT / (DXI*DXSUM)
         C (I)  = CIM1
         A (I)  = CI
         B (I)  = CIP1
         S (I)  = CIM1 * Y (I) + CI * Y (I+1) + CIP1 * Y (I+2)

         IF (NOYPPTABLE) THEN
            IF (.NOT. NOYPPCONSTRAINT (ISRF)) THEN
C ...          Flatten the curvature plot in given x interval:

               IF (X (I+1)>XLO (ISRF) .AND. X (I+1)<XHI (ISRF))
     >            S (I) = 0.5E+0 * WEIGHT * YPPCONST (ISRF)
            ELSE
C ...          Don't constrain the y" target - shoot for original.
               CONTINUE
            END IF
         ELSE
            INDEX = 1
            YPPINTERP = TABLE1 (NPTSYPPTABLE, XYPPTABLE,
     >                          YPPTABLE, INDEX, X (I+1), IER)
            IF (IER == 0) THEN
               S (I) = 0.5E+0 * WEIGHT * YPPINTERP
            ELSE IF (IER <= 4) THEN
               WRITE (LUNCRT, 1005)
               WRITE (LUNCRT, 1015) ' Table look-up error. IER:', IER,
     >                              ' NTABLE:', NPTSYPPTABLE,
     >                              ' ISRF:', ISRF
               WRITE (LUNCRT, 1020) ' Abscissa:', X (I+1)
               STOP
            ELSE
C ...          Abscissa was out of table range - use original y".
               CONTINUE
            END IF
         END IF

      END DO

C ... Adjust RHS for known leading and trailing edge values:

      Y (1)    = YSCALE (1) * Y (1)
      Y (NSRF) = YSCALE (NSRF) * Y (NSRF)
      S (1) = S (1) - C (1) * Y (1)
      S (N) = S (N) - B (N) * Y (NSRF)

      IF (DETAILS) THEN
         WRITE (LUNTAB, 1005) ' ', ' RHS vector, length 2*(N-2):', ' '
         WRITE (LUNTAB, 1040) R (1 : N)
         WRITE (LUNTAB, 1040) S (1 : N)
      END IF

C ... Solve the system:

      CALL DTDLSQ (N, A, B, C, D, R, S, SSQMIN)

C ... Extract the solution:

      DO I = 1, N
         Y (I+1) = R (I)
      END DO

      WRITE (LUNTAB, 1025) ' Sum of squares obtained  :', SSQMIN
C
  999 RETURN

C ... Formats:

 1005 FORMAT (A)
 1010 FORMAT (/, A, A, A)
 1015 FORMAT (A, I5)
 1020 FORMAT (A, 3F15.6)
 1025 FORMAT (/, A, 1P, E15.6)
 1030 FORMAT (5 (1X, F12.4, F10.5))
 1040 FORMAT (1X, 1P, 10E13.4)

      END SUBROUTINE REFSRF
C+---------------------------------------------------------------------
C
      SUBROUTINE ROTATE (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD)
C
C  PURPOSE: ROTATE applies twist to one profile about a variable center
C           of rotation.  The result may be optionally renormalized.
C
C  ARGUMENTS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   NU       -      I     I     Upper surface info.
C   XU       NU     R    I/O      "     "     "
C   YU       NU     R    I/O      "     "     "
C   NL       -      I     I     Lower surface info.
C   XL       NL     R     I       "     "     "
C   YL      NPTS    R    I/O      "     "     "
C  LUNCRT    -      I     I     Logical unit for prompting.
C  LUNKBD    -      I     I     Logical unit for responses.
C
C  PROCEDURES:
C    NRMLIZ    Normalizes X or Y for one surface of an airfoil
C    RDREALS   Prompts for more than one value
C    READER    Prompting utility
C    ROTATE2D  Rotates point(s) (X,Y) about (P,Q)
C
C  AUTHOR:     David Saunders, Sterling Software, Palo Alto, CA.
C
C  HISTORY:
C  06/29/90   DAS   Initial implementation.
C  10/21/99    "    Positive rotation should be clockwise (like AOA).
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, NL, NU

      REAL
     >   XL(NL), XU(NU), YL(NL), YU(NU)

C  *  Local constants:

      REAL, PARAMETER ::
     >   FOURTH = 0.25E+0, HALF = 0.5E+0, ONE = 1.E+0, ZERO = 0.E+0

C  *  Local variables:

      INTEGER
     >   CHOICE, NVALS

      REAL
     >   CHORD, TWIST, VALS (2), XTE, YTE

      LOGICAL
     >   DEFAULT, QUIT, RENORM

C  *  Procedures:

      EXTERNAL
     >   NRMLIZ, RDREALS, READI, READR, READY, ROTATE2D

C  *  Execution:

      WRITE (LUNCRT, 1001)
     >   ' Options for applying twist:',
     >   '    1 = Rotate about leading edge',
     >   '    2 = Rotate about trailing edge',
     >   '    3 = Rotate about quarter-chord (on center-line)',
     >   '    4 = Rotate about specified point'

      CHOICE = 1
      CALL READI (LUNCRT, 'Enter selection (<CR> = 1): ',
     >            LUNKBD, CHOICE, DEFAULT, QUIT)

      CHORD = MAX (XU (NU), XL (NL)) - XU (1)
      XTE   = (XU (NU) + XL (NL)) * HALF
      YTE   = (YU (NU) + YL (NL)) * HALF

      IF (CHOICE == 1) THEN
         VALS (1) = XU (1)
         VALS (2) = YU (1)
      ELSE IF (CHOICE == 2) THEN
         VALS (1) = XTE
         VALS (2) = YTE
      ELSE IF (CHOICE == 3) THEN
         VALS (1) = XU (1) + (XTE - XU (1)) * FOURTH
         VALS (2) = YU (1) + (YTE - YU (1)) * FOURTH
      ELSE IF (CHOICE == 4) THEN
  200    NVALS = 2
         CALL RDREALS (LUNCRT, '$Enter center of rotation (2 values): ',
     >      LUNKBD, NVALS, VALS)
         IF (NVALS /= 2) GO TO 200
      END IF

  300 CALL READR (LUNCRT,
     >   'Enter twist (degrees; positive is clockwise): ',
     >   LUNKBD, TWIST, DEFAULT, QUIT)
      IF (DEFAULT) GO TO 300

C  *  Provide for renormalizing (but in general this does not recover
C     the original abscissas):

      IF (XU (1) == ZERO .AND. YU (1) == ZERO .AND. CHORD == ONE) THEN
         RENORM = .TRUE.
         CALL READY (LUNCRT,
     >      'Do you want the result renormalized?  (Y/N; <CR> = Yes): ',
     >       LUNKBD, RENORM, DEFAULT, QUIT)
      ELSE
         RENORM = .FALSE.
      END IF

C  *  Perform the rotation, in-place:

      TWIST = -TWIST ! ROTATE2D thinks in quadrants

      CALL ROTATE2D (NU, XU, YU, TWIST, VALS (1), VALS (2))
      CALL ROTATE2D (NL, XL, YL, TWIST, VALS (1), VALS (2))

C  *  Warn user if leading edge point has been displaced:

      IF (XU (2) < XU (1) .OR. XL (2) < XL (1)) THEN
         WRITE (LUNCRT, 1001)
     >   ' *** WARNING: LE point is no longer the true leading edge.',
     >   '              PROFILE''s "rectify" option may be appropriate.'
         IF (RENORM) WRITE (LUNCRT, 1000)
     >   '              Normalization has been affected also.'
      END IF

      IF (RENORM) THEN
         CHORD = MAX (XU (NU), XL (NL)) - XU (1) ! Not if warning above applied.
         CALL NRMLIZ (NU, XU, XU, XU (1), CHORD)
         CALL NRMLIZ (NU, YU, YU, YU (1), CHORD)
         CALL NRMLIZ (NL, XL, XL, XL (1), CHORD)
         CALL NRMLIZ (NL, YL, YL, YL (1), CHORD)
         WRITE (LUNCRT, 1001)
     >    ' *** Note that the original Xs are NOT retained in general.',
     >    '     Rerun PROFILE in "redistribute" mode to recover them.'
      END IF

 1000 FORMAT (A)
 1001 FORMAT (/, (A))

      END SUBROUTINE ROTATE
C+----------------------------------------------------------------------
C
      SUBROUTINE SMOOTH (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD, LUNOUT)
C
C  PURPOSE:  SMOOTH drives smoothing of an airfoil (either surface or
C            both surfaces) by a choice of methods.  The Y coordinates
C            are smoothed in-place.
C
C  ARGUMENTS:
C     VAR   DIM   I/O/S   DESCRIPTION
C     NU     -    I       # points on upper surface
C     XU     NU   I       Upper surface coordinates, leading to trailing
C     YU     NU   I/O     edge (not necessarily normalized)
C     NL     -    I       # points on lower surface
C     XL     NL   I       Lower surface coordinates as for upper
C     YL     NL   I/O
C     LUNCRT -    I       Logical units for screen, keyboard, and
C     LUNKBD -    I       printed coefficients
C     LUNOUT -    I
C
C  PROCEDURES:
C     WAGSMOOTH  Smooths an airfoil surface using Wagner functions
C     IMPSMOOTH  Smooths an airfoil surface via implicit + explicit methods
C
C  HISTORY:
C     08/11/86  DAS  Initial Wagner fn. form (from program SMOOTH).
C     12/20/96  DAS  Added implicit/explicit smoothing option.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C-------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER  LUNCRT, LUNKBD, LUNOUT, NU, NL
      REAL     XU (NU), XL (NL), YU (NU), YL (NL)

C     Local variables:

      INTEGER  METHOD
      LOGICAL  CR, EOF

C     Procedures:

      EXTERNAL IMPSMOOTH, READI, WAGSMOOTH

C     Execution:

      WRITE (LUNCRT, '(/, (1X, A))')
     >   'Smoothing methods:',
     >   '   1:  Least squares fitting of Wagner functions 1:N',
     >   '   2:  Implicit + explicit smoothing, variable along chord'

  10  METHOD = 1
      CALL READI (LUNCRT, 'Pick one. Default = 1: ',
     >            LUNKBD, METHOD, CR, EOF)
      IF (EOF) GO TO 99

      IF (METHOD == 1) THEN  ! Wagner function smoothing

         CALL WAGSMOOTH (NU, XU, YU, 1, LUNCRT, LUNKBD, LUNOUT)
         CALL WAGSMOOTH (NL, XL, YL, 2, LUNCRT, LUNKBD, LUNOUT)

      ELSE IF (METHOD == 2) THEN  ! Implicit (+ explicit) smoothing

         CALL IMPSMOOTH (NU, XU, YU, 1, LUNCRT, LUNKBD, LUNOUT)
         CALL IMPSMOOTH (NL, XL, YL, 2, LUNCRT, LUNKBD, LUNOUT)

      ELSE
         GO TO 10
      END IF

   99 RETURN

      END SUBROUTINE SMOOTH
C+---------------------------------------------------------------------
C
      SUBROUTINE TRANSFORM (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD,
     >                      ICASE)
C
C  PURPOSE: TRANSFORM  drives transformation of an airfoil from  upper/
C           lower surface representation to thickness/camber represent-
C           ation, or vice versa.   It removes from the calling program
C           such necessities as prompting for which way to go, and dis-
C           allowing it if upper/lower abscissas are not the same, etc.
C           This version provides for zeroing the camber of a section.
C
C  METHOD:  See XFORM for the essentials - separated out in case  it is
C           reusable elsewhere.  Note that ICASE is an OUTput, and that
C           the transformation is done in-place.
C
C           Checking for dissimilar abscissas might save somebody  some
C           grief - quitting here is better than not bothering to check.
C           
C  ARGUMENTS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   NU       -      I     I     ICASE = 0 or 1: no. of upper surface pts.
C   XU       NU     R     I     and corresp. abscissas, else same as NL,XL.
C   YU       NU     R    I/O    ICASE = 0: upper Ys in/decambered Ys out;
C                               ICASE = 1: upper Ys in/camber out;
C                               ICASE = 2: camber in/upper Ys out.
C   NL       -      I     I     ICASE = 2: no. of values in camber/thickness
C   XL       NL     R     I                distribns. and corresp. abscissas
C                                          else # lower surface pts. and Xs.
C   YL      NPTS    R    I/O    ICASE = 0: lower Ys in/decambered Ys out;
C                               ICASE = 1: lower Ys in/semi-thickness out;
C                               ICASE = 2: semi-thickness in; lower Ys out.
C  LUNCRT    -      I     I     Logical unit for prompting.
C  LUNKBD    -      I     I     Logical unit for responses.
C  ICASE     -      I     O     Lets user say which way to go, and tells the
C                               calling program for tabulation purposes, etc.
C                               ICASE = 0 means zero out the camber;
C                                     = 1 means Y-coordinates in, and
C                                         camber/semi-thickness out;
C                                     = 2 means transform the other way.
C
C  PROCEDURES:
C    READER   Prompting utility.
C    XFORM    Does the actual transformation (in same source module as
C             TRANSFORM because the two go hand in hand, but XFORM may
C             be useful on its own some day).
C
C  HISTORY:
C  02/12/85   DAS   Initial design and code.
C  08/31/87   DAS   'A' format instead of list-directed.
C  11/21/95   DAS   Installed a decambering option.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C----------------------------------------------------------------------
C
      IMPLICIT NONE

C  *  Arguments:

      INTEGER  LUNCRT, LUNKBD, NL, NU, ICASE
      REAL     XL (NL), XU (NU), YL (NL), YU (NU)

C  *  Local variables:

      INTEGER  I
      REAL     YLE
      LOGICAL  DEFAULT, QUIT, SAMEXS, YIN

C  *  Execution:

      WRITE (LUNCRT, 1001)
     >   ' Options for TRANSFORM mode:',
     >   '  0 = Decamber upper/lower surface coordinates',
     >   '  1 = Upper/lower surfaces to camber/thickness distributions',
     >   '  2 = Camber/thickness distributions to upper/lower surfaces',
     >   ' '

      ICASE = 1
      CALL READI (LUNCRT, 'Enter selection (<CR> = 1): ',
     >            LUNKBD, ICASE, DEFAULT, QUIT)

      IF (NU /= NL) THEN
         WRITE (LUNCRT, 1001)
     >      ' Cannot transform - different-sized distributions.',
     >      ' Stopping here.'
         GO TO 900
      END IF

      IF (ICASE /= 2) THEN  ! Allow round-off X differences for ICASE = 2
         SAMEXS = .TRUE.
         DO I = 1, NU
            IF (XU (I) /= XL (I)) SAMEXS = .FALSE.
         END DO

         IF (.NOT. SAMEXS) THEN
            WRITE (LUNCRT, 1001)
     >         ' Cannot transform - different abscissa distributions.',
     >         ' Use REDISTRIBUTE mode and try again.  Stopping here.'
            GO TO 900
         END IF

C        Ensure zero camber at the nose:

         YLE = YU (1)
         IF (YLE /= 0.) THEN
            DO I = 1, NU
               YU (I) = YU (I) - YLE
               YL (I) = YL (I) - YLE
            END DO
         END IF
      END IF

      YIN = ICASE /= 2

      CALL XFORM (YIN, NU, YU, YL)

      IF (ICASE == 0) THEN  ! Zero the camber ...
         DO I = 1, NU
            YU (I) = 0.
         END DO

         CALL XFORM (.FALSE., NU, YU, YL)  ! ... and reverse the transformation
      END IF

      GO TO 999

  900 STOP

  999 RETURN

 1001 FORMAT (/, (A))

      END SUBROUTINE TRANSFORM
C+---------------------------------------------------------------------
C
      SUBROUTINE XFORM (YIN, NPTS, YU, YL)
C
C  PURPOSE: XFORM  transforms one  airfoil  from  upper/lower  surface
C           representation to thickness/camber representation, or vice
C           versa.    Actually, semi-thickness is used, not thickness,
C           while "camber" is more accurately the mean line here:
C
C                 C = (YU + YL)/2               YU = C + T
C                 T = (YU - YL)/2               YL = C - T
C
C  METHOD:  The required output overwrites the input arrays.  Note the
C           assumption of equal numbers of points on both surfaces and
C           matching abscissas.
C
C  ARGUMENTS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   YIN      -      L     I     YIN=T means y-coordinates in; camber/semi-
C                                     thickness out;
C                               YIN=F means transform the other way.
C   NPTS     -      I     I     No. of points in each of the distributions
C                               involved.
C   YU      NPTS    R    I/O    If YIN=T: Upper surface in; camber out.
C                               If YIN=F: Camber in; upper surface out.
C   YL      NPTS    R    I/O    If YIN=T: Lower surface in; semi-thickness out.
C                               If YIN=F: Semi-thickness in; lower surface out.
C
C  HISTORY:
C  02/12/85   DAS   Initial design and code.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      LOGICAL  YIN
      INTEGER  NPTS
      REAL     YU (NPTS), YL (NPTS)

C  *  Local variables:

      INTEGER  I
      REAL     FACTOR, TEMP

C  *  Execution:

      IF (YIN) THEN
         FACTOR = 0.5E+0
      ELSE
         FACTOR = 1.E+0
      END IF

      DO I = 1, NPTS
         TEMP   = (YU (I) + YL (I)) * FACTOR
         YL (I) = (YU (I) - YL (I)) * FACTOR
         YU (I) = TEMP
      END DO

      END SUBROUTINE XFORM
C+-----------------------------------------------------------------------
C
      SUBROUTINE WAGSMOOTH (NPTS, X, Y, ISRF, LUNCRT, LUNKBD, LUNOUT)
C
C  PURPOSE:  WAGSMOOTH smooths the given airfoil surface by fitting a
C            linear combination of the first N Wagner shape functions.
C            The airfoil need not be normalized.
C
C  METHOD:   Prompt for whether to smooth this surface. If so, ...
C
C         >  Prompt for number of Wagner functions to use.
C         >  Normalize coordinates: [0.,1.] for X; corresp. scaling for Y.
C         >  If YTE = Y(NPTS) is not zero, use BEVAL's 'RAMP' option to
C            subtract Y = YTE * X at every abscissa, in place.  The
C            Wagner functions will be fitted to the adjusted ordinates.
C         >  Perform the least squares fit.  (See WAGFIT.)
C         >  Evaluate the result, including adding back the ramp term
C            if necessary, all in-place.
C         >  Denormalize.
C            
C            Note that local work-space is allocated for WAGFIT -
C            PROFILE's main program does not have enough available,
C            and WAGSMOOTH is somewhat self-contained as a result.
C
C  ARGUMENTS:
C  VAR   DIM   I/O/S   DESCRIPTION
C  NPTS   -    I       Number of points on current surface
C  X     NPTS  I       Abscissas of current surface   (not necessarily
C  Y     NPTS  I/O     Ordinates of current surface        normalized)
C  ISRF   -    I       1 means uppper surface;
C                      2 means lower surface
C  LUNCRT -    I       Logical units for screen, keyboard, and
C  LUNKBD -    I       printed coefficients
C  LUNOUT -    I
C
C  PROCEDURES:
C  READER   Prompting utility
C  BEVAL    Evaluates numerous shape functions at given Xs
C  NRMLIZ   Normalizes/denormalizes coordinates
C  WAGFIT   Sets up and solves the linear least squares problem
C
C  HISTORY:
C  08/11/86   DAS   Initial implementation (originally installed in
C                   program SMOOTH, but PROFILE is the logical place;
C                   expect X in [0.,1.], because PROFILE can normalize).
C  10/22/88   DAS   Introduced normalizing/denormalizing here - too many
C                   steps otherwise for arbitrary coordinates.
C  12/13/96   DAS   Added option to replace points near the trailing edge.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C-------------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER,INTENT(IN):: 
     >   ISRF, LUNCRT, LUNKBD, LUNOUT, NPTS

       REAL,INTENT(IN OUT):: X(NPTS)
       REAL,INTENT(IN OUT):: Y(NPTS)

C  *  Local constants:

      INTEGER, PARAMETER ::
     >   MXPTS  = 200,   ! Arbitrary - more pts./surface OK if NWAG < MXWAG
     >   MXWAG  = 24,
     >   MXWORK = (MXWAG + 2) * MXPTS

      REAL, PARAMETER ::
     >   ONE    = 1.E+0,
     >   ZERO   = 0.E+0

      CHARACTER, PARAMETER ::
     >   RAMP * 4 = 'RAMP'

C  *  Local variables:

      INTEGER
     >   I, IER, J, NFIX, NLEFT, NWAG

      REAL
     >   B (2), C (MXWAG), CHORD, RMSDEV, WORK (MXWORK), XLE, XTE, YLE,
     >   YN(1), YTE

      LOGICAL 
     >   ADD, CR, EOF

      CHARACTER
     >   SURFCE (2) * 5

C  *  Procedures:

      EXTERNAL
     >   BEVAL, NRMLIZ, READI, WAGFIT

C  *  Storage:

      SAVE
     >   NWAG
      DATA
     >   SURFCE / 'UPPER', 'LOWER' /

C  *  Execution:

      IF (ISRF == 1) THEN
         NWAG = 10
         WRITE (LUNCRT, '(/, A)')
     >      ' Enter # Wagner functions to use on the UPPER surface.'
         CALL READI (LUNCRT,
     >      '(<CR> = 10; 0 or EOF (^D) = leave surface alone): ',
     >      LUNKBD, NWAG, CR, EOF)
      ELSE  ! ISRF = 2
         IF (NWAG == 0) THEN
            NWAG = 10
            CALL READI (LUNCRT,
     >         ' # Wagner functions to use on the LOWER surface? [10] ',
     >         LUNKBD, NWAG, CR, EOF)
         ELSE
            CALL READI (LUNCRT,
     >   '# Wagner fns. for the LOWER surface? (<CR> = same as upper) ',
     >         LUNKBD, NWAG, CR, EOF)
         END IF
      END IF

      IF (EOF .OR. NWAG <= 0) GO TO 899

      IF (NWAG > MXWAG) THEN
         NWAG = MXWAG
         WRITE (LUNCRT, '(A, I3)') ' Limiting N to maximum of', MXWAG
      END IF

  200 NFIX = 3
      CALL READI (LUNCRT,
     > '# pts. to fix (after smoothing) forward of the TE? (>=0; [3]) ',
     >   LUNKBD, NFIX, CR, EOF)
      IF (EOF) GO TO 899
      IF (NFIX < 0 .OR. NFIX > NPTS / 4) GO TO 200

C  *  Normalize coordinates, since Wagner functions are defined on [0,1]:

      XLE = X (1)
      YLE = Y (1)
      XTE = X (NPTS)
      YTE = Y (NPTS)
      CHORD = XTE - XLE

      CALL NRMLIZ (NPTS, X, X, XLE, CHORD)
      CALL NRMLIZ (NPTS, Y, Y, YLE, CHORD)

C  *  To preserve a thick trailing edge exactly, a "ramp" function
C     must be subtracted first then added back in:

      YN(1) = -Y (NPTS)

      IF (YN(1) /= ZERO) THEN

C  *     Subtract "ramp" function in-place:

         ADD = .TRUE.
         CALL BEVAL (RAMP, 1, YN, ADD, NPTS, X, Y)
      END IF
                              
C  *  Set up and solve the linear least squares problem:

      CALL WAGFIT (NPTS, X, Y, NWAG, MXWORK, WORK, C (1),
     >             RMSDEV, IER)
      IF (IER /= 0) GO TO 910

      WRITE (LUNOUT, 1120)
     >   SURFCE (ISRF), NWAG, NFIX, RMSDEV, (C (I), I = 1, NWAG)

C  *  Evaluate the fit (in-place):

      ADD = .FALSE.
      DO J = 1, NWAG
         B (1) = J
         B (2) = C (J)

         CALL BEVAL ('WAGNER', 2, B, ADD, NPTS, X, Y)
         ADD = .TRUE.
      END DO

      IF (YN(1) /= ZERO) THEN

C  *     Add back the "ramp" function:

         YN(1) = -YN(1)
         CALL BEVAL (RAMP, 1, YN, ADD, NPTS, X, Y)
      END IF
                              
C  *  Patch bad points near the trailing edge?

      IF (NFIX > 0) THEN

         DO I = 1, NPTS
            WORK (I) = X (I)
            WORK (I + NPTS) = Y (I)
         END DO

         NLEFT = NPTS - NFIX
         WORK (NLEFT) = WORK (NPTS)  ! Remove NFIX pts
         WORK (NPTS + NLEFT) = WORK (NPTS + NPTS)

         CALL LCSFIT (NLEFT, WORK (1), WORK (1 + NPTS), .TRUE., 'M',
     >                NFIX, X (NLEFT), Y (NLEFT), Y (NLEFT))
      END IF

C  *  Denormalize:

      CALL NRMLIZ (NPTS, X, X, XLE, -CHORD)
      CALL NRMLIZ (NPTS, Y, Y, YLE, -CHORD)

      X (1)    = XLE
      X (NPTS) = XTE
      Y (1)    = YLE
      Y (NPTS) = YTE

  899 RETURN

C  *  Error handling:
      
  910 WRITE (LUNCRT, 1110) IER
      STOP

C  *  Formats:

 1110 FORMAT (' WAGSMOOTH: Bad return from WAGFIT - abort. IER: ', I1,
     >        /)
 1120 FORMAT (//, ' Details of smoothing for ', A5, ' surface:',
     >        //, ' # Wagner functions:', I3,
     >        /,  ' # patched points forward of TE:', I3,
     >        /,  ' RMS Deviation:', 1P, E14.6,
     >        /,  ' Coefficients:',
     >        /,  (5E14.6))

      END SUBROUTINE WAGSMOOTH


C+---------------------------------------------------------------------
C
      SUBROUTINE BOUNDS ( M, N, MDIM, F, FMIN, FMAX )
C
C  PURPOSE:
C     BOUNDS determines the maximum and minimum values contained in the
C     input array F, which is treated as 2-dimensional.  This version
C     provides for finding extremes across multiple arrays (multiple
C     calls, one per array).
C
C  ARGUMENTS:
C     ARG   TYPE I/O/S DIM    DESCRIPTION
C     M,N    I     I    -     This routine checks elements F(1:M,1:N).
C                             N=1 for a 1-D array.
C     MDIM   I     I    -     Effective row dimension of F in the
C                             routine that sets up this array.
C     F      R     I  MDIM,N  Array to be scanned for min and max.
C     FMIN,  R    I/O   -     Minimum and maximum values found (so far).
C     FMAX
C            WARNING:  Do not initialize FMIN or FMAX with extreme
C                      values, since they may never be reset because
C                      the code takes advantage of the fact that a new
C                      maximum cannot also be a new minimum.
C
C            Sample usage (1-D case):
C
C            FMIN = F1(1)
C            FMAX = FMIN
C            CALL BOUNDS ( NPTS, 1, MXPTS, F1, FMIN, FMAX )
C            CALL BOUNDS ( NPTS, 1, MXPTS, F2, FMIN, FMAX )
C
C            Then FMIN, FMAX are the extreme values across 1-D arrays F1 & F2.
C
C  ENVIRONMENT:  VAX, CRAY -- FORTRAN 77
C
C  HISTORY:
C   01/04/82    PJT    Original design and coding.
C   11/23/82    RGL    Extended to handle 2D arrays as well as vectors.
C   02/11/87    DAS    Provided for usage across multiple arrays.
C   10/27/88    DAS    Clarified usage for the 1-D array case.  Revised
C                      code to take advantage of fact that a new maximum
C                      cannot also be a new minimum.  (But use of MIN/MAX
C                      intrinsics may be better on vector machines.)
C   06/12/89    DAS    Aaargh!  Minimum is never set if inputs are
C                      monotonically increasing.  Therefore, it must be
C                      input as a legitimate value and not some extreme.
C                      Chose to warn user in header rather than change
C                      the code, but is the efficiency worth the risk?
C
C  AUTHOR:   Jeff Trosin, Informatics Inc.
C
C----------------------------------------------------------------------

      INTEGER   M, N, MDIM
      REAL      F(MDIM,N), FMIN, FMAX
      INTEGER   I, J

      DO 60, J = 1, N
         DO 50, I = 1, M
            IF ( F(I,J) .GT. FMAX ) THEN
               FMAX = F(I,J)
            ELSE IF ( F(I,J) .LT. FMIN ) THEN
               FMIN = F(I,J)
            END IF
   50    CONTINUE
   60 CONTINUE

      RETURN
      END
*DECK CENDIF
C+----------------------------------------------------------------------
C
      SUBROUTINE CENDIF (N, X, FX, EPSOBJ, H, GRAD, DIAG, NFUNCT,
     >   LWRITE, FUNCT)
C
C
C     Description:
C
C           CENDIF calculates central difference approximations to the
C        gradient and Hessian of an optimization objective function.  If
C        requested, improved finite difference stepsizes are calculated.
C        CENDIF/FDSTEP/DIFFER were developed for use with QNMDIF (a
C        quasi-Newton optimizer), but the stepsizes may be used by any
C        forward-difference method.
C
C
C     Parameters:
C
C        Name   Dimension   Type  I/O/S   Description
C        N                   I    I       Number of optimization variables.
C        X       N           R    I       Optimization variables.
C        FX                  R    I       Initial function value.
C        EPSOBJ              R    I       Estimate of absolute error in
C                                         objective.
C        H       N           R    I/O     Initial stepsizes used for
C                                         differencing if >= 0.  If H(I)
C                                         is negative, indicates that an
C                                         improved value is to be estimated
C                                         by FDSTEP and returned.
C        GRAD    N           R      O     Central difference approximation to
C                                         gradient.
C        DIAG    N           R      O     Central difference approximation
C                                         to diagonal elements of Hessian.
C        NFUNCT              I      O     Number of function evaluations
C                                         (increment).
C        LWRITE              I      I     Logical unit number for output
C                                         diagnostics, provided >= 0.
C                                         Use a negative unit number to
C                                         suppress all output.
C        FUNCT                            Objective function subroutine
C                                         (EXTERNAL).
C
C     External References:
C
C        Name     Description
C        DIFFER   Computes first and second derivative approximations
C                 by central differences.
C        FDSTEP   Produces estimates of optimal one sided finite
C                 difference stepsizes.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Author:  Robert Kennelly, Sterling Software.
C
C
C     Development History:
C
C        16 July 1982    RAK    Complete revision of earlier version by
C                               Trosin - now use FDSTEP and DIFFER.
C         2 Sep. 1982    RAK    Reordered parameters, redefined NFUNCT.
C        19 Jan. 1983    RAK    Updated NFDEL following DIFFER.
C        10 Mar. 1983    RAK    Repaired call to DIFFER (... DIAG(I) ...),
C                               added IF-THEN-ELSE blocks.
C         9 Aug. 1986    RAK    Added LWRITE input parameter (had been
C                               in COMMON); reordered parameters again
C                               (sorry!).  If LWRITE is negative,
C                               suppress all output.  Use PARAMETER
C                               for local array dimensions, and test
C                               input N for consistency.  Increase
C                               MAXN from 30 to 40.  Add IMPLICIT NONE.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      INTEGER
     >   MAXN
      PARAMETER
     >   (MAXN = 40)

      REAL
     >   ZERO, HALF
      PARAMETER
     >   (ZERO = 0.0E+0,
     >    HALF = 0.5E+0)

C     Arguments.

      INTEGER
     >   LWRITE, N, NFUNCT
      REAL
     >   DIAG (N), FX, EPSOBJ, GRAD (N), H (N), X (N)

C     Local variables.

      INTEGER
     >   I, IERROR (MAXN), NFDEL (MAXN)
      REAL
     >   CFPBAC, CFPFOR, CFPP, ERRBND, FPBAC, FPFOR, HI, URDIAG,
     >   URH (MAXN)

C     Procedures.

      EXTERNAL
     >   DIFFER, FDSTEP, FUNCT


C     Execution.
C     ----------

      IF (LWRITE .GE. 0) WRITE (LWRITE, 1000)
      IF (N .GT. MAXN) STOP 'CENDIF:  ERROR - too many variables!'

C     Loop over components of H and compute derivatives or improved
C     stepsizes (as well as first and second derivatives).

      NFUNCT = 0
      DO 10 I = 1, N
         HI = H (I)
         URH (I) = HI
         IF (HI .GT. ZERO) THEN

C           Use the stepsize passed in for differencing.

            CALL DIFFER (I, N, X, EPSOBJ, HI, FX, FPFOR, CFPFOR,
     >         FPBAC, CFPBAC, DIAG (I), CFPP, FUNCT)
            IERROR (I) = 0
            NFDEL (I) = 2
            GRAD (I) = HALF * (FPFOR + FPBAC)
         ELSE

C           Generate H from scratch, with GRAD and DIAG as by-products.

            CALL FDSTEP (I, N, X, EPSOBJ, FX, GRAD (I), DIAG (I), H (I),
     >         ERRBND, IERROR (I), NFDEL (I), FUNCT)
         END IF
         NFUNCT = NFUNCT + NFDEL (I)

C        Safeguard Hessian approx. by bounding diagonals away from zero
C        (to preserve positive definitness).

         IF (DIAG (I) .LT. EPSOBJ) THEN
            URDIAG = DIAG (I)
            DIAG (I) = MAX (ABS (DIAG (I)), EPSOBJ)
            IF (LWRITE .GE. 0) WRITE (LWRITE, 1010) I, URDIAG, DIAG (I)
         END IF
   10 CONTINUE

C     Recap results.

      IF (LWRITE .GE. 0) THEN
         WRITE (LWRITE, 1020)
         WRITE (LWRITE, 1030) (I, X (I), GRAD (I), DIAG (I), URH (I),
     >      H (I), NFDEL (I), IERROR (I), I = 1, N)
         WRITE (LWRITE, 1040)
      END IF


C     Termination.
C     ------------

      RETURN


C     Formats.
C     --------

 1000 FORMAT (37H0CENDIF:  Begin gradient calculation.)
 1010 FORMAT (25H0CENDIF:  WARNING - DIAG(, I3, 15H) changed from ,
     >   E12.6, 4H to , E12.6)
 1020 FORMAT (8H1CENDIF:/ 1H0, 4X, 1HI, 4X, 1HV, 20X, 8HGradient,
     >   13X, 15HDiag of Hessian, 6X, 5HOld H, 16X, 5HNew H, 14X,
     >   6HNo Fun, 2X, 6HIERROR)
 1030 FORMAT (1H , 3X, I2, 3X, E18.12, 3X, E18.12, 3X, E18.12, 3X,
     >   E18.12, 3X, E18.12, 4X, I2, 6X, I1)
 1040 FORMAT (1H1)

      END
*DECK FDSTEP
C+----------------------------------------------------------------------
C
      SUBROUTINE FDSTEP (I, N, X, EPSOBJ, FX, FP, FPP, HOPT, ERRBND,
     >   IERROR, NFUNCT, FUNCT)
C
C
C     Description:
C
C           This routine is intended for automatic estimation of optimal
C        finite difference intervals for numerical optimization of functions
C        whose derivatives cannot be calculated analytically.  Each call to
C        FDSTEP will handle one component of the stepsize vector.
C
C
C     Parameters:
C
C        Name    Dimension   Type  I/O/S  Description
C        I                    I    I      Index of the component of X to be
C                                         varied.
C        N                    I    I      Dimension of X.
C        X          N         R    I      Vector of optimization variables.
C        EPSOBJ               R    I      Smallest significant absolute change
C                                         in the objective function FUNCT.
C        FX                   R    I      Initial value of objective function.
C        FP                   R      O    Best central difference approximation
C                                         to first partial derivative of
C                                         objective with respect to the I-th
C                                         component of X.
C        FPP                  R      O    Best central difference approximation
C                                         to second partial derivative of
C                                         objective with respect to the I-th
C                                         component of X.
C        HOPT                 R      O    Estimate of optimal finite difference
C                                         stepsize for I-th component of X.
C        ERRBND               R      O    Estimated bound on (absolute) error
C                                         in first derivatives calculated
C                                         using one sided differencing with
C                                         stepsize HOPT.
C        IERROR               I      O    Error status upon termination. Zero
C                                         means OK; for other values see
C                                         error handling section below.
C        NFUNCT               I      O    Number of calls to FUNCT (increment).
C        FUNCT                     I      Routine for calculating objective
C                                         function (EXTERNAL).
C
C
C     External References:
C
C        DIFFER     Utility for calculating derivatives by finite differences,
C                   along with cancellation error estimates.
C        FUNCT      User-supplied objective function subroutine.
C
C
C     Notes:
C
C        (1)  The goal is to find the smallest stepsize H(KPHI) which yields
C             acceptable relative cancellation errors in the forward and
C             backward difference first derivatives and the central difference
C             second derivative of FUNCT.  The second derivative computed with
C             H(KPHI) is then used to obtain HOPT, an interval which should
C             result in good forward difference approximations to the first
C             derivative.  The "optimality" of the interval chosen is based
C             on several assumptions - consult reference (1) for details.
C
C        (2)  The algorithm also returns estimates FP, the first derivative,
C             and FPP, the second derivative, which are both computed by
C             central differences.
C
C        (3)  The stepsize and derivatives calculated each iteration are
C             saved in arrays for possible later use to save function
C             evaluations.  These are relevant to the step associated with
C             optimization variable I only, and should not be confused with
C             the stepsize and gradient arrays in the calling routine.
C
C        (4)  Constant TINY, set below, is a machine dependent quantity
C             intended to provide some protection against division by zero.
C             The arithmetic expressions involved could still overflow for
C             badly scaled problems in which EPSOBJ is much larger than one.
C
C        (5)  FDSTEP and DIFFER are based on Algorithm FD in Reference (1).
C             The parenthesized labels in the code (e.g. FD1) refer to steps
C             of the algorithm as originally published.
C
C
C     Bibliography:
C
C        (1)  Gill, P., Murray, W., and Wright, M.  Practical Optimization,
C                pp. 341-344.  London:  Academic Press, 1981.
C
C        (2)  Gill, P., Murray, W., Saunders, M., and Wright, M.  "A Procedure
C                for Computing Forward Difference Intervals for Numerical
C                Optimization."  Tech. Rep. SOL 81-25, Dept. of Operations
C                Research, Stanford University, Dec. 1981.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development History:
C
C        16 June 1982   RAK   Original coding ("literal" transcription).
C        23 June 1982   RAK   Revised and extensively restructured.
C        13 July 1982   RAK   Section FD3 revised (small CFPP now OK).
C         2 Sep. 1982   RAK   Redefined function counter NFUNCT.
C        24 Sep. 1982   RAK   Some minor changes adapted from Ref. (2).
C        11 Aug. 1986   RAK   Add IMPLICIT NONE.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      INTEGER
     >   KBOUND
      PARAMETER
     >   (KBOUND = 6)

      REAL
     >   ZERO, ONE, TWO, HALF, BNDFP, BNDLO, BNDUP, ETA, OMEGA, RHO,
     >   TINY
      PARAMETER
     >   (ZERO  = 0.0E+0,
     >    ONE   = 1.0E+0,
     >    TWO   = 2.0E+0,
     >    HALF  = ONE / TWO,
     >    BNDFP = 1.0E-1,
     >    BNDLO = 1.0E-3,
     >    BNDUP = 1.0E-1, 
     >    ETA   = ONE,
     >    OMEGA = ONE,
     >    RHO   = 1.0E+1,
     >    TINY  = 1.0E-32)

C     Arguments.

      INTEGER
     >   I, IERROR, N, NFUNCT
      REAL
     >   EPSOBJ, ERRBND, FP, FPP, FX, HOPT, X (N)

C     Local variables.

      LOGICAL
     >   DONE
      INTEGER
     >   KOUNT, KPHI, KSAVE
      REAL
     >   CFPP, CFPBAC, CFPFOR, H (0:KBOUND), HBAR, FPFOR (0:KBOUND),
     >   FPBAC (0:KBOUND), FPPCEN (0:KBOUND)

C     Procedures.

      EXTERNAL
     >   DIFFER, FUNCT

C     Statement functions.

      LOGICAL
     >   INSIDE, LARGE, MAXLE, SMALL
      REAL
     >   A, B

      INSIDE (A)   = (A .GE. BNDLO) .AND. (A .LE. BNDUP)
      LARGE (A)    = A .GT. BNDUP
      MAXLE (A, B) = MAX (A, B) .LE. BNDFP
      SMALL (A)    = A .LT. BNDLO


C     Execution.
C     ----------

C     Initialization (FD1).

      KOUNT = 0
      NFUNCT = 0
      IERROR = 0
      DONE = .FALSE.

      HBAR = TWO * (ETA + ABS (X (I))) *
     >   SQRT (EPSOBJ / (OMEGA + ABS (FX)))
      H (KOUNT) = HBAR * RHO

      CALL DIFFER (I, N, X, EPSOBJ, H (KOUNT), FX, FPFOR (KOUNT),
     >   CFPFOR, FPBAC (KOUNT), CFPBAC, FPPCEN (KOUNT), CFPP, FUNCT)
      NFUNCT = NFUNCT + 2

C     Accept, increase, or decrease H?  (FD2).

      IF (MAXLE (CFPFOR, CFPBAC) .AND. .NOT.LARGE (CFPP)) THEN
         IF (INSIDE (CFPP)) THEN

C           Accept H and fall through to estimate of optimal H.

            KPHI = KOUNT
         ELSE

C           Decrease H to reduce truncation error (FD4).

C           The cancellation errors in the first derivatives are OK, while
C           that for the second derivative is smaller than necessary.
C           Try to reduce H without letting the errors get too big.

   10       CONTINUE
               KOUNT = KOUNT + 1
               H (KOUNT) = H (KOUNT - 1) / RHO
               CALL DIFFER (I, N, X, EPSOBJ, H (KOUNT), FX,
     >            FPFOR (KOUNT), CFPFOR, FPBAC (KOUNT), CFPBAC,
     >            FPPCEN (KOUNT), CFPP, FUNCT)
               NFUNCT = NFUNCT + 2
               IF (.NOT.MAXLE (CFPFOR, CFPBAC) .OR. LARGE (CFPP)) THEN

C                 We've gone too far - back up one iteration and quit.

                  DONE = .TRUE.
                  KPHI = KOUNT - 1
               ELSE IF (INSIDE (CFPP)) THEN

C                 The current stepsize H is acceptable.

                  DONE = .TRUE.
                  KPHI = KOUNT
               ELSE

C                 Second derivative cancellation error is still smaller
C                 than necessary, but check iteration counter before
C                 continuing.

                  IF (KOUNT .GE. KBOUND) THEN

C                    The iteration limit has been reached.  The error
C                    flag is set as a warning that the stepsize may not
C                    be as small as possible.  Note:  this can also
C                    indicate trouble - a slope discontinuity can also
C                    give this error, so the objective function should be
C                    double checked in this case.

                     IERROR = 1
                     KPHI = KOUNT
                  END IF
               END IF
               IF (.NOT.DONE .AND. IERROR .EQ. 0) GO TO 10

         END IF
      ELSE

C        Increase H to reduce cancellation error bounds in first or second
C        derivatives (FD3).

         KSAVE = -1
   20    CONTINUE

C           Keep track of index of smallest H (if any) with sufficiently
C           small cancellation errors in the one sided first derivatives.

            IF (MAXLE (CFPFOR, CFPBAC) .AND. KSAVE .LT. 0)
     >         KSAVE = KOUNT
            KOUNT = KOUNT + 1
            H (KOUNT) = H (KOUNT - 1) * RHO
            CALL DIFFER (I, N, X, EPSOBJ, H (KOUNT), FX, FPFOR (KOUNT),
     >         CFPFOR, FPBAC (KOUNT), CFPBAC, FPPCEN (KOUNT), CFPP,
     >         FUNCT)
            NFUNCT = NFUNCT + 2
            IF (MAXLE (CFPFOR, CFPBAC) .AND. .NOT.LARGE (CFPP)) THEN

C              Current H is acceptable.

               DONE = .TRUE.
               KPHI = KOUNT
            ELSE IF (KOUNT .GE. KBOUND) THEN

C              No satisfactory interval has been found (FD6).

               IF (.NOT.MAXLE (CFPFOR, CFPBAC)) THEN

C                 Cancellation errors in first derivatives are still too
C                 big - the function is nearly constant.  Use the step
C                 appropriate for the "well scaled" case.

                  IERROR = 2
                  HOPT = HBAR
                  FP = ZERO
                  FPP = ZERO
               ELSE

C                 Second derivative cancellation error is excessive - the
C                 first derivative appears to be nearly constant.  Use the
C                 smallest H with acceptable cancellation error in first
C                 derivatives.

                  IERROR = 3
                  HOPT = H (KSAVE)
                  FP = FPFOR (KSAVE)
                  FPP = ZERO
               END IF
            END IF
            IF (.NOT.DONE .AND. IERROR .EQ. 0) GO TO 20

      END IF
      IF (IERROR .LE. 1) THEN

C        Set output values of derivatives and estimate of optimal H (FD5).

         FP = HALF * (FPFOR (KPHI) + FPBAC (KPHI))
         FPP = FPPCEN (KPHI)
         HOPT = MAX (TWO * SQRT (EPSOBJ / MAX (ABS (FPP), TINY)), TINY)
      END IF

C     Calculate estimate of error bound on one sided derivatives based
C     on finite difference interval HOPT.  If interval should be OK but
C     bound on relative error is large, set warning flag.

      ERRBND = HALF * HOPT * ABS (FPP) + TWO * EPSOBJ / HOPT
      IF ((IERROR .EQ. 0) .AND. (ERRBND .GT. HALF * ABS (FP)))
     >   IERROR = 4


C     Termination.
C     ------------

      RETURN
      END
*DECK DIFFER
C+----------------------------------------------------------------------
C
      SUBROUTINE DIFFER (I, N, X, EPSOBJ, H, FX, FPFOR, CFPFOR,
     >   FPBAC, CFPBAC, FPPCEN, CFPP, FUNCT)
C
C
C     Description:
C
C           Utility for calculating derivatives and cancellation error
C        estimates by finite differences.  Called repeatedly by subroutine
C        FDSTEP during estimation of optimal stepsizes to be passed to a
C        nonlinear optimization routine.  See CENDIF/FDSTEP headers for
C        more details.
C
C
C     Parameters:
C
C        Name    Dimension   Type  I/O/S  Description
C        I                    I    I      Index of the component of X to be
C                                         varied.
C        N                    I    I      Dimension of X in calling routine.
C        X          N         R    I      Vector of optimization variables.
C        EPSOBJ               R    I      Smallest significant absolute change
C                                         in the objective function FUNCT.
C        H                    R    I      Finite difference stepsize for I-th
C                                         component of X.  Assumed > 0.
C        FX                   R    I      Initial value of objective function.
C        FPFOR                R      O    Forward difference approximation
C                                         to first partial derivative of
C                                         objective with respect to the I-th
C                                         component of X.
C        CFPFOR               R      O    Relative cancellation error in FPFOR.
C        FPBAC                R      O    Backward difference approximation
C                                         to first partial derivative of
C                                         objective with respect to the I-th
C                                         component of X.
C        CFPBAC               R      O    Relative cancellation error in FPBAC.
C        FPPCEN               R      O    Central difference approximation
C                                         to second partial derivative of
C                                         objective with respect to the I-th
C                                         component of X.
C        CFPP                 R      O    Relative cancellation error in
C                                         FPPCEN.
C        FUNCT                     I      Routine for calculating objective
C                                         function (EXTERNAL).
C
C
C     Notes:
C
C        (1)  Constant TINY, set below, is a machine dependent quantity
C             intended to provide some protection against division by zero.
C             The arithmetic expressions involved could still overflow for
C             badly scaled problems in which EPSOBJ is much larger than one.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development History:
C
C        23 June 1982    RAK    Original coding.
C         2 Sep. 1982    RAK    Deleted local function counter.
C        11 Aug. 1986    RAK    Add IMPLICIT NONE.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   TWO, FOUR, TINY
      PARAMETER
     >   (TWO  = 2.0E+0,
     >    FOUR = 4.0E+0,
     >    TINY = 1.E-32)

C     Arguments.

      INTEGER
     >   I, N
      REAL
     >   CFPBAC, CFPFOR, CFPP, EPSOBJ, FPBAC, FPFOR, FPPCEN, FX, H,
     >   X (N)

C     Local variables.

      REAL
     >   FXMH, FXPH, XINIT

C     Procedures.

      EXTERNAL
     >   FUNCT


C     Execution.
C     ----------

      XINIT = X (I)

C     Compute first derivative by forward difference...

      X (I) = XINIT + H
      CALL FUNCT (N, X, FXPH)
      FPFOR = (FXPH - FX) / H
      CFPFOR = TWO * EPSOBJ / MAX (H * ABS (FPFOR), TINY)

C     ... and by backward difference.

      X (I) = XINIT - H
      CALL FUNCT (N, X, FXMH)
      FPBAC = (FX - FXMH) / H
      CFPBAC = TWO * EPSOBJ / MAX (H * ABS (FPBAC), TINY)

C     Calculate second derivative by central differences.

      FPPCEN = (FXPH - TWO * FX + FXMH) / MAX (H ** 2, TINY)
      CFPP = FOUR * EPSOBJ / MAX (H ** 2 * ABS (FPPCEN), TINY)

C     Restore X (I) to original value.

      X (I) = XINIT


C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION CHORD (X, Y, I1, I2)
C
C     One-liner: Summed chord-lengths for X-Y curve over range of indices
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes the sum of the Euclidean distance between adjacent
C     points in a curve represented by two arrays. The calculation is
C     rearranged so as to avoid (almost) all chance of overflow, and
C     any unnecessary loss of precision when one component of the
C     distance from one point to the next is small relative to the other.
C     The calling routine must supply beginning and ending indices for
C     the summation (this is intended to facilitate operations with
C     packed data arrays). The result does not depend on the order of
C     I1 and I2.
C
C        CHORD was originally written for use with PLSFIT, which performs
C     parametric cubic interpolation with cumulative chord length as the
C     curve parameter. In use, it is a good idea to try to use all
C     available information to avoid (expensively) re-calculating the
C     lengths of the same intervals over and over; CHORD should be
C     thought of as providing length increments.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     X         R (*)         I      Array of abscissas.
C
C     Y         R (*)         I      Array of ordinates.
C
C     I1,I2                   I      Indices for summation. The loop
C                                    runs from MIN(I1,I2) to MAX(I1,I2)
C                                    so the result is independent of
C                                    order.
C
C     CHORD   R                 O    Function value is the sum of the
C                                    chord lengths along the curve
C                                    defined by the X and Y arrays
C                                    between indices I1 and I2.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C      8 Apr. 1988  RAK/DAS  Reformulated to reduce chance of overflow
C                            or unnecessary underflow. Result is not
C                            dependent on order of I1, I2.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE
      PARAMETER
     &  (ZERO = 0.0E+0,
     &   ONE  = 1.0E+0)

C     Arguments.

      INTEGER
     &   I1, I2
      REAL
     &   CHORD, X (*), Y (*)

C     Local variables.

      INTEGER
     &   I
      REAL
     &   DENOM, DX, DY

C     Execution.
C     ----------

      CHORD = ZERO

      DO 10, I = MIN (I1, I2), MAX (I1, I2) - 1
         DX    = ABS (X (I + 1) - X (I))
         DY    = ABS (Y (I + 1) - Y (I))
         DENOM = MAX (DX, DY)
         IF (DENOM .GT. ZERO) CHORD = CHORD +
     &      DENOM * SQRT (ONE + (MIN (DX, DY) / DENOM) ** 2)
   10 CONTINUE

C     Termination.
C     ------------

      RETURN
      END

C+------------------------------------------------------------------------------
C
      SUBROUTINE CHORDS2D (N, X, Y, NORMALIZ, TOTAL, CHORD)
C
C  ONE-LINER: Cumulative [relative] chord-lengths for 2-space geometric curve
C
C  DESCRIPTION:
C
C        CHORDS2D computes the cumulative Euclidean distances for two or
C     more points on a 2-space curve represented by two arrays, with an
C     option to normalize all values by the TOTAL distance.  Thus CHORD(1)
C     is returned as 0. and CHORD(N) may be either 1. or TOTAL depending
C     on whether NORMALIZ is true or false.
C
C        CHORDS2D was introduced in spite of the existing CHORD function
C     because use of CHORD with the double precision DTNURBS library in a
C     single precision application such as SMOOTH would clash with prior
C     use of the single precision version of CHORD.  The functionality
C     is a little different (multiple values per call, and an option to
C     normalize), and since NURBS are intended for geometric data (i.e.,
C     X and Y expected to have similar units), the careful safeguarding
C     of CHORD is eschewed.
C
C  ARGUMENTS:
C
C     Name    Type/Dimension  I/O/S  Description
C
C     N         I             I      Number of points on the curve. N >= 2.
C
C     X         R (N)         I      Array of abscissas.
C
C     Y         R (N)         I      Array of ordinates.
C
C     NORMALIZ  L             I      .TRUE. means normalize results
C                                    to the interval [0, 1].
C
C     TOTAL     R               O    Total chord length, returned
C                                    because it is otherwise lost
C                                    if results are normalized.
C
C     CHORD     R (N)           O    Cumulative chord lengths as
C                                    described above.
C
C  ENVIRONMENT:  FORTRAN 77 with minor extensions
C
C  HISTORY:
C
C     13 Mar. 1992   DAS   Adapted from CHORD for use with DTNURBS.
C     08 Dec. 1993    "    Added TOTAL argument.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   N
      REAL
     &   CHORD (N), TOTAL, X (N), Y (N)
      LOGICAL
     &   NORMALIZ

C     Local constants.

      REAL
     &   ZERO, ONE
      PARAMETER
     &  (ZERO = 0.0E+0,
     &   ONE  = 1.0E+0)

C     Local variables.

      INTEGER
     &   I
      REAL
     &   DINV

C     Execution.

      CHORD (1) = ZERO

      DO 10, I = 2, N
         CHORD (I) = CHORD (I - 1) +
     &      SQRT ((X (I) - X (I - 1)) ** 2 + (Y (I) - Y (I - 1)) ** 2)
   10 CONTINUE

      TOTAL = CHORD (N)

      IF (NORMALIZ) THEN
         DINV = ONE / TOTAL
         DO 20, I = 2, N
            CHORD (I) = CHORD (I) * DINV
   20    CONTINUE
         CHORD (N) = ONE
      END IF

      RETURN
      END
C+**********************************************************************
C
      SUBROUTINE COPY ( NX, X1, X2 )
C
C PURPOSE:  COPY copies the elements of X1(*) into X2(*).
C
C PARAMETERS:
C    ARG   DIM   TYPE I/O/S DESCRIPTION 
C     NX    -      I    I   No. of elements in input and output arrays
C     X1   NX      R    I   Array to be copied
C     X2   NX      R    O   Copy of X1(*)
C    
C ENVIRONMENT: FORTRAN IV
C
C AUTHOR: David Saunders, Informatics, Palo Alto, CA.  (07/30/82)
C
C-**********************************************************************
C
      DIMENSION  X1(NX), X2(NX)
C
      DO 20 I=1,NX
         X2(I) = X1(I)
 20   CONTINUE
C
      RETURN
      END
C+---------------------------------------------------------------------
C
      SUBROUTINE CSDVAL (NX, X, Y, NU, U, B, C, D, S, SP, SPP)
C
C     One-liner: Cubic spline evaluation on an array, with derivatives.
C     ----------
C
C     Purpose:
C     --------
C
C        CSDVAL evaluates the 1-dimensional cubic spline function
C
C     S(U) = Y(I) + B(I)*(U - X(I)) + C(I)*(U - X(I))**2 + D(I)*(U - X(I))**3
C
C     along with its first and second derivatives
C
C     SP(U)  =      B(I)          + 2*C(I)*(U - X(I))  + 3*D(I)*(U - X(I))**2
C     SPP(U) =    2*C(I)          + 6*D(I)*(U - X(I))
C
C     using Horner's Rule for polynomials.  If U < X(1), then I = 1 is used;
C     if U > X(NX) then I = NX-1 is used. The data must be monotonic, but
C     may be increasing or decreasing.  (The coefficients associated with
C     the NXth point are not required.)  Most of the effort is in deciding
C     which set of spline coefficients to use - an "interpolation search"
C     is employed which should be highly efficient if the knots are fairly
C     uniformly spaced (commonly the case).
C
C        CSDVAL was developed as a companion routine to CSFIT, but it should
C     work with any method of generating cubic spline coefficients.  CSEVAL
C     should be used if derivatives are not required.  See also NOTES.
C
C     Arguments:
C     ----------
C
C     Name     Dimension  Type  I/O/S  Description
C     NX                   I    I      Number of data points defining the
C                                      spline ("knots"); >= 2.
C
C     X        NX          R    I      Abscissas of knots. Must be monotone
C                                      increasing or decreasing (consistent
C                                      with the data used to fit the spline).
C
C     Y        NX          R    I      Ordinates of knots.
C
C     NU                   I    I      Number of points at which to evaluate
C                                      the spline.
C
C     U        NU          R    I      Abscissas at which to evaluate spline.
C
C     B,C,D    NX          R    I      Spline coefficients, e.g. as computed
C                                      by CSFIT.
C
C     S       NU           R      O    Spline values at points U(I).
C
C     SP      NU           R      O    1st derivative values at points U(I).
C
C     SPP     NU           R      O    2nd derivative values at points U(I).
C
C     External modules:
C     -----------------
C
C     INTERVAL  Interpolation search for interval containing a point.
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN 77 V4.7
C     ------------
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and eight character symbols are not (yet) standard.
C
C     (2)  If derivatives are not desired, use CSEVAL instead. Or, take
C          advantage of the fact that the returned values are computed in
C          the order  SPP, SP, S so that if only first derivative is needed,
C          for example, the same array can be passed for SPP as for SP.
C
C     Bibliography:
C     -------------
C
C     (1)  Forsythe, G., M. Malcolm, and C. Moler.  Computer Methods for
C             Mathematical Computations.  Englewood Cliffs: Prentice-Hall,
C             1977.  (Chapter 4)
C
C     Author:  Robert Kennelly, Informatics General Corp.
C     -------
C
C     History:
C     --------
C
C      8/27/84    RGL    Adapted from CSEVAL to provide derivatives.
C     10/20/87    RAK    Abstract the interval search (with revisions) as
C                        a separate module.
C     06/09/91    DAS    Ensured S (U) = Y (NX) if U = X (NX).  (This
C                        is not guaranteed by CSDVAL's use of the cubic
C                        in the (NX-1)th interval for such a U, except
C                        if the arithmetic is exact, which it isn't.)
C
C----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE, TWO, THREE, SIX
      PARAMETER
     &  (ONE   = 1.0E+0,
     &   TWO   = 2.0E+0,
     &   THREE = 3.0E+0,
     &   SIX   = 6.0E+0)

C     Arguments.

      INTEGER
     &   NU, NX
      REAL
     &   B (NX), C (NX), D (NX), S (NU), SP (NU), SPP (NU), U (NU),
     &   X (NX), Y (NX)

C     Local variables.

      INTEGER
     &   IU, LEFT
      REAL
     &   ARROW, DX, XRIGHT

C     Execution.
C     ----------

      ARROW = SIGN (ONE, X (2) - X (1))
      XRIGHT = X (NX)
      LEFT = 1

      DO 10, IU = 1, NU

C        Search for the best "left-hand" endpoint to use for interpolation.

         CALL INTERVAL (NX, X, U (IU), ARROW, LEFT)

C        Evaluate the spline.  Note that if U is off-scale on the left,
C        the coefficients of the first interval are used for extrapolation,
C        while if U is greater than X(NX), the (NX-1)th set is used.
C        This means when U = X (NX), roundoff could mean S (U) is not
C        exactly Y (NX).  Hence the test for this special case.
C        Note also the order, in case y" is not really needed - see NOTES.

         DX = U (IU) - X (LEFT)

         SPP (IU) = TWO * C (LEFT) +
     &       (DX * (SIX * D (LEFT)))

         SP  (IU) = B (LEFT) +
     &       (DX * (TWO * C (LEFT) +
     &       (DX * (THREE * D (LEFT)))))

         S   (IU) = Y (LEFT) +
     &       (DX * (B (LEFT) +
     &       (DX * (C (LEFT) +
     &       (DX * (D (LEFT)))))))

         IF (U (IU) .EQ. XRIGHT) S (IU) = Y (NX)
   10 CONTINUE

C     Termination.
C     ------------

      RETURN
      END
C+---------------------------------------------------------------------
C
      SUBROUTINE CSEVAL (NDATA, X, Y, NU, U, B, C, D, S)
C
C     One-liner: Cubic spline evaluation on an array.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        CSEVAL evaluates the 1-dimensional cubic spline function
C
C     S(U) = Y(I) + B(I)*(U - X(I)) + C(I)*(U - X(I))**2 + D(I)*(U - X(I))**3
C
C     using Horner's Rule for polynomials. Normally, if U < X(1), then I = 1
C     is used, and if U > X(NX) then I = NX - 1 is used.  However, this version
C     treats periodic data properly - see the usage of NDATA.
C
C        The data must be monotonic, but may be increasing or decreasing.
C     (The coefficients associated with the NXth point are not required.)
C     Most of the effort is in deciding which set of spline coefficients to
C     use.  An "interpolation" search is employed which should be highly
C     efficient if the knots are fairly uniformly spaced (commonly the case).
C
C        CSEVAL was developed as a companion routine to CSFIT, but it should
C     work with any method of generating cubic spline coefficients.  CSDVAL
C     should be used if derivatives of the spline are required.
C
C     Arguments:
C     ----------
C
C     Name     Dimension  Type  I/O/S  Description
C     NDATA                I    I      NX = |NDATA| is the number of data points
C                                      defining the spline ("knots");
C                                      use NDATA = -NX to signify the periodic
C                                      data case; NX >= 2.
C
C     X        NX          R    I      Abscissas of knots. Must be monotone
C                                      increasing or decreasing (consistent
C                                      with the data used to fit the spline).
C
C     Y        NX          R    I      Ordinates of knots.
C
C     NU                   I    I      Number of points at which to evaluate
C                                      the spline.
C
C     U        NU          R    I      Abscissas at which to evaluate spline.
C
C     B,C,D    NX          R    I      Spline coefficients, e.g. as computed
C                                      by CSFIT.
C
C     S       NU           R      O    Spline values at points U(I).
C
C     External modules:
C     -----------------
C
C        INTERVAL   Interpolation search for interval containing a point.
C
C     Environment:  FORTRAN 77 + IMPLICIT NONE
C     ------------
C
C     Bibliography:
C     -------------
C
C     (1)  Forsythe, G., M. Malcolm, and C. Moler.  Computer Methods for
C             Mathematical Computations.  Englewood Cliffs: Prentice-Hall,
C             1977.  (Chapter 4)
C
C     Author:  Robert Kennelly, Informatics General Corp.
C     -------
C
C     History:
C     --------
C
C      7/27/84    RAK    Original design and coding, based in part on
C                        SEVAL from Forsythe, Malcolm, and Moler, and
C                        using a search method adapted from Sedgewick.
C                        An earlier rewrite of SEVAL by Trosin included
C                        conversion from a FUNCTION to a SUBROUTINE to
C                        allow more than one evaluation per call.
C     10/19/87    RAK    Abstracted the interval search (with revisions)
C                        as a separate module.
C     06/09/91    DAS    Ensured S (U) = Y (NX) if U = X (NX).  (This
C                        is not guaranteed by CSEVAL's use of the cubic
C                        in the (NX-1)th interval for such a U, except
C                        if the arithmetic is exact, which it isn't.)
C     08/16/97     "     Handled the periodic case via NX < 0 kludge.
C
C----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   NDATA, NU
      REAL
     &   B (*), C (*), D (*), S (NU), U (NU), X (*), Y (*)

C     Local variables.

      INTEGER
     &   IU, LEFT, NX
      REAL
     &   ARROW, DX, PERIOD, XEVAL, XLEFT, XRIGHT

C     Execution.
C     ----------

      ARROW  = SIGN (ONE, X (2) - X (1))
      NX     = ABS (NDATA)
      XRIGHT = X (NX)
      LEFT   = 1

      IF (NDATA .GT. 0) THEN  ! Non-periodic case

         DO IU = 1, NU

            XEVAL = U (IU)

C           Search for the best "left-hand" endpoint to use for interpolation.

            CALL INTERVAL (NX, X, XEVAL, ARROW, LEFT)

C           Evaluate the spline.  Note that if U is off-scale on the left,
C           the coefficients of the first interval are used for extrapolation,
C           while if U is greater than X(NX), the (NX-1)th set is used.
C           This means when U = X (NX), roundoff could mean S (U) is not
C           exactly Y (NX).  Hence the test for this special case.

            DX = XEVAL - X (LEFT)

            S (IU) = Y (LEFT) +
     &        (DX * (B (LEFT) +
     &        (DX * (C (LEFT) +
     &        (DX * (D (LEFT)))))))

            IF (XEVAL .EQ. XRIGHT) S (IU) = Y (NX)

         END DO

      ELSE  ! Periodic case

         IF (ARROW .EQ. ONE) THEN
            XLEFT = X (1)
         ELSE
            XLEFT = XRIGHT
            XRIGHT = X (1)
         END IF

         PERIOD = XRIGHT - XLEFT  ! Definitely positive

         DO IU = 1, NU

            XEVAL = U (IU)

            IF (XEVAL .LT. XLEFT) THEN

               XEVAL = XRIGHT - MOD (XLEFT - XEVAL, PERIOD)

            ELSE IF (XEVAL .GT. XRIGHT) THEN

               XEVAL = XLEFT + MOD (XEVAL - XRIGHT, PERIOD)

            END IF

            CALL INTERVAL (NX, X, XEVAL, ARROW, LEFT)

            DX = XEVAL - X (LEFT)

            S (IU) = Y (LEFT) +
     &        (DX * (B (LEFT) +
     &        (DX * (C (LEFT) +
     &        (DX * (D (LEFT)))))))

            IF (XEVAL .EQ. X (NX)) S (IU) = Y (NX)

         END DO

      END IF

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE CSFIT ( N, X, Y, IENDL, DERIVL, IENDR, DERIVR,
     +                   B, C, D, IER )
C
C  ACRONYM:  Cubic Spline FIT
C            -     -      ---
C  PURPOSE:
C    CSFIT computes the coefficients B(I), C(I), and D(I), I=1,2,...,N
C    for the conventional cubic spline defined by
C
C    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
C
C    for X(I) <= X <= X(I+1).  The spline interpolates the data points
C    (X(I), Y(I)) (also called the knots) and has continuous first and
C    second derivatives.
C
C    This implementation contains a number of enhancements to the sub-
C    routine SPLINE of Forsythe, Malcolm, and Moler (ref. below):
C
C       1. Several boundary conditions are provided, including period-
C          icity and the familiar "natural" spline (S"(X)=0 at the end
C          points). Different conditions may be specified at each end.
C          
C       2. Degenerate cases (N=2 and 3) are handled as far as possible
C          depending on the end-conditions specified - all of them are
C          applicable if N is at least 4. (A lot of algebra was needed
C          to derive the code for these cases.)
C
C       3. An Nth set of coefficients is computed,  corresponding to X
C          beyond X(N).  (This was needed by one version of CSEVAL and
C          is retained because some applications expect it.)
C
C       4. The search in CSEVAL now permits decreasing abscissas;  the
C          algebra in CSFIT applies without change in this case.
C
C  METHOD:
C    (0)  Represent the spline in the ith interval by
C
C         s(x) = w*y(i+1) + (1 - w)*y(i) + h(i)**2 *
C                [(w**3 - w)*sigma(i+1) + ((1 - w)**3 - w)*sigma(i)],
C
C         where h(i) = x(i+1) - x(i)  and  w = (x - x(i)) / h(i).
C
C         (as on p.71 of Forsythe, Malcolm, and Moler).
C    (1)  Check for and process any degenerate case of the cubic spline.
C         (See description for N under ARGUMENTS.)
C    (2)  Construct the bulk of the symmetric tridiagonal system
C         involved (common to all cases).
C    (3)  Apply the indicated end conditions for this system.
C    (4)  Solve the system for the sigmas.
C    (5)  Derive the corresponding polynomial coefficients.
C
C  ARGUMENTS:
C    NAME   DIM TYPE I/O/S DESCRIPTION
C    N       -   I     I   Number of data points, or knots.  Minimum N
C                          depends on the spline to be constructed--
C                          > 3: All possible end conditions are feasible;
C                          = 3: If IENDL=2, then IENDR must also be 2;
C                          = 2: The same order derivative must be used
C                          at both ends in all cases; if third deriva-
C                          tive is used, the same value also must be
C                          applied at each end.
C    X       N   R     I   Abscissas of the knots, strictly increasing
C                          or strictly decreasing.
C    Y       N   R     I   Ordinates of the knots.  Y(N)=Y(1) for the
C                          cyclic case (IENDL=4).
C    IENDL   -   I     I   = 0: The 3rd derivative of the spline at the
C                               left endpoint is to match the 3rd deriv-
C                               ative of the cubic passing through the
C                               first 4 data points; if N=3 and IENDL=
C                               IENDR=0, both parts of the spline become
C                               the same quadratic defined by the pts.;
C                          = 1: The 1st derivative of the spline at the
C                               left endpoint is to be the given DERIVL;
C                          = 2: The 2nd derivative of the spline at the
C                               left endpoint is to be the given DERIVL;
C                          = 3: The 3rd derivative of the spline at the
C                               left endpoint is to be the given DERIVL;
C                          = 4: The spline and its first 3 derivatives
C                               at X(N) are to match the values at X(1).
C                               This is the cyclic, or periodic, case.
C    DERIVL  -   R     I   Value of derivative used if IENDL=1, 2, or 3;
C                          ignored if IENDL=0 or 4.
C    IENDR,  -             As for IENDL, DERIVL, but pertaining to the
C    DERIVR                right endpoint; ignored if IENDL=4.
C    B,C,D   N   R     O   Spline coefficients (see PURPOSE and NOTES).
C    IER     -   I     O   =0: No errors were detected.
C                          =1: Too few data points; N < 2.
C                          =2: Degenerate cases of N=2 for all end con-
C                              ditions and of N=3 with second derivative
C                              end condition require that the same order
C                              of derivative is applied at both ends.
C                          =3: Degenerate case where N=2 with third der-
C                              ivative applied must have the same value
C                              for the derivative at each endpoint.
C                          =4: Cyclic mode -- Y(N) not equal to Y(1).
C                          =5: Non-cyclic mode -- IENDL or IENDR is
C                              beyond the range [0,4].
C
C  EXTERNAL REFERENCES:
C    NAME     DESCRIPTION
C    TRICPS   Solves a cyclic, positive-definite, symmetric tridiagonal
C             system of equations.  Also used here for true tridiagonal
C             cases.
C
C  SYSTEM DEPENDENCIES:
C    (1)  IMPLICIT NONE is an extension of FORTRAN 77.
C
C  NOTES:
C    (1)  Y(I) = S(X(I))
C         B(I) = S'(X(I))
C         C(I) = S''(X(I))/2
C         D(I) = S'''(X(I))/6
C
C    (2)  To evaluate the spline, use subroutine CSEVAL.
C
C  ENVIRONMENT:  VAX, CRAY -- FORTRAN 77
C
C  BIBLIOGRAPHY:
C    Forsythe, Malcolm, and Moler, Computer Methods for Mathematical
C       Computations. Englewood Cliffs, NJ:  Prentice-Hall, 1977.
C
C  DEVELOPMENT HISTORY:
C    DATE   INITIALS   DESCRIPTION 
C   c. 1977  F,M,M     Original implementation (IENDL=IENDR=0 case only)
C   06/25/84 RGL/RCL/  Now provides for multiple end conditions, includ-
C            RAK/DAS   ing cyclic functions.
C   10/23/87 DAS       Revised PURPOSE description (prompted by CSEVAL's
C                      new ability to deal with descending abscissas).
C
C  AUTHORS:  Forsythe, Malcolm, and Moler, c. 1977
C            Informatics, Palo Alto, CA., June 1984
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER  N, IENDL, IENDR, IER
      REAL     X(N), Y(N), B(N), C(N), D(N)

C ... Local constants:

      REAL      ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX
      PARAMETER ( ZERO=0.E+0, ONE=1.E+0, TWO=2.E+0, THREE=3.E+0,
     +            FOUR=4.E+0, FIVE=5.E+0, SIX=6.E+0 )

C ... Local variables:

      INTEGER  I, ISYS, NSYS
      REAL     DERIVL, DERIVR, DYDX1, DYDX2, DYDX3, H1, H2, H3

C ... Procedures:

      EXTERNAL  TRICPS

C ... Execution:

C ... Check for too few data points:

      IER = 0
      IF ( N.LT.2 ) THEN
         IER = 1
         GO TO 999
      END IF

C ... N will commonly be more than enough, so:

      IF ( N.GT.4 ) GO TO 10

C ... Handle degenerate cases:

      IF ( N.EQ.2 ) THEN
         IF ( IENDL.NE.IENDR ) THEN
            IER = 2
            GO TO 999
         END IF

         H1    = X(2) - X(1)
         DYDX1 = (Y(2) - Y(1)) / H1

         IF ( IENDL.EQ.0 ) THEN
            B(1) = DYDX1
            C(1) = ZERO
            D(1) = ZERO
         ELSE IF ( IENDL.EQ.1 ) THEN
            B(1) = DERIVL
            C(1) = (DYDX1 * THREE - DERIVL * TWO - DERIVR) / H1
            D(1) = (DERIVL + DERIVR - DYDX1 * TWO) / H1**2
         ELSE IF ( IENDL.EQ.2 ) THEN
            B(1) = DYDX1 - (DERIVL * FOUR - DERIVR) * (H1 / SIX)
            C(1) = DERIVL / TWO
            D(1) = (DERIVR - DERIVL) / (H1 * SIX)
         ELSE IF ( IENDL.EQ.3 ) THEN
            IF ( DERIVL.EQ.DERIVR ) THEN
               B(1) = (Y(2) - Y(1)) / H1 - DERIVL / SIX * H1**2
               C(1) = ZERO
               D(1) = DERIVL / SIX
            ELSE
               IER = 3
               GO TO 999
            END IF
         ELSE IF ( IENDL.EQ.4 ) THEN
            IF ( Y(N).NE.Y(1) ) THEN
               IER = 4
               GO TO 999
            END IF
            B(1) = ZERO
            C(1) = ZERO
            D(1) = ZERO
         END IF
         GO TO 50

      ELSE IF ( N.EQ.3 ) THEN
         H1  = X(2) - X(1)
         H2  = X(3) - X(2)
         DYDX1 = (Y(2) - Y(1)) / H1
         DYDX2 = (Y(3) - Y(2)) / H2

         IF ( IENDL.EQ.2 .OR. IENDR.EQ.2 ) THEN
            IF ( IENDL.NE.IENDR ) THEN
               IER = 2
               GO TO 999
            END IF

            C(1) = DERIVL / SIX
            C(2) = (DYDX2 - DYDX1 - (DERIVL * H1 + DERIVR * H2) /
     +             SIX) / ((H1 + H2) * TWO)
            C(3) = DERIVR / SIX
            GO TO 30
         ELSE IF ( IENDL.EQ.4 ) THEN
            IF ( Y(N).NE.Y(1) ) THEN
               IER = 4
               GO TO 999
            END IF
            C(1) = (DYDX1 - DYDX2) / (H1 + H2)
            C(2) =-C(1)
            C(3) = C(1)
            GO TO 30
         END IF

      ELSE IF ( N.EQ.4 ) THEN
         IF ( IENDL.EQ.2 .AND. IENDR.EQ.2 ) THEN
            H1  = X(2) - X(1)
            H2  = X(3) - X(2)
            H3  = X(4) - X(3)
            DYDX1 = (Y(2) - Y(1)) / H1
            DYDX2 = (Y(3) - Y(2)) / H2
            DYDX3 = (Y(4) - Y(3)) / H3

            C(1)  = DERIVL / SIX
            C(2)  = (DYDX3 - DYDX2 - DERIVR * H3 / SIX -
     +              (H3 / H2 + ONE) * (DYDX2 - DYDX1 -
     +              DERIVL * H1 / SIX) * TWO) /
     +              (H2 - (H2 + H3) *
     +              (H1 / H2 + ONE) * FOUR)
            C(3)  = (DYDX2 - DYDX1 - DERIVL * H1 / SIX -
     +              (H1 / H2 + ONE) * (DYDX3 - DYDX2 -
     +              DERIVR * H3 / SIX) * TWO) /
     +              (H2 - (H2 + H3) *
     +              (H1 / H2 + ONE) * FOUR)
            C(4)  = DERIVR / SIX
            GO TO 30
         END IF
      END IF

   10 CONTINUE

C ... Set up the bulk of the tridiagonal system (common to all cases).
C     B = diagonal, D = offdiagonal, C = right hand side.

      D(1) = X(2) - X(1)

      DO 20 I = 2, N-1
         D(I) = X(I+1) - X(I)
         B(I) = ( X(I+1) - X(I-1) ) * TWO
         C(I) = ( Y(I+1) - Y(I) ) / D(I) -
     +          ( Y(I) - Y(I-1) ) / ( X(I) - X(I-1) )
   20 CONTINUE

C ... Now for the boundary conditions.  Cyclic case ignores IENDR.

      ISYS = 1
      IF ( IENDL.NE.4 ) THEN

C ...    The two ends of the spline are independent.

         NSYS = N
         D(N) = ZERO

C ...    Set the left end condition:

         IF ( IENDL.EQ.0 ) THEN

C ...       Use divided differences for the 3rd derivative of the
C           cubic through the first 4 points:

            B(1) = -D(1)
            C(1) = ZERO
            IF ( N.GT.3 )
     +         C(1) = (C(3)/(X(4)-X(2)) - C(2)/(X(3)-X(1)))*
     +                D(1)**2/(X(4)-X(1))
         ELSE IF ( IENDL.EQ.1 ) THEN
            B(1) = D(1) * TWO
            C(1) = - ( DERIVL - ( Y(2) - Y(1) )/D(1) )
         ELSE IF ( IENDL.EQ.2 ) THEN
            NSYS = NSYS - 1
            ISYS = 2
            C(1) = DERIVL/SIX
            C(2) = C(2) - D(1) * C(1)
         ELSE IF ( IENDL.EQ.3 ) THEN
            B(1) = -D(1)
            C(1) = ( DERIVL/SIX ) * D(1)**2
         ELSE
            IER  = 5
            GO TO 999
         END IF

C ...    Set the right end condition similarly:

         IF ( IENDR.EQ.0 ) THEN
            B(N) = -D(N-1)
            C(N) = ZERO
            IF ( N.GT.3 )
     +         C(N) = (C(N-2)/(X(N-1)-X(N-3)) - C(N-1)/(X(N)-X(N-2)))*
     +                D(N-1)**2/(X(N)-X(N-3))
         ELSE IF ( IENDR.EQ.1 ) THEN
            B(N) = D(N-1) * TWO
            C(N) = DERIVR - (Y(N)-Y(N-1))/D(N-1)
         ELSE IF ( IENDR.EQ.2 ) THEN
            NSYS   = NSYS - 1
            C(N)   = DERIVR/SIX
            C(N-1) = C(N-1) - D(N-1) * C(N)
            D(N-1) = ZERO
         ELSE IF ( IENDR.EQ.3 ) THEN
            B(N) = -D(N-1)
            C(N) = -(DERIVR/SIX) * D(N-1)**2
         ELSE
            IER  = 5
            GO TO 999
         END IF
      ELSE

C ...    Cyclic boundary conditions, IENDL = 4.  Sigma(N) = sigma(1) is
C        used to keep the system symmetric, giving one equation less:

         IF ( Y(N).NE.Y(1) ) THEN
            IER = 4
            GO TO 999
         END IF

         NSYS = N-1
         B(1) = (D(1) + D(NSYS)) * TWO
         C(1) = (Y(2) - Y(1)) / D(1) - (Y(N) - Y(NSYS)) / D(NSYS)

      END IF

C ... Solve the tridiagonal system for the sigmas (returned in C(*)):

      CALL TRICPS ( NSYS, B(ISYS), D(ISYS), C(ISYS), C(ISYS) )

      IF ( IENDL.EQ.4 ) C(N) = C(1)

C ... Derive the spline coefficients (noting that the DXs are lost):

   30 DO 40 I = 1, N-1
         H1   = X(I+1) - X(I)
         B(I) = (Y(I+1) - Y(I)) / H1 - (C(I+1) + C(I)*TWO) * H1
         D(I) = (C(I+1) - C(I)) / H1
         C(I) = C(I) * THREE
   40 CONTINUE

C ... Set coefficients for rightmost knot (and beyond).

   50 IF ( IENDL.NE.4 ) THEN

C ...    1st deriv. at X(N) uses right-sided deriv. for X in X(N-1):X(N)

         H1   = X(N) - X(N-1)
         B(N) = B(N-1) + C(N-1) * H1 * TWO + D(N-1) * H1**2 * THREE
         C(N) = C(N-1) + D(N-1) * H1 * THREE
         D(N) = D(N-1)
      ELSE

C ...    Use periodicity for IENDL=4 case:

         B(N) = B(1)
         C(N) = C(1)
         D(N) = D(1)
      END IF

  999 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE CSQUAD (N, X, F, OFFSET, C, FINT)
C
C
C     Description and usage:
C
C           CSQUAD (Cubic Spline QUADrature) is a specialized routine for
C        integrating a function using its cubic spline representation.  The
C        (definite) integral up to each knot is returned, with an optional
C        constant added to the value of the integral array at the first point.
C        Output array FINT may share storage with spline coefficient array C.
C
C           This version uses a compact, but tricky, form of the integral
C        taken from Forsythe, et al., which is faster than the naive
C        integration formula used in the original.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C         N                  I    I      Number of knots.
C         X        N         R    I      Abscissas (spline knots Xi).
C         F        N         R    I      Function values.
C         OFFSET             R    I      Value assigned to the first element
C                                        of FINT (will be added to all).
C         C        N         R    I      Array of coefficients of (X - Xi)**2
C                                        from spline fit (= F"/2 at X = Xi).
C         FINT     N         R      O    Integral of F up to each knot.
C                                        FINT (1) = OFFSET; FINT (N) =
C                                        OFFSET + integral from X (1) to X (N).
C
C
C     Standards violations:  IMPLICIT NONE is non-standard.
C
C
C     Bibliography:
C
C        (1)  Forsythe, G., M. Malcolm, and C. Moler.  Computer
C                Methods for Mathematical Computations (Englewood Cliffs:
C                Prentice-Hall, 1977).  Chapter 5.
C
C
C     Development environments:  Digital VAX-11/780  VMS/V4.1   FORTRAN
C                                Cray X-MP/48        COS/V1.14  CFT/V1.11
C
C
C     Author:  Robert Kennelly, Sterling Software, Palo Alto, CA.
C
C
C     History:
C
C         6 July 1984    RAK    Initial design/coding as INTEG in FLO6QNM.
C        10 July 1984    RAK    Modified to permit shared storage,
C                               added OFFSET parameter.
C        11 Feb. 1986    RAK    Use IMPLICIT NONE, and declare all
C                               constants locally.
C        20 Mar. 1986    RAK    Switch to cute FM&M form of integral.
C         9 Sep. 1986    RAK    Renamed CSQUAD for general use.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   HALF, TWELTH
      PARAMETER
     >   (HALF   = 1.0E+0 / 2.0E+0,
     >    TWELTH = 1.0E+0 / 12.0E+0)

C     Local variables.

      INTEGER
     >   I, N
      REAL
     >   C (N), F (N), FINT (N), H, OFFSET, X (N)


C     Execution.
C     ----------

C     Compute integral on each subinterval.  Work backwards to avoid
C     overwriting the C array if it was also passed in as FINT.

      DO 10, I = N - 1, 1, -1
         H = X (I + 1) - X (I)
         FINT (I + 1) = H      * ((F (I) + F (I + 1)) * HALF -
     >                  H ** 2 *  (C (I) + C (I + 1)) * TWELTH)
   10 CONTINUE

C     Sum the contributions from each subinterval.

      FINT (1) = OFFSET
      DO 20, I = 2, N
         FINT (I) = FINT (I) + FINT (I - 1)
   20 CONTINUE


C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE CURV2D (N, XP, XPP, YP, YPP, KAPPA)
C
C  PURPOSE: CURV2D returns signed values for the curvature at N points
C           on a curve in (X, Y) space given the 1st and 2nd derivatives
C           with respect to arc-length at those points.  The formulas for
C           the magnitude and sign of curvature are:
C
C                   |KAPPA| = SQRT (XPP ** 2 + YPP ** 2)
C
C              sign (KAPPA) = sign (XP * YPP - XPP * YP)
C
C           See module FDCURV for the non-parametric case.
C
C  ARGUMENTS:  Obvious from the description: N >= 1; all others are REAL (N).
C
C  HISTORY: 08/27/91  Parametric form for one point at a time, analogous
C                     to FDCURV's Y = Y (X) form, which was adapted from
C                     FD12K when the latter found not to lend itself to
C                     application to one point at a time.
C           10/25/91  One point only was a bad decision.  Introduced N.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    N
      REAL       XP (N), XPP (N), YP (N), YPP (N), KAPPA (N)

C     Local variables:

      INTEGER    I

C     Execution:

      DO 10, I = 1, N
         KAPPA (I) = SIGN (SQRT (XPP (I) ** 2 + YPP (I) ** 2),
     >                     XP (I) * YPP (I) - XPP (I) * YP (I))
   10 CONTINUE

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE DECOMP(NDIM,N,A,COND,IPVT,WORK)
C
      INTEGER NDIM,N
      REAL A(NDIM,N),COND,WORK(N)
      INTEGER IPVT(N)
C
C     DECOMPOSES A SINGLE PRECISION MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     USE SOLVE TO COMPUTE SOLUTIONS TO LINEAR SYSTEMS.
C
C     INPUT..
C
C        NDIM = DECLARED ROW DIMENSION OF THE ARRAY CONTAINING  A.
C
C        N = ORDER OF THE MATRIX.
C
C        A = MATRIX TO BE TRIANGULARIZED.
C
C     OUTPUT..
C
C        A  CONTAINS AN UPPER TRIANGULAR MATRIX  U  AND A PERMUTED
C          VERSION OF A LOWER TRIANGULAR MATRIX  I-L  SO THAT
C          (PERMUTATION MATRIX)*A = L*U .
C
C        COND = AN ESTIMATE OF THE CONDITION OF  A .
C           FOR THE LINEAR SYSTEM  A*X = B, CHANGES IN  A  AND  B
C           MAY CAUSE CHANGES  COND  TIMES AS LARGE IN  X .
C           IF COND+1 = COND, A IS SINGULAR TO WORKING PRECISION
C           COND = 1.0E+32  IF EXACT SINGULARITY IS DETECTED.
C
C        IPVT = THE PIVOT VECTOR.
C           IPVT(K) = THE INDEX OF THE K-TH PIVOT ROW
C           IPVT(N) = (-1)**(NUMBER OF INTERCHANGES)
C
C     WORK SPACE..  THE VECTOR  WORK  MUST BE DECLARED AND INCLUDED
C                   IN THE CALL.  ITS INPUT CONTENTS ARE IGNORED.
C                   ITS OUTPUT CONTENTS ARE USUALLY UNIMPORTANT.
C
C     THE DETERMINANT OF A CAN BE OBTAINED ON OUTPUT BY
C        DET(A) = IPVT(N) * A(1,1) * A(2,2) * ... * A(N,N).
C
C-----------------------------------------------------------------------
C
      REAL    EK, T, ANORM, YNORM, ZNORM
      INTEGER NM1, I, J, K, KP1, KB, KM1, M
      INTRINSIC ABS, SIGN
C
      IPVT(N) = 1
      IF (N .EQ. 1) GO TO 80
      NM1 = N - 1
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0E+0
      DO 10 J = 1, N
         T = 0.0E+0
         DO 5 I = 1, N
            T = T + ABS(A(I,J))
    5    CONTINUE
         IF (T .GT. ANORM) ANORM = T
   10 CONTINUE
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      DO 35 K = 1,NM1
         KP1= K+1
C
C        FIND PIVOT
C
         M = K
         DO 15 I = KP1,N
            IF (ABS(A(I,K)) .GT. ABS(A(M,K))) M = I
   15    CONTINUE
         IPVT(K) = M
         IF (M .NE. K) IPVT(N) = -IPVT(N)
         T = A(M,K)
         A(M,K) = A(K,K)
         A(K,K) = T
C
C        SKIP STEP IF PIVOT IS ZERO
C
         IF (T .EQ. 0.0E+0) GO TO 35
C
C        COMPUTE MULTIPLIERS
C
         DO 20 I = KP1,N
             A(I,K) = -A(I,K)/T
   20    CONTINUE
C
C        INTERCHANGE AND ELIMINATE BY COLUMNS
C
         DO 30 J = KP1,N
             T = A(M,J)
             A(M,J) = A(K,J)
             A(K,J) = T
             IF (T .EQ. 0.0E+0) GO TO 30
             DO 25 I = KP1,N
                A(I,J) = A(I,J) + A(I,K)*T
   25        CONTINUE
   30    CONTINUE
   35 CONTINUE
C
C     COND = (1-NORM OF A)*(AN ESTIMATE OF 1-NORM OF A-INVERSE)
C     ESTIMATE OBTAINED BY ONE STEP OF INVERSE ITERATION FOR THE
C     SMALL SINGULAR VECTOR.  THIS INVOLVES SOLVING TWO SYSTEMS
C     OF EQUATIONS, (A-TRANSPOSE)*Y = E  AND  A*Z = Y  WHERE  E
C     IS A VECTOR OF +1 OR -1 CHOSEN TO CAUSE GROWTH IN Y.
C     ESTIMATE = (1-NORM OF Z)/(1-NORM OF Y)
C
C     SOLVE (A-TRANSPOSE)*Y = E
C
      DO 50 K = 1, N
         T = 0.0E+0
         IF (K .EQ. 1) GO TO 45
         KM1 = K-1
         DO 40 I = 1, KM1
            T = T + A(I,K)*WORK(I)
   40    CONTINUE
   45    EK = 1.0E+0
         IF (T .LT. 0.0E+0) EK = -1.0E+0
         IF (A(K,K) .EQ. 0.0E+0) GO TO 90
         WORK(K) = -(EK + T)/A(K,K)
   50 CONTINUE
      DO 60 KB = 1, NM1
         K = N - KB
         T = 0.0E+0
         KP1 = K+1
         DO 55 I = KP1, N
            T = T + A(I,K)*WORK(K)
   55    CONTINUE
         WORK(K) = T
         M = IPVT(K)
         IF (M .EQ. K) GO TO 60
         T = WORK(M)
         WORK(M) = WORK(K)
         WORK(K) = T
   60 CONTINUE
C
      YNORM = 0.0E+0
      DO 65 I = 1, N
         YNORM = YNORM + ABS(WORK(I))
   65 CONTINUE
C
C     SOLVE A*Z = Y
C
      CALL SOLVE(NDIM, N, A, WORK, IPVT)
C
      ZNORM = 0.0E+0
      DO 70 I = 1, N
         ZNORM = ZNORM + ABS(WORK(I))
   70 CONTINUE
C
C     ESTIMATE CONDITION
C
      COND = ANORM*ZNORM/YNORM
      IF (COND .LT. 1.0E+0) COND = 1.0E+0
      RETURN
C
C     1-BY-1
C
   80 COND = 1.0E+0
      IF (A(1,1) .NE. 0.0E+0) RETURN
C
C     EXACT SINGULARITY
C
   90 COND = 1.0E+32
      RETURN
      END
C+---------------------------------------------------------------------
C
      SUBROUTINE DSTRIB (MODE, NP, P, NPTS, XMIN, XMAX, X)
C
C  PURPOSE:  DSTRIB generates a distribution of points in the range
C            [XMIN,XMAX].  A choice of distributions is offered via
C            argument MODE.  Most of the choices require additional
C            parameters,  unlike those offered by subroutine XGRID,
C            from which DSTRIB was derived.
C
C            (Two reasons to leave XGRID intact: (1) It is not easy
C            to track down all existing usages to change the calls;
C            (2) XGRID avoids fractional powers (i.e., logarithms);
C            DSTRIB does fractional exponentiation even if the  ex-
C            ponent is 1.)
C            
C  METHOD:   The uniform distribution is retained for completeness.
C            MODES 1 and 2 are generalizations of XGRID's sinusoid-
C            al distributions, providing for powers of cosine other
C            than 1.   MODE 3 provides for bunching at an arbitrary
C            internal point only,  or at an internal point and each
C            end as well.
C 
C  NOTES:    The sinusoidal distributions (MODE = -1,1,2,3)  have a
C            common feature of giving greater bunching in the  same
C            regions as "WIDTH" parameter P(1)=1.0 gives if P(1) is
C            LESS than 1.0 (and greater than zero - i.e., fraction-
C            al power as opposed to squaring, cubing, etc.).   When
C            P(1) > 1. is used,  results are less easily described,
C            and likely to be other than what was intended.  Brief-
C            ly, the higher this exponent, the less the bunching at
C            the regions given by P(1)=1.,  and the more the bunch-
C            ing towards either XMAX (MODE=1) or the midpoint (MODE
C            =2) or the "CENTER" (P(2), MODE=3).
C
C  ARGUMENTS:
C    ARG      DIM  TYPE  I/O/S  DESCRIPTION
C    MODE      -     I     I    MODE controls type of distribution,
C                               as described below.
C    NP        -     I     I    Number of input parameters P(*). In
C                               the case of MODE=3, NP should be 3,
C                               not 2, to allow for returning info.
C    P        NP     R     I    Parameters depending on MODE below.
C    NPTS      -     I     I    No. of pts in desired distribution.
C  XMIN,XMAX   -     R     I    First and last values to include.
C    X       NPTS    R     O    Desired distribution of points.
C
C  USAGE:
C       MODE = -1 means same as for MODE = 1 except the bunching is
C                 is towards XMAX instead of XMIN.
C       MODE = 0  means uniform distribution. Use NP=1; P not used.
C       MODE = 1  means sinusoidal bunching towards XMIN if P(1) is
C                 in the range (0.,1.].   P(1) > 1. gives bunchings
C                 at both XMIN and XMAX that differ from each other
C                 in general, and are probably not what you want.
C       MODE = 2  means sinusoidal bunchings of the  same  type  at
C                 both XMIN and XMAX if P(1) is in (0.,1.].   These
C                 bunchings are reduced, and additional bunching is
C                 given at the mid-point, if P(1) > 1. is used. The
C                 distribution is symmetric about the center point.
C       MODE = 3  means sinusoidal bunching around the internal pt.
C                 defined by P(2), where XMIN < P(2) < XMAX.   P(1)
C                 controls bunching as for MODE=2,  with additional
C                 bunching at the ends of the interval if P(1) > 1.
C
C                 MODE 3 NOTES:
C
C              1. P(3) is RETURNED with  the corresponding interior
C                 index M (as a REAL) such that X(M) = P(2).
C
C              2. If X = P(2) is to be included exactly (as assumed
C                 here),  ensuring smoothness in  dX is nontrivial.
C                 DSTRIB does NOT attempt to ensure this.   It does
C                 no more than determine the  proportion  of points
C                 to be placed on either side of the interior point
C                 then generate two similar distributions.  Routine
C                 SMOOTHX can be used to smooth the result if reqd.
C
C  ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C  HISTORY:
C  05/24/85    DAS    Added NP, P(*) arguments to handle generation
C                     of internally-bunched distributions (MODE=3);
C                     generalized XGRID's options (MODE=1,2) in the
C                     process.  Left XGRID intact.
C  06/25/85    DAS    Revised MODE=3: straightforward sine function
C                     used now, where the modified sine function of
C                     subroutine BEVAL had been used originally. It
C                     turned out that the latter gives a  different
C                     type of distribution at the end  points  that
C                     is hard to predict and is therefore unwise.
C  01/22/87    DAS    Patch in MODE=3 to handle symmetric case like
C                     P(2)=0.5 on [0.,1.] with odd N properly. (For
C                     even N, there is no symmetric solution here.)
C                     Also: had to return subscript of interior pt.
C                     Also: introduced MODE=-1 for completeness.
C  08/14/89    DAS    Return MODE=3's interior point index in P(3),
C                     not NP.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   MODE, NP, NPTS
      REAL
     >   P (NP), X (NPTS), XMAX, XMIN

C     Local constants.

      REAL
     >   EPS, HALF, ONE
      PARAMETER
     >  (EPS = 1.E-6, HALF = 5.E-1, ONE = 1.E+0)

C     Local variables.

      INTEGER
     >   I, MIDDLE, NPTS1, NPTS2
      REAL
     >   DTHETA, DX, PIBY2, RANGE, YNORM

C     Execution.


      PIBY2  = ASIN (ONE)
      DTHETA = PIBY2 / REAL (NPTS - 1)
      RANGE  = XMAX - XMIN


      IF (MODE .EQ. 0) THEN

C  *     Uniform distribution:

         DX = RANGE / REAL (NPTS - 1)

         DO 50 I = 2, NPTS - 1
            X (I) = REAL (I - 1) * DX + XMIN
   50    CONTINUE


      ELSE IF (MODE .EQ. 1) THEN

C  *     Sinusoidal bunching at the low end if 0.0 < P (1) <= 1.0;
C        unsymmetric bunching at both ends if P (1) > 1.0:

         DO 100 I = 2, NPTS - 1
            X (I) = XMIN +
     +         RANGE * (ONE - (COS (REAL (I-1) * DTHETA)) ** P (1))
  100    CONTINUE


      ELSE IF (MODE .EQ. -1) THEN

C  *     Sinusoidal bunching at the high end if 0.0 < P (1) <= 1.0;
C        unsymmetric bunching at both ends if P (1) > 1.0:

         DO 150 I = 2, NPTS - 1
            X (I) = XMIN + RANGE * (SIN (REAL (I-1) * DTHETA)) ** P (1)
  150    CONTINUE


      ELSE IF (MODE .EQ. 2) THEN

C  *     Symmetric sinusoidal bunching at both ends if 0. < P (1) <= 1.;
C        additional bunching in the middle if P (1) > 1.0:

         DTHETA = DTHETA + DTHETA
         RANGE  = RANGE * HALF
         MIDDLE = NPTS / 2

         DO 200 I = 2, MIDDLE
            YNORM = ONE - (COS (REAL (I-1) * DTHETA)) ** P (1)
            X (I)        = XMIN + RANGE * YNORM
            X (NPTS+1-I) = XMAX - RANGE * YNORM
  200    CONTINUE

C  *     The following avoids difficulties with cosine at what is
C        supposed to be exactly PI/2 but may be slightly greater,
C        giving a negative cosine, when NPTS is odd:

         IF (MIDDLE + MIDDLE .LT. NPTS) X (MIDDLE+1) = XMIN + RANGE


      ELSE IF (MODE .EQ. 3) THEN

C  *     Modified sinusoidal bunching around internal point P(2) if
C        0.0 < P(1) <= 1.0; bunching also towards the end-points if
C        P(1) > 1.0.    The following gives proportional numbers of
C        points in the two intervals involved.    For instance,  if
C        NPTS = 100 and CENTER = XMIN + 0.7 * RANGE, then 70 of the
C        points should be in [XMIN,CENTER] and 30 in [CENTER,XMAX],
C        with bunching towards CENTER (and possibly towards the end
C        points XMIN and XMAX if P(1) > 1.).
C
C        See NOTE above about lack of continuity in dX across P(2).
C
C        Patch for symmetric case with odd N:  we need to round up,
C        but e.g. for N = 101 on [0.,1.] with roundoff we might get
C        NINT (50.49999...) = 50, not 51 - hence the EPS.
C
C        Also:  we need to return the subscript of the interior pt.
C        NP was the only available argument - caller beware.

         NPTS1 = NINT (REAL (NPTS) * (P (2) - XMIN) / RANGE + EPS)
         DTHETA = PIBY2 / (NPTS1 - 1)
         RANGE = P (2) - XMIN

         DO 300 I = 2, NPTS1 - 1
            X (I) = XMIN + RANGE * (SIN (REAL (I-1) * DTHETA)) ** P (1)
  300    CONTINUE

         P (3) = REAL (NPTS1)
         X (NPTS1) = P (2)
         NPTS2 = NPTS - NPTS1 + 1
         DTHETA = PIBY2 / REAL (NPTS2 - 1)
         RANGE = XMAX - P (2)

         DO 400 I = 2, NPTS2 - 1
            X (I+NPTS1-1) = P (2) + RANGE *
     +         (ONE - (SIN (PIBY2 + REAL (I-1) * DTHETA)) ** P (1))
  400    CONTINUE
         
      END IF


C  *  Ensure precise values at end-points:

      X (1) = XMIN
      X (NPTS) = XMAX

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE DTDLSQ ( N, A, B, C, D, R, S, SSQ )
C
C ACRONYM: Diagonal + TriDiagonal system; Least SQuares solution
C          -          -  -                -     --
C
C PURPOSE: DTDLSQ solves one overdetermined linear system of the form
C
C                     |       |  |   |       |   |
C                     |   D   |  | x |       | r |
C                     |       |  |   |       |   |
C                     | - - - |          =   | - |
C                     |       |              |   |
C                     |   T   |              | s |
C                     |       |              |   |
C
C          where  D  is diagonal, and  T  is tridiagonal.
C
C METHOD:  The QR decomposition of the left-hand-side matrix is used,
C          where  R  is readily shown to be upper tridiagonal.  R  is
C          produced by a sequence of symmetric plane rotations, which
C          are also applied to the right-hand-side vector.  The least
C          squares solution follows in the usual way by  solving  the
C          (square upper tridiagonal) portion of the transformed sys-
C          tem.
C
C          Annihilation of the tridiagonal lower portion is  done  in
C          3-part steps involving a(i), c(i+1), b(i) in that order as
C          indicated for i=1:
C
C                           | d         |     | r r r     |
C                           |   d       |     |   x       |
C                           |     d     |     |     d     |
C                           |       d   |     |       d   |
C                           |         d |     |         d |
C                  Q3 Q2 Q1 | - - - - - |  =  | - - - - - |
C                           | a b       |     | 0 0       |
C                           | c a b     |     | 0 x x     |
C                           |   c a b   |     |   c a b   |
C                           |     c a b |     |     c a b |
C                           |       c a |     |       c a |
C
C          Note that all elements of a, b, c, d  are numbered by row.
C          Vectors a, b, c are overwritten with the  upper triangular
C          factor R.  The right-hand-side vectors are transformed in-
C          place, and the least squares solution overwrites the upper
C          half, r.   The minimum sum of squares is derived from  the
C          rest of the transformed right-hand-side vector s.
C
C ARGUMENTS:
C   ARG  DIM TYPE I/O/S DESCRIPTION
C    N    -    I    I   Order of the system (N >= 3).
C    A    N    R   I/S  Input with main diagonal of  tridiagonal part
C                       of the matrix; destroyed upon return.
C    B    N    R   I/S  Input with upper diagonal of tridiagonal part
C                       of the matrix in B(1):B(N-1);  B(N) is needed
C                       to avoid special handling of i=N-1 iteration.
C    C    N    R   I/S  Input with lower diagonal of tridiagonal part
C                       of the matrix in  C(2):C(N);  destroyed  upon
C                       return.
C    D    N    R   I/S  Input with upper diagonal of the matrix;  de-
C                       stroyed upon return.
C    R    N    R  I/S/O Input with upper half of right-hand side vec-
C                       tor; returned with the least squares solution
C                       of the overdetermined system.
C    S    N    R   I/S  Input with lower half of right-hand-side vec-
C                       tor; destroyed upon return.
C    SSQ  -    R    O   Minimum sum of squares corresponding  to  the
C                       calculated solution.
C
C ERROR HANDLING:  None.  Divide by zero can mean the matrix is sing-
C                  ular but is more likely to mean  unwitting  errors
C                  in the call to DTDLSQ.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   11/01/84   MAS    Algorithm description.
C   Nov. '84   DAS    Initial implementation,  for efficient solution
C                     of the thickness/curvature  refinement  problem
C                     in airfoil design.
C
C AUTHORS: Michael Saunders, Stanford University, Palo Alto, CA.
C          David Saunders,   Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C ... Arguments:

      INTEGER   N
      REAL      A(N), B(N), C(N), D(N), R(N), S(N), SSQ

C ... Local variables:

      INTEGER   I
      REAL      A2, B1, B2, CS, RI, R1, R11, R12, R13, R2, R22,
     +          SI, SN, S1, T1, T2

C ... Intrinsics:

      INTRINSIC SQRT

C ... Execution:

      DO 20 I = 1, N-1

C ...    Eliminate a(i):

         T1 = SQRT ( D(I) ** 2 + A(I) ** 2 )
         CS = D(I) / T1
         SN = A(I) / T1
         T2 = SN * B(I)
         B1 =-CS * B(I)

C ...    Modify affected elements of RHS:

         R1 = R(I)
         S1 = S(I)
         RI = CS * R1 + SN * S1
         SI = SN * R1 - CS * S1

C ...    Eliminate c(i+1):

         R11 = SQRT ( T1 ** 2 + C(I+1) ** 2 )
         CS  = T1 / R11
         SN  = C(I+1) / R11
         R12 = CS * T2 + SN * A(I+1)
         A2  = SN * T2 - CS * A(I+1)
         R13 =           SN * B(I+1)
         B2  =          -CS * B(I+1)
         R(I)   = CS * RI + SN * S(I+1)
         S(I+1) = SN * RI - CS * S(I+1)

C ...    Eliminate b(i):

         R22 = SQRT ( D(I+1) ** 2 + B1 ** 2 )
         CS  = D(I+1) / R22
         SN  = B1 / R22
         R2  = R(I+1)
         R(I+1) = CS * R2 + SN * SI
         S(I)   = SN * R2 - CS * SI

C ...    Here is the ith row of upper triangular factor R ...

         A(I) = R11
         B(I) = R12
         C(I) = R13

C ...    ... and the modified (i+1)th row of T and D:

         A(I+1) = A2
         B(I+1) = B2
         D(I+1) = R22

   20 CONTINUE

C ... Nth iteration eliminates just one element of T:

      T1 = SQRT ( D(N) ** 2 + A(N) ** 2 )
      CS = D(N) / T1
      SN = A(N) / T1
      A(N) = T1
      RI   = R(N)
      R(N) = CS * RI + SN * S(N)
      S(N) = SN * RI - CS * S(N)

C ... Solve the upper portion of the triangularized system:

      R(N) = R(N) / A(N)
      R(N-1) = ( R(N-1) - B(N-1) * R(N) ) / A(N-1)

      DO 30 I = N-2, 1, -1
         R(I) = ( R(I) - B(I) * R(I+1) - C(I) * R(I+2) ) / A(I)
   30 CONTINUE

C ... Minimum residual is contained in lower part of modified RHS:

      SSQ = S(1)
      DO 40 I = 2, N
         SSQ = S(I) ** 2 + SSQ
   40 CONTINUE

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE FD12K (N, X, F, FP, FPP, FK)
C
C  PURPOSE: FD12K returns estimates of 1st and 2nd derivatives and of
C           curvature, by finite differencing, for each of the points
C           (X(I), F(I)), I = 1 : N.  The abscissas are assumed to be
C           nonuniform, and they must be monotonic.
C
C           This routine combines calls to FDCNTR, FD1SID, FDCURV for
C           the common case of needing results for N >= 2 points.
C
C           If (say) curvature is wanted at a single point only, call
C           either FDCNTR or FD1SID and FDCURV directly.
C
C  INPUTS:  X(*) & F(*) are N-vectors defining some curve in 2-space.
C           For N > 2, the 3-pt formulas are used for all I (with the
C                      one-sided forms used at the end-points).
C           For N = 2, the 2-pt formulas are used.
C
C  OUTPUTS: FP, FPP, FK are N-vectors representing 1st and 2nd deriv-
C           atives and curvature respectively.  These are assigned in
C           reverse order (FK, FPP, FP) so that a call such as
C
C                     CALL FD12K (N, X, Y, YP, YP, YP)
C
C           can be used if just 1st derivatives are desired, to avoid
C           declaring storage for FPP and FK. (Similarly for the case
C           when 1st and 2nd derivatives are desired but curvature is
C           not. The unnecessary arithmetic in these cases is consid-
C           ered preferable to another argument and extra logic.)
C
C  METHOD:  Central differencing is used at all interior points, with
C           one-sided 3-point formulas used at each end-point.
C
C           The curvature formula is safeguarded against overflow  in
C           the presence of large 1st derivatives.  The basic formula
C           used here is:
C
C               FK (I) = FPP (I) / (1 + FP(I) ** 2) ** 3/2
C
C           Note that if X is not necessarily monotonic, curvature is
C           defined as
C
C               CURVATURE = (X" ** 2  +  Y" ** 2) ** 1/2
C
C           where " means 2nd derivative with respect to  arc-length.
C           See modules CURV2D and CURV3D for these parametric cases.
C
C  NOTES:   1. Finite differencing errors can be large if the delta-X
C              values are too small,  especially if the precision  in
C              the function values is less than full.
C           2. Nevertheless, finite differences have been observed to
C              behave better than the analytic derivatives of splines
C              in airfoil geometry applications.
C
C  EXTERNALS:
C           FDCNTR modularizes the central 3-point formulas for first
C                  and second derivatives by finite differences.
C           FDCURV modularizes the curvature formula (safe-guarded).
C           FD1SID modularizes the 1-sided forward and backward 3-pt.
C                  formulas for first and second derivatives.
C
C  HISTORY:
C           09/15/83   DAS   Initial implementation (interior pts only).
C           12/27/85   DAS   End points are handled by FD1SID now.
C           09/18/87   DAS   The N=2 case is handled now.
C           08/21/89   DAS   Formulation revised to use separate dF/dX
C                            terms instead of a common denominator.
C           08/17/91   DAS   Introduced FDCNTR and FDCURV when it was
C                            found that FD12K did not lend itself to
C                            application to one point at a time.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    N
      REAL       X (N), F (N), FP (N), FPP (N), FK (N)

C     Local constants:

      REAL       ZERO
      PARAMETER (ZERO = 0.E+0)

C     Local variables:

      INTEGER    I
      REAL       FPI, FPPI

C     Procedures:

      EXTERNAL   FDCNTR, FDCURV, FD1SID

C     Execution:

C     Assign values in the order  curvature, f", f'  so that the
C     user can pass the same array for, say, curvature and f" if
C     curvature is not wanted:

      IF (N .EQ. 2) THEN

         FK (1)  = ZERO
         FPP (1) = ZERO
         FP (1)  = (F (2) - F (1)) / (X (2) - X (1))
         FK (2)  = ZERO
         FPP (2) = ZERO
         FP (2)  = FP (1)

      ELSE

C        Forward 3-pt. differencing for the first point:

         CALL FD1SID (1, 1, X, F, FPI, FPPI)
         CALL FDCURV (FPI, FPPI, FK (1))
         FPP (1) = FPPI
         FP (1)  = FPI
      
C        Central 3-pt. differencing for the bulk of the points:

         DO 20 I = 2, N-1
            CALL FDCNTR (I, X, F, FPI, FPPI)
            CALL FDCURV (FPI, FPPI, FK (I))
            FPP (I) = FPPI
            FP (I)  = FPI
   20    CONTINUE

C        Backward differencing for the last point:

         CALL FD1SID (N, -1, X, F, FPI, FPPI)
         CALL FDCURV (FPI, FPPI, FK (N))
         FPP (N) = FPPI
         FP (N)  = FPI
      END IF

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE FD1SID (I, INC, X, F, FP, FPP)
C
C  PURPOSE: FD1SID returns one-sided 3-point finite-difference estimates
C           of the 1st and 2nd derivatives at the point  ( X(I), F(I) ).
C           If INC = 1, points I, I+1, I+2 are used,  while if INC = -1,
C           points  I-2, I-1, I are used. The abscissas need not be uni-
C           formly spaced.  
C
C  ARGS:    Obvious from PURPOSE.
C
C  METHOD:  FPP is computed first,  in case only FP is desired,  so that
C           the same item may be passed for both arguments. The formula-
C           tion is similar to that of central differencing - see FD12K.
C
C  HISTORY: 12/27/85  DAS  Initial implementation  (prompted by the need
C                          to approximate an airfoil leading edge with a
C                          cubic having specified slope at one end).
C           08/21/89  DAS  Formulation revised as for centrals (FD12K).
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    INC, I
      REAL       X (*), F (*), FP, FPP

C     Local constants:

      REAL       ONE
      PARAMETER (ONE = 1.E+0)

C     Local variables:

      INTEGER  I1, I2
      REAL     DEL1, DEL2, DIV, DX1, DX2, W

C     Execution:

C     The minus signs take care of themselves for the backward case.

      I1   = I  + INC
      I2   = I1 + INC
      DX1  =  X (I1) - X (I)
      DEL1 = (F (I1) - F (I)) / DX1
      DX2  =  X (I2) - X (I1)
      DEL2 = (F (I2) - F (I1)) / DX2
      DIV  = ONE / (DX1 + DX2)
      W    = -DX1 * DIV
      FPP  = (DEL2 - DEL1) * (DIV + DIV)
      FP   = W * DEL2 + (ONE - W) * DEL1

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE FDCNTR (I, X, F, FP, FPP)
C
C  PURPOSE: FDCNTR returns central 3-point finite-difference estimates
C           of the 1st and 2nd derivatives at the point  (X(I), F(I)).
C           Use FD1SID for the end-point cases.
C
C  ARGS:    Obvious from PURPOSE.
C
C  METHOD:  FPP is computed first, in case only FP is desired, so that
C           the same item may be passed for both arguments.
C
C  HISTORY: 08/17/91  DAS  FDCNTR adapted from FD12K's in-line code,
C                          as for FD1SID, for the case of a single I
C                          at a time (which FD12K can't do).
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    I
      REAL       X (*), F (*), FP, FPP

C     Local constants:

      REAL       ONE
      PARAMETER (ONE = 1.E+0)

C     Local variables:

      REAL     DEL1, DEL2, DIV, DX1, DX2, W

C     Execution:

      DX1  =  X (I) - X (I-1)
      DEL1 = (F (I) - F (I-1)) / DX1
      DX2  =  X (I+1) - X (I)
      DEL2 = (F (I+1) - F (I)) / DX2
      DIV  = ONE / (DX1 + DX2)
      W    = DX2 * DIV
      FPP  = (DEL2 - DEL1) * (DIV + DIV)
      FP   = W * DEL1 + (ONE - W) * DEL2

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE FDCURV (DYDX, D2YDX2, KAPPA)
C
C  PURPOSE: FDCURV (Finite-Difference CURVature estimate) returns a
C           safe-guarded value for curvature at one point on the
C           curve Y = Y (X) given the 1st and 2nd derivatives at
C           that point, using the formula
C
C               KAPPA = D2YDX2 / (1 + DYDX ** 2) ** 3/2
C
C           The sign of KAPPA is clearly that of the 2nd derivative.
C           The derivatives could well be obtained from a spline, but
C           experience shows finite differencing can be preferable.
C
C           See modules CURV2D and CURV3D for the parametric cases.
C
C  ARGUMENTS:  Obvious from the description.  KAPPA is REAL.
C
C  HISTORY: 08/17/91  Derived FDCURV from FD12K, along with FDCNTR
C                     when it was found that FD12K did not lend
C                     itself to application to one point at a time.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      REAL       DYDX, D2YDX2, KAPPA

C     Local constants:

      REAL       DYDXMAX, ONE
      PARAMETER (DYDXMAX = 1.E+10, ONE = 1.E+0)

C     Local variables:

      REAL       TERM

C     Execution:

      TERM = ONE + (MIN (ABS (DYDX), DYDXMAX)) ** 2
      KAPPA = D2YDX2 / (TERM * SQRT (TERM))

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE FOILGRD (N, XMIN, XMAX, WTLIN, WTQUAD, WTSIN, WTCOS, X)
C
C  ONE-LINER:  1-D grid distribution tailored for airfoils
C
C  DESCRIPTION:
C
C        FOILGRD generates N abscissas in the range [XMIN, XMAX] using a
C     combination of linear, quadratic, sine, and cosine shape functions.
C     The four weights should add to 1. and each should be in [0., 1.].
C     Values suggested for a typical airfoil:  0.04, 0.0, 0.3, 0.66
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S  DIM   DESCRIPTION
C    N        I     I     -    |N| = desired number of abscissas;
C                              N < 0 can be used to reverse the bunching
C                              of the N > 0 case, to suit the lower
C                              surface of a wrap-around airfoil
C  XMIN,XMAX  R     I     -    First and last values to include in X (*)
C   WTLIN     R     I     -    Weight applied to the linear component
C   WTQUAD    R     I     -    Weight applied to the quadratic term
C   WTSIN     R     I     -    Weight applied to the sinusoidal term
C   WTCOS     R     I     -    Weight applied to the cosine term
C    X        R     O     N    Desired distribution of abscissas;
C                              X (1) = XMIN and X (|N|) = XMAX
C
C  08/23/95  JJR  Adaptation of FOILGRID (sine + quadratic).
C  08/29/95  DAS  Tidied up the nomenclature.
C  07/21/97   "   Implemented N < 0 option to simplify usage.
C
C  AUTHOR:  James Reuther, Ames Research Center, Mt. View, CA
C
C----------------------------------------------------------------------

C     Arguments:

      INTEGER    N
      REAL       XMIN, XMAX, WTLIN, WTQUAD, WTSIN, WTCOS, X (*)

C     Local constants:

      REAL       HALF, ONE, PI

      PARAMETER (HALF = 0.5,
     >           ONE  = 1.0,
     >           PI   = 180.)

C     Local variables:

      INTEGER    I, J, NPTS
      REAL       RANGE, RN, SNORM, WL, WQ, WS, WC, XNORM
      REAL       TOTAL

C     Intrinsics:

      REAL       COSD

C     Execution.

      TOTAL = ONE / (WTLIN + WTQUAD + WTSIN + WTCOS)
      WL    = WTLIN * TOTAL
      WQ    = WTQUAD* TOTAL
      WS    = WTSIN * TOTAL
      WC    = WTCOS * TOTAL
      RANGE = XMAX - XMIN
      NPTS  = ABS (N)
      RN    = ONE / REAL (NPTS - 1)
      X (1) = XMIN

      IF (N .GT. 0) THEN  ! Ascending order

         DO I = 2, NPTS - 1
            XNORM = REAL (I - 1) * RN
            SNORM = WL * XNORM +
     >              WQ * XNORM ** 2 +
     >              WS * (ONE - COSD (XNORM * PI * HALF)) +
     >              WC * (ONE - COSD (XNORM * PI)) * HALF
            X (I) = XMIN + SNORM * RANGE
         END DO

      ELSE  ! N < 0 - reverse the order

         J = NPTS - 1
         DO I = 2, NPTS - 1
            J = J - 1
            XNORM = ONE - REAL (J) * RN
            SNORM = WL * XNORM +
     >              WQ * XNORM ** 2 +
     >              WS * (ONE - COSD (XNORM * PI * HALF)) +
     >              WC * (ONE - COSD (XNORM * PI)) * HALF
            X (I) = XMIN + SNORM * RANGE
         END DO

      END IF

      X (NPTS) = XMAX

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE GETLINE (LUN, COMMENT, LINE, LAST, IOS)
C
C  One-liner:  Low-level data reader; suppresses trailing comments/blanks
C
C  Description and usage:
C
C        GETLINE is intended to be a standard low level text input utility,
C     providing a uniform input format which permits free use of blank
C     lines, comment lines, trailing comments, and "commented-out" lines.
C     It reads a record from logical unit LUN and returns the "significant"
C     portion in LINE, with only trailing COMMENTs, blanks, or tabs after
C     LAST.  It returns LAST = 0 if the line is effectively empty.
C
C        Double-COMMENTs are replaced with single ones so that COMMENT may
C     still be used in a string if required.  SPECIAL CASE: if COMMENT is
C     the FIRST significant character, the line is considered empty.  This
C     covers the common case of "commenting-out" lines with more than one
C     COMMENT character, as in  !!!! 0.800000   1.
C
C  Arguments:
C
C     Name    Type/Dimension  I/O/S  Description
C
C     LUN     I               I      Logical unit for reading text.
C
C     COMMENT C*1             I      Character signalling the end of the
C                                    "significant" part of a text string
C                                    and the beginning of comments. If
C                                    blank, this feature is disabled and
C                                    much less scanning is needed.
C
C     LINE    C*(*)             O/S  Buffer for reading one record;
C                                    returned with "significant" text in
C                                    LINE (1 : LAST) unless LAST is zero.
C
C     LAST    I                 O    Index of last significant character
C                                    in LINE.  If LINE is null or all
C                                    blanks, tabs, and/or comment, LAST
C                                    is returned as zero.
C
C     IOS     I                 O    Error status of read;  = 0 means no
C                                    read error, negative is EOF, and the
C                                    meaning of positive IOS is system-
C                                    dependent.
C
C  Environment:  Fortran 90
C
C  Notes:
C
C     (1)  Note that GETLINE does NOT keep reading if it encounters an
C          empty line.  Part of the purpose of GETLINE is to avoid this
C          very drawback of FORTRAN's list-directed I/O.
C
C  Authors:  Ronald Langhi/Robert Kennelly,  Sterling Software, Palo Alto, CA
C
C  History:
C
C     10 Mar. 1987  RGL/RAK  Initial design and code as GETSTRING.
C      4 Apr. 1987    RAK    Cosmetics, documentation revised.
C      8 May  1987    DAS    Name changed to GETLINE.
C     22 Apr. 1988  DAS/RAK  Trim trailing tabs as well as the blanks.
C      7 July 1990  DAS/RAK  Handle special case of commenting-out lines
C                            with more than one leading COMMENT.
C     19 May  1999    DAS    Fortran 90's SIZE=N replaces nonstandard
C                            '(Q, A)' read for counting # characters.
C     20 May  1999    DAS    NO!  SIZE=... raises EOR/EOF issues that
C                            appear insoluble in a portable way.
C                            Use LEN_TRIM (LINE) and avoid SIZE=...
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Constants.

      CHARACTER
     >   BLANK * 1, HTAB * 1
      PARAMETER
     >   (BLANK = ' ',
     >    HTAB  = CHAR (9))   ! HTAB = 9 for Absoft FORTRAN on 68000 machines

C     Arguments.

      INTEGER
     >   LUN, LAST, IOS
      CHARACTER
     >   COMMENT * 1, LINE * (*)

C     Local variables.

      INTEGER
     >   I, J, N              ! N can be LAST throughout, but may be inefficient
      LOGICAL
     >   SUPPRESD

C     Execution.
C     ----------

C     Read the next line.  The ADVANCE = 'NO' must be present with SIZE = N.
C     NO:  Even with PAD = 'YES' on the OPEN, we get an end-of-record value
C          for IOS (-2 for DEC, -4006 for SGI) which cannot be portably
C          distinguished from EOF (IOS < 0).  It seems SIZE = ... is unusable.

C***  READ (LUN, '(A)', SIZE = N, IOSTAT = IOS, ADVANCE = 'NO') LINE
C***  IF (IOS /= 0) GO TO 99

C***  IF (N == 0) THEN
C***     GO TO 99
C***  ELSE
C***     N = MIN (N, LEN (LINE))
C***  END IF

      N = 0
      READ (LUN, '(A)', IOSTAT = IOS) LINE
      IF (IOS /= 0) GO TO 99

      N = LEN_TRIM (LINE)

C     Perform search for comments?
C     ----------------------------

      IF (COMMENT /= BLANK) THEN

C        Examine one character at a time in LINE (1 : N) for a
C        transition between significant characters and comments.

         I = 0
   10    CONTINUE
            I = I + 1
            IF (LINE (I : I) == COMMENT) THEN

C              Handle a special case here rather than impact all lines.
C              If the FIRST significant character is a COMMENT, consider
C              the line suppressed.  This covers the common case where
C              more than one COMMENT character is used to "comment-out"
C              an input line.
C
C              Search for preceding significant text:

               SUPPRESD = .TRUE.
               DO 20, J = 1, I-1
                  IF (LINE (J : J) /= BLANK .AND.
     >                LINE (J : J) /= HTAB) SUPPRESD = .FALSE.
   20          CONTINUE

               IF (SUPPRESD) THEN
                  N = 0
                  GO TO 40       ! For consistency; GO TO 99 is more direct.
               ELSE IF (I == N) THEN

C                 Single comment at the very end of the text - done.

                  N = I - 1
                  GO TO 40
               ELSE IF (LINE (I+1 : I+1) == COMMENT) THEN

C                 Double comment - strip the first one from the text,
C                 decrement N, and keep searching.  The DO-loop is
C                 clunky but required by the standard.

                  DO 30, J = I, N - 1
                     LINE(J : J) = LINE(J+1 : J+1)
   30             CONTINUE
                  N = N - 1
               ELSE

C                 Single comment embedded in the text - done.

                  N = I - 1
                  GO TO 40
               END IF
            END IF
            IF (I < N) GO TO 10

   40    CONTINUE
      END IF

C     Remove trailing blanks and tabs.
C     --------------------------------

C     We're taking advantage of the fact that if the DO-loop runs to
C     completion, then N = 0 on exit.

      DO 50, N = N, 1, -1
         IF (LINE (N : N) /= BLANK .AND.
     >       LINE (N : N) /= HTAB) GO TO 99
   50 CONTINUE

C     Termination.
C     ------------

   99 LAST = N
   
      write(*,*) 'DEBUG: end of getline ',IOS, TRIM(LINE) 
   
      RETURN
      END SUBROUTINE GETLINE 
C+----------------------------------------------------------------------
C
      SUBROUTINE GETXFORM (A, B, P, Q, SCALE, SHIFT)
C
C PURPOSE:
C     GETXFORM calculates the coefficients of the linear transformation
C
C                        X' <-- SCALE * X + SHIFT
C
C     such that X in interval [A, B] is transformed to X' in [P, Q].
C
C     USESCALE can be used to apply (or reverse) the transformation.
C
C     GETXFORM is a variation of GETSCALE for just the 1-D case.
C     It was introduced when the design of XFORMX was found to be
C     improper for the case of applying the same transformation to
C     more than one dataset.
C
C ARGUMENTS:
C   ARG   DIM  TYPE  I/O/S  DESCRIPTION 
C   A      -    R    I      Bounds of interval [A, B] to be transformed
C   B      -    R    I
C   P      -    R    I      Bounds of desired interval [P, Q]
C   Q      -    R    I
C   SCALE  -    R      O    Coefficients of the transformation
C   SHIFT  -    R      O
C
C HISTORY: 01/17/91  DAS  Ideas from XFORMX and GETSCALE combined for
C                         greater flexibility than XFORMX provided.
C
C AUTHOR:  David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C     Arguments:

      REAL      A, B, P, Q, SCALE, SHIFT

      SCALE = (Q - P) / (B - A)
      SHIFT = (B * P - A * Q) / (B - A)

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE HDESOL (MDIM, M, N1, A1, X, SSQMIN)
C
C  Purpose:
C
C     HDESOL (Householder DEcomposition and SOLution) solves one over-
C     determined system  Ax ~ b  where A is m x n with m >= n.  This is
C     the familiar linear least squares problem:  minimize the 2-norm of
C     Ax - b  with respect to x.
C
C  Notes:
C
C     1. If more than one right-hand-side b is to be used with the same
C        matrix A, use subroutines HDECOM, HSOLVE.
C     2. For the common polynomial case, see PNFIT and PNEVAL.
C     3. The direct factorization of A used here avoids the squaring of
C        the condition number that results from the alternative "normal
C        equations" approach, which solves the square system A'Ax = A'b.
C     4. No attempt is made here to handle rank-deficient cases.
C
C  Method:
C
C     A sequence of (orthogonal) Householder transformations (product =
C     Q) is used to triangularize matrix A one column at a time:
C                                                                T
C        Q   Q   ... Q Q A = R            so that           A = Q R
C         n-1 n-2     2 1
C
C     where R is a right (upper) triangular m x n matrix.
C
C     Inputting the right-hand-side b as the (n+1)th column makes for
C     convenient transformation of b during the factorization.  Then
C     the optimal solution is obtained from the upper n x n portion of
C     the transformed triangular system.
C
C  Arguments:
C
C     MDIM      is the declared row dimension of the array containing
C               matrix A.
C     M         is the number of equations in the system.  M >= N.
C     N1        is N+1 where N is the number of linear parameters X(1:N)
C               being computed.  N >= 1, so N1 >= 2.
C     A1(M,N1)  is the augmented matrix (A b), where A is the matrix
C               of coefficients in columns 1:N, and B is the RHS vector
C               in column N1 = N+1.  A is assumed to be full-rank.
C               A1 is destroyed upon return.
C     X(M)      is used for work space before being output with the
C               required linear parameters in X(1) through X(N).
C     SSQMIN    is output with the square of the  2-norm of the minimum
C               residual (minimum sum of squares).
C
C  Reference:
C
C     "Linear Least Squares Solutions by Householder Transformations,"
C     by Businger and Golub (1965).  Numer. Math. 7, pp 269-276.
C
C  Acknowledgement:
C
C     G.Golub, C.Moler, M.Saunders, Stanford University, 1972.
C
C     Jan. 1989 D.Saunders  Description clarified.
C     July 2000   "   "     Fortran 90 translation.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER, INTENT (IN)    :: MDIM, M, N1
      REAL,    INTENT (INOUT) :: A1(MDIM,N1), X(M)
      REAL,    INTENT (OUT)   :: SSQMIN

C  *  Local constants:

      REAL, PARAMETER :: ZERO = 0.

C  *  Local variables:

      INTEGER  I, J, K, N
      REAL     ALPHA, BETA, GAMMA, T1

C     Execution:

      N = N1 - 1

      DO K = 1, N1

C  *     Find the reflection which zeros A1(I,K), I = K+1, ... M.

         SSQMIN = ZERO

         IF (K < N1) THEN

            DO I = K, M
               X(I) = A1(I,K)
               SSQMIN = SSQMIN + X(I) ** 2
            END DO

         ELSE ! K = N1

            IF (M > N) THEN

               DO I = K, M
                  X(I) = A1(I,K)
                  SSQMIN = SSQMIN + X(I) ** 2
               END DO

            END IF ! SSQMIN is thus left with its required value.

            EXIT

         END IF

         ALPHA = SQRT (SSQMIN)
         T1 = X(K)
         IF (T1 < ZERO) ALPHA = -ALPHA
         T1 = T1 + ALPHA
         BETA = ALPHA * T1
         X(K) = T1

C  *     Apply the reflection to the remaining columns of A1.

         DO J = K + 1, N1

            GAMMA = ZERO
            DO I = K, M
               GAMMA = GAMMA + A1(I,J) * X(I)
            END DO
            GAMMA = GAMMA / BETA

            DO I = K, M
               A1(I,J) = A1(I,J) - GAMMA * X(I)
            END DO

         END DO

         X(K) = -ALPHA

      END DO

C  *  Solve the triangular system  Rx = <upper N part of transformed b>.

      X(N) = -A1(N,N1) / ALPHA

      IF (N > 1) THEN

         DO I = N - 1, 1, -1
            T1 = A1(I,N1)
            DO J = I + 1, N
               T1 = T1 - A1(I,J) * X(J)
            END DO
            X(I) = T1 / X(I)
         END DO

      END IF

      END SUBROUTINE HDESOL
C+----------------------------------------------------------------------
C
      SUBROUTINE HTDIS2 (DSINPUT, XA, XB, D1, D2, N, X, LUNOUT, IER)
C
C  ACRONYM:  Hyperbolic Tangent-type DIStribution, 2-sided, 2nd version
C            -          -            ---                    -
C  PURPOSE
C  -------
C
C     HTDIS2 generates abscissas on the interval [XA, XB] using the
C     method of Marcel Vinokur to achieve asymmetric bunching.  The
C     input D1 and D2 may be either the desired initial and final
C     INCREMENTS ("deltas") if DSINPUT is .TRUE., or the desired
C     end-point SLOPES of the (normalized) stretching function if
C     DSINPUT is .FALSE.
C
C  METHOD
C  ------
C
C     The mathematics appear in the paper
C
C        "On One-Dimensional Stretching Functions Functions for Finite-
C        Difference Calculations" by Marcel Vinokur,
C        Journal of Computational Physics, Vol. 50, No. 2, May 1983.
C
C     The paper develops criteria for stretching functions that optimize
C     in some sense the truncation errors inherent in finite-difference
C     approximations to derivatives.  The simplest "universal function"
C     (of which a scaled portion is employed) satisfying the criteria is
C     w = tan z, where  z = x + iy is complex.
C
C     The analysis uses the quantity B = SQRT (S0 * S1), where S0 and S1
C     are dimensionless end-point SLOPES of the stretching function.
C     I.e.,  S0 = dXI / dT  at  T = 0,  and  S1 = dXI / dT  at  T = 1,
C     where XI and T are normalized variables and XI = XI (T) is the
C     stretching function.
C
C     For the cases of usual interest, B > 1 and z is purely imaginary,
C     leading to relationships involving sinh and tanh - hence the name
C     HTDIS2.  However, for B < 1, z is real and the analogous solution
C     involves sine and tangent.  Furthermore, for B ~ 1, both formula-
C     tions break down, so a third relation is used involving B - 1.
C
C     In this implementation, a Newton iteration is used to solve
C
C        SINH (DEL Y) / (DEL Y) = B   if  B > 1, or
C
C        SIN (DEL X) / (DEL X)  = B   if  B < 1.
C
C     Then for each XI = (I - 1) / (N - 1) an ANTISYMMETRIC function is
C     calculated:
C
C        U = 0.5 + 0.5 TANH [DEL Y * (XI - 0.5)] / TANH (0.5 * DEL Y) or
C
C        U = 0.5 + 0.5 TAN  [DEL X * (XI - 0.5)] / TAN  (0.5 * DEL X)
C
C     from which the unsymmetric general case is derived as
C
C        T = U / [A + (1 - A) * U]    where  A = SQRT (S0 / S1).
C
C     For B ~ 1. (actually | B - 1. | < 1.E-5), the approximation is
C
C        U = XI * [1. + 2. * (B - 1.) * (XI - 0.5) * (1. - XI)]
C
C     with no nonlinear equation to solve.  B = 1 if the distribution
C     is uniform, and this is handled as a special case if D1 = D2.
C
C  CLARIFICATION OF WHAT D1 AND D2 REALLY MEAN
C  -------------------------------------------
C
C     Confusion has arisen from the tendency to specify first and last
C     INTERVALS (D1 and D2), while the method calls for specifying SLOPES
C     at the stretching function's end-points.  In fact, one is related
C     to the reciprocal of the other.  We have, for the end-point slope,
C
C        Slope ~ [XI(2) - XI(1)] / [T(2) - T(1)]
C              = [(2 - 1) / (N - 1)] / dT(1)
C              = [1 / (N - 1)] * [1 / dT(1)]
C
C     Vinokur's use of the end slopes was intended to enable blending
C     of two distributions end-on-end smoothly.  However, there may be
C     situations where specifying and achieving precise increments is
C     preferred (for instance if it is desired to prepare a table
C     showing the behavior of the distribution in terms of largest and
C     smallest increments for varying N or for fixed N and varying D1, D2).
C
C     Therefore, this version has an option to perform an outer iteration
C     (via the secant method) in order to achieve specified increments
C     more precisely.  DSINPUT = .TRUE. invokes the outer iteration.
C
C     [This outer iteration really has to solve a nonlinear equation
C     in TWO variables - the D1, D2 arguments to lower-level routine
C     HTDIS3.  But independent parallel iterations with the secant
C     (non-derivative) method - one for adjusting D1 and one for
C     adjusting D2 - appear to be reliable in this peculiar circumstance.
C     Extreme cases can cause these iterations to converge before
C     the normal tolerances are met, but the results should be close
C     anyway.  IER = 3 in these cases.]
C
C  ARGUMENTS
C  ---------
C
C   NAME   DIM   TYPE I/O/S DESCRIPTION
C
C   DSINPUT -     L     I   .TRUE. means D1 and D2 are the desired first
C                                  and last increments;
C                           .FALSE. means they are the desired end-point
C                                  slopes of the normalized stretching
C                                  function.
C   XA      -     R     I   Desired X(1); passing X(1) here is safe.
C   XB      -     R     I   Desired X(N); passing X(N) here is safe.
C                           Note that XA = 0. and XB = 1. can provide one
C                           normalized distribution for multiple reuse,
C                           in some cases.  XA > XB is required.
C   D1,     -     R     I   See the DSINPUT description.  D1 and D2 should
C   D2                      be positive in either mode.  If they represent
C                           slopes, they refer to the end-slopes of the
C                           curve formed by plotting (I - 1) / (N - 1)
C                           (vertical axis) vs. (X (I) - XA) / (XB - XA)
C                           (horizontal axis).  This curve is independent
C                           of N for given slopes.  A larger slope means
C                           a higher density of points at the corresponding
C                           end, in normalized space.
C   N       -     I     I   Number of points; N > 3.
C   X       N     R     O   Reqd. distribn., with X(1) = XA, X(N) = XB.
C   LUNOUT  -     I     I   LUNOUT > 0 shows the iteration history;
C                           LUNOUT < 0 suppresses it (but any error
C                                      message goes to |LUNOUT|).
C   IER     -     I     O   IER = 0 means no problems;
C                                 1 means bad combination of D1 and D2,
C                                   such as D1 + D2 > XB - XA;
C                                 2 means the inner Newton iteration
C                                   failed - must have been bad inputs.
C                                   X (*) is not usable in this case;
C                                 3 means the outer iteration by the
C                                   secant method did not meet the
C                                   normal tolerances, presumably
C                                   because the case was extreme,
C                                   but things stopped changing so
C                                   it did the best it could.
C                                   X (*) should still be usable;
C                                 4 means the outer iteration did not
C                                   converge after 20 passes.  X (*)
C                                   is not usable in this case.
C  ERROR HANDLING:
C
C     See LUNOUT and IER.  Failure handling is left to the calling program,
C     which can proceed or reprompt for parameters and try again.
C
C  INTERNAL PROCEDURE:  HTDIS3 (the original HTDIS2 with additional arguments
C                       for performing an efficient outer iteration).
C
C  ENVIRONMENT:  FORTRAN 90
C
C  HISTORY:
C  c. 1980   M.Vinokur   Original analysis.
C  04/01/89  J.E.Melton  Initial implementation of HTDIS2 (bisection
C                        used for the inner iteration).
C  05/01/89  R.G.Langhi  Introduced error handling; patterned after
C                        one-sided routine HTDIS.
C  08/10/89  D.Saunders  Internal arithmetic is now in DOUBLE PRECISION
C                        to help application to Navier-Stokes grids.
C                        Results prove to be very similar, with the
C                        bisection just less likely to fail.  It seems
C                        great precision in DEL does not help:  X(2) and
C                        X(N-1) are still strangely imprecise except for
C                        the uniform case (where DEL in SINGLE and DOUBLE
C                        are quite different but give similar X(I)s).
C  11/30/90  D.Saunders  Safeguarded the B <= 1 cases, which have no soln.
C  04/23/91   "    "     Introduced an outer iteration for precise D1, D2:
C                        pushed the main algorithm down a level with added
C                        arguments for more efficient solution of the
C                        nonlinear equation.  6 steps and a fixed starting
C                        guess for DEL seem to suffice for likely cases.
C  04/28/91   "    "     Introduced a Newton inner iteration.  (Bisection
C                        hardly benefits from good starting guesses.)
C  06/04/91   "    "     Laborious explanation of the misconception that
C                        had arisen from the apparently imprecise results.
C                        Retained the more-precise option by setting the
C                        number of outer iterations to 1 or 6 according
C                        to the signs of the input D1, D2.
C  08/20/91   "    "     Incorporated the B < 1 and B ~ 1 cases.
C  08/28/91   "    "     The small-correction outer iteration was found to
C                        diverge on extreme cases with small N.  Replaced
C                        it with the secant method (up to 20 iterations
C                        for each of D1 and D2 in parallel).  IER = 3 or
C                        4 are new possibilities.
C  09/10/91   "    "     Introduced the DSINPUT argument to clarify the
C                        "deltas" or "slopes" options.
C  03/03/95   "    "     Encountered a near-singular case (D1 * D2 ~ Du**2)
C                        during iteration, which improperly terminated.
C                        EPS = 1.E-5 instead of 1.E-3 in HTDIS3 helps, and
C                        the outer iteration may now continue for CASE = 3.
C                        Improved the second estimates derived from the
C                        first call to HTDIS3.
C  08/07/99   "    "     Fortran 90 upgrade: HTDIS3 is internal now, using
C                        just 2 arguments.  The other 10 become global.
C
C  AUTHORS: Marcel Vinokur, NASA/Ames Research Center, Moffett Field, CA
C           and John Melton, Ron Langhi, David Saunders
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IER, LUNOUT, N
      REAL
     >   D1, D2, XA, XB, X (N)
      LOGICAL
     >   DSINPUT

C     Local constants:

      INTEGER, PARAMETER ::
     >   MXITER = 20      ! Max. number of secant method iterations

      REAL, PARAMETER ::
     >   DLOW   = 1.E-2,  ! Smallest reduction in trial D1, D2 inputs to HTDIS3.
     >   ONE    = 1.E+0,
     >   THREE  = 3.E+0,
     >   TOL    = 1.E-4,  ! Relative tolerance on D1, D2
     >   ZERO   = 0.E+0

C     Local variables:

      INTEGER
     >   I

      REAL
     >   DEL0, DELSOLN, DMIN1, DMIN2, DX1 (3), DX2 (3), F1 (2), F2 (2),
     >   R, TOL1, TOL2, X1TARG, X2TARG

      LOGICAL
     >   CLOSE, TRIAL

C     Execution:
C     ----------

      TRIAL   = DSINPUT        ! Forces iteration to match Ds precisely
      DX1 (1) = D1
      DX2 (1) = D2

      IF (TRIAL) THEN
         TOL1   = DX1 (1) * TOL
         DMIN1  = DX1 (1) * DLOW
         X1TARG = XA + DX1 (1)      ! Target X (2)

         TOL2   = DX2 (1) * TOL
         DMIN2  = DX2 (1) * DLOW
         X2TARG = XB - DX2 (1)      ! Target X (N-1)

         IF (X1TARG > X2TARG) THEN  ! Bad D1 and D2
            IER = 1
            GO TO 99
         END IF
      END IF

      DEL0 = THREE  ! B > 1 gives DEL ~ 0.1 to 10 (increasingly non-uniform);
                    ! B < 1 gives DEL in the range [0+, PI-].
      CLOSE = .FALSE.     ! Used to detect "converged but above tolerance"

C     The first call to the basic Vinokur method is all we need if we are
C     matching slopes, else we need it to initialize the outer iteration:

C*****CALL HTDIS3 (DSINPUT, XA, XB, DX1 (1), DX2 (1), N, X, TRIAL, DEL0,
C    >             DELSOLN, LUNOUT, IER)

      CALL HTDIS3 (DX1 (1), DX2 (1)) ! Internal procedure now

      IF (IER /= 0) GO TO 99

      IF (.NOT. DSINPUT) GO TO 99    ! Slopes are matched without iteration.

      IF (DELSOLN == ZERO) GO TO 99  ! DEL = 0. signals uniform distribution
                                     ! which is specially handled.


C     Otherwise apply two secant iterations in parallel - one for D1 and
C     one for D2 - to solve what is really a nonlinear equation in TWO
C     variables.  Finish the first function evaluation:

      F1 (1) = X (2) - X1TARG      
      F2 (1) = X2TARG - X (N - 1)

C     Now for the second evaluation.  Attempt a better estimate.

      DX1 (2) = DX1 (1) ** 2 / (X (2) - XA)
      DX2 (2) = DX2 (1) ** 2 / (XB - X (N - 1))
      DEL0 = DELSOLN

      CALL HTDIS3 (DX1 (2), DX2 (2))

      IF (IER /= 0) GO TO 99

      F1 (2) = X (2) - X1TARG      
      F2 (2) = X2TARG - X (N - 1)

C     Here is the rest of the two secant iterations in parallel:

      DO I = 1, MXITER

         DEL0 = DELSOLN

         R = F1 (1) - F1 (2)
         IF (R /= ZERO) THEN
            R = F1 (1) / R
            DX1 (3) = R * DX1 (2) + (ONE - R) * DX1 (1)
            DX1 (3) = MAX (DX1 (3), DMIN1)
         ELSE
            DX1 (3) = DX1 (2)
         END IF

         R = F2 (1) - F2 (2)
         IF (R /= ZERO) THEN
            R = F2 (1) / R
            DX2 (3) = R * DX2 (2) + (ONE - R) * DX2 (1)
            DX2 (3) = MAX (DX2 (3), DMIN2)
         ELSE
            DX2 (3) = DX2 (2)
         END IF

         CALL HTDIS3 (DX1 (3), DX2 (3))

         IF (IER /= 0) GO TO 99

         IF (.NOT. TRIAL) THEN 
            IF (CLOSE) IER = 3
            GO TO 99
         END IF

C        Normal convergence test:
C        Reset TRIAL to get the full distribution on a final pass.

         F1 (1) = F1 (2)
         F1 (2) = X (2) - X1TARG      
         F2 (1) = F2 (2)
         F2 (2) = X2TARG - X (N - 1)

         IF (ABS (F1 (2)) < TOL1 .AND. ABS (F2 (2)) < TOL2) THEN
            TRIAL = .FALSE.
         ELSE IF (DX1 (3) == DX1 (2) .AND. DX2 (3) == DX2 (2)) THEN
            CLOSE = .TRUE.
            TRIAL = .FALSE.
         END IF

         DX1 (1) = DX1 (2)
         DX1 (2) = DX1 (3)
         DX2 (1) = DX2 (2)
         DX2 (2) = DX2 (3)

      END DO

      IER = 4  ! No convergence - X (*) is not usable.

   99 RETURN   ! From HTDIS2


C     Internal procedure for HTDIS2:
C     ------------------------------

      CONTAINS

!        -----------------------------------------------------------------------
!
         SUBROUTINE HTDIS3 (D1, D2)  ! The essence of the Vinokur method
!
!        -----------------------------------------------------------------------
!
!        ACRONYM:  Hyperbolic Tangent-type DIStribution, 2-sided (3rd version)
!                  -          -            ---                    -
!
!        HTDIS3 is the essence of the original HTDIS2.  It allows the specified
!        initial and final increments to be obtained with minimal work during
!        the intermediate iterations.
!
!        ADDITIONAL VARIABLES (See HTDIS2 for the others):
!
!        VAR      TYPE I/O/S DESCRIPTION
!
!        TRIAL     L     I   TRIAL = .TRUE. suppresses calculation of the
!                            interior X (*), since only X (2) and X (N-1)
!                            are being made use of at the higher level.
!        DEL0      R     I   Initial guess for the DELSOLN which solves
!                            the relevant nonlinear equation.  If B > 1,
!                            DEL is somewhere around 1. to 10. (larger
!                            for more nonuniform distributions); 1.E-10
!                            and 25. are used as bounds if necessary.
!                            If B < 1, the bounds are 1.E-10 and PI - eps.
!        DELSOLN   R     O   This solution is returned for reuse in
!                            a possible next call.
!
!        -----------------------------------------------------------------------

!        Arguments:

         REAL
     >      D1, D2

!        Local constants:

         INTEGER, PARAMETER ::
     >      MXITER = 30

         DOUBLE PRECISION, PARAMETER ::
     >      EPS = 1.D-5, HALF = 0.5D+0, ONE = 1.D+0,
     >      PIMINUS = 3.14159D+0, TOL = 1.D-12, TWO = 2.D+0,
     >      ZERO = 0.D+0

!        Local variables:

         INTEGER
     >      CASE, I, INC, ITER, LUNMSG

         DOUBLE PRECISION
     >      A, B, DEL, DELHI, DELLO, DELINV, DUNIFM, EXPDEL, FACTOR,
     >      FPRIME, FRACT, FUNDEL, RANGE, RFNM1, SLOPE0, SLOPE1, U, XI


!        Execution:
!        ----------

         IER = 0

!        Normalize the input control parameters:

         RANGE  = DBLE (XB - XA)
         RFNM1  = ONE / DBLE (N - 1)
         DUNIFM = RANGE * RFNM1
         SLOPE0 = DBLE (D1)
         SLOPE1 = DBLE (D2)

         IF (DSINPUT) THEN
            SLOPE0 = DUNIFM / SLOPE0
            SLOPE1 = DUNIFM / SLOPE1
         END IF

         A = SQRT (SLOPE0 / SLOPE1)
         B = SQRT (SLOPE0 * SLOPE1)

!        B distinguishes the cases.

         IF (B == ONE) THEN

            IF (SLOPE0 == SLOPE1) THEN  ! D1 = D2 = uniform DX
               DEL = ZERO               ! Signal no need to iterate
               DO I = 1, N - 1
                  X (I) = XA + REAL (I-1) * DUNIFM
               END DO
               GO TO 40
            ELSE
               CASE = 3                 ! Singular case, but not uniform
               DEL = -ONE               ! Something other than zero
            END IF

         ELSE IF (B > ONE + EPS) THEN   ! Need to solve SINH (DEL) / DEL = B

            CASE = 1
            DELHI = 25.D+0
            DELLO = 1.D-10

         ELSE IF (B < ONE - EPS) THEN   ! Need to solve SIN (DEL) / DEL = B

            CASE = 2
            DELHI = PIMINUS
            DELLO = EPS

         ELSE                           ! Simpler approximation

            CASE = 3
            DEL = -ONE

         END IF


         IF (CASE /= 3) THEN

!           Newton iteration for the two main cases:

            IF (LUNOUT > 0) WRITE (LUNOUT, 1001) CASE

            DEL = DBLE (DEL0)

            DO ITER = 1, MXITER

               IF (CASE == 1) THEN  ! B > 1 + EPS

!                 Solve   SINH (DEL) / DEL = B  in the form
!
!                 F (DEL) = SINH (DEL) / (DEL * B) - 1 = 0  (better scaled?)
!                 The exponential form is used since EXP (DEL) and its
!                 reciprocal should be more efficient than SINH and COSH.

                  EXPDEL = EXP (DEL)
                  DELINV = ONE / DEL
                  FRACT  = HALF * DELINV / B
                  FUNDEL = FRACT * (EXPDEL - ONE / EXPDEL) - ONE
                  FPRIME = FRACT * ((ONE - DELINV) * EXPDEL +
     >                              (ONE + DELINV) / EXPDEL)

               ELSE  ! CASE = 2; B < 1 - EPS
         
!                 Solve   SIN (DEL) / DEL = B  in the form
!
!                 F (DEL) = SIN (DEL) / (DEL * B) - 1 = 0  (better scaled?)

                  FUNDEL = SIN (DEL) / (DEL * B) - ONE
                  FPRIME = (COS (DEL) / B - FUNDEL - ONE) / DEL

               END IF

               IF (LUNOUT > 0)
     >            WRITE (LUNOUT, 1002) ITER, DEL, FUNDEL, FPRIME, TOL

               IF (ABS (FUNDEL) < TOL) GO TO 20

               DEL = MIN (MAX (DEL - FUNDEL / FPRIME, DELLO), DELHI)

            END DO

            LUNMSG = ABS (LUNOUT)
            WRITE (LUNMSG, 1003) MXITER, DEL, FUNDEL, TOL
            IER = 2
            GO TO 99

   20       CONTINUE

         END IF

!        Calculate the antisymmetric function, from which the desired
!        stretching function is derived by another simple transformation.
!        The three cases are becoming cumbersome...

         IF (CASE == 1) THEN

            FACTOR = ONE / TANH (HALF * DEL)

         ELSE IF (CASE == 2) THEN

            FACTOR = ONE / TAN (HALF * DEL)

         ELSE

            FACTOR = TWO * (B - ONE)

         END IF

         X (1) = XA

         IF (TRIAL) THEN       ! Only the 1st & last increments are looked at
            INC = MAX (1, N - 3)
         ELSE
            INC = 1
         END IF

         DO I = 2, N - 1, INC

            XI = DBLE (I-1) * RFNM1 - HALF  ! Subtracting 0.5 helps a bit

            IF (CASE == 1) THEN

               U = HALF * (ONE + TANH (DEL * XI) * FACTOR)

            ELSE IF (CASE == 2) THEN

               U = HALF * (ONE + TAN (DEL * XI) * FACTOR)

            ELSE

               U = (XI + HALF) * (ONE + FACTOR * XI * (HALF - XI))

            END IF

            FRACT = U / (A + (ONE - A) * U)
            X (I) = XA + REAL (FRACT * RANGE)

         END DO

   40    X (N) = XB
         DELSOLN = REAL (DEL)

         IF (LUNOUT > 0)
     >      WRITE (LUNOUT, 1004) X (2) - X (1), X (N) - X (N - 1)

   99    RETURN

 1001    FORMAT (/, ' ITER (Case = ', I1, ')    DEL', 7X, 'F(DEL)', 11X,
     >           'F''', 10X, 'TOL')
 1002    FORMAT (1X, I2, 1P, E20.12, 3E13.5)
 1003    FORMAT (/' *** HTDIS2:  Newton iteration failed ***',
     >           /' Number of iterations:', I4,
     >           /' DEL, F(DEL), TOL:', 1P, 3E14.6,
     >           /' No distribution computed.'/)
 1004    FORMAT (/, ' dX(1): ', 1P, E14.6, '   dX(N-1): ', E14.6)

         END SUBROUTINE HTDIS3

      END SUBROUTINE HTDIS2
C+----------------------------------------------------------------------
C
      SUBROUTINE INTERVAL (NX, X, XFIND, ARROW, LEFT)
C
C     One-liner: Interpolation search for interval containing a point.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Written primarily for interval-based interpolations such as
C     piecewise linear or cubic spline, INTERVAL performs a search to
C     locate the best interval for evaluating the interpolant at a
C     given point. The normal case returns the "left-hand" endpoint of
C     the interval bracketing the point, but for the out-of-range cases
C     below or above the range of the knots, the interval to be used is
C     the first or last. The array of knots must be monotonic, either
C     increasing or decreasing. Diagrammatically, LEFT is returned as
C     shown below for the normal case (no extrapolation):
C
C          X (1)  ...   X (LEFT)   X (LEFT+1)   ...      X (NX)
C                                ^
C                              XFIND
C
C     And for extrapolation:
C
C                     X (LEFT = 1)  ...   X (NX)
C             ^
C           XFIND
C
C     or,
C                X (1)  ...   X (LEFT = NX-1)    X (NX)
C                                                           ^
C                                                         XFIND
C
C     If the point to be bracketed (XFIND) matches one of the knots, the
C     index of that knot is returned as LEFT, i.e., the condition for a
C     bracket of an interior point is:
C
C        X (LEFT) <= XFIND < X (LEFT+1)  if  ARROW = +1.0,  or
C        X (LEFT) >= XFIND > X (LEFT+1)  if  ARROW = -1.0.
C
C        This is a low-level routine with minimal error checking. The
C     calling program is assumed to have verified the following:
C
C     (1)  NX >= 2
C     (2)  X strictly monotonic
C     (3)  ARROW = +1.0 or -1.0
C
C     Subroutine PROTECT is available from the author for easily checking
C     conditions (2) and (3). LEFT is verified on input, but efficiency in
C     loops will benefit from passing the best estimate available, usually
C     just the result of the last call.
C
C        INTERVAL was originally written for use with CSEVAL and TABLE1.
C     The interpolation search was adapted from ideas in Sedgewick's book
C     referenced below.
C
C     Arguments:
C     ----------
C
C     Name  Dimension  Type  I/O/S  Description
C     NX                I    I      Number of points in array X; must
C                                   be >= 2 (no check performed).
C
C     X        NX       R    I      Array of points defining the set
C                                   of intervals to be examined. Only
C                                   the first NX-1 points are required.
C
C     XFIND             R    I      The point for which a bracketing
C                                   interval is sought.
C
C     ARROW             R    I      Monotonicity indicator for input
C                                   array X:
C                                     -1.0  strictly decreasing
C                                      0.0  NOT ALLOWED!
C                                     +1.0  strictly increasing
C                                   Supplied by the calling routine for
C                                   reasons of speed (not checked).
C
C     LEFT              I    I/O    Input: guessed index of left-hand
C                                   endpoint of the interval containing
C                                   the specified point.
C
C                                   Output: index of the largest array
C                                   value <= specified point (if ARROW=+1.0).
C                                   Special case for data out of range:
C                                   return left endpoint of closest interval.
C                                   Thus, LEFT = 1 for XFIND < X (2), and
C                                   LEFT = NX-1 for XFIND >= X (NX-1).
C                                   (If ARROW=-1.0, reverse the inequalities.)
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and eight character symbols are not (yet) standard.
C
C     (2)  In speed-critical applications, it might be a good idea to build
C          this algorithm in-line since it is typically called many times
C          from within a loop. Another potential speed-up is removal of the
C          ARROW multiplies, which restricts the method to increasing data.
C          So far, the simplicity of separating out the messy search details
C          and the generality of bi-directional searching have outweighed
C          the modest speed penalty incurred.
C
C     Bibliography:
C     -------------
C
C     (1) Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983.
C            (Chap. 14)
C
C     Author:  Robert Kennelly and David Saunders, Sterling Federal Systems
C     -------
C
C     Development history:
C     --------------------
C
C     20 Oct. 1987    RAK    Interpolation search adapted (with mods.
C                            for bidirectional search and some minor
C                            repair) from CSEVAL (RAK) and TABLE1 (DAS).
C     08 Aug. 1988    DAS    Clarified descriptions of bracketing, where
C                            the inequalities depend upon ARROW.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   LEFT, NX
      REAL
     &   ARROW, X (NX), XFIND

C     Local variables.

      INTEGER
     &   LENGTH, NXLESS1, RIGHT, TRIAL
      REAL
     &   XBYARROW

C     Execution.
C     ----------

      XBYARROW = XFIND * ARROW

C     Simplify things by disposing of two important special cases so that
C     X (LEFT) and X (RIGHT) can really bracket XFIND. As a by-product,
C     this also takes care of the NX = 2, 3 cases.

      NXLESS1 = NX - 1

      IF (XBYARROW .GE. X (NXLESS1) * ARROW) THEN
         LEFT = NXLESS1
         GO TO 990
      ELSE IF (XBYARROW .LT. X (2) * ARROW) THEN
         LEFT = 1
         GO TO 990
      END IF

C      ---------------------------------
C     |                                 |
C     |   X (2) <= XFIND < X (NX - 1)   |
C     |            - or -               |
C     |   X (2) > XFIND >= X (NX - 1)   |
C     |                                 |
C     |   NX > 3                        |
C     |                                 |
C      ---------------------------------

C     Adjust the pointers. We hope that the calling routine has provided
C     a reasonable guess (since it's probably working on an ordered array
C     of points to evaluate), but check anyway.

      LEFT = MIN (MAX (2, LEFT), NX - 2)

      IF (XBYARROW .GE. X (LEFT) * ARROW) THEN
         IF (XBYARROW .LT. X (LEFT + 1) * ARROW) THEN

C           XFIND is in the original guessed-at interval.

            GO TO 990
         ELSE

C           We'll look farther to the right. Evidently LEFT was < NX - 2.

            RIGHT = NXLESS1
            LEFT  = LEFT + 1
         END IF
      ELSE

C        Look to the left of the guess. Evidently LEFT was > 2.

         RIGHT = LEFT
         LEFT  = 2
      END IF

C      ----------------------------------
C     |                                  |
C     |   2 <= LEFT < RIGHT <= NX - 1    |
C     |                                  |
C      ----------------------------------

C     The interval length must decrease each time through - terminate
C     when the correct interval is found or when the interval length
C     cannot be decreased.

   10 CONTINUE
         LENGTH = RIGHT - LEFT
         IF (LENGTH .GT. 1) THEN

C           The trial value is a "linear" estimate of the left-hand endpoint
C           of the interval bracketing the target XFIND, with protection
C           against round-off (which can affect convergence).

            TRIAL = MIN (RIGHT - 1, LEFT + MAX (0, INT (REAL (LENGTH) *
     &         (XFIND - X (LEFT)) / (X (RIGHT) - X (LEFT)))))

C            ------------------------------------------
C           |                                          |
C           |   2 <= LEFT <= TRIAL < RIGHT <= NX - 1   |
C           |                                          |
C            ------------------------------------------

C           Adjust pointers. Increase LEFT or decrease RIGHT until done.

            IF (XBYARROW .GE. X (TRIAL + 1) * ARROW) THEN
               LEFT  = TRIAL + 1
            ELSE IF (XBYARROW .LT. X (TRIAL) * ARROW) THEN
               RIGHT = TRIAL
            ELSE

C              We're done: XFIND is in the interval [X (TRIAL), X (TRIAL+1)).

               LEFT  = TRIAL
               GO TO 990
            END IF
            GO TO 10

         END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE INTSEC2 (N1, X1, Y1, T1, I1, CALCT1,
     >                    N2, X2, Y2, T2, I2, CALCT2, TOL,
     >                    XINT, YINT, LUNOUT, IER)
C
C ONE-LINER: INTerSECtion of two curves in 2-space, allowing extrapolation
C
C PURPOSE:
C
C        INTSEC2 determines where two curves in 2-space meet.  The curves
C     are defined as discrete points, and extrapolation of either curve or
C     both is permitted.  It was prompted by a need to sharpen the nose
C     of a rounded airfoil.  The earlier INTSEC does not allow extrapolation,
C     and is less precise in its use of linear interpolation.
C
C        This version allows the calling program to perform the arc-length
C     computations if that is more efficient (as it would be for multiple
C     calls where one of the curves is not changing).
C
C METHOD:
C
C        Each curve is represented parametrically by "local" splines
C
C           x1 = x1 (t1), y1 = y1 (t1) and x2 = x2 (t2), y2 = y2 (t2).
C
C     (Local spline coefficients may be calculated as needed from local data,
C     thus avoiding storage problems.  They are normally quite adequate for
C     interpolations.)
C
C        The problem then becomes that of solving a nonlinear system of
C     two equations in two variables:
C
C           f1 (t1, t2) = x1 (t1) - x2 (t2) = 0    and
C           f2 (t1, t2) = y1 (t1) - y2 (t2) = 0
C
C     which is solved by a safeguarded Newton iteration.  Note that the
C     partial derivatives of f1, f2, are functions of the derivatives of
C     x1, y1, x2, y2 which are readily determined by the spline utility.
C     Safeguarding involves halving each pure Newton step until || f ||
C     is reduced.
C
C        The starting guess is determined from the input values of I1 and
C     I2 (pointers to the data point on each curve as near as is known to
C     the point of intersection).  If 0 is entered for I1 or I2, then the
C     corresponding middle point is used for a starting guess.
C
C        The parametric technique avoids difficulties with vertical lines
C     which would affect INTSEC.  Parallel lines may still cause failure,
C     and this should be the only situation that could cause the maximum
C     iteration count to be reached.  Tests for convergence, however, are
C     subject to the choice of tolerance, which is dependent upon the data
C     scaling.  The input TOL is used relative to the average of the total
C     arc lengths of the two curves.
C
C ENVIRONMENT:  FORTRAN 77, with minor extensions
C
C HISTORY:
C     10/08/91  D.A.Saunders  Initial implementation, for sharpening airfoils.
C     07/31/93       "        Substituted LUSOLVE for DECOMP & SOLVE, and
C                             added the CALCT1, CALCT2 arguments.
C     12/08/93       "        Substituted CHORDS2D for CHORD; made TOL an
C                             argument.
C     02/03/98       "        Changed FNORM .GT. FNORM0 test to .GE.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER N1, N2           ! (I) Number of points in each curve, >= 2.

      REAL    X1 (N1), Y1 (N1) ! (I) Coordinates defining the two curves.
      REAL    X2 (N2), Y2 (N2)

      REAL    T1 (N1), T2 (N2) ! (I/S) Cumulative chord lengths, either input
                               !       or computed here - see CALCT1, CALCT2.

      INTEGER I1, I2           ! (I) Estimates of indices near the point of
                               !     intersection.  If 0, middle point is used.
                               ! (O) Actual indices nearest to the point of
                               !     intersection in INTERVAL's sense (that
                               !     is, between 1 and N - 1 inclusive).

      LOGICAL CALCT1, CALCT2   ! (I) .FALSE. means the calling program has
                               !     supplied T1 (*) or T2 (*).

      REAL    TOL              ! (I) Tolerance on || f || used relative to
                               !     the average curve length.  Try 1.E-6 or
                               !     1.E-14 for 32- or 64-bit arithmetic.

      REAL    XINT, YINT       ! (O) Estimated point of intersection.

      INTEGER LUNOUT           ! (I) Logical unit for displaying iterations,
                               !     which are suppressed if LUNOUT < 0.
                               !     |LUNOUT| is used for error messages.

      INTEGER IER              ! (O) 0 means no problem was encountered;
                               !     1 means the safeguarded iteration failed
                               !       to satisfy TOL but did satisfy 10*TOL;
                               !     2 means the iteration limit was reached
                               !       without (close) convergence;
                               !     3 means the step-halving failed - somehow
                               !       the iteration diverged;
                               !     4 means the linearized system was singular:
                               !       the lines were probably parallel.
C     Procedures.

      EXTERNAL  CHORDS2D       ! Chord-length utility.
      EXTERNAL  INTERVAL       ! Search utility, also used by LCSFIT.
      EXTERNAL  LCSFIT         ! Local spline utility.  "Loose" fits are used.
      EXTERNAL  LUSOLVE        ! Square system solver (LU decomposition).

C-----------------------------------------------------------------------

C     Local constants.

      INTEGER   ITMAX
      REAL      HALF, ONE, STEPMIN, ZERO
      LOGICAL   NEW
      CHARACTER METHOD * 1
      PARAMETER
     >  (ITMAX   = 10,          ! Outer iteration limit.  (Halving the step
     >   HALF    = 0.5E+0,      !   is the inner iteration.)
     >   METHOD  = 'B',         ! "Bessel" = "loose" spline fits.
     >   NEW     = .TRUE.,      ! Always alternating between X1, Y1, etc.
     >   ONE     = 1.E+0,
     >   STEPMIN = 1.E-3,
     >   ZERO    = 0.E+0)

C     Local variables.

      INTEGER   I, ITER, LUNERR
      REAL      A (2, 2), DT (2), F (2), T (2), TLAST (2),
     >          X (2), XP (2), Y (2), YP (2)
      REAL      ALPHA, FNORM, FNORM0, TOLER, TTOTAL1, TTOTAL2

C     Execution.

      LUNERR = ABS (LUNOUT)
      IF (I1 .EQ. 0) I1 = (N1 + 1) / 2
      IF (I2 .EQ. 0) I2 = (N2 + 1) / 2

C     Set up the cumulative chord lengths unless they are supplied
C     by the higher level for efficiency reasons:

      IF (CALCT1) THEN
         CALL CHORDS2D (N1, X1, Y1, .FALSE., TTOTAL1, T1)
      ELSE
         TTOTAL1 = T1 (N1)
      END IF

      IF (CALCT2) THEN
         CALL CHORDS2D (N2, X2, Y2, .FALSE., TTOTAL2, T2)
      ELSE
         TTOTAL2 = T2 (N2)
      END IF

C     Perform a safeguarded Newton iteration to solve the 2x2 nonlinear system.

      IER = 0
      ITER = 0
      T (1) = T1 (I1)
      T (2) = T2 (I2)      
      ALPHA = ZERO            ! For printout purposes only
      FNORM0 = 1.E+20         ! I.e., big to avoid iteration 0 test
      TOLER = TOL * HALF * (TTOTAL1 + TTOTAL2)  ! Scale tolerance by avg. length

  300 CONTINUE

C        Evaluate x1, dx1/dt, etc., at the current t = (t(1), t(2)).

         CALL LCSFIT (N1, T1, X1, NEW, METHOD, 1, T (1), X (1), XP (1))
         CALL LCSFIT (N1, T1, Y1, NEW, METHOD, 1, T (1), Y (1), YP (1))
         CALL LCSFIT (N2, T2, X2, NEW, METHOD, 1, T (2), X (2), XP (2))
         CALL LCSFIT (N2, T2, Y2, NEW, METHOD, 1, T (2), Y (2), YP (2))

C        Find the norm of the residual vector f, which should converge to ~0.
C        Set up the RHS vector for the system J dT = f in the process.

         F (1) = X (1) - X (2)
         FNORM = ABS (F (1))
         F (2) = Y (1) - Y (2)
         FNORM = MAX (FNORM, ABS (F (2)))

C        Halve the step until || f || is reduced (except first time through).

         IF (FNORM .GE. FNORM0) THEN
            IF (ALPHA .GT. STEPMIN) THEN
               ALPHA = HALF * ALPHA
               T (1) = TLAST (1) - ALPHA * DT (1)
               T (2) = TLAST (2) - ALPHA * DT (2)
               GO TO 300
            END IF
            GO TO 830
         END IF

         IF (LUNOUT .GT. 0) THEN
            WRITE (LUNOUT, '(A, I3, A, 1P, E13.6, A, E9.2, A, 2E14.6)')
     >         ' INTSEC2:', ITER, '  ||f||:', FNORM, '  step:', ALPHA,
     >         '  t:', T
         END IF

         IF (FNORM .LT. TOLER) GO TO 900  ! Success


         ITER = ITER + 1
         IF (ITER .EQ. ITMAX) GO TO 810

         FNORM0 = FNORM

C        Set up the LHS matrix for the next iteration.

         A (1, 1) = XP (1)
         A (2, 1) = YP (1)
         A (1, 2) = -XP (2)
         A (2, 2) = -YP (2)

C        Solve  J dt = f; dt overwrites f.

         CALL LUSOLVE (2, 2, A, F, IER)
         IF (IER .NE. 0) GO TO 840

         DO 500, I = 1, 2         
            DT (I) = F (I)
            TLAST (I) = T (I)
            T (I) = T (I) - DT (I)
  500    CONTINUE

         ALPHA = ONE

      GO TO 300                           ! Do another iteration


C     Error handling.

  810 IER = 1
      WRITE (LUNERR, 1010) 'Iteration limit reached.'
      IF (FNORM .LT. 10. * TOLER) GO TO 900

      IER = 2
      GO TO 999

  830 IER = 3
      WRITE (LUNERR, 1010) 'Step-halving iteration failed.'
      GO TO 999

  840 IER = 4
      WRITE (LUNERR, 1010) 'Singular system. Parallel lines?'
      GO TO 999
      

  900 CONTINUE
C     Wrap up, having converged or almost converged.

      XINT = HALF * (X (1) + X (2))
      YINT = HALF * (Y (1) + Y (2))

      CALL INTERVAL (N1, T1, T (1), ONE, I1)
      CALL INTERVAL (N2, T2, T (2), ONE, I2)

  999 RETURN

C     Formats.

 1010 FORMAT ('0INTSEC2: ', A)

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE LCSFIT (NDATA, X, Y, NEW, METHOD, NEVAL, XEVAL, YEVAL,
     &   YPEVAL)
C
C     Two-liner:  Storage-efficient local cubic spline fit (2-space)
C     ----------  (monotonic and piecewise linear options too)
C
C     Description and usage:
C     ----------------------
C
C        LCSFIT is the non-parametric analog of PLSFIT (parametric).
C     It is intended for spline applications which do not require the
C     spline coefficients as output.  It is efficient for repeated
C     calls with the same data, so repeated use with NEVAL = 1 may be
C     preferable to storing vectors of results.
C
C        LCSFIT offers monotonic spline and piecewise linear options
C     also.  And it returns an interpolated first derivative along
C     with the function value.  (The second derivative is omitted
C     because Y" is not guaranteed to be continuous by local methods.)
C
C        See PLSFIT for more details on local methods.  As with most
C     numerical methods, scaling of the data to the unit interval (and
C     unscaling of the result) is recommended to avoid unnecessary
C     effects attributable to the data units.  Utilities GETSCALE and
C     USESCALE from the present authors are appropriate.  The data
C     abscissas should be distinct and either ascending or descending.
C     PROTECT is available to check this.  Extrapolation is permitted
C     (mainly in case of round-off; it is normally inadvisable).
C
C        The CSFIT/CSEVAL or CSDVAL pair are probably preferable if
C     efficiency is not an issue, since CSFIT gives Y" continuity.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     NDATA   I               I      Length of X, Y input data arrays.
C
C     X,      R (NDATA)       I      Input data coordinates.  The Xs
C     Y                              must be distinct and monotonic,
C                                    either ascending or descending.
C                                    (No check here.) 
C
C     NEW     L               I      If control flag NEW is .TRUE., the
C                                    search for a bracket starts from
C                                    scratch, otherwise locally-saved
C                                    search and fit information will be
C                                    assumed to be correct. If calling
C                                    LCSFIT from within a loop, set
C                                    NEW = .FALSE. after the first call.
C
C     METHOD   C*1            I      (Uppercase) Type of fit to be used:
C                                    'M' means Monotonic piecewise cubics;
C                                    'B' means non-monotonic "Bessel"-type
C                                        piecewise cubics (looser fit);
C                                    'L' means piecewise Linear fit;
C                                    'C' means Cyclic (periodic) end
C                                        conditions: loose fit assumed.
C
C     NEVAL   I               I      Number of interpolations requested.
C                                    NEVAL >= 1.  One call per result
C                                    (NEVAL = 1) may save storage, and is
C                                    not too inefficient as long as NEW
C                                    is set to .FALSE. after the first.
C
C     XEVAL   R (NEVAL)       I      Abscissa(s) to interpolate to.  These
C                                    are normally in the data range, but
C                                    extrapolation - probably due to
C                                    round-off - is not prevented.
C
C     YEVAL   R (NEVAL)       O      Interpolated function value(s).
C
C     YPEVAL  R (NEVAL)       O      Interpolated 1st derivative value(s).
C                                    Pass the same storage as for YEVAL
C                                    if no derivatives are required.
C
C     Significant local variables:
C     ----------------------------
C
C     MEMORY         Indicates that coefficients are correct for the
C                    current point.
C
C     H, DEL         Delta X and forward difference derivative arrays.
C
C     B, C, D        Coefficients of cubic on the bracketing interval.
C
C     Procedures:
C     -----------
C
C     INTERVAL  1-D "interpolation" search.
C     BESSEL    First derivative (central 3-point formula).
C     BRODLIE   First derivative (central), adjusted for monotonicity.
C     BUTLAND   First derivative (non-central), adjusted for monotonicity.
C     THREEPT   First derivative (non-central 3-point formula).
C
C     Environment:  FORTRAN 90
C     ------------
C
C     Error handling:  None
C     ---------------
C
C     Notes:
C     ------
C
C     (1)  Since many of the calculations must be repeated at both ends
C          of an interval, the various finite difference quantities used
C          are stored as arrays. The following "map" of a typical interior
C          interval and its neighbors should help in understanding the
C          notation.  The local array indices are all numbered relative
C          to the left-hand end of the interval which brackets the point
C          to be evaluated.
C
C                                  LEFT       RIGHT
C
C          Point         -1          0         +1          +2
C
C          Data           x -------- x -------- x --------- x
C
C          Interval      -1          0         +1
C
C
C     Author: Robert Kennelly, Sterling Software/NASA Ames  (PLSFIT)
C     -------
C
C     History:
C     --------
C
C     27 Feb. 1987  R.A.Kennelly  Initial implementation of PLSFIT.
C     23 Aug. 1989  D.A.Saunders  LCSFIT adapted as non-parametric form,
C                                 for embedding in other utilities where
C                                 minimizing work-space is desirable.
C     20 June 1991    "    "      THREEPT (monotonic) renamed BUTLAND;
C                                 THREEPT (pure 3-pt. formula) now used
C                                 for nonmonotonic end-point handling;
C                                 METHOD='C' case belatedly added, as
C                                 needed by PLSINTRP for closed curves.
C     23 July 1991    "    "      The tests for being in the same interval
C                                 as before were not allowing for the
C                                 descending-Xs case.
C     06 May  1998    "    "      Minor Fortran 90 updates.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER,   INTENT (IN)  :: NDATA, NEVAL
      REAL,      INTENT (IN)  :: X (NDATA), Y (NDATA), XEVAL (NEVAL)
      REAL,      INTENT (OUT) :: YEVAL (NEVAL), YPEVAL (NEVAL)
      LOGICAL,   INTENT (IN)  :: NEW
      CHARACTER, INTENT (IN)  :: METHOD * 1

C     Local constants:

      REAL, PARAMETER :: ZERO = 0., ONE = 1., TWO = 2., THREE = 3.

C     Local variables:

      LOGICAL
     &   CYCLIC, LINEAR, MEMORY, MONO
      INTEGER
     &   IEVAL, J, K, LEFT, RIGHT
      REAL
     &   ARROW, B (0:1), C, DELY (-1:1), D, DX, H (-1:1), XBYARROW, XE

C     Procedures:

      REAL, EXTERNAL ::
     &   BESSEL, BRODLIE, BUTLAND, THREEPT

C     Storage:

      SAVE
     &   ARROW, B, C, D, LEFT, RIGHT

C     Execution:
C     ----------

      MONO   = METHOD == 'M'
      CYCLIC = METHOD == 'C'
      LINEAR = METHOD == 'L'

      IF (CYCLIC) THEN
         IF (Y (NDATA) /= Y (1)) STOP 'LCSFIT: End points must match.'
      END IF

C     Initialize search or avoid it if possible:

      IF (NEW) THEN
         MEMORY = .FALSE.
         ARROW  = SIGN (ONE, X (2) - X (1))
         LEFT   = 1
      END IF

      IEVAL = 1
      XE = XEVAL (1)
      XBYARROW = XE * ARROW

      IF (.NOT. NEW) THEN
      
C        We can save a lot of time when LCSFIT is being called from within
C        a loop by setting MEMORY if possible. The out-of-range checking
C        relies on the fact that RIGHT = LEFT + 1 after the last return.
C        Cater to the more likely case of XE in the previous, interior
C        interval.

         MEMORY = XBYARROW >= X (LEFT)  * ARROW .AND.
     &            XBYARROW <  X (RIGHT) * ARROW

         IF (.NOT. MEMORY) THEN
            MEMORY =
     &         LEFT  == 1     .AND. XBYARROW <  X (RIGHT) * ARROW
     &         .OR.
     &         RIGHT == NDATA .AND. XBYARROW >= X (LEFT)  * ARROW
         END IF

      END IF

      IF (MEMORY) GO TO 70 ! Skip the bulk of the computation

C     Loop over evaluation points requiring a new search:
C     ---------------------------------------------------

   10 CONTINUE

C        Interpolation search for bracketing interval:
C        ---------------------------------------------

         CALL INTERVAL (NDATA, X, XE, ARROW, LEFT)

         RIGHT = LEFT + 1

C         -------------------------------------------
C        |                                           |
C        |   1 <= LEFT < RIGHT = LEFT + 1 <= NDATA   |
C        |                                           |
C         -------------------------------------------

C        Compute derivatives by finite-differences:
C        ------------------------------------------

         IF (NDATA > 2 .AND. .NOT. LINEAR) THEN

C           Interval and derivative approximations:
C           ---------------------------------------

C           The following duplicates more code than PLSFIT's approach,
C           but eliminates some indirection - no need to wrap-around here.
C           Handle the end conditions first to minimize testing LEFT, RIGHT.

            IF (LEFT == 1) THEN

               H (0) = X (2) - X (1)
               DELY (0) = (Y (2) - Y (1)) / H (0)
               H (1) = X (3) - X (2)
               DELY (1) = (Y (3) - Y (2)) / H (1)

               IF (CYCLIC) THEN ! Loose fit assumed
                  H (-1) = X (NDATA) - X (NDATA - 1)
                  DELY (-1) = (Y (NDATA) - Y (NDATA - 1)) / H (-1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE
                  IF (MONO) THEN
                     B (0) = BUTLAND (0, H, DELY)
                     B (1) = BRODLIE (1, H, DELY)
                  ELSE
                     B (0) = THREEPT (0, H, DELY)
                     B (1) = BESSEL  (1, H, DELY)
                  END IF
               END IF

            ELSE IF (RIGHT == NDATA) THEN

               H(-1) = X (LEFT) - X (LEFT-1)
               DELY(-1) = (Y (LEFT) - Y (LEFT-1)) / H (-1)
               H (0) = X (RIGHT) - X (LEFT)
               DELY (0) = (Y (RIGHT) - Y (LEFT))  / H (0)

               IF (CYCLIC) THEN
                  H (1) = X (2) - X (1)
                  DELY (1) = (Y (2) - Y (1)) / H (1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE

                  IF (MONO) THEN
                     B (0) = BRODLIE (0, H, DELY)
                     B (1) = BUTLAND (1, H, DELY)
                  ELSE
                     B (0) = BESSEL  (0, H, DELY)
                     B (1) = THREEPT (1, H, DELY)
                  END IF
               END IF

            ELSE

               K = LEFT
               DO J = -1, +1
                  H (J)    =  X (K) - X (K-1)
                  DELY (J) = (Y (K) - Y (K-1)) / H (J)
                  K = K + 1
               END DO

C              Select interpolation scheme:
C              ----------------------------

C              Compute (possibly adjusted) first derivatives at both
C              left- and right-hand endpoints of the interval.

               IF (MONO) THEN

C                 Monotone - use Brodlie modification of Butland's
C                 formula to adjust the derivatives at the knots.

                  B (0) = BRODLIE (0, H, DELY)
                  B (1) = BRODLIE (1, H, DELY)

               ELSE ! METHOD = 'B'

C                 Bessel - use central difference formula at the knots.

                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)

               END IF

            END IF

C           Compute the remaining cubic coefficients.

            C = (THREE * DELY (0) - TWO * B (0) - B (1)) / H (0)
            D = ( -TWO * DELY (0) + B (0) + B (1)) / H (0) ** 2

         ELSE ! NDATA = 2 .OR. METHOD = 'L'

C           Degenerate case (linear).
C           -------------------------

            B (0) = (Y (RIGHT) - Y (LEFT)) / (X (RIGHT) - X (LEFT))
            C     = ZERO
            D     = ZERO

         END IF

C        Evaluate the cubic (derivative first in case only YEVAL is reqd.):
C        ------------------------------------------------------------------

   70    CONTINUE ! Start of same-interval loop inside new-interval loop

            DX = XE - X (LEFT)
            YPEVAL (IEVAL) = B (0) + DX * (TWO * C + DX * THREE * D)
            YEVAL  (IEVAL) = Y (LEFT) + DX * (B (0) + DX * (C + DX * D))

C           The next evaluation point may be in the same interval:
C           ------------------------------------------------------

            IF (IEVAL < NEVAL) THEN ! Skips this if NEVAL = 1

               IEVAL = IEVAL + 1
               XE = XEVAL (IEVAL)
               XBYARROW = XE * ARROW
               IF (XBYARROW >= X (LEFT)  * ARROW  .AND.
     &             XBYARROW <  X (RIGHT) * ARROW) GO TO 70
            
               GO TO 10 ! Else much more work required.

            END IF

C     Termination:
C     ------------

      RETURN

      END SUBROUTINE LCSFIT
C+----------------------------------------------------------------------
C
      SUBROUTINE LINTRP ( NY, YM, YP, XM, XP, X, Y )
C
C PURPOSE:  LINTRP  interpolates linearly for all elements of the given
C           parallel arrays YM(*), YP(*), into the third parallel array
C           Y(*),  according to the abscissas XM, XP, X that are common
C           to all elements of the corresponding array...   (Originally
C           developed for interpolating wing pressure data from comput-
C           ational span stations to arbitrary span stations.)
C
C PARAMETERS:
C   ARG   DIM   TYPE I/O/S DESCRIPTION
C   NY     -      I    I   No. of interpolated elements required (>=1).
C   YM,YP  NY     R    I   Parallel arrays between which Y(*) is reqd.
C   XM,XP  -      R    I   Abscissas applying to YM(*), YP(*) resp.
C   X      -      R    I   Abscissa applying to desired array Y(*).
C                          Normally XM < X < XP, but extrapolation is
C                          permissible with this code.
C   Y      NY     R    O   Interpolated value(s).
C    
C ENVIRONMENT: FORTRAN IV
C
C AUTHOR: David Saunders, Informatics, Palo Alto, CA.   (07/30/82)
C
C-----------------------------------------------------------------------
C
      DIMENSION  YM(NY), YP(NY), Y(NY)
C
      RATIO = (X-XM) / (XP-XM)
C
      DO 20 I = 1, NY
         Y(I) = YM(I) + (YP(I)-YM(I))*RATIO
 20   CONTINUE
C
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE LOOKUP (NDICT, DICTRY, ALPHA, KEY, ENTRY)
C
C
C     Description and usage:
C
C           Performs dictionary lookups.  A pointer is returned if a
C        match is found between the input key and the corresponding
C        initial characters of one of the elements of the dictionary.
C        If a "synonym" has been provided for an entry, the search is
C        continued until a match to a primary dictionary entry is found.
C        Cases of no match, or multiple matches, are also provided for.
C
C           Dictionary entries must be left-justified, and may be alphabetized
C        for faster searches.  Secondary entries, if any, are composed of
C        two words separated by one or more characters such as blank, tab,
C        comma, colon, or equal sign which are treated as non-significant
C        by SCANNR.  The first entry of each such pair serves as a synonym
C        for the second, more fundamental keyword.
C
C           The ordered search stops after the section of the dictionary
C        having the same first letters as the key has been checked, or
C        after a specified number of entries have been examined.  A special
C        dictionary entry, the vertical bar '|', will also terminate the
C        search.  This will speed things up if an appropriate dictionary
C        length parameter cannot be determined.  Both types of search are
C        sequential.  See "Notes" below for some suggestions if efficiency
C        is an issue.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        NDICT               I    I      Number of dictionary entries to be
C                                        examined.
C        DICTRY  NDICT       C    I      Array of dictionary entries,
C                                        left-justified in their fields.
C                                        May be alphabetized for efficiency,
C                                        in which case ALPHA should be .TRUE.
C                                        Entries with synonyms are of the form
C                                        'ENTRY:SYNONYM', where 'SYNONYM'
C                                        is a more fundamental entry in the
C                                        same dictionary.  NOTE: Don't build
C                                        "circular" dictionaries!
C        ALPHA               L    I      Indicates whether the dictionary
C                                        is in alphabetical order, in which
C                                        case the search can be terminated
C                                        sooner.
C        KEY                 C    I/O    String to be compared against the
C                                        dictionary.  Abbreviations are legal
C                                        provided they correspond to a unique
C                                        entry in the dictionary.  KEY is
C                                        replaced on termination by its most
C                                        fundamental equivalent dictionary
C                                        entry (uppercase, left-justified) if
C                                        a match was found.
C        ENTRY               I      O    Dictionary pointer.  If > 0, it
C                                        indicates which entry matched KEY.
C                                        In case of trouble, a negative value
C                                        means that a UNIQUE match was not
C                                        found - the absolute value of ENTRY
C                                        points to the second dictionary entry
C                                        which matched KEY.  Zero means that
C                                        NO match could be found.  ENTRY
C                                        always refers to the last search
C                                        performed - in searching a chain of
C                                        synonyms a non-positive value will be
C                                        returned if there is any break, even
C                                        if the original input key was found.
C
C
C     External references:
C
C        Name    Description
C        SCANNR  Finds first and last significant characters.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  We have assumed that the dictionary is not too big.  If
C             many searches are to be done or if the dictionary has more
C             than a dozen or so entries, it may be advantageous to build
C             an index array of pointers to the beginning of the section
C             of the dictionary containing each letter, then pass in the
C             portion of the dictionary beginning with DICTRY (INDEX).
C             (This won't generally work for dictionaries with synonyms.)
C             For very large problems, a completely different approach may
C             be advisable, e.g. a binary search for ordered dictionaries.
C
C        (3)  LOOKUP is case sensitive.  In most applications it will be
C             necessary to use an uppercase dictionary, and to convert the
C             input key to uppercase before calling LOOKUP.  Companion
C             routines TOKENS and PAIRS, available from the author, already
C             take care of this.
C
C        (4)  The key need not be left-justified.  Any leading (or
C             trailing) characters which are "non-significant" to SCANNR
C             will be ignored.  These include blanks, horizontal tabs,
C             commas, colons, and equal signs.  See SCANNR for details.
C
C        (5)  The ASCII collating sequence for character data is assumed.
C             (Note that this means the numerals precede the alphabet, unlike
C             common practice!)  On some machines, it may be necessary to
C             use the FORTRAN lexical library routines to force use of the
C             ASCII sequence.
C
C        (6)  Parameter NUMSIG sets a limit on the length of significant
C             dictionary entries.  Special applications may require that this
C             be increased.  (It is 16 in the present version.)
C
C        (7)  No protection against "circular" dictionaries is provided: don't
C             claim that A is B, and that B is A.  All synonym chains must
C             terminate!  Other potential errors not checked for include
C             duplicate or mis-ordered entries.
C
C        (8)  The handling of ambiguities introduces some ambiguity:
C
C                ALPHA = .TRUE.  A potential problem, when one entry
C                                looks like an abbreviation for another
C                                (eg. does 'A' match 'A' or 'AB'?) was
C                                resolved by dropping out of the search
C                                immediately when an "exact" match is found.
C
C                ALPHA = .FALSE. The programmer must ensure that the above
C                                situation does not arise: each dictionary
C                                entry must be recognizable, at least when
C                                specified to full length.  Otherwise, the
C                                result of a search will depend on the
C                                order of entries.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C        24 Feb. 1984  RAK/DAS  Initial design and coding.
C        25 Feb. 1984    RAK    Combined the two searches by suitable
C                               choice of terminator FLAG.
C        28 Feb. 1984    RAK    Optional synonyms in dictionary, no
C                               longer update KEY.
C        29 Mar. 1984    RAK    Put back replacement of KEY by its
C                               corresponding entry.
C        21 June 1984    RAK    Corrected bug in error handling for cases
C                               where no match was found.
C        23 Apr. 1985    RAK    Introduced test for exact matches, which
C                               permits use of dictionary entries which
C                               would appear to be ambiguous (for ordered
C                               case).  Return -I to point to the entry
C                               which appeared ambiguous (had been -1).
C                               Repaired loop termination - had to use
C                               equal length strings or risk quitting too
C                               soon when one entry is an abbreviation
C                               for another.  Eliminated HIT, reduced
C                               NUMSIG to 16.
C        16 May  1988    DAS    Had to use LEN in definition of LAST to
C                               suit revised SCANNR; early termination if
C                               KEY length exceeds LEN (DICTRY (1)).
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   ENTRY, NDICT
      LOGICAL
     >   ALPHA
      CHARACTER
     >   DICTRY (NDICT) * (*), KEY * (*)

C     Local constants.

      INTEGER
     >   NUMSIG
      CHARACTER
     >   BLANK, CURLY
      PARAMETER
     >   (BLANK = ' ', CURLY = '{', NUMSIG = 16)

C     Local variables.

      INTEGER
     >   FIRST, I, LAST, LENDIC, LENGTH, MARK
      CHARACTER
     >   FLAG * (NUMSIG), TARGET * (NUMSIG)

C     Procedures.

      EXTERNAL
     >   SCANNR


C     Execution.
C     ----------

      ENTRY = 0

C     Isolate the significant portion of the input key (if any).

      FIRST = 1
      LAST = LEN (KEY)
      CALL SCANNR (KEY, FIRST, LAST, MARK)
      IF (MARK .EQ. 0) GO TO 99

C     Can't hope to find a match if the key is longer than dictionary entries.

      LENGTH = MARK - FIRST + 1
      LENDIC = LEN (DICTRY (1))
      IF (LENGTH .GT. LENDIC) GO TO 99


C     The search starts with the input key, but may be repeated if that
C     target is just a synonym for a more fundamental dictionary entry.
C     NUMSIG = LEN (TARGET) is assumed to be plenty big enough.

      TARGET = KEY (FIRST:MARK)

   10 CONTINUE

C        Select search strategy by cunning choice of termination test
C        flag.  The left curly bracket follows all the alphabetic
C        characters in the ASCII collating sequence, but precedes the
C        vertical bar.

         IF (ALPHA) THEN
            FLAG = TARGET
         ELSE
            FLAG = CURLY
         END IF


C        Perform search.
C        ---------------

         I = 0
   20    CONTINUE
            I = I + 1
            IF (TARGET (1:LENGTH) .EQ. DICTRY (I) (1:LENGTH)) THEN
               IF (ENTRY .EQ. 0) THEN

C                 First "hit" - must still guard against ambiguities
C                 by searching until we've gone beyond the key (ordered
C                 dictionary), or until the end-of-dictionary mark is
C                 reached (exhaustive search).

                  ENTRY = I

C                 Special handling if match is exact - terminate search.
C                 We thus avoid confusion if one dictionary entry looks
C                 like an abbreviation of another.  This fix won't
C                 generally work for un-ordered dictionaries!

                  FIRST = 1
                  LAST = LENDIC
                  CALL SCANNR (DICTRY (ENTRY), FIRST, LAST, MARK)
                  IF (MARK .EQ. LENGTH) I = NDICT
               ELSE


C                 Oops - two hits!  Abnormal termination.
C                 ---------------------------------------

                  ENTRY = -I
                  GO TO 99
               END IF
            END IF

C           Check whether we've gone past the appropriate section of the 
C           dictionary.  The test on the index provides insurance and an
C           optional means for limiting the extent of the search.

            IF (DICTRY (I) (1:LENGTH) .LE. FLAG .AND. I .LT. NDICT)
     >   GO TO 20


C        Check for a synonym.
C        --------------------

         IF (ENTRY .GT. 0) THEN

C           Look for a second entry "behind" the first entry.  (FIRST
C           and MARK were determined above when the hit was detected.)

            FIRST = MARK + 2
            CALL SCANNR (DICTRY (ENTRY), FIRST, LAST, MARK)
            IF (MARK .GT. 0) THEN

C              Reset the target and dictionary pointer and repeat the
C              search for the synonym instead of the original key.

               TARGET = DICTRY (ENTRY) (FIRST:MARK)
               LENGTH = MARK - FIRST + 1
               ENTRY = 0
               GO TO 10

            ELSE

C              Expand the key to the full dictionary entry as a possible aid
C              to the calling program (which may prefer to avoid dealing with
C              entry numbers).

               KEY = DICTRY (ENTRY)

            END IF
         END IF


C     Normal termination.
C     -------------------

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE LSTFIL (LUNIN, LUNOUT, BUFFER)
C
C  PURPOSE: LSTFIL lists the indicated file in the sense of the TYPE
C           command of VMS.  The output may be the screen or another
C           file.  The file is assumed to be sequential/formatted.
C
C  METHOD:  Assume the file is already open and rewound.    Leave it
C           open (and not rewound) upon return.    Don't use a local
C           buffer -  the application is likely to have such scratch
C           available. Trailing blanks are suppressed - a frill per-
C           haps, but why not?
C
C  ENVIRONMENT: VAX/VMS - FORTRAN.
C
C  ARGUMENTS:
C  VAR    TYPE I/O/S    DIM    DESCRIPTION
C  LUNIN   I     I       -     Logical unit number for input file
C  LUNOUT  I     I       -     Logical unit number for output copy
C  BUFFER  C*(*) S       -     Buffer for one line of file - normally
C                              C*80, or maybe C*132.
C
C  FILES USED:  See arguments.
C
C  HISTORY:
C    April 1986  DAS/RAK  Adapted from ECHO - don't want the extra disk copy.
C    July  1989    DAS    Cosmetics in preparation for Mac translation.
C
C  AUTHORS:  David Saunders, Rob Kennelly, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNIN, LUNOUT
      CHARACTER
     >   BUFFER * (*)

C     Local constants:

      CHARACTER
     >   BLANK
      PARAMETER
     >  (BLANK = ' ')

C     Local variables:

      INTEGER
     >   LAST, MAXLEN

C     Execution:
      
C     Put out an initial blank line - saves the calling program from doing it:

      WRITE (LUNOUT, '(A)')

      MAXLEN = LEN (BUFFER)

   10 CONTINUE
         READ  (LUNIN, '(A)', END = 40) BUFFER

C        Find the last non-blank (if any).

         DO 20 LAST = MAXLEN, 1, -1
            IF (BUFFER (LAST : LAST) .NE. BLANK) GO TO 30
   20    CONTINUE
         LAST = 1

   30    WRITE (LUNOUT, '(1X, A)') BUFFER (1 : LAST)
         GO TO 10

   40 RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE LUSOLVE (N, NDIM, A, B, IER)
C
C     ONE-LINER:  Square-system solution by LU decomposition
C
C     DESCRIPTION:
C
C        LUSOLVE solves the N*N system  A x = b  by Gaussian elimination
C     with partial pivoting.  The right-hand side  b  is overwritten by the
C     solution  x.  (Matrix  A  is also overwritten.)
C
C        For more than one right-hand side, use instead DECOMP (once) and
C     SOLVE once per RHS  b.  See also DECSLV for the single-right-hand-side
C     case:  it is slightly more efficient than LUSOLVE in its processing of
C     b  as column N+1 of the matrix, but may be less convenient.
C
C        LUSOLVE was adapted from DECSLV, itself a merger of the original
C     DECOMP and SOLVE (reference below).  Note that one version of DECOMP
C     and SOLVE provides an estimate of the matrix condition number, while
C     another version avoids this extra work for applications where cond (A)
C     of no interest.  All of these variations vary the first index in the
C     inner loops.
C
C     INPUT:
C        N   : Order of matrix  A;  N >= 2
C        NDIM: Row dim. of array  A  declared in the calling program
C        A   : N*N array containing matrix  A
C        B   : Right-hand-side N-vector  b
C
C     OUTPUT:
C        B   : Solution vector  x
C        IER : 0 means no problem was detected;
C              1 means  A  was found to be singular
C
C     REFERENCE:  Forsythe, Malcolm & Moler, 1972.
C
C     HISTORY:    07/09/93  DAS  Revision of DECSLV to avoid having to
C                                store b in an extra column of A.
C                                IER = 0 now means no error (as it should
C                                have from the start).
C
C     PROGRAMMING: David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N, NDIM, IER
      REAL
     >   A (NDIM, N), B (N)

C     Local constants:

      REAL
     >   ONE, ZERO
      PARAMETER
     >  (ONE = 1.E+0, ZERO = 0.E+0)

C     Local variables:

      INTEGER
     >   I, J, K, M
      REAL
     >   T

C     Execution:

C     Perform the LU factorization, solving L y = b as we go:

      DO 60 K = 1, N - 1

C        Determine pivot element:

         M = K
         DO 20 I = K + 1, N
            IF (ABS (A (I, K)) .GT. ABS (A (M, K))) M = I
 20      CONTINUE

         T = A (M, K)
         A (M, K) = A (K, K)
         A (K, K) = T
         IF (T .EQ. ZERO) GO TO 90

C        Store the multipliers:

         T = -ONE / T
         DO 30 I = K + 1, N
            A (I, K) = A (I, K) * T
 30      CONTINUE

C        Apply the multipliers to current submatrix, including RHS:

         DO 50 J = K + 1, N
            T = A (M, J)
            A (M, J) = A (K, J)
            A (K, J) = T

            DO 40 I = K + 1, N
               A (I, J) = A (I, J) + A (I, K) * T
 40         CONTINUE

 50      CONTINUE

         T = B (M)
         B (M) = B (K)
         B (K) = T

         DO 55 I = K + 1, N
            B (I) = B (I) + A (I, K) * T
 55      CONTINUE

 60   CONTINUE

      IF (A (N, N) .EQ. ZERO) GO TO 90

C     Back substitution (solution of U x = y):

      DO 80 K = N, 2, -1
         T = B (K) / A (K, K)
         B (K) = T
         DO 70 I = 1, K - 1
            B (I) = B (I) - A (I, K) * T
 70      CONTINUE
 80   CONTINUE

      B (1) = B (1) / A (1, 1)
      IER = 0

      GO TO 99


 90   IER = 1  ! Matrix was singular

 99   RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE MSFIT (N, X, Y, CYCLIC, B, C, D, IER)
C
C  PURPOSE:
C     MSFIT computes the coefficients B(I), C(I), and D(I), I=1,2,...,N
C     for the MONOTONE cubic interpolating spline defined by
C
C     S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
C
C     for  X(I) <= X <= X(I+1).  It is derived from CSPLIN by J. Cordova
C     (Sterling Software, Oct. '86), which implements the algorithm of
C     Fritsch and Butland in SIAM J. SCI. STAT. COMPUT., VOL. 5, NO. 2,
C     June, 1984. MSFIT separates fitting the spline from evaluating it.
C     Use CSEVAL to evaluate the fit.
C
C     This class of spline ensures continuity of the first derivative at
C     the knots, but not of the second derivative.
C
C     This version provides a periodic (cyclic) option, needed to close
C     a curve smoothly as part of a parametric fit - see PSFIT.  (But
C     PLSFIT or PLSINTRP may be preferable for parametric fits.)
C
C     See also the LCSFIT implementation (Y vs. X derivation of PLSFIT's
C     X vs. T, Y vs. T).  These both calculate the spline coefficients
C     as they are needed, to avoid storing them.
C
C  ARGUMENTS:
C     NAME   DIM TYPE I/O/S DESCRIPTION
C     N       -   I     I   Number of data points, or knots.  N >= 2.
C     X       N   R     I   Abscissas of the knots, strictly increasing
C                           or strictly decreasing.
C     Y       N   R     I   Ordinates of the knots; Y(N)=Y(1) if CYCLIC.
C     CYCLIC  -   L     I   .TRUE. means periodic case; Y(N)=Y(1) reqd.;
C                           .FALSE. otherwise.
C     B,C,D   N   R     O   Spline coefficients (see PURPOSE and NOTES).
C     IER     -   I     O   0: No errors were detected;
C                           1: Too few data points; N < 2;
C                           4: Cyclic mode requested but Y(N) .NE. Y(1).
C  PROCEDURES:
C     PCHST   Function returning +1, -1, or 0. (No obvious way to make it
C             a statement function.  Source is in same module as MSFIT.)
C
C  NOTES:
C         Y(I) = S (X(I))
C         B(I) = S' (X(I))
C         C(I) = S'' (X(I))/2
C         D(I) = S''' (X(I))/6
C
C     Note that C (*) and D (*) are used for intermediate results before
C     their final definitions, to avoid local storage.
C
C  ENVIRONMENT:  VAX/VMS; FORTRAN 77
C                IMPLICIT NONE is an extension of FORTRAN 77.
C
C  HISTORY:
C  10/26/86  J. Cordova   Fit + evaluation, as routine CSPLIN
C  12/23/86  D. Saunders  Fit only, as MSFIT, for use with CSEVAL; avoid
C                         local storage; CSFIT-type nomenclature
C  01/14/87  D. Saunders  Introduced CYCLIC argument, for use by PSFIT.
C  02/27/87    "    "     Revised nomenclature following Kennelly's
C                         PLSFIT and the PCHIP package of Fritsch, et al.
C  10/23/87    "    "     Descending abscissas permissible now in light
C                         of revised CSEVAL.
C  06/21/91    "    "     References made to later related subroutines.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N, IER
      REAL
     >   X (N), Y (N), B (N), C (N), D (N)
      LOGICAL
     >   CYCLIC

C     Local constants:

      REAL
     >   ZERO, ONE, TWO, THREE, THIRD
      PARAMETER
     >  (ZERO  = 0.E+0,
     >   ONE   = 1.E+0,
     >   TWO   = 2.E+0,
     >   THREE = 3.E+0,
     >   THIRD = ONE/THREE)

C     Local variables:

      INTEGER
     >   I, ILAST
      REAL
     >   ALPHA, BMAX, H, TOP, W1, W2

C     Procedures:

      REAL
     >   PCHST
      EXTERNAL
     >   PCHST

C     Execution:

      IER = 0
      IF (N .LT. 2) THEN
         IER = 1
         GO TO 999
      END IF

      IF (CYCLIC) THEN
         IF (Y (N) .NE. Y (1)) THEN
            IER = 4
            GO TO 999
         END IF         
      END IF         

      DO 100, I = 1, N-1
         C (I) = X (I+1) - X (I)
         D (I) = (Y (I+1) - Y (I)) / C (I)
  100 CONTINUE

C     Linear interpolation for the N = 2 case.

      IF (N .EQ. 2) THEN
         B (1) = D (1)
         B (2) = D (1)
         C (1) = ZERO
         C (2) = ZERO
         D (1) = ZERO
         D (2) = ZERO
         GO TO 999
      END IF

C     Monotone cubic spline for N >= 3 case.

      IF (CYCLIC) THEN

C        Periodic case.  B (N) is like other B (*) if ...

         ILAST = N
         C (N) = C (1)
         D (N) = D (1)
      ELSE

C        First boundary point: 3-pt formula modified to preserve shape.

         ILAST = N - 1
         H = C (1) + C (2)
         W2 = -C (1) / H
         W1 = ONE - W2
         B (1) = W1 * D (1) + W2 * D (2)

         IF (PCHST (B (1), D (1)) .LE. ZERO) THEN
            B (1) = ZERO
         ELSE IF (PCHST (D (1), D (2)) .LT. ZERO) THEN
            BMAX = THREE * D (1)
            IF (ABS (B (1)) .GT. ABS (BMAX)) B (1) = BMAX
         END IF

C        Last boundary point, as for first.

         H = C (N-2) + C (N-1)
         W1 = -C (N-1) / H
         W2 = ONE - W2
         B (N) = W1 * D (N-2) + W2 * D (N-1)

         IF (PCHST (B (N), D (N-1)) .LE. ZERO) THEN
            B (N) = ZERO
         ELSE IF (PCHST (D (N-2), D (N-1)) .LT. ZERO) THEN
            BMAX = THREE * D (N-1)
            IF (ABS (B (N)) .GT. ABS (BMAX)) B (N) = BMAX
         END IF
      END IF

C     Interior points (Brodlie modification of Butland formula).
C     (This loop is vectorizable using SIGN for TOP and adding an
C     epsilon to the denominator, but at a cost in clarity.)

      DO 200, I = 2, ILAST
         TOP = D (I-1) * D (I)
         IF (TOP .GT. ZERO ) THEN          ! Slopes are same sign
            ALPHA = THIRD * (ONE + C (I) / (C (I-1) + C (I)))
            B (I) = TOP / (ALPHA * D (I) + (ONE - ALPHA) * D (I-1))
         ELSE
            B (I) = ZERO                   ! Make this point an extremum
         END IF
  200 CONTINUE

      IF (CYCLIC) THEN
         B (1) = B (N)
      END IF         

C     Quadratic and cubic coefficients for intervals 1:N-1.

      DO 300, I = 1, N-1
         H = ONE / C (I)
         C (I) = (THREE * D (I) - TWO * B (I) - B (I+1)) * H
         D (I) = (B (I) + B (I+1) - TWO * D (I)) * H * H
  300 CONTINUE

      IF (CYCLIC) THEN

C        Coefs. for interval N are same as coefs. for interval 1,
C        which in turn use interval 2.

         H = ONE / C (N)
         C (N) = (THREE * D (N) - TWO * B (N) - B (2)) * H
         D (N) = (B (N) + B (2) - TWO * D (N)) * H * H
      ELSE

C        Coefs. for X >= X (N) are related to derivatives at X (N),
C        using coefs. for interval N-1.

         H = X (N) - X (N-1)
         B (N) = B (N-1) + C (N-1) * TWO * H + D (N-1) * THREE * H * H
         C (N) = C (N-1) + D (N-1) * THREE * H
         D (N) = D (N-1)
      END IF

  999 RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION PCHST (ARG1, ARG2)
C
C  PURPOSE:
C     PCHST returns +1., -1., or 0. according to the signs of its two
C     arguments. It is convenient for monotonic spline routine MSFIT.
C     Nomenclature is that of PCHIP (Fritsch, et al.).
C
C-----------------------------------------------------------------------

      REAL
     >   ARG1, ARG2, ONE, PCHST, ZERO

      PARAMETER
     >   (ONE = 1.E+0, ZERO = 0.E+0)

      IF (ARG1 * ARG2 .EQ. ZERO) THEN
         PCHST = ZERO
      ELSE
         PCHST = SIGN (ONE, ARG1) * SIGN (ONE, ARG2)
      END IF

      END
C+----------------------------------------------------------------------
C
      FUNCTION NUMBER( STRING )
C
C     Acronym:
C
C        if ( NUMBER( string ) ) then
C           string very likely represents a number - integer or real
C
C     Description and usage:
C
C        A simple(-minded) test for numeric data is implemented by
C        searching the input string for legitimate characters:
C                digits 0 to 9, D, E, -, + and .
C        Insurance is provided by requiring that a numeric string
C        have at least one digit, at most one D, E or .
C        and at most two -s or +s.  Note that a few ambiguities remain:
C
C           (a)  A string might have the form of numeric data but be
C                intended as text.  No general test can hope to detect
C                such cases.
C
C           (b)  There is no check for correctness of the data format.
C                For example a meaningless string such as 'E1.+2-'
C                will be accepted as numeric.
C
C        Despite these weaknesses, the method should work in the
C        majority of cases.
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        NUMBER      -       L      O    Set .TRUE. if STRING appears
C                                        to be numerical data, else
C                                        set .FALSE.
C        STRING      *       C    I      Input data to be tested,
C                                        assumed to be in upper case.
C
C     Notes:
C
C        (1)  It is assumed that STRING has been extracted by
C             a "token" utility - hence the upper case assumption.
C
C        (2)  The scan of STRING stops at the first blank.
C
C        (3)  COMPLEX data with parentheses will not look numeric.
C
C     Environment:  ANSI FORTRAN 77.
C
C     Michael Saunders, Systems Optimization Lab., Stanford University.
C     12 Nov  1985    Initial design and coding, starting from the
C                     routine ALPHA from Informatics General, Inc.
C     23 May  1986    OPNUMB name changed to NUMBER, analogous to ALPHA,
C                     for use at NASA Ames - D. Saunders, Informatics.
C
C-----------------------------------------------------------------------

C  Arguments:

      LOGICAL          NUMBER
      CHARACTER*(*)    STRING

C  Local variables:

      LOGICAL         NUM
      INTEGER         J, LENGTH, NDIGIT, NEXP, NMINUS, NPLUS, NPOINT
      CHARACTER*1     ATOM

C  Executable statements:

      NDIGIT = 0
      NEXP   = 0
      NMINUS = 0
      NPLUS  = 0
      NPOINT = 0
      NUM    = .TRUE.
      LENGTH = LEN (STRING)
      J      = 0

   10    J    = J + 1
         ATOM = STRING (J:J)
         IF      (ATOM .GE. '0'  .AND.  ATOM .LE. '9') THEN
            NDIGIT = NDIGIT + 1
         ELSE IF (ATOM .EQ. 'D'  .OR.   ATOM .EQ. 'E') THEN
            NEXP   = NEXP   + 1
         ELSE IF (ATOM .EQ. '-') THEN
            NMINUS = NMINUS + 1
         ELSE IF (ATOM .EQ. '+') THEN
            NPLUS  = NPLUS  + 1
         ELSE IF (ATOM .EQ. '.') THEN
            NPOINT = NPOINT + 1
         ELSE IF (ATOM .EQ. ' ') THEN
            J      = LENGTH
         ELSE
            NUM    = .FALSE.
         END IF

         IF (NUM  .AND.  J .LT. LENGTH) GO TO 10

      NUMBER = NUM
     $         .AND.  NDIGIT .GE. 1
     $         .AND.  NEXP   .LE. 1
     $         .AND.  NMINUS .LE. 2
     $         .AND.  NPLUS  .LE. 2
     $         .AND.  NPOINT .LE. 1

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE OPENER (LUNCRT, PROMPT, LUNKBD, FNAME, LUNFIL, FSTAT)
C
C
C  PURPOSE:
C
C        OPENER modularizes the common situation of prompting for a file
C     name, opening the file, and reprompting if the file is not found.
C     Isolating any system dependencies here (such as VAX/VMS's READONLY
C     extension) enhances the portability of typical applications.
C
C        This version is restricted to sequential files, either formatted
C     or unformatted, old or new, with a few of the other occasionally
C     desirable attributes.
C
C        This version also permits proceeding anyway if a file specified
C     as old is not found.  This option required FSTAT to be used as an
C     OUTPUT as well as an input; all other options use it as input only.
C
C        System-dependent feature:
C
C        If there is trouble opening the file, the program user has the
C     option to execute an operating-system command - probably to look in
C     some directory to see what the filename should be - then try again.
C
C
C  ARGUMENTS:
C
C   NAME     DIM    TYPE I/O/S DESCRIPTION
C  LUNCRT     -      I     I   Logical unit for screen prompts.
C  PROMPT     *      C     I   Prompt (possibly indicating default filename).
C                              May be blank to suppress the prompt (but FNAME
C                              must be non-blank in this case).
C  LUNKBD     -      I     I   Logical unit for keyboard entries.
C  FNAME      *      C    I/O  Name of file to be opened - see METHOD/NOTES.
C                              May be blank to permit termination of an
C                              indefinite loop over file names (carriage
C                              return response to prompt, no open attempted,
C                              and a check in the calling program to see if
C                              FNAME is still blank).  A non-blank FNAME on
C                              input will be treated as the default file name.
C  LUNFIL     -      I     I   Logical unit for file to be opened.
C  FSTAT      *      C   I[/O] INPUT:  String of file status attributes,
C                              separated by commas, colons, or blanks.
C                              Unique abbreviations suffice, in upper or
C                              lower case.  The more common possibilities
C                              provided for appear below.
C                              OUTPUT:  None, with one exception: if the
C                              keyword 'IfPresent' is included in the input 
C                              FSTAT string, then FSTAT = 'MISSING' (upper
C                              case) is returned as output if the file is
C                              not found.  ('Old' is implied by 'IfPresent'
C                              here.)  FSTAT must be a character VARIABLE
C                              in this one case; a CONSTANT is fine in all
C                              other cases.
C                              
C  FSTAT examples:
C
C     'NEW, LIST, 160'      New formatted file with up to 160 characters per
C                           record, and implied (not explicit) carriage control
C
C     'old:binary:write'    Existing unformatted file where READONLY access (the
C                           default) is not enough
C
C     'IfPresent'           Suppresses reprompting if the specified file is
C                           not found.
C
C
C  FSTAT token summary:
C
C     (Meaningful combinations should be concatenated.)
C
C  Attribute        Description          Corresponding OPEN keyword   Default
C
C     'OLD'        File already exists            STATUS              'UNKNOWN'
C     'NEW'        New file desired                "  "                 "   "
C     'SCRATCH'    New file; deleted upon closing  "  "                 "   "
C
C     'BINARY'     Unformatted file                FORM              'FORMATTED'
C     'UNFORMATTED'     "      "                   "  "                 "   "
C
C     'NONE'       No implied carriage ctrl     "  "  "    Default for unfrmttd.
C
C     'WRITE'      Write (and read) access        READONLY           The default
C                  is 'READONLY' if STATUS='OLD', else it is 'WRITE'=read/write.
C
C     'nnn'        Max. record length             RECL         <System default> 
C
C     'IfPresent'  See FSTAT description and examples.
C
C  Further FSTAT Notes:
C
C     Integers nnn for RECL refer to bytes (characters) for formatted files
C     (where 132 is normally the upper limit), or longwords (4-byte units)
C     for unformatted files (where the default is system-dependent).
C
C     READONLY is normally needed for reading files of another owner.
C
C     The defaults are also legitimate (if redundant) input values.
C
C
C  METHOD:
C
C     (1) Default the file attributes.  (The READONLY one will be adjusted
C         later if the file status is new or unknown.)
C
C     (2) Decode the subset of attributes indicated in FSTAT, one token at
C         a time - just itemize the finite number of cases in a way that
C         could be extended.  Any unknown attributes are considered fatal
C         programmer errors - stop with a diagnostic.
C
C     (3) Prompt for the name of the file to be opened (unless the prompt
C         is blank - handy for opening files with fixed names).
C
C     (4) If a carriage return is entered, and FNAME was passed to OPENER
C         as all blanks, return immediately without opening a file.  The
C         calling program can detect this case by checking if FNAME is
C         still blank - handy for performing an indefinite loop over
C         multiple files.
C
C     (5) If a carriage return is entered, and FNAME was NOT all blanks,
C         then FNAME is assumed to represent a default filename (presumably
C         indicated as part of the prompt).
C
C     (6) If "EOF" is entered (^Z under VMS; ^D under Unix), this is assumed
C         to mean "STOP."  Perhaps this should be indicated to the calling
C         program via a returned argument value.  However, the original
C         design had the stop occurring here, and there is no good, upwardly
C         compatible way of changing that.
C
C     (7) If the keyword 'IfPresent' has been input in FSTAT, use INQUIRE to
C         detect existence, and return with FSTAT = 'MISSING' if appropriate.
C
C     (8) Attempt to open the file:
C
C         IF the file cannot be opened THEN
C            Inform the user.
C            Prompt for either the correct file or an operating
C               system command (probably for a directory listing).
C            IF the first character found is a '$' THEN
C               Execute the associated command (system-dependent)
C               Go back to (3)
C            ELSE
C               Go back to (8)
C            END IF
C         END IF
C
C  EXTERNAL REFERENCES:
C
C     LOOKUP        Dictionary lookup: permits abbreviations
C     NUMBER        Identifies RECL parameter
C     READER        Prompting utility
C     SCANNR        Identifies tokens in FSTAT string
C     SYSTEM        IRIS utility analogous to VMS's LIB$SPAWN
C     UPCASE        Upper-casifies FSTAT prior to SCANNR/LOOKUP
C
C  ERROR HANDLING:
C
C        Normally, if a file cannot be opened, a message to that effect is
C     sent to the screen.  The user then has three options:  re-enter the
C     file, which is consequently opened; type a command to list directory
C     contents, whereupon the prompt/response cycle is reinitiated; or quit.
C
C        Alternatively, a missing file may be handled by the calling program.
C     See the FSTAT description for more.
C
C        File attributes found in FSTAT (e.g. 'NULL' for the BLANK keyword)
C     which are not (yet) handled by OPENER are fatal: OPENER will STOP.
C
C  NOTES:
C
C        The declared length of the name of the file to be opened should be
C     enough to accommodate a possibly lengthy directory specification.  A
C     generous size, such as CHARACTER * 60, is suggested.
C
C  SYSTEM DEPENDENCIES:
C
C     (1) IMPLICIT NONE is nonstandard.
C     (2) READONLY is a VAX/VMS extension - removed in this version.
C     (3) CALL SYSTEM is IRIS-specific.
C
C  ENVIRONMENT:  IRIS/IRIX, FORTRAN 77
C
C  HISTORY:
C
C  02/27/86  R.G.Langhi  Initial design and code (formatted files only).
C
C  05/24/86  D.Saunders  Provided for defaulting filename, and for
C                        returning with FNAME=' ' and no file opened
C                        to handle the indefinite-loop-over-files case;
C                        documentation clarified (sequential/formatted).
C
C  08/28/86  RGL         Added unformatted sequential capability.
C
C  05/28/87  RGL         Bug:  the default filename was lost after an
C                        erroneous response to the prompt; need to keep
C                        a copy locally.  On cancellation, FNAME is
C                        returned now with the default instead of blank.
C                        (Note:  The default filename MAY be blank.)
C
C  05/10/88  DAS         Generalized FSTAT usage to indicate more than
C                        just old/new and formatted/unformatted, giving
C                        CARRIAGECONTROL, READONLY, and RECL control too;
C                        suppressed prompt if it is blank, so that OPENER
C                        can still be used with fixed file names.
C
C  02/15/90  DAS         Generalized FSTAT further to permit the calling
C                        program to proceed without the specified (old)
C                        file (as opposed to OPENER's reprompting, which
C                        remains the usual option).
C                        Also: changed from READC to READS to accommodate
C                        the case-sensitive Unix community.
C
C  02/22/90  DAS         IRIS 4D version: suppressed LIB$SPAWN feature
C                        (is there an IRIX equivalent?); ^D instead of ^Z
C                        indicated in the prompts.
C
C  03/14/90  DLH/DAS     Dexter found "call system" as a substitute for
C                        LIB$SPAWN.
C
C  08/28/90  DAS         Unformatted files needed CC='NONE'.
C
C  06/10/91  DAS         NUMBER was declared INTEGER - should be LOGICAL.
C
C  10/15/99  DAS         Removed CARRIAGECONTROL & READONLY for SGI/IRIX f90.
C
C  AUTHOR:  Ronald Langhi, NASA Ames/Sterling Software, Palo Alto, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNFIL
      CHARACTER
     >   FNAME * (*), FSTAT * (*), PROMPT * (*)

C     Local constants:

      INTEGER
     >   MXCHAR, MXWORD
      CHARACTER
     >   BLANK * 1

      PARAMETER
     >  (BLANK  = ' ',
     >   MXCHAR = 11,
     >   MXWORD = 13)

C     Local variables:

      INTEGER
     >   ENTRY, FIRST, LAST, MARK, RECLEN
      CHARACTER
     >   ATTRIB * 80, CC * 7, DICT (MXWORD) * (MXCHAR), FORIG * 80,
     >   FORM * 11, KEY * (MXCHAR), STATUS * 7
      LOGICAL
     >   CR, DFRECL, ENQUIRE, EOF, PRESENT, RDONLY

C     Procedures:

      LOGICAL
     >   NUMBER
      EXTERNAL
     >   LOOKUP, NUMBER, SCANNR, READS, SYSTEM, UPCASE

C     Storage:

      DATA
     >   DICT
     >      /'BINARY', 'FORMATTED', 'FORTRAN', 'LIST', 'NEW', 'NONE',
     >       'OLD', 'READONLY', 'SCRATCH', 'UNFORMATTED', 'UNKNOWN',
     >       'WRITE', 'IFPRESENT'/
C     The dictionary should be in upper case.  It need not be alphabetized.


C     Execution:

C     Save the default input filename in case of an erroneous response 
C     to the prompt:

      FORIG = FNAME

C     Default the file attributes so that corresponding lookup hits can
C     be ignored (and input string FSTAT can be short):

      FORM    = DICT (2)
      CC      = DICT (3)
      STATUS  = DICT (11)
      RDONLY  = .TRUE.
      DFRECL  = .TRUE.
      ENQUIRE = .FALSE.

C     Ensure that the attributes text is upper case:

      ATTRIB = FSTAT
      CALL UPCASE (ATTRIB)

C     Start of loop over tokens in attributes text:

      FIRST = 1
      LAST = LEN (FSTAT)
   50 CONTINUE

         CALL SCANNR (ATTRIB, FIRST, LAST, MARK)
         IF (MARK .EQ. 0) GO TO 70

         KEY = ATTRIB (FIRST : MARK)
         CALL LOOKUP (MXWORD, DICT, .FALSE., KEY, ENTRY)

         IF (ENTRY .GT. 0) THEN

C           There is only a modest number of possibilities.
C           Avoid replicating the dictionary text by working with
C           subscripts rather than text.  Adding keywords at the
C           end of the dictionary will not affect this code, which
C           is why LOOKUP's non-alphabetic option is used above.

            IF (ENTRY .EQ. 1 .OR. ENTRY .EQ. 10) THEN
               FORM = DICT (10)
               CC = DICT (6)
            END IF
            IF (ENTRY .EQ. 4 .OR.
     >          ENTRY .EQ. 6)    CC = DICT (ENTRY)
            IF (ENTRY .EQ. 12)   RDONLY = .FALSE.
            IF (ENTRY .EQ. 5 .OR.
     >          ENTRY .EQ. 7 .OR.
     >          ENTRY .EQ. 9)    STATUS = DICT (ENTRY)
            IF (ENTRY .EQ. 13)   ENQUIRE = .TRUE.

         ELSE IF (NUMBER (ATTRIB (FIRST : MARK))) THEN
            DFRECL = .FALSE.
            READ (ATTRIB (FIRST : MARK), '(BN, I11)') RECLEN
         ELSE
            GO TO 810
         END IF

         FIRST = MARK + 2
         IF (FIRST .LE. LAST)
     >GO TO 50

   70 CONTINUE

C     Check for inconsistencies:

      IF (STATUS .NE. DICT (7)) RDONLY = .FALSE.


  100 CONTINUE

C     Start of loop over retries at opening the file:

      IF (PROMPT .EQ. BLANK) GO TO 300

C     The prompt should indicate any default file name here.
C     Use READS instead of the original READC, for Unix reasons.

      CALL READS (LUNCRT, PROMPT, LUNKBD, FNAME, CR, EOF)

  200 CONTINUE

      IF (CR .AND. FNAME .EQ. BLANK)  GO TO 999
      IF (EOF) GO TO 800

  300 CONTINUE

C     Try to open the file, but test for its existence first in one case:

      IF (ENQUIRE) THEN
         INQUIRE (FILE=FNAME, EXIST=PRESENT)
         IF (.NOT. PRESENT) THEN
            FSTAT = 'MISSING'
            GO TO 999
         END IF
      END IF

C     Oddball READONLY keyword, and uncertain system default for RECL
C     force four possibilities here:

      IF (RDONLY .AND. DFRECL) THEN

         OPEN (UNIT=LUNFIL, FILE=FNAME, FORM=FORM, STATUS=STATUS,
     >      ERR=400)

      ELSE IF (RDONLY .AND. .NOT.DFRECL) THEN

         OPEN (UNIT=LUNFIL, FILE=FNAME, FORM=FORM, STATUS=STATUS,
     >      RECL=RECLEN, ERR=400)

      ELSE IF (.NOT.RDONLY .AND. DFRECL) THEN

         OPEN (UNIT=LUNFIL, FILE=FNAME, FORM=FORM, STATUS=STATUS,
     >      ERR=400)

      ELSE
C****    IF (.NOT.RDONLY .AND. .NOT.DFRECL) THEN

         OPEN (UNIT=LUNFIL, FILE=FNAME, FORM=FORM, STATUS=STATUS,
     >      RECL=RECLEN, ERR=400)

      END IF

      GO TO 999


C     OPEN error handling:

  400 CONTINUE

      WRITE (LUNCRT, 1000) 'Error opening following file:', FNAME
      FNAME = FORIG

      IF (FNAME .EQ. BLANK) THEN
         CALL READS (LUNCRT,
     >      'Try again, look in directory ("$ls ..."), cancel open ' //
     >      '(CR), or stop (^D): ', LUNKBD, FNAME, CR, EOF)
      ELSE
         CALL READS (LUNCRT,
     >      'Try again, look in directory ("$ls ..."), open default ' //
     >      '(CR), or stop (^D): ', LUNKBD, FNAME, CR, EOF)
      END IF

C     Either another attempt at the file name was entered ...

      IF (FNAME (1 : 1) .NE. '$') GO TO 200

C     ... or an operating-system command was requested (system-dependent) ...

!!!!!!      CALL SYSTEM (FNAME (2 :))
      CALL EXECUTE_COMMAND_LINE(FNAME(2:))

C     ... and start over:

      FNAME = FORIG
      GO TO 100


  800 WRITE (LUNCRT, 1000) 'Stopping as requested.'
      GO TO 990

  810 WRITE (LUNCRT, 1000) 'OPENER: Bad file attribute in FSTAT:',
     >   FSTAT (FIRST : MARK), 'Aborting.'
C**** GO TO 990

  990 STOP ' '

  999 RETURN

C     Formats:

 1000 FORMAT (1X, A)

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE PAIRS (STRING, NUMBER, KEYS, VALUES)
C
C
C     Description and usage:
C
C           An aid to parsing input data in keyword format.  A character
C        string is broken up into keywords and corresponding values.  The
C        "keywords" and "values" are simply alternating fields of significant,
C        contiguous characters, converted to uppercase.  The following
C        characters are NON-significant, and hence may serve as separators:
C           (a)  blanks,
C           (b)  horizontal tabs,
C           (c)  commas,
C           (d)  colons, and
C           (e)  equal signs.
C        See SCANNR for details.  Processing continues until the requested
C        number of (KEYS, VALUES) pairs have been found or the end of the
C        input string is reached.  A trailing key without a corresponding
C        value will be recognized and included in the count of pairs returned.
C        (Entries requested, but not found, will be set to the blank string.)
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING              C    I      Input string.
C        NUMBER              I    I/O    Number of (KEYS,VALUES) pairs
C                                        requested (input) and found (output).
C        KEYS    NUMBER      C      O    Array of keywords, left-justified and
C                                        converted to uppercase.  "Leftover"
C                                        array entries are filled with blanks.
C        VALUES  NUMBER      C      O    Array of values in the form of
C                                        character strings, left-justified and
C                                        converted to uppercase.  "Leftover"
C                                        array entries are filled with blanks.
C
C
C     External references:
C
C        Name    Description
C        SCANNR  Finds positions of the first and last significant characters.
C        UPCASE  Converts a string to uppercase.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C        21 Feb. 1984  RAK/DAS    Initial design and coding based on TOKENS.
C        16 Mar. 1984    RAK      Cosmetic changes.
C         1 July 1985    RAK      Don't change input string to uppercase,
C                                 change keys/values individually instead.
C
C-----------------------------------------------------------------------


C     Variable declarations.
C     ----------------------

      IMPLICIT NONE

C     Parameters.

      CHARACTER
     >   BLANK
      PARAMETER
     >   (BLANK = ' ')

C     Variables.

      LOGICAL
     >   ODD
      INTEGER
     >   COUNT, I, FIRST, LAST, MARK, NUMBER
      CHARACTER
     >   KEYS (NUMBER) * (*), STRING * (*), VALUES (NUMBER) * (*)

C     Procedures.

      EXTERNAL
     >   SCANNR, UPCASE


C     Executable statements.
C     ----------------------

C     WHILE there are tokens left, loop UNTIL enough have been found.

      FIRST = 1
      LAST = LEN (STRING)

      ODD = .FALSE.
      COUNT = 0
   10 CONTINUE
         CALL SCANNR (STRING, FIRST, LAST, MARK)
         IF (LAST .GT. 0) THEN

C           The tokens go alternately into KEYS and VALUES.

            ODD = .NOT.ODD         
            IF (ODD) THEN
               COUNT = COUNT + 1
               KEYS (COUNT) = STRING (FIRST:MARK)
               CALL UPCASE (KEYS (COUNT))
            ELSE
               VALUES (COUNT) = STRING (FIRST:MARK)
               CALL UPCASE (VALUES (COUNT))
            END IF

            FIRST = MARK + 2
            IF (COUNT .LT. NUMBER .OR.
     >         (COUNT .EQ. NUMBER .AND. ODD)) GO TO 10

         END IF


C     Fill the rest of the output arrays with blanks and set NUMBER for output.

      IF (ODD) VALUES (COUNT) = BLANK
      DO 20 I = COUNT + 1, NUMBER
         KEYS (I) = BLANK
         VALUES (I) = BLANK
   20 CONTINUE

      NUMBER = COUNT


C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE PROTECT (NX, X, Y, ARROW, DISTINCT)
C
C     One-liner: Test data for monotonicity and distinctness
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        PROTECT scans a curve represented by arrays X and Y to check
C     strict monotonicity (increasing or decreasing) of the X data and
C     to verify that no two successive points match. One or the other
C     of these tests is required by conventional or parametric spline
C     fitting methods. Since it is often inefficient to perform these
C     tests in the fitting routine, which may be called more than once
C     with the same data or with data known to be good, the present
C     modular approach was chosen.
C
C        The initial application is program SMOOTH, where a single set
C     of data points may be fit by several different techniques - the
C     idea was to check each dataset just once as it was read in, and
C     then test flags ARROW and DISTINCT as needed.
C
C     Arguments:
C     ----------
C
C     Name     Type/Dimension  I/O/S  Description
C     NX       I               I      Dimension of X and Y arrays. If
C                                     NX <= 1, then ARROW will be
C                                     0.0 and DISTINCT will be .FALSE.
C
C     X        R (NX)          I      Array of abscissas.
C
C     Y        R (NX)          I      Array of ordinates.
C
C     ARROW    R                 O    Monotonicity indicator:
C                                       -1.0  strictly decreasing
C                                        0.0  neither
C                                       +1.0  strictly increasing
C
C     DISTINCT L                 O    Indicates whether successive points
C                                     are distinct in both X and Y.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     (2)  There is no provision for roundoff error in the monotonicity
C          test, i.e., all interval lengths are compared to zero.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     22 Feb. 1987    RAK    Initial design and coding.
C     14 Apr. 1988    RAK    Corrected value returned by DISTINCT when
C                            NX <= 1. Revised comments.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE
      PARAMETER
     &  (ZERO   = 0.0E+0,
     &   ONE    = 1.0E+0)

C     Arguments.

      LOGICAL
     &   DISTINCT
      INTEGER
     &   NX
      REAL
     &   ARROW, X (NX), Y (NX)

C     Local variables.

      LOGICAL
     &   FIRST
      INTEGER
     &   J
      REAL
     &   DX

C     Execution.
C     ----------

C     Set the output flags and bail out if the input data is trivial
C     or improper.

      ARROW    = ZERO
      DISTINCT = .FALSE.
      IF (NX .LE. 1) GO TO 990

C     Reset the direction flag according to the first interval, and reset
C     distinctness flag prior to testing.

      DX = X (2) - X (1)
      IF (DX .NE. ZERO) THEN
         FIRST = (DX .GT. ZERO)
         ARROW = SIGN (ONE, DX)
      END IF

      DISTINCT = .TRUE.

C     The approach is to try to set ARROW = ZERO and DISTINCT = .FALSE.
C     as we go, and quit as soon as possible. No harm is done if ARROW is
C     set repeatedly - it's not worth testing for.

      DO 10, J = 1, NX - 1
         DX = X (J + 1) - X (J)
         IF (DX .NE. ZERO) THEN

C           Compare the sign of the increment in this interval to that of
C           the first interval.

            IF ((DX .GT. ZERO) .NEQV. FIRST) ARROW = ZERO
         ELSE

C           If a pair of X's match, the data is not strictly monotonic. We
C           must still check whether the Y's are distinct.

            ARROW = ZERO
            IF (Y (J + 1) - Y (J) .EQ. ZERO) THEN

C              The data is neither monotonic nor distinct - time to die.

               DISTINCT = .FALSE.
               GO TO 990

            END IF
         END  IF
   10 CONTINUE

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+---------------------------------------------------------------------
C
      SUBROUTINE PSEVAL (NX, X, Y, YP, YPP, T1, T2, NGEOM, XGEOM,
     >                   YGEOM, TOL, IER)
C
C  ACRONYM:    Parametric Spline EVALuation (Y for given X, not T)
C              -          -      ----
C  PURPOSE:
C                 PSEVAL will normally be called after subroutine PSFIT
C              has fit a parametric cubic spline to a given  curve,  to
C              evaluate the monotonic sub-curve defined by  T1 and  T2,
C              at arbitrary abscissas which are assumed to be  in-range
C              for this sub-curve.  PSEVAL computes ordinates and first
C              and second derivatives.  (See NOTES below if derivatives
C              are not needed, and for tips on choosing T1, T2.)
C
C                 Use PSTVAL if evaluations for given T, not X, are all
C              that is needed.
C
C  METHOD:
C                 A zero finder is used to approximate the value of the
C              parametric variable corresponding to the given X  on the
C              subcurve defined by T1 and T2.   The corresponding Y  is
C              theneasily evaluated, as are the derivatives.   Outline:
C
C              FOR each abscissa X=X(J), J=1:NX, DO
C               * Use a zero finder to find the zero of the function
C                 represented by routine GETT in the range  [T1,T2].
C                 This value, TBEST,  is such that evaluation of the
C                 spline from PSFIT on  TGEOM, XGEOM  at TBEST gives
C                 an XTEST arbitrarily close to the given X.
C               * Evaluate the spline from PSFIT on TGEOM, YGEOM for
C                 the corresponding ordinate, Y(J).
C               * Compute 1st and 2nd derivatives by the chain rule.
C
C NOTES:
C         (1)  PSEVAL calculates derivatives whether they are needed or
C              not, for reasons of simplicity.  Undesirable storage may
C              be avoided, however, since the output values are assign-
C              ed in the order YPP(*), YP(*), Y(*).  For example, if no
C              derivatives are required, use the same array for each of
C              these arguments; if just ordinates and first derivatives
C              are desired  (but no second derivatives),  pass the same
C              array for YPP as for YP.
C
C         (2)  Derivatives with respect to  X  require derivatives with
C              respect to the parametric variable,  so CSDVAL is needed
C              rather than CSEVAL.  GETT does not need derivatives, but
C              it was considered cleaner to use CSDVAL in both places.
C
C         (3)  FUNCTION TSUBJ may be used to determine T1 and T2. TSUBJ
C              is included in the same source module as  PSEVAL & GETT.
C              Note that the end-points of sub-curves are not necessar-
C              ily at original data points.   (Consider, say,  a modest
C              number of points roughly in the form of a  circle.   The
C              fitted curve may "bulge" at the sides...)
C
C  ARGUMENTS:
C   ARG    TYPE  I/O/S   DIM   DESCRIPTION
C   NX      I      I      -    Number of evaluations required.  NX>=1.
C   X       R      I      NX   Arbitrary abscissas (except all must be
C                              valid for the indicated sub-curve).
C   Y       R      O      NX   Ordinates corresponding to X(*).
C   YP      R      O      NX   First derivatives, DY/DX (see NOTES).
C   YPP     R      O      NX   Second derivatives w.r.t. X (see NOTES).
C   T1,     R      I      -    Values of the parametric variable that
C    T2     R      I      -    delimit the desired monotonic sub-curve.
C                              T1 < T2 is required.  See NOTES.  Using
C                              TSUBJ ( J ) to define the value of T at
C                              data point J will commonly suffice, but
C                              the X-extent of sub-curves is not nece-
C                              ssarily bounded by the original points.
C   NGEOM   R      I      -    Number of data points fitted by PSFIT.
C   XGEOM   R      I    NGEOM  Geometry abscissas fitted by PSFIT.
C   YGEOM   R      I    NGEOM  Geometry ordinates fitted by PSFIT.
C   TOL     R      I      -    Convergence criterion for zero finder, which
C                              returns TBEST as its best approximation to
C                              the zero of function GETT(T).  TOL is applied
C                              relative to the data range as given by the
C                              total arc length.  For normal applications
C                              and 32-bit arithmetic, TOL=1.E-6 is suggested;
C                              for 64-bit arithmetic, TOL~1.E-14.
C   IER     I      O      -    Error return code from ZEROIN:
C                                0 if no error;
C                                1 if F(T1), F(T2) are not of opposite sign;
C                                2 if T1 >= T2.
C
C  ERROR HANDLING: Nonzero IER from ZEROIN will force early return.
C
C  PROCEDURES:
C   CSDVAL  Evaluates a 1-dimensional cubic spline along with its first
C           and second derivatives.  (Also used by GETT.)
C   GETT    Real function passed as argument to ZEROIN.  GETT is in
C           the same source module as PSEVAL, and should not concern
C           the user unless the following COMMON block is changed.
C   ZEROIN  Finds zero of a function which changes sign in a given interval.
C
C  INTERNAL COMMON BLOCK:
C   /PSWORK/ Work arrays, set up by PSFIT.
C            Dimensions are set to MAXPS by PARAMETER statement.
C   VAR    TYPE I/O/S DIM  DESCRIPTION
C   TGEOM   R     I MAXPS  Values of the paramateric variable at the
C                          data points, set up by PSFIT.
C   PSCOFS  R     I MAXPS  Cubic spline coefficients from CSFIT.
C                   *3*2   An array of size MAXPS*3 is required for
C                          XGEOM and another for YGEOM.
C   XLOC    R     S   -    Copy of X, abscissa at which PSEVAL is seeking
C                          an ordinate.
C   NCURV   I     O   -    Copy of NGEOM.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C  HISTORY:
C  Jan.  83  R.Lefkowitz  Original design and coding.
C  Sept. 84    "    "     CSDVAL used in place of IMSL's ICSEVU, DCSEVU.
C  Sept. 85  D.Saunders   Description of T1, T2 changed to reflect move
C                         to cumulative chord length from point index as
C                         the parametric variable. No code changes needed.
C  April 86    "    "     Calling sequence to ZEROIN changed (no LUN now).
C  Oct.  87    "    "     EPS=1.E-12 was too extreme - try 1.E-8. (What IS
C                         the proper fix here for zero DX/DT?)
C  Oct.  91    "    "     TOL is applied relative to total arc-length now.
C                         So is EPS (back to 1.E-10, with EPS**3 avoided).
C                         Other cosmetics.
C
C  AUTHOR: Rosalie Lefkowitz, Sterling Software/NASA Ames, Moffett Field, CA
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   NX, NGEOM, IER
      REAL
     >   T1, T2, TOL, X (NX), Y (NX), YP (NX), YPP (NX), XGEOM (NGEOM),
     >   YGEOM (NGEOM)

C  *  Internal COMMON:

      INTEGER
     >   MAXPS
      PARAMETER
     >  (MAXPS = 1001)
      INTEGER
     >   NCURV
      REAL
     >   PSCURV, TGEOM, PSCOFS, XLOC
      COMMON /PSWORK/
     >   PSCURV (MAXPS),TGEOM (MAXPS),PSCOFS(MAXPS,3,2),XLOC(1),NCURV    ! RLC 20211024

C  *  Local constants:

      REAL
     >   EPS
      PARAMETER
     >  (EPS = 1.E-10)  ! EPS is a limit on |DX/DT| to guard against
                        ! division by zero.

C  *  Local variables:

      INTEGER
     >   J
      REAL
     >   EPSIL, TBEST(1), TOLER, XD1(1), XD2(1), YD1(1), YD2(1), YJ(1)  ! RLC 20211024 

C  *  Procedures:

      REAL
     >   GETT, ZEROIN
      EXTERNAL
     >   CSDVAL, GETT, ZEROIN

C  *  Execution:

      TOLER = TOL * TGEOM (NGEOM)
      EPSIL = EPS * TGEOM (NGEOM)

      DO 100, J = 1, NX

C  *     Estimate parametric variable value TBEST corresponding to X (J),
C        using spline for TGEOM, XGEOM.

         XLOC = X(J)  ! not sure what this is for?? XLOC is set by CSDVAL

         TBEST(1) = ZEROIN (T1, T2, GETT, TOLER, IER)   ! RLC 20211024 
         IF (IER .NE. 0) GO TO 999

C  *     Evaluate spline TGEOM, XGEOM for derivatives at TBEST:

         CALL CSDVAL (NGEOM, TGEOM, XGEOM, 1, TBEST, PSCOFS (1:,1,1),    ! RLC 20211024 
     >                PSCOFS (1:,2,1), PSCOFS (1:,3,1), XLOC, XD1, XD2)

C  *     Evaluate spline TGEOM, YGEOM for Y and its derivatives at TBEST:

         CALL CSDVAL (NGEOM, TGEOM, YGEOM, 1, TBEST, PSCOFS (1:,1,2),   ! RLC 20211024 
     >                PSCOFS (1:,2,2), PSCOFS (1:,3,2), YJ, YD1, YD2)

C  *     Evaluate derivatives using chain rule.  (Do it in the order
C        y", y', y because derivatives may not be needed - see NOTES.)

         IF (ABS (XD1(1)) .LT. EPSIL) XD1(1) = SIGN (EPSIL, XD1(1))   ! RLC 20211024 

C****    YPP(J) = (XD1 * YD2 - YD1 * XD2) / XD1 ** 3
         YPP(J) = (YD2(1) - YD1(1) * XD2(1) / XD1(1)) / XD1(1) ** 2   ! RLC 20211024 
         YP(J)  = YD1(1) / XD1(1)                                     ! RLC 20211024 
         Y(J)   = YJ(1)                                               ! RLC 20211024 
  100 CONTINUE

  999 RETURN
      END
C+---------------------------------------------------------------------
C
      FUNCTION GETT (T)
C
C  ACRONYM:  GET T corresponding to given X
C            --- -
C  PURPOSE:
C             GETT is used by PSEVAL as part of evaluating a parametric
C          spline as fit by PSFIT, for some X.  It evaluates
C                              F = XLOC - XTEST
C          where XLOC is an arbitrary X, and XTEST is the ordinate cor-
C          responding to the abscissa T of the spline for TGEOM & XGEOM
C          (TGEOM & PSCURV actually, since these must be passed through
C          COMMON).
C
C  METHOD:
C             * Call CSDVAL to evaluate spline for TGEOM & PSCURV at T,
C               with XTEST as output.
C             * Set GETT = XLOC - XTEST.
C
C             Eventually the zero-finder ZEROIN will find a T such that
C          this difference is arbitrarily small, meaning a value of the
C          parametric variable for the given abscissa has been found.
C
C  ARGUMENTS:
C   ARG    TYPE  I/O/S   DIM   DESCRIPTION
C    T      R    I        -    Used as test abscissa of spline for TGEOM
C                              & PSCURV.  Desired ordinate corresponding
C                              to T is XLOC.
C
C  PROCEDURES:
C   CSDVAL  Spline evaluation with optional derivatives
C
C  INTERNAL COMMON BLOCK:
C           /PSWORK/ Work space set up by PSFIT/PSEVAL.
C           Row dimensions are set to MAXPS by PARAMETER statement.
C   VAR    TYPE I/O/S DIM  DESCRIPTION
C   PSCURV  R     I MAXPS  Copy of original data abscissas, XGEOM as
C                          passed to PSFIT, for FUNCTION GETT, which
C                          is an argument of function  ZEROIN and is
C                          itself limited to a single argument.
C   TGEOM   R     I MAXPS  Values of the parametric variable at the
C                          original data points.
C   PSCOFS  R     I MAXPS  Cubic spline coefficients from CSFIT.
C                   *3*2   An array of size MAXPS*3 is required for
C                          XGEOM and another for YGEOM.
C   XLOC    R     I   -    Copy of abscissa at which PSEVAL is seeking
C                          an ordinate.
C   NCURV   I     I   -    Number of elements defining the curve.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C  HISTORY:
C     Jan.  83    RCL    Original design and coding.
C     Sept. 84    RCL    Now uses CSDVAL instead of IMSL's ICSEVU.
C
C  AUTHOR: Rosalie Lefkowitz, Sterling Software/NASA Ames, Moffett Field, CA
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      REAL
     >   GETT, T

C  *  Internal COMMON:

      INTEGER
     >   MAXPS
      PARAMETER
     >  (MAXPS = 1001)
      INTEGER
     >   NCURV
      REAL
     >   PSCURV, TGEOM, PSCOFS, XLOC(1)
      COMMON /PSWORK/
     >   PSCURV (MAXPS), TGEOM (MAXPS), PSCOFS (MAXPS,3,2), XLOC,    ! RLC 20211024 
     >   NCURV

C  *  Local variables:

      REAL
     >   DUM(1), XTEST(1), tlocal(1)         ! RLC 20211024 

C  *  Execution:

C  *  Evaluate spline TGEOM, PSCURV at T to get trial X, XTEST:

      tlocal(1)=t                                                     ! RLC 20211024 
      CALL CSDVAL (NCURV, TGEOM, PSCURV, 1, tlocal, PSCOFS (1,1,1),   ! RLC 20211024 
     >             PSCOFS (1,2,1), PSCOFS (1,3,1), XTEST, DUM, DUM)

C  *  Compare XTEST with the target XLOC:

      GETT = XLOC(1) - XTEST(1)                  ! RLC 20211024 

      RETURN
      END
C+---------------------------------------------------------------------
C
      FUNCTION TSUBJ (J)
C
C  TWO-LINER: Find T (J), where T is the parametric variable defined at
C             points (X (J), Y (J)) (probably cumulative chord length).
C
C  PURPOSE:
C             TSUBJ returns the value of the parametric variable set up
C          by PSFIT for the given point J.   It may be used to find the
C          interval of this variable needed by PSEVAL.
C
C  METHOD:
C             Simply use the internal COMMON set up by PSFIT.  The user
C          is thereby isolated from this COMMON.
C
C  ARGUMENTS:
C   ARG  TYPE  I/O/S  DIM  DESCRIPTION
C    J    I      I     -   Subscript of a point known to PSFIT at which
C                          the value of the parametric variable set  up
C                          by PSFIT is required.
C
C  INTERNAL COMMON BLOCK:
C           /PSWORK/ Work space set up by PSFIT/PSEVAL, q.v.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C  HISTORY:
C     Sept. 85  DAS  Introduced when cumulative chord length was
C                    substituted for index number as the parametric
C                    variable used by PSFIT/PSEVAL.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   J
      REAL
     >   TSUBJ

C  *  Internal COMMON:

      INTEGER
     >   MAXPS
      PARAMETER
     >  (MAXPS = 1001)
      INTEGER
     >   NCURV
      REAL
     >   PSCURV, TGEOM, PSCOFS, XLOC
      COMMON /PSWORK/
     >   PSCURV (MAXPS), TGEOM (MAXPS), PSCOFS (MAXPS,3,2), XLOC(1), 
     >   NCURV

C  *  Execution:

      TSUBJ = TGEOM (J)

      RETURN
      END
C+---------------------------------------------------------------------
C
      SUBROUTINE PSFIT (NGEOM, XGEOM, YGEOM, METHOD, CLOSED, IER)
C
C  ACRONYM: Parametric Spline FIT (for evaluating Y vs. X or X/Y vs. T)
C           -          -      ---
C  PURPOSE:
C                PSFIT fits an interpolating parametric cubic spline to
C             the given open or closed curve represented  by  abscissas
C             which are not necessarily monotonic, and by corresponding
C             ordinates.    (The parametric spline actually consists of
C             two splines, one for the abscissas and one for the ordin-
C             ates, each fitted against the same values of a parametric
C             variable - cumulative chord length - generated here.)
C
C                Evaluation for given abscissas of a specified subcurve
C             of the fitted parametric spline can be done using PSEVAL,
C             and probably TSUBJ, which returns the value of the  para-
C             metric variable T at a given data point J.  (PSEVAL needs
C             an interval in which to estimate the T which  corresponds
C             to the given X, then Y is computed from this T.)
C
C                Evaluation for a given value or values of T (a simpler
C             problem) is also provided for by subroutine PSTVAL.
C
C                This package was originally written for  working  with
C             airfoil sections with either sharp or rounded leading and
C             trailing edges.
C
C  METHOD:
C                The parametric variable is the usual approximation  to
C             arc length: cumulative chord length from the first point.
C
C                The splines for X vs. T and Y vs. T are interpolatory,
C             with a choice of conventional or monotonic.  (A monotonic
C             spline follows the data better in difficult cases such as
C             adjacent points with equal Y where flatness is  required.
C             However, monotonic splines are not appropriate for curva-
C             ture calculations, because the second derivatives are not
C             necessarily continuous.)
C
C                In either case,  the  curve  may be closed smoothly by
C             invoking periodic boundary conditions in the spline fits.
C
C                This implementation is intended to permit  evaluations
C             of Y for given values of either X or T.   It  necessarily
C             uses an internal COMMON area for the spline coefficients,
C             since this is the only way that the zero-finder  employed
C             by PSEVAL can get at the information it needs.  
C
C                One consequence is that PSFIT and PSEVAL or PSTVAL can
C             be used on only one curve at a time (unless copies of the
C             internal COMMON data  are kept - not recommended - appli-
C             cation programs should normally NOT reference /PSWORK/).
C
C  ARGUMENTS:
C   ARG       TYPE I/O/S  DIM    DESCRIPTION
C  NGEOM        I    I     -     Number of data points.
C  XGEOM,YGEOM  R    I   NGEOM   Coordinates of the data points.
C  METHOD     C*(*)  I     -     'C' or 'CC': conventional spline fits
C                                             using IENDL=IENDR=0;
C                                'M' or 'MM': monotonic spline fits.
C                                'MC': monotonic for X vs. T, but
C                                      conventional for Y vs. T;
C                                'CM': converse of 'MC'.
C  CLOSED       L    I     -     .TRUE. means the data points form a
C                                closed curve. In this case, the calling
C                                program must duplicate the first data
C                                point (XGEOM (1), YGEOM (1)) as
C                                (XGEOM (NGEOM), YGEOM (NGEOM)), and
C                                smoothness in the spline will be sought
C                                here as at all the other data points;
C                                .FALSE. means the end-points are treated
C                                as distinct, and are not joined (though
C                                in fact they MAY be the same point, as
C                                in the case of an airfoil with a sharp
C                                trailing edge).
C  IER          I    O     -     Error return code:
C                                IER = 0: No errors were detected;
C                                IER = 1: Too few data points:
C                                         (NGEOM < 2 for open curve or
C                                          NGEOM < 3 for closed curve);
C                                IER = 4: Closed curve case only:
C                                         First and last data points do
C                                         not match;
C                                IER = 5: Too many data points: NGEOM
C                                         exceeds MAXPS used in the
C                                         internal COMMON block.
C
C  INTERNAL COMMON BLOCK:
C   /PSWORK/ Work space required for parametric spline interpolation.
C            Row dimensions are set to MAXPS by PARAMETER statement.
C            (/PSWORK/ is normally NOT referenced by the calling program.)
C
C   VAR    TYPE I/O/S DIM  DESCRIPTION
C   PSCURV  R     O MAXPS  Copy of XGEOM for FUNCTION GETT, which is
C                          an argument of a zero-finder and is itself
C                          limited to a single argument.
C   TGEOM   R     O MAXPS  Values of parametric variable for the data
C                          points, generated here; available to the
C                          user via FUNCTION TSUBJ (to avoid COMMON).
C   PSCOFS  R     O MAXPS  Cubic spline coefficients from CSFIT, one
C                   *3*2   set for XGEOM vs. TGEOM, one for YGEOM vs.
C                          TGEOM.
C   NCURV   I     O   -    Copy of NGEOM.
C
C  PROCEDURES:
C   CHORD   Safeguarded calculation of arc lengths
C   CSFIT   Fits conventional interpolating spline
C   MSFIT   Fits monotonic interpolating spline
C
C  RELATED MODULES:
C   PSEVAL  Evaluates parametric spline set up by PSFIT for 
C           ordinates and derivatives at arbitrary abscissas (Xs)
C           on a specified monotonic sub-curve.  Uses FUNCTION
C           GETT internally, which should not concern the user.
C   TSUBJ   FUNCTION routine which may be used to determine intervals
C           defining subcurves for evaluation purposes.  Both GETT
C           and TSUBJ are in the same source module as PSEVAL.
C   PSTVAL  Evaluates parametric spline set up by PSFIT for given
C           value(s) of T (as opposed to X).  Returns X, Y, Y', Y''.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C  HISTORY:
C     Jan. 83  R.Lefkowitz  Original implementation using conventional
C                           splines from IMSL, which did not permit
C                           proper handling of closed curves.
C     Sep. 84    "    "     Handled closed curves properly via CSFIT.
C     Sep. 85  D.Saunders   Changed from index number to cumulative chord
C                           length for the parametric variable.
C     Jan. 87    "    "     Introduced monotonic spline option (prompted
C                           by a mesh generation requirement).
C     Sep. 87    "    "     Made it clear that /PSWORK/ is internal.
C     Oct. 91    "    "     Introduced CHORD for more careful arc-length
C                           calculations; added METHOD argument in place
C                           of NGEOM < 0 for indicating monotonic option;
C                           introduced IMPLICIT NONE; other cosmetics.
C     Aug. 93    "    "     Provided METHOD='MC' (and 'CM') option to
C                           accommodate redistribution of blunt-nosed
C                           airfoils where 'C' can allow the curve fit
C                           to exceed the data range at the leading edge.
C
C  AUTHOR: Rosalie C. Lefkowitz, Sterling Software/NASA Ames, Moffett Field, CA
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   NGEOM, IER
      REAL
     >   XGEOM (NGEOM), YGEOM (NGEOM)
      LOGICAL
     >   CLOSED
      CHARACTER
     >   METHOD * (*)

C     Internal COMMON.

      INTEGER
     >   MAXPS
      PARAMETER
     >  (MAXPS = 1001)
      INTEGER
     >   NCURV
      REAL
     >   PSCURV, TGEOM, PSCOFS, XLOC
      COMMON /PSWORK/
     >   PSCURV (MAXPS), TGEOM (MAXPS), PSCOFS (MAXPS,3,2), XLOC, NCURV

C     Local constants.

      REAL
     >   ZERO
      PARAMETER
     >  (ZERO = 0.E+0)

C     Local variables.

      INTEGER
     >   IENDL, J
      CHARACTER
     >   FIT * 2

C     Procedures.

      REAL
     >   CHORD
      EXTERNAL
     >   CHORD, CSFIT, MSFIT


C     Execution.

C     Too few data points?

      IF (NGEOM .LT. 2 .OR. (CLOSED .AND. NGEOM .LT. 3)) THEN
         IER = 1
         GO TO 999
      END IF

C     Too many data points?

      IF (NGEOM .GT. MAXPS) THEN
         IER = 5
         GO TO 999
      END IF


C     Before METHOD = 'MC' or 'CM' were allowed, applications are assumed
C     to have passed single-character strings:

      FIT = METHOD
      IF (LEN (METHOD) .EQ. 1) FIT (2 : 2) = FIT (1 : 1)


C     Set up parametric variable values, and transfer original abscissas
C     to internal COMMON for later use by zero-finder called by PSEVAL.

      NCURV = NGEOM
      TGEOM (1) = ZERO
      PSCURV (1) = XGEOM (1)

      DO 20, J = 2, NGEOM
         TGEOM (J) = CHORD (XGEOM, YGEOM, J-1, J) + TGEOM (J-1) 
         PSCURV (J) = XGEOM (J)
   20 CONTINUE


C     For conventional splines, smooth closure means periodic end conditions.

      IF (CLOSED) THEN
         IENDL = 4
      ELSE
         IENDL = 0
      END IF


C     Fit spline for XGEOM vs. TGEOM.

      IF (FIT (1 : 1) .EQ. 'C') THEN

C        Conventional spline:

         CALL CSFIT (NGEOM, TGEOM, XGEOM, IENDL, ZERO, 0, ZERO,
     >      PSCOFS (1,1,1), PSCOFS (1,2,1), PSCOFS (1,3,1), IER)

      ELSE

C        Monotonic spline: CLOSED used for "cyclic" is misleading, but valid.

         CALL MSFIT (NGEOM, TGEOM, XGEOM, CLOSED,
     >      PSCOFS (1,1,1), PSCOFS (1,2,1), PSCOFS (1,3,1), IER)
      END IF

      IF (IER .NE. 0) GO TO 999


C     Likewise for YGEOM vs. TGEOM.

      IF (FIT (2 : 2) .EQ. 'C') THEN

         CALL CSFIT (NGEOM, TGEOM, YGEOM, IENDL, ZERO, 0, ZERO,
     >      PSCOFS (1,1,2), PSCOFS (1,2,2), PSCOFS (1,3,2), IER)

      ELSE

         CALL MSFIT (NGEOM, TGEOM, YGEOM, CLOSED,
     >      PSCOFS (1,1,2), PSCOFS (1,2,2), PSCOFS (1,3,2), IER)

      END IF


  999 RETURN
      END
C+---------------------------------------------------------------------
C
      SUBROUTINE PSTVAL (NT, T, X, Y, XP, YP, XPP, YPP, NGEOM, XGEOM,
     >                   YGEOM)
C
C  ACRONYM:    Parametric Spline eVALuation for given T(s)
C              -          -      ----                 -
C  PURPOSE:    PSTVAL is an alternative to PSEVAL for the case when the
C              curve fitted by PSFIT is to be evaluated (X and Y) for a
C              given value or values of parametric variable T.  This is
C              in contrast to evaluation for a given X.
C
C              Evaluation for given T is much more straightforward than
C              for given X,  and may be appropriate for graphical work.
C
C              Provision is made for convenient evaluations at the data
C              points (no need to pass in the corresponding Ts).
C
C              The optional derivatives with respect to T will probably
C              be used for curvature calculations.  See also CURV2D and
C              XDERIVS.
C
C  METHOD:     The splines XGEOM vs. TGEOM and YGEOM vs. TGEOM, already
C              calculated by PSFIT, are evaluated at each of the values
C              of T supplied, or at TGEOM (1 : -NT) if NT < 0 is input.
C              First and second derivatives with respect to  T  (not X)
C              are also calculated.  As with PSEVAL, these are computed
C              in the order x", x', x and y", y', y.   If x" and y" are
C              not needed, pass the same arrays as for x' and y'. If x'
C              and y' are not needed either, pass x and y appropriately
C              to save storage.
C
C  ARGUMENTS:
C   ARG    TYPE  I/O/S   DIM   DESCRIPTION
C   NT      I      I      -    Number of evaluations required.  NT < 0
C                              is taken to mean the evaluations are to
C                              be at the data points 1:-NT.  |NT|>=1.
C   T       R      I      NT   Arbitrary values of parametric variable
C                              T for which evaluations are desired.  If
C                              NT < 0, the existing values of T for the
C                              data points (TGEOM(*)) stored internally
C                              by PSFIT are used - no need to rederive
C                              them in the application.
C   X       R      O      NT   Corresponding values of X on the curve.
C                              May be same locations as T(*).
C   Y       R      O      NT   Corresponding values of Y on the curve.
C   XP      R      O      NT   First derivatives, DX/DT.
C   YP      R      O      NT   First derivatives, DY/DT.
C   XPP     R      O      NT   Second derivatives of X and Y w.r.t. T.
C   YPP     R      O      NT   (See above if these are not needed.)
C   NGEOM   R      I      -    Number of data points fitted by PSFIT.
C   XGEOM   R      I    NGEOM  Geometry abscissas fitted by PSFIT.
C   YGEOM   R      I    NGEOM  Geometry ordinates fitted by PSFIT.
C
C  INTERNAL COMMON BLOCK:
C   /PSWORK/ Work arrays, set up by PSFIT.
C            Dimensions are set to MAXPS by PARAMETER statement.
C   VAR    TYPE I/O/S DIM  DESCRIPTION
C   TGEOM   R     I MAXPS  Values of the paramateric variable at the
C                          data points, set up by PSFIT.
C   PSCOFS  R     I MAXPS  Cubic spline coefficients from CSFIT.
C                   *3*2   An array of size MAXPS*3 is required for
C                          XGEOM and another for YGEOM.
C   XLOC    R     S   -    Copy of X, abscissa at which PSEVAL is seeking
C                          an ordinate.
C   NCURV   I     O   -    Copy of NGEOM.
C
C  ERROR HANDLING: None.
C
C  EXTERNAL REFERENCES:
C   CSDVAL  Evaluates a 1-dimensional cubic spline along with its first
C           and second derivatives.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C  HISTORY:
C  Sept. 85  DAS  Adapted from PSEVAL for the simpler situation,
C                 in case it proves handy some time.
C  April 86  DAS  Allowed for caller to overwrite T(*) with X(*).
C  10/16/91  DAS  Derivatives returned are with respect to T, not X,
C                 now (to facilitate curvature calculations).  Gave
C                 up on allowing T to be overwritten with X, to
C                 avoid 2*NT calls to CSDVAL where two calls will do.
C  10/25/91  DAS  Derivatives at the data points required the internal
C                 TGEOM values.  Use of TSUBJ for all of them would
C                 be ugly.  Thus NT < 0 now avoids having to input Ts.
C  04/01/95  DAS  MAXPS had not been raised to 1001 to match PSFIT.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   NT, NGEOM
      REAL
     >   T (*), X (*), XP (*), XPP (*), Y (*), YP (*), YPP (*),
     >   XGEOM (NGEOM), YGEOM (NGEOM)

C  *  Internal COMMON:

      INTEGER
     >   MAXPS
      PARAMETER
     >  (MAXPS = 1001)
      INTEGER
     >   NCURV
      REAL
     >   PSCURV, TGEOM, PSCOFS, XLOC
      COMMON /PSWORK/
     >   PSCURV (MAXPS), TGEOM (MAXPS), PSCOFS (MAXPS,3,2), XLOC, NCURV

C  *  Local variables:

      INTEGER
     >   NPTS

C  *  Execution:

      NPTS = ABS (NT)

      IF (NT .GT. 0) THEN

C  *     Evaluate the XGEOM vs. TGEOM spline at each input T:

         CALL CSDVAL (NGEOM, TGEOM, XGEOM, NPTS, T,
     >                PSCOFS (1,1,1), PSCOFS (1,2,1), PSCOFS (1,3,1),
     >                X, XP, XPP)

C  *     ... and the YGEOM vs. TGEOM spline at each input T:

         CALL CSDVAL (NGEOM, TGEOM, YGEOM, NPTS, T,
     >                PSCOFS (1,1,2), PSCOFS (1,2,2), PSCOFS (1,3,2),
     >                Y, YP, YPP)

      ELSE  ! NT < 0 means evaluate at the original data points:

C        First, the X derivatives:

         CALL CSDVAL (NGEOM, TGEOM, XGEOM, NPTS, TGEOM,
     >                PSCOFS (1,1,1), PSCOFS (1,2,1), PSCOFS (1,3,1),
     >                X, XP, XPP)

C  *     ... and now the Y derivatives:

         CALL CSDVAL (NGEOM, TGEOM, YGEOM, NPTS, TGEOM,
     >                PSCOFS (1,1,2), PSCOFS (1,2,2), PSCOFS (1,3,2),
     >                Y, YP, YPP)

      END IF

      RETURN
      END
*DECK QNMDIF
C+----------------------------------------------------------------------
C
      SUBROUTINE QNMDIF (N, NLDIM, NFTOTL, NITER, NTYPE, LWRITE,
     >   X, F, FM, G, H, EL, D, ETA, TOL, STEPMX, EPSMCH, EPSOBJ,
     >   UNITL, LOCAL, CONV, CONTIN, RESCUE, PRINT, FUN, USER)
C
C
C     Description and usage
C     ---------------------
C
C        QNMDIF uses a quasi-Newton method with difference
C     approximations to derivatives to minimize a function F(X)
C     of the N independent variables X1,X2,...,XN.  The method is
C     iterative, and so requires an initial estimate of the position
C     of the minimum.  It is intended for functions which are 
C     continuous and have continuous first and second derivatives in
C     the region between the initial estimate and the minimum, although
C     it may work if the derivatives have occasional discontinuities.
C
C        The user is required to provide a subroutine to evaluate
C     the function to be minimized.  The first derivatives are
C     calculated by differencing function values.
C
C        QNMDIF is a revised version of the program given in the NPL
C     report by Gill, Murray, and Pitfield, "The Implementation of
C     Two Revised Quasi-Newton Algorithms for Unconstrained
C     Optimization", National Physical Laboratory, Div. of Numerical
C     Analysis and Computing, Report No. DNAC 11, 1972.
C
C        From an initial point it generates a sequence of points which
C     is intended to converge to a local minimum of F(X).  A typical
C     iteration starts at the point X, with an estimate G of the
C     gradient vector, and a lower triangular matrix EL and a diagonal
C     matrix D such that EL * D * ELT is a positive definite approximation
C     to the current matrix of second derivatives (Hessian matrix).
C     (note that ELT denotes the transpose of EL).  The linear equations
C
C              EL * D * ELT  *  P = -G
C
C     are solved to give an estimate P of the direction of search.
C     A scalar ALPHA is found so that X + ALPHA * P is approximately a
C     minimum with respect to ALPHA, and then X is replaced by
C     X + ALPHA * P.  The matrices EL and D are then updated according
C     to the complementary Davidon-Fletcher-Powell update formula.
C     Most iterations use forward difference approximations to compute
C     G, but central differences are used when they seem necessary.
C
C        QNMDIF is the main subroutine in a set of routines that
C     includes the following six subroutines:  QNMDIF, LINSCH, OUTPUT,
C     APROXG, UPCHOL, and C1MDCH.  In addition, the user must supply
C     a subroutine called USER which may be a dummy - see below.
CRAK  A second group of subroutines comprised of CENDIF, FDSTEP, and
CRAK  DIFFER may be used to estimate good finite difference step sizes.
C
C        The subroutines QNMDIF, LINSCH, and UPCHOL contain locally
C     dimensioned arrays, which must all be of length at least N.  In
C     the current version, these arrays are of length MAXN, a PARAMETER.
C     QNMDIF and companion routine CENDIF both check that N <= MAXN.
C
CRAK     NOTE: This particular version of QNMDIF has been revised for
CRAK  expensive objective functions which cannot be evaluated to full
CRAK  machine precision.
C
C
C
C     Parameters
C     ----------
C
C     N - The number of independent variables, an integer constant.
C        Present maximum value is 40, set by parameter MAXN, below.
C
C     NLDIM - The value N * (N-1)/ 2, which is the number of elements
C        in the matrix EL, stored as a vector.  Note, however, that
C        NLDIM must always be at least 1, so that if QNMDIF is called
C        with N = 1, NLDIM should be set to 1.
C
C     NFTOTL - An integer, which is set by QNMDIF on exit to the total
C        number of calls to the objective function executed during the
C        location of the minimum.  It need not be initialized.
CRAK
CRAK     MODIFICATION NOTE:  NFTOTL must now be initialized by the
CRAK     user, even if only to zero.  When F, G, EL, D are supplied
CRAK     (i.e. UNITL = .FALSE.), NFTOTL may be initialized to reflect
CRAK     the work done prior to calling QNMDIF.
CRAK
C
C     NITER - An integer, which is set by QNMDIF to the number of
C        iterations executed during the location of the minimum,
C        where an iteration is the procedure described in the comments
C        at the beginning of QNMDIF.  It need not be initialized.
CRAK
CRAK     MODIFICATION NOTE:  NITER may now be used to limit the number
CRAK     of opimization iterations to be performed by QNMDIF.  Set
CRAK     NITER to a large positive value in calling routine to disable
CRAK     the test (and hence revert to original mode).
CRAK
CRAK
CRAK  NTYPE - An integer flag set by QNMDIF to indicate whether the
CRAK     gradient evaluated before returning used forward (NTYPE = 0)
CRAK     or central differences (NTYPE = 1).  On input, must be set to
CRAK     correspond to the user-supplied gradient if UNITL = .FALSE.
CRAK     or to the desired initial gradient type if starting from
CRAK     scratch (UNITL = .TRUE.).
CRAK
CRAK
CRAK  LWRITE - Required integer input denoting the unit number of
CRAK     the primary output device (often device 6).
CRAK
C
C     X - The array of independent variables, a vector of length N.
C        On entry, it should be set to the initial estimate of the
C        location of the minimum.  On output it will contain the best
C        estimate found of the minimum.
C
CRAK
CRAK  F - A real scalar.  On entry, it must be set to initial function
CRAK     value.  On exit it will be the lowest value found of the
CRAK     objective function.
CRAK
CRAK
CRAK  FM - A real scalar.  On entry, it should be set to a lower
CRAK     bound on the value of the function at the minimum, or to
CRAK     a large negative value if no estimate is available.
CRAK
CRAK
CRAK  G - A real array of length N, containing the gradient, which
CRAK     must be initialized if UNITL = .FALSE., and otherwise not.
CRAK     On exit, equals the last gradient evaluated, which may be
CRAK     useful.
CRAK
C
C     H - A difference array, of length N, where H (I) is the interval
C        used to compute the partial derivative of F with respect to
C        X (I).  Great care should be used in choosing the H array, and
C        the user should consider that the smaller the value of H, the
C        less error is incurred in computing the derivative when exact
C        arithmetic is used, since the error in the forward difference
C        approximation is of order H (using the Taylor series expansion);
C        on the other hand, the smaller H is, the more cancellation
C        error will occur, since values of the function computed at close
C        points will have many matching figures, which results in a loss
C        of accuracy in the computed derivative when these close values
C        are subtracted.  Each H (I) should be chosen small enough for the
C        difference estimates to be close to the true derivatives, but
C        not so small that cancellation error destroys the accuracy.  It
C        is suggested that the user experiment with different values of
C        H (I) before using this subroutine, to see how the function
C        values behave for various sizes of H (I).
C
C
C     EL - A vector of length NLDIM (at least 1), used to hold the
C        lower triangular factor of the approximate Hessian.  If the
C        user wishes to provide an initial estimate of the Hessian, EL
C        should be filled with the elements of the lower triangle of the
C        Cholesky factorization EL * D * ELT, where ELT is EL tranpose,
C        stored by rows.  Otherwise, set UNITL to .TRUE., and QNMDIF will
C        use the identity matrix as the initial approximation (i.e., the
C        first iteration will be steepest descent).  On exit, EL contains
C        the lower triangle of the factors of the final approximate Hessian.
C
C
C     D - A vector of length N, used to hold the diagonal matrix of the
C        factors of the approximate Hessian.  If the user wishes to
C        provide an initial estimate of the Hessian, D should be set in
C        the same way as EL (described above).  Otherwise, when UNITL is
C        set to .TRUE., the identity is used for the first iteration, and
C        the user need not set EL and D. On exit, D contains the diagonal
C        factor of the final approximate Hessian.
C
C
C     ETA - A scalar constant, set by the user to a number strictly
C        between 0 and 1.  It controls the accuracy of the search for the
C        minimum of the function along the current direction of search,
C        in the following way:  if the exact minimum along the line were
C        found at each step, the projected gradient (gradient of the
C        function along the line) would be zero.  However, it has been
C        found that finding the exact minimum at each step is not
C        computationally efficient, in that it requires many extra
C        function evaluations, and does not in general assist the speed
C        of convergence dramatically.  The parameter ETA is multipled by
C        the projected gradient at the start of each minimization along a
C        line.  When the projected gradient at a point is less than ETA
C        times the initial value, the linear search is terminated.  Thus,
C        a small value of ETA means that the projected gradient must be
C        very small, i.e., close to the minimum along the line.  A value
C        of ETA close to 1 means that the projected gradient need
C        decrease only slightly from the starting point, so that the
C        resulting point is only a rough approximation to the minimum
C        along the line.  A guide to its value can be found in the NPL
C        report mentioned above, but the user may wish to try different
C        values on experimental runs to find the best value for his problem.
C        A value of 0.2E+0 is generally acceptable for QNMDIF.
C
C
C     TOL - A real positive number, which specifies the accuracy to
C        which the minimum is to be found.  QNMDIF terminates when the
C        norm of the approximate gradient is less than TOL ** (1/ 3) and
C        the distance moved during the last iteration is less than a
C        number close to TOL.  TOL should be chosen as approximately
C        the Euclidean distance which can be allowed between the point
C        returned by QNMDIF and the actual minimum.  Its value should
C        be chosen with care, keeping in mind that the error in the
C        approximate gradient is of order MAXH (I) ** 2 at best (when
C        central differences are used), where MAXH is the maximum
C        component of H in absolute value, that the function may be
C        known only to a given number of figures, and that the
C        parameters may be known only to a given number of figures
C        (for example, when the problems arises from a physical
C        situation, or a problem involving measured data).  A lower
C        bound for an acceptable value of TOL for most problems
C        is SQRT (EPSOBJ), but a larger value may be required.
C
C
C     STEPMX - The maximum allowed step from any point to the minimum
C        along a line from that point.  this parameter should be used
C        when the user is confident that the parameter vector should
C        not move more than the distance STEPMX from any particular
C        location.  If the linear search fails to find a minimum closer
C        to the current point than STEPMX, then the distance moved is
C        STEPMX.  It should be set to an upper bound on the distance
C        between the starting point and the minimum.  If the user does
c        not have a reasonable guess for this value, it can simply be
C        set to a very large number (say 1.0E+10), and then will have
C        essentially no effect on the linear search.
C
C
C     EPSMCH -  The machine precision of the computer to be used, i.e.,
C        the smallest value such that, in floating point arithmetic,
C        1 + EPSMCH > 1.  on the IBM 360, EPSMCH is 16.0 ** (-13).
C
CRAK
CRAK  EPSOBJ - Must be set by calling routine to value of smallest
CRAK     significant (absolute) change in objective function.  QNMDIF
CRAK     resets EPSOBJ to EPSMCH if input value is less than EPSMCH.
CRAK
C
C     UNITL - A logical parameter, set by the user.  In general, the
C        user will not know an approximation to the Hessian matrix
C        at the starting point, and in this event UNITL should be set
C        to .TRUE., so that the identity matrix is used as the starting
C        approximation to the Hessian.  If, however, the user does have
C        an idea of the approximate Hessian, UNITL should be set to
C        .FALSE., and the EL and D arrays should be set to the Cholesky
C        factors of the initial approximation (see comments under EL
C        and D, above).
CRAK
CRAK     MODIFICATION NOTE:  UNITL = .FALSE now implies that all of the
CRAK     information needed to begin an iteration (F, X, G, EL, and D)
CRAK     has been supplied (not merely EL and D).
CRAK
C
C     LOCAL - A logical parameter, to be set to .TRUE. or .FALSE. by the
C        user.  Gill and Murray have developed a sophisticated (?) "local
C        search" which follows the normal minimization, and which checks
C        through a random search whether a point can be found lower than
C        the estimated minimum.  The purpose of this search is to avoid
C        convergence to a saddle point, where the gradient is zero, but
C        which is not a minimum.  This search requires a fair number of
C        extra function evaluations, but will essentially guarantee
C        convergence to a true minimum.  If the user feels that there is
C        no danger of convergence to a saddle point, and in particular,
C        if the function is extremely expensive to evaluate, then LOCAL
C        should be set to .FALSE., and the local search will not be
C        executed.
C
C
C     CONV - A logical flag set by QNMDIF to .TRUE. if there has been
C        satisfactory convergence to a minimum, and to .FALSE. if not.
C        If LOCAL is .TRUE., then CONV = .TRUE. means that the local
C        search failed to find a significantly lower point.  if LOCAL
C        was set to .FALSE., then CONV = .TRUE. means that the gradient
C        at the final point was essentially zero, and that the steps to
C        the final point had converged.  Note that if LOCAL is set to
C        .FALSE., CONV = .TRUE. can happen if there is convergence to
C        a saddle point, since there is no way to check the second order
C        conditions to distinguish a minimum from a saddle point.
C
CRAK
CRAK  CONTIN - A logical input flag which must be .TRUE. if QNNDIF is
CRAK     to update the gradient and Hessian before returning due to
CRAK     iteration count.  Efficiency is not compromised if the final
CRAK     values are saved and later used by the calling routine for
CRAK     restarting.
CRAK
CRAK
CRAK  RESCUE - A logical input flag which may be set .TRUE. if the user
CRAK     wishes to begin the run with the local search procedure, which
CRAK     counts as two iterations.  Normally, use RESCUE = .FALSE. and set
CRAK     LOCAL = .TRUE. if some assurance against premature convergence
CRAK     is desired.  The RESCUE option may be useful if the objective
CRAK     function routine has blundered into a bad region.
CRAK
C
C     PRINT - A logical flag, set by the user to control printout.
C        for PRINT = .FALSE., most printing is suppressed.
C
C
C     FUN - A subroutine provided by the user, which must be declared
C        as EXTERNAL in the calling subprogram.  Its call must be of
C        the form
C
C        CALL FUN (N, X, FVAL)
C
C        where N is the number of variables, X is the vector of
C        variables and FVAL is the computed value of the function.
C
CRAK
CRAK  USER - A subroutine provided by the user, which must be declared
CRAK     as EXTERNAL in the calling subprogram.  Its call must be of
CRAK     the form
CRAK
CRAK     CALL USER (IENTRY, N, NLDIM, X, F, G, H, EL, D, NFTOTL,
CRAK    >   NITER, NTYPE, CONTIN)
CRAK
CRAK     where the variables are the current values of those described
CRAK     above.  QNMDIF calls USER at the end of each successful line
CRAK     search with IENTRY = 1, at which point X, F, H, NFTOTL, and
CRAK     NITER are current, but G, EL, and D have not yet been updated.
CRAK     This permits certain housekeeping chores to be performed before
CRAK     the calculations for the next iteration are begun.  QNMDIF sets
CRAK     IENTRY = 2 and calls USER again following calculation of gradient
CRAK     and Hessian updates, as well as after the local search if
CRAK     requested.  Disk files for saving restart information should be
CRAK     written when IENTRY = 2, with an initial REWIND command.
CRAK     Different applications may dictate other choices for parameter
CRAK     list and calling points.  Note that either or both of these
CRAK     options may be simply turned off by providing subroutine USER
CRAK     with RETURN statements for the unneeded cases.
CRAK
CRAK
C
C     Development history
C     -------------------
C
C        Date           Initials   Description
C            ?   1972              Initial coding in ALGOL at the NPL
C                                     (Gill, Murray, and Pitfield)
C            ?                     FORTRAN translation, Stanford Univ.
C                                     (Margaret Wright)
C        20 Sep. 1982   RAK        Version 1.1, modified as noted for
C                                     NASA release
C        13 Aug. 1986   RAK        Version 1.2 - mostly cosmetic mods.
C                                     Included CENDIF/FDSTEP/DIFFER for
C                                     stepsize selection (though they must
C                                     be called separately).  Added PARAMETER
C                                     statements for local array dimensions,
C                                     with test that N <= MAXN on entry.
C                                     Increase MAXN from 30 to 40.  Some
C                                     protection added in LINSCH.
C
C
C     Author
C     ------
C
C        Robert Kennelly
C        NASA - Ames Research Center
C        Mail Stop 227-6
C        Moffett Field, CA  94035
C
C        Telephone:  (415) 694-5944
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT REAL (A - H, O - Z)

C     Constants.

      INTEGER
     >   MAXN
      PARAMETER
     >   (MAXN = 40)

C     Arguments.

      LOGICAL
     >   CONTIN, CONV, COUNT, FINAL, LOCAL, PRINT, RESCUE, SUCCES, UNITL
      REAL
     >   D(N), EL(NLDIM), G(N), H(N), X(N)

C     Local variables.

      LOGICAL
     >   DONE
      REAL
     >   OLDG(MAXN), P(MAXN), PP(MAXN), Y(MAXN)

C     Procedures.

      EXTERNAL
     >   FUN, USER


C     Execution.
C     ----------

      IF (N .GT. MAXN) STOP 'QNMDIF:  ERROR - too many variables!'


C     Initialization.
C     ---------------

CRAK  Save input value of NITER for use as limit.

      NITMAX = NITER
      NITER = 0

C     If UNITL is .TRUE., use identity matrix as approximation to Hessian
C     for the first iteration (steepest descent).

      FNEW = F
      OLDF = F
      IF (.NOT.UNITL) GO TO 30
         DO 10 I = 1, N
            G(I) = 0.0E+0
            D(I) = 1.0E+0
   10    CONTINUE
         DO 20 I = 1, NLDIM
            EL(I) = 0.0E+0
   20    CONTINUE
         IF (.NOT.RESCUE) CALL APROXG (N, NTYPE, X, G, H, PP, FUN,
     >      FNEW, NFTOTL, PRINT, LWRITE)
   30 CONTINUE
      IF (EPSOBJ .LT. EPSMCH) EPSOBJ = EPSMCH
      RTEPS = SQRT (EPSOBJ)
      TINY = EPSMCH ** 2
      CONV = .TRUE.
      FINAL = .FALSE.
      COUNT = .FALSE.
      BOUNDK = 0.01E+0/ (SQRT (REAL (N)) * EPSMCH)

C     Calculate maximum component of vector H and ALPMIN, which
C     is used as a measure of the smallest step that should be
C     taken along any direction.

      HMAX = 0.0E+0
      DO 40 I = 1, N
         U = ABS (H(I))
         IF (U .GT. HMAX) HMAX = U
         P(I) = 0.0E+0
   40 CONTINUE
      ALPMIN = HMAX ** (2.0E+0/ 3.0E+0)

CRAK  Start with local search if rescue option selected.

      IF (RESCUE) GO TO 290


C     Main iteration loop.
C     --------------------

   50 CONTINUE

CRAK  Save current variable and gradient vectors.

      OLDF = FNEW
      DO 60 I = 1, N
         OLDG(I) = G(I)
         Y(I) = X(I)
   60 CONTINUE


C     Search direction.
C     -----------------

C     Calculate the direction of search, P, by solving the system of
C     linear equations EL * D * ELT * P = -G, where EL is a unit lower
C     triangular matrix , ELT is its transpose, and D is a diagonal
C     matrix represented in the program by the vector D.

C     Forward solution.

      P(1) = - OLDG(1)
      IF (N .EQ. 1) GO TO 90
         IR = 1
         DO 80 I = 2, N
            SUM = - OLDG(I)
            IT  = I - 1
            DO 70 K = 1, IT
               SUM = SUM - P(K) * EL(IR)
               IR  = IR + 1
   70       CONTINUE
            P(I) = SUM
   80    CONTINUE
   90 CONTINUE

C     Back substitution.

      P(N) = P(N)/ D(N)
      IF (N .EQ. 1) GO TO 120
         DO 110 I = 2, N
            IN = N + 1 - I
            IR = IR - 1
            IS = IR
            IT = IN + 1
            SUM = P(IN)/ D(IN)
            DO 100 K = IT, N
               KN = N + IT - K
               SUM = SUM - P(KN) * EL(IS)
               IS = IS + 2 - KN
  100       CONTINUE
            P(IN) = SUM
  110    CONTINUE
  120 CONTINUE

C     Compute the norm of P (the direction of search), the norm of
C     G (the approximate gradient vector), and the dot product
C     of P and G (which should be negative if P is a descent
C     direction).

      SUM = 0.0E+0
      GTP = 0.0E+0
      GNM = 0.0E+0
      DO 130 I = 1, N
         SUM = SUM + P(I) ** 2
         GTP = GTP + G(I) * P(I)
         GNM = GNM + G(I) ** 2
  130 CONTINUE
      PNORM = SQRT (SUM) + TINY
      GNM = SQRT (GNM) + TINY

CRAK  Protect GTP from being unreasonably small.

      IF (ABS (GTP) .LT. TINY) GTP = SIGN (TINY, GTP)


C     Line search.
C     ------------

C     Initialize ALPHA, the step to be taken along the direction of
C     search P.  The initial value is input to the linear search
C     subroutine, which estimates the approximate minimum of the
C     objective function along P.

CRAK  EPSOBJ has been substituted for SQRT (EPSMCH) as a tolerance
CRAK  here.  In all other appearances of EPSOBJ, it was substituted
CRAK  directly for EPSMCH.

      DELTAF = MAX (ABS (FNEW - FM), EPSOBJ)
      ALPHA = MIN (- 2.0E+0  *  DELTAF/  GTP, 1.E+0)

CRAK  Line search tolerance (T) is scale dependent.  For problems with
CRAK  second derivatives very different from order one (especially
CRAK  where desired accuracy is limited by available objective function
CRAK  precision), it may be desirable to try other values.  The original
CRAK  version of QNMDIF used T = RTEPS/ PNORM; we use EPSOBJ/ PNORM to
CRAK  to force LINSCH to work a bit harder.

      T = EPSOBJ/ PNORM
      SFTBND = 0.E+0
      IF (NTYPE .EQ. 0) SFTBND = ALPMIN/ PNORM

C     Find suitably lower point along direction of search.  If the linear
C     search is successful, the vector X and the function value FNEW are
C     modified within the subroutine LINSCH.  If the search fails, they
C     remain unchanged.

      CALL LINSCH (N, FUN, RTEPS, T, ETA, SFTBND, STEPMX/ PNORM, P,
     >    X, FNEW, ALPHA, GTP, NFTOTL, SUCCES, PRINT, LWRITE)
      IF (SUCCES .OR. NTYPE.EQ.1) NITER = NITER + 1
      COUNT = NITER .GE. NITMAX
      IF (SUCCES) GO TO 150
         IF (PRINT) WRITE (LWRITE, 1000)
 1000    FORMAT (34H0QNMDIF:  UNSUCCESSFUL LINE SEARCH)

C        If central differences were used to approximate the gradient,
C        branch to failure exit.

         IF (NTYPE .EQ. 1) GO TO 290

CRAK     Start over using central differences if first line search based
CRAK     on forward difference input derivative data failed (since array
CRAK     PP not initialized).

         IF (NITER.GT.0 .OR. UNITL) GO TO 140
            NTYPE = 1
            CALL APROXG (N, NTYPE, X, G, H, PP, FUN, FNEW, NFTOTL,
     >         PRINT, LWRITE)
            GO TO 50
  140    CONTINUE

C        Otherwise, switch to central difference approximation (using
C        intermediate values stored in array PP) and try again.

         NTYPE = 1
         CALL APROXG (N, 2, X, G, H, PP, FUN, FNEW, NFTOTL,
     >      PRINT, LWRITE)
         GO TO 50
  150 CONTINUE

C     The linear search was successful - print if required.

      IF (.NOT.PRINT) GO TO 160
         WRITE (LWRITE, 1010)
 1010    FORMAT (1H1)
         CALL OUTPUT (N, NLDIM, Y, G, H, P, OLDF, ALPHA, D,
     >      EL, NITER, NFTOTL, FINAL, LWRITE)
         WRITE (LWRITE, 1010)
  160 CONTINUE

CRAK  Invoke user-supplied routine for housekeeping, then jump out if
CRAK  maximum number of iterations have been performed and final
CRAK  derivative information was not requested.  Note that the values
CRAK  printed out in "final information" will not be up-to-date.

      IENTRY = 1
      CALL USER (IENTRY, N, NLDIM, X, FNEW, G, H, EL, D, NFTOTL,
     >    NITER, NTYPE, CONTIN)
      IF (COUNT .AND. .NOT.CONTIN) GO TO 290


C     Calculate gradient.
C     -------------------

C     Check whether central differences were used for the gradient.
C     If so, and if the step taken was sufficiently large, switch
C     back to forward differences.

      IF (NTYPE .EQ. 1 .AND. ALPHA .GT. ALPMIN/ PNORM)  NTYPE = 0

C     Calculate the new approximate gradient.

      CALL APROXG (N, NTYPE, X, G, H, PP, FUN, FNEW, NFTOTL,
     >   PRINT, LWRITE)
      GNMSQ = 0.0E+0
      DO 170 I = 1, N
         GNMSQ = GNMSQ + G(I) ** 2
  170 CONTINUE

C     If we are using forward differences, and the norm of the gradient
C     is small, recalculate the gradient using central differences.

      IF (GNMSQ .GT. HMAX .OR. NTYPE .NE. 0) GO TO 190
         NTYPE = 1
         CALL APROXG (N, 2, X, G, H, PP, FUN, FNEW, NFTOTL,
     >      PRINT, LWRITE)
         GNMSQ = 0.0E+0
         DO 180 I = 1, N
            GNMSQ = GNMSQ + G(I) ** 2
  180    CONTINUE
  190 CONTINUE


C     Update Hessian.
C     ---------------

C     Update the Cholesky factors of the Hessian approximation using
C     the rank-two COMDFP (Complementary Davidon-Fletcher-Powell)
C     update formula, providing that new GTP > old GTP.

      GTPFIN = 0.E+0
      DO 200 I = 1, N
         GTPFIN = GTPFIN + G(I) * P(I)
  200 CONTINUE
      V = ALPHA * (GTPFIN - GTP)
      IF (V .GT. 0.E+0) GO TO 210
         IF (PRINT) WRITE (LWRITE, 1020) NITER, GTP, GTPFIN
 1020    FORMAT (37H0QNMDIF:  WARNING - ON ITERATION NO. , I6,
     >      29H, ABS (GTP) DID NOT DECREASE./
     >      10X, 15HGTP(INITIAL) = , E25.15/
     >      10X, 15HGTP(FINAL)   = , E25.15/
     >      10X, 37HUPDATES TO THE HESSIAN WERE BYPASSED.)
         GO TO 270
  210 CONTINUE
         V = SQRT (V)
         DO 220 I = 1, N
            P(I) = (G(I) - OLDG(I))/ V
  220    CONTINUE

C        Update the Cholesky factors after a positive rank-one change.

         CALL C1MDCH (N, NLDIM, D, EL, 1.0E+0, P, IFAIL)
         IF (IFAIL .NE. 0) GO TO 290
         V = SQRT (ABS (GTP))
         DO 230 I = 1, N
            OLDG(I) = - OLDG(I)/ V
  230    CONTINUE

C        Update Cholesky factors following negative rank-one correction.

         CALL UPCHOL (N, NLDIM, EPSMCH, D, EL, OLDG)

C        Check whether the diagonal elements of the updated factors are
C        sufficiently large.  If not, modify them appropriately to ensure
C        that the approximate Hessian is sufficiently positive definite.

         DMAX = D(1)
         IF (N .EQ. 1) GO TO 250
            DO 240 I = 2, N
               IF (D(I) .GT. DMAX)  DMAX = D(I)
  240       CONTINUE
  250    CONTINUE
         BL = DMAX/ BOUNDK
         DO 260 I = 1, N
            IF (D(I) .LT. BL)  D(I) = BL
  260    CONTINUE
  270 CONTINUE

CRAK  Invoke user-supplied routine for housekeeping.

      IENTRY = 2
      CALL USER (IENTRY, N, NLDIM, X, FNEW, G, H, EL, D, NFTOTL,
     >   NITER, NTYPE, CONTIN)


C     Check convergence.
C     ------------------

      SUM = 0.E+0
      DO 280 I = 1, N
         SUM = SUM + X(I) ** 2
  280 CONTINUE

C     Check whether norm of gradient and last step taken were sufficiently
C     small.  If not, verify that objective decreased and check iteration
C     count before continuing.

      IF (GNMSQ.LT.TOL ** .667E+0 .AND.
     >   ALPHA * PNORM.LT.SQRT (EPSMCH * SUM) + TOL) GO TO 300
      IF (OLDF.GT.FNEW .AND. .NOT.COUNT) GO TO 50

CRAK  Exits from main loop:
CRAK
CRAK     290 - Failure, reset CONV
CRAK     300 - Success, leave CONV set to .TRUE.

  290 CONTINUE
      CONV = .FALSE.
  300 CONTINUE

CRAK  Perform local search if requested, but check iteration count.

      IF (.NOT.(LOCAL.OR.RESCUE) .OR. COUNT) GO TO 500


C     Local search.
C     -------------

C     Carry out local search, which is a procedure of searching along
C     random orthogonal directions to see if a lower point can be found.

      IF (PRINT) WRITE (LWRITE, 1030)
 1030 FORMAT (38H0QNMDIF:  BEGIN LOCAL SEARCH PROCEDURE)
      U = ALPMIN
      EPS10 = 10.E+0 * EPSOBJ * ABS (FNEW)

C     Modify independent variable vector with a small change.

  310 CONTINUE
      DO 320 I = 1, N
         Y(I) = X(I) + U
  320 CONTINUE
      CALL FUN (N, Y, V)
      NFTOTL = NFTOTL + 1

C     Check whether a sufficiently large change in the function occurred.
C     If not, increase step size and try again.

      IF (ABS (V - OLDF) .GE. EPS10) GO TO 330
         U = 5.0E+0 * U
         GO TO 310
  330 CONTINUE

C     Calculate orthogonal direction at the modified point, W.

      P(1) = U
      IF (N .EQ. 1) GO TO 350
         DO 340 I = 2, N
            P(I) = - P(I - 1)
  340    CONTINUE
         IF (MOD(N, 2) .EQ. 1)  P(N) = 0.0E+0
  350 CONTINUE
      U = SQRT (ALPMIN)
      EPS10 = EPS10 * ABS (V)

CRAK  Begin first "iteration" (rewritten without loop - RAK, 16 June 81)

      DO 360 I = 1, N
         OLDG(I) = Y(I) + U * P(I)
  360 CONTINUE
      CALL FUN (N, OLDG, OLDF)
      NFTOTL = NFTOTL + 1
      GTP = OLDF - V

C     Ensure that the step has made a significant change in the function
C     value - if not, increase step size.

      IF (ABS (GTP) .GE. EPS10) GO TO 370
         U = 5.0E+0 * U
         GO TO 310
  370 CONTINUE
      GTP = - ABS (GTP/ U) - 1.0E+0
      ALPHA = U
      IF (OLDF .LE. V) GO TO 390
         V = OLDF
         DO 380 I = 1, N
            Y(I) = OLDG(I)
            P(I) = - P(I)
  380    CONTINUE
  390 CONTINUE
      SUM = 0.0E+0
      DO 400 I = 1, N
         SUM = SUM + P(I) ** 2
  400 CONTINUE
      PNORM = SQRT (SUM) + TINY
      SFTBND = 0.0E+0
      T = RTEPS/ PNORM
      CALL LINSCH (N, FUN, RTEPS, T, EPSMCH, SFTBND, STEPMX/ PNORM,
     >   P, Y, V, ALPHA, GTP, NFTOTL, SUCCES, PRINT, LWRITE)

CRAK  Begin second iteration.

      DO 410 I = 1, N
         P(I) = X(I) - Y(I)
         Y(I) = X(I)
  410 CONTINUE
      IF (V .LT. FNEW) GO TO 440
         SUM = 0.0E+0
         DO 420 I = 1, N
            SUM = SUM + P(I) ** 2
  420    CONTINUE
         PNORM = SQRT (SUM) + TINY
         U = HMAX/ PNORM
         DO 430 I = 1, N
            OLDG(I) = X(I) + U * P(I)
  430    CONTINUE
         CALL FUN (N, OLDG, OLDF)
         NFTOTL = NFTOTL + 1
         GTP = - ABS ((OLDF - FNEW)/ U)
         GO TO 450
  440 CONTINUE
         GTP = (V - FNEW) * 4.0E+0/ 3.0E+0
  450 CONTINUE
      ALPHA = 1.0E+0
      IF (V.GE.FNEW .AND. OLDF.LE.FNEW) GO TO 470
         DO 460 I = 1, N
            P(I) = -P(I)
  460    CONTINUE
  470 CONTINUE
      V = FNEW
      SFTBND = 0.0E+0
      T = RTEPS/ PNORM
      CALL LINSCH (N, FUN, RTEPS, T, EPSMCH, SFTBND, STEPMX/ PNORM,
     >   P, Y, V, ALPHA, GTP, NFTOTL, SUCCES, PRINT, LWRITE)
      NITER = NITER + 2
      COUNT = NITER .GE. NITMAX

CRAK  Invoke user-supplied routine for housekeeping.  The second call
CRAK  here (with IENTRY = 2) will usually be repeated below with more
CRAK  complete information, if requested.

      IENTRY = 1
      CALL USER (IENTRY, N, NLDIM, Y, V, G, H, EL, D, NFTOTL,
     >    NITER, NTYPE, CONTIN)
      IENTRY = 2
      CALL USER (IENTRY, N, NLDIM, Y, V, G, H, EL, D, NFTOTL,
     >    NITER, NTYPE, CONTIN)
      IF (PRINT) WRITE (LWRITE, 1040)
 1040 FORMAT (26H0QNMDIF:  END LOCAL SEARCH)

C     Terminate if a lower point was not found.  No new gradient will
C     be calculated, even if CONTIN = .TRUE., but previous values are
C     still current.  In special RESCUE = .TRUE. case, the gradient
C     will be either the input value or set to zero, depending on UNITL.

      IF (V .GE. FNEW) GO TO 500

C     If a lower point was found, the norm of the step to the lower
C     point is small, and the function decreased during the last
C     iteration, or iteration limit has been reached, and no final
C     gradient is requested, terminate.

      FNEW = V
      SUM = 0.0E+0
      DO 490 I = 1, N
         X(I) = Y(I)
         SUM = SUM + Y(I) ** 2
  490 CONTINUE
      DONE = (ALPHA * PNORM .LT. SQRT (EPSMCH * SUM) + TOL) .AND. CONV
      IF ((DONE .OR. COUNT) .AND. .NOT.CONTIN) GO TO 500

CRAK  Calculate gradient, call user interface to permit saving new
CRAK  information, and quit if through.  Note that stopping here leaves
CRAK  the Hessian as it was, i.e. no update is used.  This is a
CRAK  temporary expedient for this special case - the matrix update
CRAK  procedure will be modularized and called following APROXG in
CRAK  a future version.

      NTYPE = 1
      IF (RESCUE) NTYPE = 0
      CALL APROXG (N, NTYPE, X, G, H, PP, FUN, FNEW, NFTOTL,
     >   PRINT, LWRITE)
      IENTRY = 2
      CALL USER (IENTRY, N, NLDIM, X, V, G, H, EL, D, NFTOTL,
     >    NITER, NTYPE, CONTIN)
      IF (DONE .OR. COUNT) GO TO 500
      RESCUE = .FALSE.
      CONV = .TRUE.
      GO TO 50


C     Termination.
C     ------------

  500 CONTINUE
      F = FNEW

CRAK  This version always prints final results - restore the test if desired.
CRAK  IF (.NOT.PRINT) GO TO 510

C        Print final optimization results.

         WRITE (LWRITE, 1010)
         IF (COUNT) WRITE (LWRITE, 1050)
 1050    FORMAT (37H0MAXIMUM NUMBER OF ITERATIONS REACHED)
         IF (CONV) WRITE (LWRITE, 1060)
 1060    FORMAT (44H0THE TEST FOR CONVERGENCE HAS BEEN SATISFIED)
         IF (.NOT. CONV) WRITE (LWRITE, 1070)
 1070    FORMAT (50H0QNMDIF HAS FAILED TO SATISFY THE CONVERGENCE TEST)
         FINAL = .TRUE.
         CALL OUTPUT (N, NLDIM, X, G, H, P, FNEW, ALPHA, D, EL, NITER,
     >      NFTOTL, FINAL, LWRITE)
  510 CONTINUE

      RETURN
      END
*DECK APROXG
C+----------------------------------------------------------------------
C
      SUBROUTINE APROXG (N, NTYPE, X, G, H, WORK, FUN, FNEW, NFTOTL,
     >   PRINT, LWRITE)
C
C
C        This routine is called by subroutine QNMDIF.  It assigns a
C     finite-difference approximation to the gradient of the function
C     F (X) to the 1 * N array G (I). The intervals for differencing F
C     along each of the coordinate directions are given in array H (I).
C     The integer variable NTYPE determines which type of approximation
C     is given.
C
C        NTYPE=0:  Forward differences are obtained and the function
C           values are calculated at N forward points and stored in
C           the 1 * N array WORK.
C
C        NTYPE=1:  Central differences requiring the full 2N function
C           evaluations are computed.
C
C        NTYPE=2:  Central differences are evaluated using the forward
C           points previously stored in array WORK.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT REAL (A - H, O - Z)

C     Arguments.

      LOGICAL
     >   PRINT
      REAL
     >   G(N), H(N), WORK(N), X(N)
      EXTERNAL
     >   FUN

C     Execution.
C     ----------

      IF (PRINT .AND. NTYPE .EQ. 0) WRITE (LWRITE, 1000) NTYPE
 1000 FORMAT (51H0APROXG:  CALCULATE GRADIENT BY FORWARD DIFFERENCES,
     >   10H (NTYPE = , I1, 1H))
      IF (PRINT .AND. NTYPE .GT. 0) WRITE (LWRITE, 1010) NTYPE
 1010 FORMAT (51H0APROXG:  CALCULATE GRADIENT BY CENTRAL DIFFERENCES,
     >   10H (NTYPE = , I1, 1H))

      IF (NTYPE .EQ. 2) GO TO 20

C        Compute forward points.

         DO 10 I = 1, N
            XI = X(I)
            X(I) = XI + H(I)
            CALL FUN (N, X, FUNX)
            WORK(I) = FUNX
            X(I) = XI
   10    CONTINUE
         NFTOTL = NFTOTL + N
   20 CONTINUE
      IF (NTYPE .GT. 0) GO TO 40

C        Compute forward difference gradient.

         DO 30 I = 1, N
            G(I) = (WORK(I) - FNEW)/ H(I)
   30    CONTINUE
         GO TO 60
   40 CONTINUE

C        Compute backward points and central difference gradient.

         DO 50 I = 1, N
            XI = X(I)
            HI = H(I)
            X(I) = XI - HI
            CALL FUN (N, X, FUNX)
            G(I) = 0.5E+0 * (WORK(I) - FUNX)/ HI
            X(I) = XI
   50    CONTINUE
         NFTOTL = NFTOTL + N
   60 CONTINUE


C     Termination.
C     ------------

      RETURN
      END
*DECK C1MDCH
C+----------------------------------------------------------------------
C
      SUBROUTINE C1MDCH (N, NLDIM, D, EL, ALPHA, Z, IFAIL)
C
C
C        The subroutine C1MDCH updates the Cholesky factorization of the
C     matrix EL * D * ELT  + ALPHA * Z * ZT, where EL is a unit lower
C     triangular matrix, D is a diagonal matrix with positive elements,
C     ELT is the transpose of EL, ALPHA is a positive scalar, Z is a
C     vector, and ZT is the transpose of Z, so that ALPHA * Z * ZT is
C     a positive rank-one change.
C
C        The algorithm used is described in Gill, Golub, Murray, and
C     Saunders, "Methods for Modifying Matrix factorizations", Stanford
C     Computer Science Report CS-72-322, Stanford University, Dec. 1972.
C
C        N is the dimension of the matrices and vectors, NLDIM is
C     N * (N-1)/ 2, the length of the vector EL in which the strict
C     lower triangle of the matrix EL is stored by rows, D is a
C     vector of length N containing the elements of the matrix D,
C     EL is a vector of length NLDIM containing the elements of EL,
C     ALPHA is the positive scalar multiple of the rank-one change,
C     Z is the vector involved, and IFAIL is an integer flag set
C     to zero if there are no errors, and to 1 if overflow occurred.
C
C        Both EL and D are overwritten with the modified factors of the
C     altered matrix, and the values of ALPHA and Z are not retained.
C
C        The subroutine sets IFAIL to zero if there are no problems.
C     IFAIL is set to 1 if any element of the diagonal formed is zero
C     or negative, or when the ratio of the new diagonal element to
C     the old is too large and would cause overflow.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT REAL (A - H, O - Z)

C     Constants.

C     RMAX should be set to the largest positive floating point value such
C     that +RMAX and -RMAX can both be represented in the computer.  For
C     convenience, the smallest useful value (VAX) has been used - it is
C     unlikely to matter, but could be increased on another system.

      REAL
     >   RMAX
      PARAMETER
     >   (RMAX = 1.E+38)

C     Arguments.

      REAL
     >   D(N), EL(NLDIM), Z(N)


C     Execution.
C     ----------

      A = ALPHA
      K = 0
      IFAIL = 1
      DO 70 I = 1, N
         P1 = Z(I)
         DI = D(I)
         T = A * P1
         D(I) = DI + T * P1
         DB = D(I)

C        Exit with IFAIL = 1 if any of the new diagonal elements is
C        non-positive, or if the ratio of diagonals will overflow.

         IF (DB .GE. 1.0E+0) GO TO 10
            IF (DB .LT. 0.0E+0 .OR. DI .GT. (RMAX * DB))  RETURN
   10    CONTINUE
         GAMMA = DI/ DB
         BETA = T/ DB
         K = K + I
         J = K
         A = A * GAMMA
         IF (I .EQ. N) GO TO 60
            IP1 = I + 1
            IF (GAMMA .GT. 0.25E+0) GO TO 30

C              If ratio of diagonals is less than 4.0, proceed with
C              normal update.

               DO 20 IB = IP1, N
                  T = EL(J)
                  EL(J) = T * GAMMA + BETA * Z(IB)
                  Z(IB) = Z(IB) - P1 * T
                  J = J + IB - 1
   20          CONTINUE
               GO TO 50
   30       CONTINUE

C              Use alternative update if ratio of diagonals exceeds 4.0.

               DO 40 IB = IP1, N
                  T = EL(J)
                  Z(IB) = Z(IB) - P1 * T
                  W = Z(IB)
                  EL(J) = BETA * W + T
                  J = J + IB - 1
   40          CONTINUE
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
      IFAIL = 0


C     Termination.
C     ------------

      RETURN
      END
*DECK LINSCH
C+----------------------------------------------------------------------
C
      SUBROUTINE LINSCH (N, FUN, EPS, T, ETA, SFTBND, STEPMX,
     >   P, X, F, ALPHA, GTP, NFTOTL, SUCCES, PRINT, LWRITE)
C
C
C        Called by QNMDIF to execute a linear search along a given
C     direction P, starting from a given point X, to locate an
C     approximation ALPHA to the point at which the objective function
C     attains its minimum value along the direction P.  The method used
C     is that of successive quadratic interpolation with safeguards.
C
C
C     History:
C
C        12 Aug. 1986    RAK    Print and test T on entry to monitor
C                               conflicts with SFTBND and STEPMX.  Check
C                               D and GTP for before main loop.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT REAL (A - H, O - Z)

C     Constants.

      INTEGER
     >   MAXN
      PARAMETER
     >   (MAXN = 40)

C     Arguments.

      LOGICAL
     >   PRINT, SUCCES
      REAL
     >   P(N), X(N)

C     Local variables.

      REAL
     >   Z(MAXN)

C     Procedures.

      EXTERNAL
     >   FUN

C     Execution.
C     ----------

      IF (PRINT) WRITE (LWRITE, 1000) SFTBND, STEPMX, T
 1000 FORMAT (36H0LINSCH:  BEGIN SEARCH, MIN. STEP = , E15.8/
     >   24X, 12HMAX. STEP = , E15.8/
     >   24X, 12HTOLERANCE = , E15.8)

      IF (T .GT. STEPMX) THEN
         WRITE (LWRITE, 1010)
 1010    FORMAT (51H0LINSCH:  MAX. STEP MUST BE GREATER THAN TOLERANCE!)
         SUCCES = .FALSE.
         GO TO 330
      END IF
         

C     Initialization.
C     ---------------

C     NLIN is the number of function evaluations during the linear
C     search.  It is used locally as a flag to compute the parabolic
C     step the first time using GTP, the directional derivative along
C     the search vector.

CRAK  MAXLIN is a (fairly loose) upper limit on the number of function
CRAK  evaluations permitted before declaring the search a failure.

      NLIN = 0
      MAXLIN = 20
      V = 0.0E+0
      W = 0.0E+0
      XMIN = 0.0E+0
      FA = F
      FV = F
      FW = F
      FMIN = F
      FOLD = F
      D = ALPHA
      TOL = T
      T2 = 2.0E+0 * TOL

C     All points in the linear search are scaled (shifted) so that
C     the currently lowest point (XMIN) is the origin.  A and B define
C     the interval of uncertainty, W is the last value of XMIN, and V
C     is the highest of the three points through which a parabola may
C     be fitted.

      A = 0.0E+0
      B = STEPMX + EPS * ABS (STEPMX) + T
      E = B
      SCXBD = STEPMX

C     NPCNT is a local count of the number of non-parabolic interpolation
C     steps used in LINSCH.  It may be useful for diagnostic purposes if
C     there are problems locating the minimum of a particular function.
C     GTEST1 and GTEST2 are used to test convergence.

      NPCNT = 0
      B1 = B
      A1 = A
      GTEST1 = - 1.0E-4 * GTP
      GTEST2 = - ETA * GTP
      ALPHA = 0.0E+0

C     D is the estimate of the next point at which the function is to be
C     evaluated, with respect to the current origin.  Make sure that it,
C     and GTP, are sensible before proceding.

      IF (D .LE. 0.0E+0 .OR. GTP .GE. 0.0E+0) THEN
         WRITE (LWRITE, 1020) D, GTP
 1020    FORMAT (32H0LINSCH:  BAD INITIAL DATA, D = , E25.15,
     >      8H, GTP = , E25.15)
         SUCCES = .FALSE.
         GO TO 330
      END IF


C     Begin main loop.
C     ----------------

   10 CONTINUE
      IF (D .LT. SCXBD) GO TO 20

C        D exceeds the shifted bound on the minimum, so adjust D and
C        the bound.

         D = SCXBD
         SCXBD = (SCXBD - TOL)/ (1.0E+0 + EPS)
   20 CONTINUE

C     U is the point at which the function will actually be evaluated,
C     scaled with respect to the shifted origin.  Thus XMIN, the actual
C     step from the original starting point, must be added to obtain
C     the true step.  The estimate D will be used as the value for U
C     only if it is sufficiently far from zero (the current minimum),
C     since the function is not allowed to be evaluated at points closer
C     together than TOL.

      U = D
      IF (ABS (D) .LT. TOL) U = SIGN (TOL, D)
      R = XMIN + U
      DO 30 I = 1, N
         Z(I) = X(I) + R * P(I)
   30 CONTINUE


C     Evaluate function at new point.
C     -------------------------------

      IF (PRINT) WRITE (LWRITE, 1030) R
 1030 FORMAT (17H0LINSCH:  STEP = , E25.15)
      CALL FUN (N, Z, FU)
      NLIN = NLIN + 1

C     Update A, B, V, W, and XMIN.  Check whether new function value is
C     lower than previous minimum.

      IF (FU .GT. FMIN) GO TO 60

C     The following code treats the case where FU is the new lowest
C     point, so that the new point, U, becomes the next origin, and the
C     other points are shifted accordingly.

C     Shift left or right endpoint depending on which side of origin
C     the new point is.

      IF (U .LT. 0.0E+0) GO TO 40
         A = 0.0E+0
         FA = FMIN
         GO TO 50
   40 CONTINUE
         B = 0.0E+0
   50 CONTINUE


C     Shift all points with respect to new minimum.
C     ---------------------------------------------

      V = W
      FV = FW
      W = 0.0E+0
      FW = FMIN
      FMIN = FU
      XMIN = XMIN + U
      A = A - U
      B = B - U
      V = V - U
      W = W - U
      SCXBD = SCXBD - U
      TOL = EPS * ABS (XMIN) + T
      T2 = 2.0E+0 * TOL
      GO TO 110

C     In the following code, the new function value was greater than
C     the previously lowest.  Check the relationship of the new point
C     to the other values (next lowest, etc.).  The origin remains
C     unchanged, but other points may be interchanged.

C     Shift either the left or right endpoint of the interval of
C     uncertainty, depending on whether the new point, U, was to the
C     left or right of the origin.

   60 IF (U .GT. 0.0E+0) GO TO 70
         A = U
         FA = FU
         GO TO 80
   70 CONTINUE
         B = U
   80 CONTINUE
      IF (FU .GT. FW .AND. W .NE. 0.0E+0) GO TO 90

C     If FU is less than or equal previous second best point, or W = 0,
C     i.e., this is the first time through this section of code.

      V = W
      FV = FW
      W = U
      FW = FU
      GO TO 100
   90 IF (FU .GT. FV .AND. V .NE. 0.0E+0 .AND. V .NE. W) GO TO 100

C     FU .LE. FW, or V = 0, or V = W (first or second time through).

      V = U
      FV = FU
  100 U = 0.0E+0

C     Compute midpoint of interval of uncertainty.

  110 XM = 0.5E+0 * (A + B)


C     Check stopping criteria.
C     ------------------------

C     Stop if interval of uncertainty is sufficiently small, or
C     if the best point is less than the required bound, or
C     if function value has decreased sufficiently, or
C     if MAXLIN function evaluations have already been performed.

      IF ((ABS (XM) .LE. (T2 - 0.5E+0 * (B - A))) .OR.
     >   ((XMIN + B) .LE. SFTBND) .OR.
     >   ((FA - FMIN).LE.(ABS (A) * GTEST2) .AND. FMIN.LT.FOLD) .OR.
     >   (NLIN .GE. MAXLIN)) GO TO 280

C     If stopping criteria are not met, continue.

      S = 0.0E+0
      Q = 0.0E+0
      R = 0.0E+0
      IF (ABS (E) .LE. TOL) GO TO 170

C     Otherwise, fit parabola through XMIN, V, and W.

      IF (NLIN .GT. 1) GO TO 130

C     If NLIN = 1, this is the first parabolic fit and the (known)
C     approximate gradient at the initial point can be used.

      Q = 2.0E+0 * (FW - FMIN - W * GTP)
      IF (XMIN .EQ. 0.0E+0) GO TO 120
      S = (2.0E+0 * (FMIN - FW) + W * GTP) * W
      GO TO 140
  120 S = GTP * W * W
      GO TO 140

C     NLIN is greater than 1, so the fit uses function values only.

  130 R = W * (FV - FMIN)
      Q = V * (FW - FMIN)
      S = R * W - Q * V
      Q = 2.0E+0 * (Q - R)
  140 IF (Q .LE. 0.0E+0) GO TO 150
      S = - S
      GO TO 160
  150 Q = - Q
  160 R = E
      IF (D .NE. B1 .OR. B .LE. SCXBD)  E = D

C     Construct an artificial bound on the estimated step length.

  170 A1 = A
      B1 = B
      IF (XMIN .NE. 0.0E+0) GO TO 180
      D = XM
      GO TO 230
  180 IF (B .LE. SCXBD) GO TO 190

C     Expand interval by 4 if minimum is still not bracketed.

      D = - 4.0E+0 * A
      GO TO 230

C     B is less than or equal to SCXBD.

  190 D1 = A

C     Determine interval of length D2 in which to set artificial bound.

      D2 = B
      IF (ABS (D2) .GT. TOL .AND.
     >  ((W .LE. 0.0E+0) .OR. (ABS (D1) .LE. TOL))) GO TO 200

C     ABS (D2) is less than or equal to TOL, or W is positive and
C     ABS (D1) is greater than TOL.  In either case, interchange D1, D2.

      U = D1
      D1 = D2
      D2 = U

C     Use parabolic interpolation only if new point lies in (A1,B1).

  200 U = - D1/ D2
      IF (U .LT. 1.0E+0) GO TO 210

C     Extrapolation - U exceeds 1.

      FACT = 5.0E+0 * (0.1E+0 + 1.0E+0/ U)/ 11.0E+0
      GO TO 220

C     Interpolation step - U less than 1.

  210 FACT = 0.5E+0 * SQRT (U)
  220 D = D2 * FACT

C     If D > 0 then B1 = D else A1 = D.

  230 IF (D .LE. 0.0E+0) GO TO 240
      B1 = D
      GO TO 250
  240 A1 = D
  250 IF (ABS (S) .GE. ABS (0.5E+0 * Q * R) .OR. (S .LE. (Q * A1))
     >  .OR. (S .GE. (Q * B1))) GO TO 260

C     A parabolic interpolation step.

      D = S/ Q

C     F must not be evaluated too close to A or B.

      IF ((D - A) .GE. T2 .AND. (B - D) .GE. T2) GO TO 10

C     Otherwise, set D to plus or minus TOL.

      D = SIGN (TOL, XM)
      GO TO 10

C     A non-interpolation step.

  260 NPCNT = NPCNT + 1
      IF (XM .LE. 0.0E+0) GO TO 270
      E = B
      GO TO 10
  270 E = A
      GO TO 10


C     Check safeguards.
C     -----------------

  280 CONTINUE
      SUCCES = .FALSE.

C     Check that new point satisfies safeguard conditions.  If the
C     function did not decrease, or step to the minimum was less than
C     SFTBND, LINSCH has failed to locate an acceptable minimum.

  290 CONTINUE
         IF (XMIN + B .LE. SFTBND) GO TO 330
         IF (FOLD - FMIN .GT. GTEST1 * XMIN) GO TO 310
         IF (NLIN .GE. MAXLIN) GO TO 330

C        A sufficiently lower point was not found - try halving step.

         XMIN = XMIN * 0.5E+0
         IF (XMIN .LE. T) GO TO 330
         DO 300 I = 1, N
            Z(I) = X(I) + XMIN * P(I)
  300    CONTINUE
         IF (PRINT) WRITE (LWRITE, 1030) XMIN
         CALL FUN (N, Z, FMIN)
         NLIN = NLIN + 1
         GO TO 290

C     A sufficiently lower point was found - set output values.

  310 CONTINUE
      SUCCES = .TRUE.
      ALPHA = XMIN
      IF (SCXBD .LE. 0.0E+0) ALPHA = STEPMX
      DO 320 I = 1, N
         X(I) = X(I) + ALPHA * P(I)
  320 CONTINUE
      F = FMIN


C     Termination.
C     ------------

  330 CONTINUE
      NFTOTL = NFTOTL + NLIN

      RETURN
      END
*DECK OUTPUT
C+----------------------------------------------------------------------
C
      SUBROUTINE OUTPUT (N, NLDIM, X, G, H, P, FNEW, ALPHA, D,
     >   EL, NITER, NFTOTL, FINAL, LWRITE)
C
C
C        This routine prints details of an iteration by QNMDIF.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT REAL (A - H, O - Z)

C     Arguments.

      LOGICAL
     >   FINAL
      REAL
     >   D(N), EL(NLDIM), G(N), H(N), P(N), X(N)


C     Execution.
C     ----------

      IF (FINAL) GO TO 30
         WRITE (LWRITE, 1000) NITER, FNEW
 1000    FORMAT (13H0AT ITERATION, I5, 22H THE FUNCTION VALUE IS,
     >      E25.15)
         WRITE (LWRITE, 1010)
 1010    FORMAT (44H0          CURRENT SOLUTION         GRADIENT,
     >      50H                 SEARCH DIRECTION         STEPSIZE)
         DO 10 I = 1, N
            WRITE (LWRITE, 1020) I, X(I), G(I), P(I), H(I)
 1020       FORMAT (1H , I5, 4E25.15)
   10    CONTINUE
         WRITE (LWRITE, 1030) ALPHA, NFTOTL
 1030    FORMAT (19H0LINEAR SEARCH STEP, E25.15, 5X,
     >      50HNUMBER OF FUNCTION EVALUATIONS AT END OF ITERATION, I9)
         WRITE (LWRITE, 1040)
 1040    FORMAT (40H0CHOLESKY FACTORS OF APPROXIMATE HESSIAN/
     >      31H0   ELEMENTS OF DIAGONAL MATRIX)
         WRITE (LWRITE, 1050) (D(I), I = 1, N)
 1050    FORMAT (1H , 5E25.15)
         IF (N .EQ. 1) GO TO 60
         WRITE (LWRITE, 1060)
 1060    FORMAT (27H0   LOWER TRIANGULAR FACTOR)
         K2 = 0
         NM1 = N - 1
         DO 20 I = 1, NM1
            K1 = K2 + 1
            K2 = K2 + I
            WRITE (LWRITE, 1050) (EL(J), J = K1, K2)
   20    CONTINUE
         GO TO 60
   30 CONTINUE

CRAK     Final iteration gets special treatment.

         WRITE (LWRITE, 1070) NITER, FNEW
 1070    FORMAT (/16H0FINAL ITERATION, I6, 5X, 14HFUNCTION VALUE,
     >      E25.15)
         WRITE (LWRITE, 1080)
 1080    FORMAT (44H0          FINAL SOLUTION           GRADIENT)
         DO 40 I = 1, N
            WRITE (LWRITE, 1020) I, X(I), G(I)
   40    CONTINUE
         WRITE (LWRITE, 1090) NFTOTL
 1090    FORMAT (31H0NUMBER OF FUNCTION EVALUATIONS, I9)
         WRITE (LWRITE, 1040)
         WRITE (LWRITE, 1050) (D(I), I = 1, N)
         IF (N .EQ. 1) GO TO 60
         WRITE (LWRITE, 1060)
         K2 = 0
         NM1 = N - 1
         DO 50 I = 1, NM1
            K1 = K2 + 1
            K2 = K2 + I
            WRITE (LWRITE, 1050) (EL(J), J = K1, K2)
   50    CONTINUE
   60 CONTINUE

CRAK  Compute and print estimate of condition number of Hessian and
CRAK  norm of gradient (all iterations).

      DMAX = D(1)
      DMIN = D(1)
      SUM = 0.E+0
      DO 70 I = 1, N
         IF (D(I) .GT. DMAX) DMAX = D(I)
         IF (D(I) .LT. DMIN) DMIN = D(I)
         SUM = SUM + G(I) ** 2
   70 CONTINUE
      BOUNDK = DMAX/ DMIN
      GNORM = SQRT (SUM)
      WRITE (LWRITE, 1100) BOUNDK, GNORM
 1100 FORMAT (43H0LOWER BOUND ON CONDITION NUMBER OF HESSIAN, E25.15/
     >   17H0NORM OF GRADIENT, E25.15)


C     Termination.
C     ------------

      RETURN
      END
*DECK UPCHOL
C+----------------------------------------------------------------------
C
      SUBROUTINE UPCHOL (N, NLDIM, EPSMCH, D, EL, Z)
C
C
C        The subroutine UPCHOL forms the updated Cholesky factorization
C     of the matrix EL * D * ELT - Z * ZT, where EL is a unit lower
C     triangular matrix, D is diagonal matrix with positive elements,
C     ELT is the transpose of EL, Z is a vector, and ZT is its transpose
C     (so that -Z * ZT is a negative rank-one correction).
C
C        The algorithm used is described in Gill, Golub, Murray, and
C     Saunders, "Methods for Modifying Matrix Factorizations", Stanford
C     Computer Science Report CS-72-322, Stanford University, Dec. 1972.
C
C        N is the dimension of the matrices and vectors, NLDIM is
C     N * (N-1)/ 2, the length of the vector EL in which the strict
C     lower triangle of the matrix EL is stored by rows, EPSMCH is
C     "machine epsilon", used to ensure that the resulting matrix
C     is sufficiently positive definite, D is a vector of length N
C     containing the elements of the matrix D, EL is a vector of length
C     NLDIM containing the elements of EL, and Z is the vector to be
C     used in the update.
C
C       The vectors EL and D are overwritten with the updated Cholesky
C    factors, and the elements of Z are overwritten.  UPCHOL modifies
C    the Cholesky factors of a matrix that could be indefinite, but it
C    ensures that the new matrix is positive definite.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT REAL (A - H, O - Z)

C     Constants.

      INTEGER
     >   MAXN
      PARAMETER
     >   (MAXN = 40)

C     Arguments.

      REAL
     >   D(N), EL(NLDIM), Z(N)

C     Local variables.

      REAL
     >   P(MAXN)


C     Execution.
C     ----------

C     Solve EL * P = Z.

      J = 1
      P(1) = Z(1)
      GAMMA = P(1) ** 2/ D(1)
      IF (N .EQ. 1) GO TO 30
         DO 20 I = 2, N
            K = I - 1
            T = Z(I)
            DO 10 IB = 1, K
               T = T - P(IB) * EL(J)
               J = J + 1
   10       CONTINUE
            P(I) = T
            GAMMA = GAMMA + T * T/ D(I)
   20    CONTINUE
   30 CONTINUE

C     If 1.0 - GAMMA < EPSMCH, then the modified matrix is not sufficiently
C     positive definite.  GAMMA is replaced by a quantity which ensures
C     that the modified matrix is positive definite regardless of subsequent 
C     rounding error.

      GAMMA = 1.0E+0 - GAMMA
      IF (GAMMA .GT. EPSMCH) GO TO 60
         IF (- GAMMA .GT. EPSMCH) GO TO 40
            GAMMA = EPSMCH
            GO TO 50
   40    CONTINUE
            GAMMA = - GAMMA
   50    CONTINUE
   60 CONTINUE
      K = J + N + N

C     Solve D * ELTRANSPOSE * Z = P.

      DO 90 JJ = 1, N
         J = N + 1 - JJ
         PJ = P(J)
         DJ = D(J)
         T = PJ/ DJ
         Z(J) = PJ
         BETA = - T/ GAMMA
         G = GAMMA + PJ * T
         D(J) = DJ * GAMMA/ G
         GAMMA = G
         K = K - J - 1
         IQ = K
         IF (J .EQ. N) GO TO 80
            JP1 = J + 1
            DO 70 IB = JP1, N
               T = EL(IQ)
               EL(IQ) = T + BETA * Z(IB)
               Z(IB) = Z(IB) + PJ * T
               IQ = IQ + IB - 1
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE


C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE QPLDAT
C
C PURPOSE: QPLDAT writes QPLOTable data to the indicated file, including
C          titles, axis labels and most of the optional namelist inputs
C          accepted by QPLOT.
C
C METHOD:  Separate entry points are used to distinguish between new
C          frame information (plot titles, etc. - entry QPLFRM), new
C          curve information (namelist options  - entry QPLCRV), and
C          strictly numerical values (possibly points to be inserted
C          in an existing curve - entry QPLNXY).
C
C          Floating point arguments are checked against a flag (also
C          an argument) to suppress any unwanted variables in the
C          namelist.  Blanks should be passed for any of the character
C          variables not wanted in the namelist.
C
C ARGUMENTS FOR ENTRY QPLFRM:
C
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C
C  LUN          -        I    I   Logical unit number for QPLOT file.
C  TITLE        -      C*(*)  I   Title for a plot frame; QPLOT has a
C                                 limit of 80 characters.
C  SBTITL       -      C*(*)  I   Subtitle for plot frame.
C  XLABEL,      -      C*(*)  I   Axis labels (also up to 80 chars.).
C  YLABEL
C
C ARGUMENTS FOR ENTRY QPLCRV:
C
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C
C  LUN          -        I    I   Logical unit number for QPLOT file.
C  NPTS         -        I    I   Length of data arrays. NPTS = 0 permits
C                                 just the namelist options to be written.
C  X,Y        NPTS       R    I   Data representing one curve.
C  ENDSTR       -        C    I   Either 'END CURVE', 'END FRAME', or ' '.
C  OMIT         -        R    I   Flag for omitting floating point
C                                 variables from the namelist. Input
C                                 arguments equal to OMIT will not be
C                                 written to the output file.
C
C  Note: The remaining arguments are variables for QPLOT's optional namelist.
C
C  PLOT         -      C*(*)  I   Plot type. Options include 'DEFAULT' =
C                                 'LINEAR', 'SCALE', 'LOGX', 'LOGY',
C                                 'LOGLOG', and 'POLAR'.
C  YMIN,YMAX,   -        R    I   Bounds used for windowing or suppress-
C  XMIN,XMAX                      ing points out of desired range.
C  WIDTH        -        R    I   Desired length of X axis in inches.
C  HEIGHT       -        R    I   Desired length of Y axis in inches.
C  LEGEND       -      C*(*)  I   Text for legend entry for this curve.
C  LINE         -      C*(*)  I   Line type for this curve. Use 'DEFAULT'
C                                 or 'CONNECT' for symbols and lines. Other
C                                 choices are 'SYMBOLS', 'SOLID', 'DASH',
C                                 'DOT', 'THICK', 'CHAINDASH', 'CHAINDOT',
C                                 OR 'LONGDASH'.
C
C ARGUMENTS FOR ENTRY QPLNXY:
C
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C
C  LUN          -        I    I   Logical unit number for QPLOT file
C  NPTS         -        I    I   Length of data arrays. Must be > 0.
C  X,Y        NPTS       R    I   Data points to be included.
C  ENDSTR       -        C    I   Either 'END CURVE', 'END FRAME', or ' '.
C
C ERROR HANDLING:  None.
C
C STANDARDS VIOLATIONS:  Multiple entry points make sense here.
C
C ENVIRONMENT:  Fortran 90
C
C HISTORY:
C
C  10/30/83  DAS  Initial implementation (point suppression done here).
C  01/31/84  LJC  Added writing of optional namelist, leaving QPLOT to
C                 suppress points out of desired range.
C  04/12/84  DAS  Made use of QPLOT's legend capability.
C  10/23/84  LJC  Added more namelist options and QPLNXY entry point.
C  06/20/86  DAS  ENDSTR*(*) instead of ENDSTR*9.
C  12/10/86  DAS  X(*), Y(*) instead of X(NPTS), ... for the NPTS=0
C                 case of QPLNXY.
C  06/18/91  DAS  Kludge for missing ORIENT option: if WIDTH is too
C                 big for portrait but small enough for landscape
C                 (8.5 x 11 in both cases), insert ORIENT='LANDSCAPE'.
C  10/18/99  DAS  SGI list-directed I/O started putting out things like
C                 2*0.E+0, which cannot be parsed by QPLOT.
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C ... Arguments:

      INTEGER
     >   LUN, NPTS

      REAL
     >   HEIGHT, OMIT, WIDTH, X (*), XMAX, XMIN, Y (*), YAXIS, YMAX,
     >   YMIN

      CHARACTER
     >   ENDSTR * (*), LEGEND * (*), LINE * (*), PLOT * (*),
     >   SBTITL * (*), TITLE * (*), XLABEL * (*), YLABEL * (*)

C ... Local constants:

      CHARACTER * 1, PARAMETER ::
     >   BLANK = ' '

C ... Local variables:

      INTEGER
     >   I

      REAL
     >   EPS

C-----------------------------------------------------------------------

      ENTRY QPLFRM (LUN, TITLE, SBTITL, XLABEL, YLABEL)

C-----------------------------------------------------------------------

C ... Write the title and axis labels for a plot frame.

      WRITE (LUN, 100) TITLE
      WRITE (LUN, 100) SBTITL
      WRITE (LUN, 100) XLABEL
      WRITE (LUN, 100) YLABEL
      GO TO 99


C-----------------------------------------------------------------------

      ENTRY QPLCRV (LUN, NPTS, X, Y, ENDSTR, OMIT, PLOT, XMIN, XMAX,
     >              YMIN, YMAX, WIDTH, HEIGHT, LEGEND, LINE)

C-----------------------------------------------------------------------

C ... Write options in namelist format:

      WRITE (LUN, 100) ' $OPTIONS'

      IF (PLOT  /= BLANK) WRITE (LUN, 120) 'PLOT', PLOT
      IF (XMIN  /=  OMIT) WRITE (LUN, 110) 'XMIN', XMIN
      IF (XMAX  /=  OMIT) WRITE (LUN, 110) 'XMAX', XMAX
      IF (YMIN  /=  OMIT) WRITE (LUN, 110) 'YMIN', YMIN
      IF (YMAX  /=  OMIT) WRITE (LUN, 110) 'YMAX', YMAX
      IF (WIDTH /=  OMIT) THEN
         WRITE (LUN, 110) 'WIDTH',  WIDTH
         IF (WIDTH > 8.0 .AND. WIDTH < 11.0) THEN
            WRITE (LUN, 120) 'ORIENT', 'LANDSCAPE'
         END IF
      END IF
      IF (HEIGHT/=  OMIT) WRITE (LUN, 110) 'HEIGHT', HEIGHT
      IF (LEGEND/= BLANK) WRITE (LUN, 120) 'LEGEND', LEGEND
      IF (LINE  /= BLANK) WRITE (LUN, 120) 'LINE',   LINE

      WRITE (LUN, 100) ' $END'
 
C-----------------------------------------------------------------------

      ENTRY QPLNXY (LUN, NPTS, X, Y, ENDSTR)

C-----------------------------------------------------------------------

C ... Save data points:

      IF (NPTS > 0) THEN

         EPS = EPSILON (EPS) ! Avoid 2*0.E+0, etc., from SGI list-directed write

         IF (EPS > 1.E-10) THEN ! Assume 32-bit compile

            WRITE (LUN, '(1P, 2E15.7)') (X (I), Y (I), I = 1, NPTS)

         ELSE

            WRITE (LUN, '(1P, 2E23.15)') (X (I), Y (I), I = 1, NPTS)

         END IF

      END IF

      IF (ENDSTR /= BLANK) WRITE (LUN, 100) ENDSTR

   99 RETURN

C ... Formats:

  100 FORMAT (A)
  110 FORMAT (1X, A, ' = ', 1P, E13.6, ',')
  120 FORMAT (1X, A, ' = ''', A, ''',')

      END SUBROUTINE QPLDAT
C+------------------------------------------------------------------------------
C
      SUBROUTINE QUANC8RC (XEVAL, FEVAL, A, B, ABSERR, RELERR, RESULT,
     >                     ERREST, NOFUN, FLAG, ISTAT)
C
C  ACRONYM: QUadrature/Adaptive/Newton-Cotes/8-panel/Reverse-Communication
C           --         -        -      -     -       -       -
C
C  PURPOSE:
C
C        QUANC8RC is a reverse-communication version of QUANC8 for estimating
C     the integral of a function of one variable on the interval [A, B] to a
C     specified tolerance using an adaptive scheme based upon the 8-panel
C     Newton-Cotes rule.
C
C        This version avoids the data-communication problems encountered with
C     functions which require more than a single argument X to be sufficiently
C     defined.  (QUANC8 forces the additional information to be passed via
C     a COMMON block, which is usually inconvenient.)
C
C
C  USAGE:
C
C        :   :
C        X = A       ! Initialize the integration
C        ISTAT = 0
C
C    100 CONTINUE
C
C           CALL FUN (X, F, ..., IER)      ! Or whatever it takes to evaluate
C                                          ! the function as F at X
C           IF (IER .NE. 0) THEN
C              <Handle an error in the function evaluation>
C           END IF
C
C           CALL QUANC8RC (X, F, A, B, ABSERR, RELERR, RESULT,
C       >                  ERREST, NOFUN, FLAG, ISTAT)
C
C           IF (ISTAT .GT. 0)
C       >GO TO 100
C
C        IF (ISTAT .NE. 0) THEN   ! Trouble
C           <Display a warning with FLAG = XXX.YYY
C            and handle the failure>
C        END IF
C        :    :
C
C
C  INPUT:
C
C     XEVAL   is not really an input, but it needs to be initialized to A
C             as shown above for all evaluations to be done in the same loop.
C     FEVAL   is input with the value of the function corresponding to XEVAL.
C     A       is the lower limit of integration.
C     B       is the upper limit of integration.  (B may be less than A.)
C     ABSERR  is an absolute error tolerance which should be non-negative.
C     RELERR  is a relative error tolerance which should be non-negative.
C     ISTAT   is input as 0 on the FIRST call and SHOULD NOT BE CHANGED by
C             the calling program before the integration is complete, since
C             a positive value output by the previous call is meaningful as
C             input on the present call if it is positive.
C
C  OUTPUT:
C
C     ISTAT   is output as 0 upon successful completion;
C             ISTAT < 0 indicates a failure - see FLAG;
C             ISTAT > 0 is convenient as an index into this routine's working
C             arrays of abscissas and ordinates: the index output is used
C             as input on the next call to indicate where to store FEVAL.
C     RESULT  is an approximation to the integral which should satisfy the
C             least stringent of the two error tolerances.
C     ERREST  is an estimate of the magnitude of the actual error.
C     NOFUN   is the number of function values used in calculation of RESULT.
C     FLAG    is a reliability indicator.  If FLAG is 0., then RESULT
C             probably satisfies the error tolerance.  If FLAG is XXX.YYY,
C             then XXX = the number of intervals which have not converged
C             and 0.YYY = the fraction of the interval left to do when the
C             (internal) limit on NOFUN was approached.
C
C
C  REFERENCE:
C
C     Forsythe, G., M. Malcolm, and C. Moler.  Computer Methods for
C     Mathematical Computations (Englewood Cliffs: Prentice-Hall, 1977).
C
C
C  HISTORY:
C
C     29 Apr. 1984   R.A.Kennelly,   Partial conversion of QUANC8 to
C                    NASA Ames       FORTRAN 77 using generic functions.
C     30 Apr. 1986   R.A.K.          Added IMPLICIT NONE.
C     21 May  1992   D.A.Saunders,   Adapted QUANC8 as QUANC8RC to avoid
C                    Sterling/       COMMON blocks.  Subscript 0 for X (*)
C                    NASA Ames       and F (*) simplifies the code.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   ISTAT, NOFUN
      REAL
     >   A, B, ABSERR, ERREST, FEVAL, FLAG, RELERR, RESULT, XEVAL

C     Local constants (not all of them):

      REAL
     >   HALF, ONE, ZERO
      PARAMETER
     >  (HALF = 0.5E+0, ONE = 1.E+0, ZERO = 0.E+0)

C     Local variables:

      INTEGER
     >   LEVMIN, LEVMAX, LEVOUT, NOMAX, NOFIN, LEV, NIM, I, J
      REAL
     >   W0, W1, W2, W3, W4, AREA, STONE, STEP, COR11, TEMP, QPREV,
     >   QNOW, QDIFF, QLEFT, ESTERR, TOLERR, QRIGHT (31),
     >   F (0 : 16), X (0 : 16), FSAVE (8, 30), XSAVE (8, 30)

C     Local storage:

      SAVE

C     Execution:


C     ***   STAGE 1 ***   GENERAL INITIALIZATION

      IF (ISTAT .EQ. 0) THEN

         LEVMIN = 1
         LEVMAX = 30
         LEVOUT = 6
         NOFUN  = 0
         NOMAX  = 5000
         NOFIN  = NOMAX - 8 * (LEVMAX - LEVOUT + 2 ** (LEVOUT + 1)) !NOFUN limit

         W0 =   3956.0E+0 / 14175.0E+0
         W1 =  23552.0E+0 / 14175.0E+0
         W2 =  -3712.0E+0 / 14175.0E+0
         W3 =  41984.0E+0 / 14175.0E+0
         W4 = -18160.0E+0 / 14175.0E+0

         FLAG   = ZERO
         RESULT = ZERO
         COR11  = ZERO
         ERREST = ZERO
         AREA   = ZERO

         IF (A .EQ. B) GO TO 99


C        ***   STAGE 2 ***   INITIALIZATION FOR FIRST INTERVAL

         LEV = 0
         NIM = 1
         STONE = (B - A) / 16.0E+0

         DO 20, J = 0, 16, 2
            X (J) = A + STONE * J
   20    CONTINUE
         X (16) = B     ! To avoid round-off

         QPREV  = ZERO
      END IF

      IF (MOD (ISTAT, 2) .EQ. 0) THEN   ! Fill in an even-numbered F (*)

         NOFUN = NOFUN + 1
         F (ISTAT) = FEVAL

         IF (ISTAT .LT. 16) THEN
            ISTAT = ISTAT + 2
            XEVAL = X (ISTAT)
            GO TO 99                    ! Go back for another function value
         ELSE
            ISTAT = -1                  ! Initialize the odd function evals.
         END IF

      END IF


C     ***   STAGE 3 ***   CENTRAL CALCULATION
C     Requires QPREV, X0, X2, X4, ..., X16, F0, F2, F4, ..., F16.
C     Calculates X1, X3, ..., X15, F1, F3, ..., F15, QLEFT, QRIGHT, QNOW,
C     QDIFF, and AREA.


   30 IF (ISTAT .EQ. -1) THEN

         DO 35, J = 1, 15, 2
            X (J) = (X (J - 1) + X (J + 1)) * HALF
   35    CONTINUE

         XEVAL = X (1)
         ISTAT = 1
         GO TO 99                    ! Go get F (1)

      END IF

      NOFUN = NOFUN + 1
      F (ISTAT) = FEVAL

      IF (ISTAT .LT. 15) THEN
         ISTAT = ISTAT + 2
         XEVAL = X (ISTAT)
         GO TO 99                    ! Go get another odd function value
      END IF

      STEP = (X (16) - X (0)) / 16.0E+0
      QLEFT = (W0 * (F (0) + F (8)) + W1 * (F (1) + F (7)) +
     >         W2 * (F (2) + F (6)) + W3 * (F (3) + F (5)) +
     >         W4 * F (4)) * STEP
      QRIGHT (LEV + 1) = (W0 * (F (8) + F (16)) + W1 * (F (9) + F (15))
     >                 +  W2 * (F (10)+ F (14)) + W3 * (F (11)+ F (13))
     >                 +  W4 * F (12)) * STEP
      QNOW = QLEFT + QRIGHT (LEV + 1)
      QDIFF = QNOW - QPREV
      AREA = AREA + QDIFF

C     ***   STAGE 4 *** INTERVAL CONVERGENCE TEST

      ESTERR = ABS (QDIFF) / 1023.0E+0
      TOLERR = MAX (ABSERR, RELERR * ABS (AREA)) * (STEP / STONE)
      IF (LEV .LT. LEVMIN) GO TO 50
      IF (LEV .GE. LEVMAX) GO TO 65
      IF (NOFUN .GT. NOFIN) GO TO 60
      IF (ESTERR .LE. TOLERR) GO TO 70

C     ***   STAGE 5   ***   NO CONVERGENCE
C     Locate next interval.

   50 NIM = 2 * NIM
      LEV = LEV + 1

C     Store right hand elements for future use:

      DO 52, I = 1, 8
         FSAVE (I, LEV) = F (I + 8)
         XSAVE (I, LEV) = X (I + 8)
   52 CONTINUE

C     Assemble left hand elements for immediate use:

      QPREV = QLEFT
      DO 55, I = 1, 8
         J = -I
         F (2 * J + 18) = F (J + 9)
         X (2 * J + 18) = X (J + 9)
   55 CONTINUE

      ISTAT = -1
      GO TO 30


C     ***   STAGE 6   ***   TROUBLE HANDLING
C     The number of function values is about to exceed limit.

   60 NOFIN = NOFIN + NOFIN
      LEVMAX = LEVOUT
      FLAG = FLAG + (B - X (0)) / (B - A)
      GO TO 70

C     CURRENT LEVEL IS LEVMAX.

   65 FLAG = FLAG + ONE

C     ***   STAGE 7   ***   INTERVAL CONVERGED
C     Add contributions into the running sums:

   70 RESULT = RESULT + QNOW
      ERREST = ERREST + ESTERR
      COR11  = COR11  + QDIFF / 1023.0E+0

C     Locate the next interval:

   75 IF (NIM .NE. 2 * (NIM / 2)) THEN
         NIM = NIM / 2
         LEV = LEV - 1
         GO TO 75
      END IF

      NIM = NIM + 1
      IF (LEV .LE. 0) GO TO 80

C     Assemble elements required for the next interval:

      QPREV = QRIGHT (LEV)
      X (0) = X (16)
      F (0) = F (16)
      DO 78, I = 1, 8
         F (2 * I) = FSAVE (I, LEV)
         X (2 * I) = XSAVE (I, LEV)
   78 CONTINUE

      ISTAT = -1
      GO TO 30

C     ***   STAGE 8   ***   FINALIZE AND RETURN

   80 RESULT = RESULT + COR11

C     Make sure ERREST is not less than roundoff level:

      IF (ERREST .NE. ZERO) THEN
   90    TEMP = ABS (RESULT) + ERREST
         IF (TEMP .EQ. ABS (RESULT)) THEN
            ERREST = ERREST + ERREST
            GO TO 90
         END IF
      END IF

C     Set final ISTAT according to original error flag:

      IF (FLAG .EQ. ZERO) THEN
         ISTAT = 0
      ELSE
         ISTAT = -2
      END IF

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE RDQPL ( MXPTS, TITLE, NPTS, X, Y, LUNRD )
C
C PURPOSE: RDQPL is the complement of SAVQPL (both of which are really
C          obsolete now).  RDQPL reads one curve (with a title) from a
C          file in QPLOT format.  No attempt is made to handle QPLOT's
C          optional namelist. In fact, any text records other than the
C          first record are ignored.  Use RDQPL with discretion, since
C          any curve other than the first in a file is not  likely  to
C          have a useful title associated with  it.
C
C METHOD:  *  Read a title as the first error.  Any error is fatal.
C
C          *  Read (x,y) pairs in a counting loop (list-directed I/O).
C             Use ERR=... to skip unwanted labels and legend info etc.
C             until no error is encountered,  then  keep reading until
C             either an EOF or another error is encountered.  No BACK-
C             SPACE is done at such an error in case it was due to  an
C             END FRAME record (meaning valid title probably follows).
C
C ARGUMENTS:
C
C    ARG     DIM     TYPE I/O/S DESCRIPTION
C  MXPTS      -        I    I   Max. no. of pts. allowed in one curve
C  TITLE      -      C*(*)  O   Title found for curve
C  NPTS       -        I    O   No. of data pts. found on this call
C  X, Y      NPTS      R    O   Coordinates found as one curve.
C  LUNRD      -        I    I   Logical unit number for QPLOT file.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.  12/05/83
C         <Resurrected 08/25/88 - missed being moved to OBSOLETE lib.>
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C     Arguments:

      INTEGER   LUNRD, MXPTS, NPTS
      REAL      X(MXPTS), Y(MXPTS)
      CHARACTER TITLE*(*)

C     Execution.

C     Look for a title record:

      READ ( LUNRD, '(A)', ERR=70 ) TITLE

   30 READ ( LUNRD, *, ERR=30, END=80 ) X(1), Y(1)

      NPTS = 1

   40 CONTINUE
         NPTS = NPTS + 1
         IF ( NPTS .GT. MXPTS ) GO TO 90

         READ ( LUNRD, *, ERR=50, END=50 ) X(NPTS), Y(NPTS)
      GO TO 40

   50 NPTS = NPTS - 1
      IF ( NPTS .GE. 1 ) RETURN

      STOP 'RDQPL: No data points found.'
   70 STOP 'RDQPL: Unable to read first record as text.'
   80 STOP 'RDQPL: Unexpected EOF.'
   90 STOP 'RDQPL: Too many data points found.'

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE RDREALS (SCREEN, PROMPT, KEYBRD, NUMBER, REALS)
C
C
C     Description and usage:
C
C        RDREALS reads an indefinite number of real values from one line or
C     record, with an optional prompt.  Standard delimiters (comma, blank,
C     or tab) are assumed between values.
C
C        One motivation here is to overcome an inherent limitation in the
C     earlier READER utility (which can accept just ONE value per prompt)
C     while still offering READER's error handling (reprompting) and its
C     flexibility with respect to null inputs (<CR> only).  [Standard list-
C     directed reads cannot handle indefinite lists.]  RDREALS also offers
C     more flexible cursor control.
C
C        RDREALS is appropriate for entering multiple values into a single
C     array.  Example: prompting for an arbitrary number of phase angles in
C     an oscillation.  Higher level routine RDTUPLES should be used for
C     entering (possibly many) pairs or triples.
C
C        RDLIST is available for the [in]definite-list-of-integers case.
C     There is no integer analog of RDTUPLES.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        SCREEN              I    I      Logical unit number to which the
C                                        prompt and diagnostics are written.
C        PROMPT              C    I      Prompt string with carriage control
C                                        in PROMPT (1:1).  (Use ' ', '0', '1',
C                                        or '$' in the usual way.)
C                                        If PROMPT is blank it is suppressed.
C        KEYBRD              I    I      Logical unit number from which data
C                                        is to be read.  (May be a disk file if
C                                        called by RDTUPLES in indirect mode.)
C        NUMBER              I    I/O    Input with the maximum number of
C                                        values provided for.  On output,
C                                        NUMBER > 0 is the no. of values read;
C                                        NUMBER = 0 corresponds to null (CR);
C                                        NUMBER =-1 corresponds to end-of-data;
C                                        NUMBER =-2 corresponds to bad data if
C                                                   PROMPT is blank (else a
C                                                   reprompt is issued here);
C                                        NUMBER =-3 means a system-dependent
C                                                   input error was encountered.
C        REALS   NUMBER      R      O    Array of reals entered (if any).
C
C
C     Implementation Notes:
C
C        The prompting part is done in-line because READER's only suitable
C     entry point, READS, does not offer the desired cursor control.
C
C        Handling of more than one line's worth of data (132 characters) is not
C     attempted in this version, but a trailing hyphen (say) could serve as a
C     flag.
C
C        RDTUPLES makes use of RDREALS, and this application required three
C     minor features here that are not expected to be made use of directly:
C
C        >  the prompt may be suppressed;
C        >  indirect input may be REPORTed (but not DONE except upon further
C           calls from RDTUPLES).  An internal COMMON was reluctantly used to
C           return the @filename string (preferable to another argument);
C        >  bad data leads to termination in indirect mode (whereas reprompting
C           is appropriate interactively).  A blank PROMPT is taken to mean
C           indirect and hence an error return if an invalid value is read.
C
C
C     Procedures:
C
C        GETLINE  Reads one record, with option to strip off trailing comments.
C        SCAN2    Finds positions of the first and last significant characters.
C
C
C     Environment:  DEC VAX/VMS FORTRAN, including:
C
C        >  IMPLICIT NONE
C        >  Trailing ! comments
C        >  Some 8-character variable names
C
C
C     Author:  David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C     02/20/89  DAS  Initial implementation, derived from RDLIST for the
C                    case of real values, with improved cursor control.
C     03/01/89  DAS  (With Robert Kennelly:)  Refinements indicated above
C                    to facilitate the RDTUPLES application.
C     02/21/90  DAS  IRIS 4D requires more cumbersome carriage control:
C                    print 'prompt' from PROMPT = '$prompt' with FORMAT
C                    (1X, A, $), else just use FORMAT (A) for entire PROMPT.
C     10/20/99  DAS  Fortran 90 version (avoiding $ carriage control).
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      CHARACTER
     >   PROMPT * (*)
      INTEGER
     >   KEYBRD, NUMBER, SCREEN
      REAL
     >   REALS (*)

C     Internal COMMON (required by the RDTUPLES application of RDREALS). (Ugh!)

      INTEGER, PARAMETER ::
     >   MAXBUF = 132       ! Arbitrary limit on input record length
      CHARACTER
     >   BUFFER * (MAXBUF)
      COMMON /TUPLES/
     >   BUFFER

C     Local constants:

      CHARACTER, PARAMETER ::
     >   BLANK    * 1 = ' ',
     >   COMMENT  * 1 = '!',
     >   FORMR    * 7 = '(F25.0)',  ! Arbitrary limit on numerical token length
     >   INDIRECT * 1 = '@'

C     Local variables:

      INTEGER
     >   COUNT, FIRST, LAST, MARK, STATUS
      LOGICAL
     >   INTERACT
      CHARACTER
     >   SEPS * 3

C     Procedures:

      EXTERNAL
     >   GETLINE, SCAN2

C     Execution:

      SEPS = ' ,' // CHAR (9)   ! Token separators are blank, comma, or tab
      INTERACT = PROMPT /= BLANK

  200 CONTINUE
      COUNT = 0
      IF (INTERACT) THEN
         IF (PROMPT (1 : 1) == '$') THEN
            WRITE (SCREEN, 1020, ADVANCE = 'NO') PROMPT (2 : )
         ELSE
            WRITE (SCREEN, 1020) PROMPT (2 : )
         END IF
      END IF

      CALL GETLINE (KEYBRD, COMMENT, BUFFER, LAST, STATUS)

      IF (STATUS > 0) THEN                     ! Read error
         IF (INTERACT) GO TO 200
         GO TO 800
      END IF
      IF (STATUS < 0) GO TO 900                ! ^Z = EOF = "quit"
      IF (LAST  == 0) GO TO 910                ! CR (or insignificant line)

C     Convert character list to reals one token at a time:

      FIRST = 1
  300 CONTINUE

C        Delimit a token, if any:

         CALL SCAN2 (BUFFER (1 : LAST), SEPS, FIRST, LAST, MARK)

         IF (MARK > 0) THEN

C           The following is relevant to RDTUPLES only.

            IF (COUNT == 0) THEN                              ! Indirection?
               IF (BUFFER (FIRST : FIRST) .EQ. INDIRECT) THEN ! E.g. @XY.DAT
                  BUFFER (FIRST : FIRST) = BLANK         ! Lets RDTUPLES look
                  BUFFER (1 : 1) = INDIRECT              ! in a fixed place.
                  NUMBER = 1           ! Helps termination test in RDTUPLES.
                  GO TO 999
               END IF
            END IF

            COUNT = COUNT + 1

C           Convert the token to a real value (internal read).

            READ (BUFFER (FIRST : MARK), FORMR, IOSTAT = STATUS)
     >         REALS (COUNT)
            IF (STATUS > 0) GO TO 810

            FIRST = MARK + 2
            IF (COUNT .LT. NUMBER) GO TO 300    ! Look for more if there's room,
                                                ! else ignore further entries.
         END IF

      GO TO 910

C     Error handling:

  800 WRITE (SCREEN, 1000) ' *** System error.  IOSTAT: ', STATUS
      IF (INTERACT) GO TO 200
      NUMBER = -3
      GO TO 999

  810 WRITE (SCREEN, 1010) ' *** Bad value found: ',
     >   BUFFER (FIRST : MARK)
      IF (INTERACT) GO TO 200
      NUMBER = -2
      GO TO 999

C     Termination:

  900 NUMBER = -1       ! EOF or ^Z with COUNT = 0
      GO TO 999

  910 NUMBER = COUNT    ! CR, insignificant line, or '@...' with COUNT = 0,
                        ! else valid list with COUNT > 0

  999 RETURN

C     Formats:

 1000 FORMAT (A, I6)
 1010 FORMAT (A, A)
 1020 FORMAT (1X, A)

      END SUBROUTINE RDREALS
C+----------------------------------------------------------------------
C
      SUBROUTINE RDXYZ (NDIM, LUNDAT, LUNERR, COLUMNS, MAXPTS,
     >                  NPTS, X, Y, Z, FINIS, IER)
C
C  PURPOSE:
C
C        RDXYZ reads one set of (X), (X,Y), or (X,Y,Z) points per call
C     from the indicated columns of a file which is assumed to be open
C     already on unit LUNDAT.
C
C        For historical reasons, the data may take one of two forms (with
C     or without integer counts, basically, with the latter now preferred):
C
C     (1) "SMOOTH" or "PROFILE" format: one or more "curves" with explicit
C         integer count preceding the real values for each "curve."  "NPTS"
C         is identified as a single leading numeric token, with no decimal
C         point.  The value of "NPTS" should match the number of "good"
C         points.
C
C     (2) "Indefinite" format: no explicit counts of the points in each
C         "curve."  A line with "END" for the first signficant characters
C         serves to separate "curves." 
C
C        In both cases, any title line (normally the first line of the
C     file) should be handled externally by the application.  See programs
C     SMOOTH and PROFILE for sample usages.  Other points to note:
C
C        > The file is assumed to be positioned ready to read the next "curve."
C          A value of NPTS = 0 is valid for data form (1).
C
C        > Blank lines and any strings starting with '!' are ignored, meaning
C          data points may be annotated with trailing comments or commented out
C          altogether.
C
C        > EOF may be encountered on the first (significant) read - it is up
C          to the application to know whether this is normal or not.
C
C        > An EOF encountered after the first significant read (causing FINIS
C          to be set .TRUE.) is normal for data form (2) but will be abnormal
C          for data form (1) - this is the only case of a mismatched "NPTS"
C          trapped.
C
C  ARGUMENTS:
C     ARG      DIM   TYPE  I/O/S   DESCRIPTION
C     NDIM      -      I     I     NDIM = 1, 2, or 3 depending on whether
C                                  (X), (X,Y), or (X,Y,Z) are to be read.
C                                  Pass X, X, X or X, Y, Y if NDIM = 1 or 2.
C     LUNDAT    -      I     I     Logical unit number for file being read.
C     LUNERR    -      I     I     Logical unit number for error messages.
C     COLUMNS  NDIM    I     I     Column number(s) for X (,Y (,Z)).
C     MAXPTS    -      I     I     Maximum number of points in any one
C                                  dataset expected by calling program.
C     NPTS      -      I     O     Number of data points found; NPTS >= 0.
C     X,Y,Z   MAXPTS   R     O     Data points found.  See NDIM.
C     FINIS     -      L     O     .TRUE. means EOF was encountered after
C                                  the first significant read - see above.
C     IER       -      I     O     IER=0 means one dataset found normally;
C                                        NPTS may still be 0.
C                                  IER=1 means EOF encountered on the first
C                                        significant read - may be normal.
C                                        FINIS is set too (but redundant).
C                                  IER=2 means an abnormal error, for which
C                                        RDXYZ will have issued a diagnostic:
C                                        out-of-range NPTS; conversion error;
C                                        or missing Y or Z.
C  PROCEDURES:
C     GETLINE    Reads a line as text; handles suppressed pts. & comments
C     TOKENS     Tokenizes a string
C     UPCASE     Converts a string to upper case
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77 with ! comments, IMPLICIT NONE, and
C                some variable names up to 8 characters.
C
C  HISTORY:
C     07/31/90  DAS  Original implementation, when PROFILE's PRREAD was
C                    found inappropriate for reading 2nd derivatives.
C                    Program SMOOTH can use the same capability, since
C                    the <original> data format handled here is that of
C                    <the original> SMOOTH.
C     12/10/91  DAS  Added multi-column and indefinite-number-of-points
C                    capabilities, mainly for SMOOTH purposes.  Upward
C                    compatibility is awkward!  (Handling all cases of
C                    mismatches between "NPTS" and actual number found
C                    for the original data format would be more trouble
C                    than it is worth, so only EOF with too few points
C                    is trapped and warned of properly.)  Maybe the line
C                    buffer should be in/out (as in READCOLS now), but
C                    the present compromise still avoids backspacing
C                    and should be simpler to reuse.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NDIM, COLUMNS (NDIM), IER, LUNDAT, LUNERR, MAXPTS, NPTS
      REAL
     >   X (MAXPTS), Y (MAXPTS), Z (MAXPTS)
      LOGICAL
     >   FINIS

C     Local constants:

      INTEGER
     >   LENTOK, MAXCOL
      CHARACTER
     >   COMMENT * 1, IFMT * 8, RFMT * 10
      PARAMETER
     >  (COMMENT = '!',
     >   IFMT    = '(BN,I25)',   ! Should match LENTOK
     >   LENTOK  = 25,           ! Largest # characters in any numeric value
     >   MAXCOL  = 20,           ! Largest column number handled
     >   RFMT    = '(BN,F25.0)') ! Should match LENTOK

C  *  Local variables:

      INTEGER
     >   I, IOS, J, LAST, LASTCOL, NGIVEN, NUMBER
      LOGICAL
     >   HAVELINE
      CHARACTER
     >   LINE * 132, LIST (MAXCOL) * (LENTOK)

C  *  Procedures:

      EXTERNAL
     >   GETLINE, TOKENS, UPCASE

C     Execution:

      IER = 0
      NPTS = 0
      NGIVEN = 0
      FINIS = .FALSE.

C     Look for the start of a "curve."  The original "NPTS" dataset form
C     forces this special treatment at the beginning.

  100 CALL GETLINE (LUNDAT, COMMENT, LINE, LAST, IOS)
      IF (IOS .LT. 0) THEN                      ! Probably a normal EOF
         IER = 1
         FINIS = .TRUE.                         ! Probably redundant
         GO TO 999
      END IF
      IF (IOS .NE. 0) GO TO 900                 ! Read error
      IF (LAST .EQ. 0) GO TO 100                ! Insignificant line

      LASTCOL = COLUMNS (1)
      IF (NDIM .GE. 2) LASTCOL = MAX (LASTCOL, COLUMNS (2))
      IF (NDIM .EQ. 3) LASTCOL = MAX (LASTCOL, COLUMNS (3))

      NUMBER = LASTCOL
      CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)

      HAVELINE = .TRUE.
      IF (NUMBER .EQ. 1) THEN                   ! Must be NPTS unless NDIM = 1

         IF (INDEX (LIST (1), '.') .EQ. 0) THEN    ! No decimal point
            READ (LIST (1), IFMT, ERR=905) NGIVEN     ! Try to read NPTS

            IF (NGIVEN .LT. 0 .OR. NGIVEN .GT. MAXPTS) GO TO 915  ! Bad NPTS
            IF (NGIVEN .EQ. 0) GO TO 999                          ! NPTS = 0 OK
            HAVELINE = .FALSE.            ! Need to get the first data proper
         END IF

      ELSE              ! More than one token found - assumed to be data proper
      END IF

C     Start of loop over the data proper.

      I = 0

  200 CONTINUE
         IF (HAVELINE) THEN    ! Annoyance caused by original "NPTS" case
            HAVELINE = .FALSE.
         ELSE
            CALL GETLINE (LUNDAT, COMMENT, LINE, LAST, IOS)

            IF (IOS .LT. 0) THEN                   ! EOF - probably OK
               FINIS = .TRUE.                      ! Leave IER = 0
               IF (I .LT. NGIVEN) GO TO 925        ! "NPTS" was too big
               GO TO 999
            END IF

            IF (IOS .NE. 0) GO TO 900              ! Read error
            IF (LAST .EQ. 0) GO TO 200             ! Insignificant line

            NUMBER = LASTCOL
            CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)
         END IF

C        Look for end of "curve" (indefinite data format case):

         CALL UPCASE (LIST (1) (1 : 3))
         IF (LIST (1) (1 : 3) .EQ. 'END') GO TO 999

         I = I + 1
         IF (I .GT. MAXPTS) GO TO 915              ! No more room in X/Y/Z
         IF (NUMBER .LT. LASTCOL) GO TO 920        ! Missing value

C        At last we can do the actual (internal) read:

         J = COLUMNS (1)
         READ (LIST (J), RFMT, ERR=910) X (I)
         IF (NDIM .GE. 2) THEN
            J = COLUMNS (2)
            READ (LIST (J), RFMT, ERR=910) Y (I)
            IF (NDIM .EQ. 3) THEN
               J = COLUMNS (3)
               READ (LIST (J), RFMT, ERR=910) Z (I)
            END IF
         END IF
         NPTS = I
         IF (I .NE. NGIVEN)         ! This works for both data formats
     >GO TO 200      

      GO TO 999


C     Error handling:

  900 WRITE (LUNERR, 1001) 'Bad return from GETLINE.'
      GO TO 990

  905 WRITE (LUNERR, 1002) 'Error reading an integer from ', LIST (1)
      GO TO 990

  910 WRITE (LUNERR, 1002) 'Error reading a real value from ', LIST (J)
      GO TO 990

  915 WRITE (LUNERR, 1001) 'Bad number of points found:', NPTS,
     >                     'Maximum provided for:', MAXPTS
      GO TO 990

  920 WRITE (LUNERR, 1001)
     >   'Real value missing from current dataset.  Line #:', I
      GO TO 990

  925 WRITE (LUNERR, 1001)
     >   'Warning: EOF before indicated # pts. found:', NGIVEN,
     >   '         Proceeding with # pts. read:      ', I
      GO TO 999

  990 IER = 2
C*****GO TO 999

  999 RETURN

C     Formats:

 1001 FORMAT (/, ' *** RDXYZ: ', A, I8)
 1002 FORMAT (/, ' *** RDXYZ: ', A, A)

      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE READER (CASE, SCREEN, PROMPT, KEYBRD, CR, EOF,
     >                   INT4, REAL4, REAL8, CHARS, STRING, YES)
C
C
C     Description and usage:
C
C           READER provides a set of utilities offering some of the
C        flexibility of list-directed input to interactive FORTRAN programs,
C        while permitting use of <CR> (i.e. Carriage Return) as a legal,
C        identifiable response.  This is impossible with a list-directed
C        read, which ignores <CR> and continues to wait for input.  In
C        case of bad input (presumably due to user error), a warning is
C        displayed and the prompt repeated.  "Quit" (as opposed to "default")
C        is provided for by the "End-of-file" argument EOF.  (See Notes.)
C        If CR or EOF is .TRUE., the argument value returned is unchanged,
C        so these flags should be checked by the calling routine.  (However,
C        assigning a default value before calling READER can save checking
C        the CR argument.)
C
C           This version avoids the original's ENTRY points for the different
C        data types because of f90 compiler difficulties on an SGI system.
C
C           READER itself should not be called directly. Rather, six ancillary
C        routines provide for reading INTEGER, REAL, DOUBLE PRECISION, STRING,
C        CHARACTER, or YES/NO data, and each makes the appropriate call to
C        READER.  The names are READx, where x = I, R, D, S, C, or Y.  These
C        ancillary routines appear at the end of the READER source module.
C
C           READC returns all the non-blank characters found, converted to
C        upper case if alphabetic, and packed from the left into the output
C        string with blanks removed.  The READS option merely returns the
C        string as entered from the keyboard, without modification.  READC
C        is normally appropriate for entering a single item or token (such as
C        an identifier or a file name), while READS is appropriate for literal
C        strings such as plot titles.  HOWEVER: since Unix file names are case
C        sensitive, READS is actually the better choice for file name prompts
C        on Unix systems.
C
C           All calling sequences are identical except for the type of the
C        value to be returned.
C
C           Sample usage:
C
C           :      :                              (^D under Unix)
C       210 NPTS = 100                              |
C           CALL READI (LUNCRT,                     |
C          >   'Enter number of points.  <CR>=100; ^Z=quit: ',
C          >   LUNKBD, NPTS, CR, EOF)
C           IF (EOF) GO TO 999
C           IF (NPTS .LT. 1 .OR. NPTS .GT. MXPTS) GO TO 210
C           :      :
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        CASE                C    I      'I', 'R', 'D', 'S', 'C', or 'Y'
C                                        (but see use of READI, etc., above)
C        SCREEN              I    I      Logical unit number to which
C                                        the prompt is written (often 6).
C        PROMPT    *         C    I      Prompt string to be written to
C                                        unit SCREEN.
C        KEYBRD              I    I      Logical unit number from which
C                                        data is to be read (often 5).
C
C
C   ---> Only ONE of the following six output choices is to be used:
C
C        INT4                I      O    INTEGER quantity to be returned.
C        REAL4               R      O    REAL quantity to be returned.
C        REAL8               D      O    DOUBLE PRECISION quantity to be
C                                        returned.
C        STRING              C      O    A CHARACTER string, as entered.
C        CHARS               C      O    A CHARACTER string, converted to
C                                        upper case if alphabetic, packed
C                                        and left-justified.
C        YES                 L      O    Logical flag set .TRUE. for a 'YES'
C                                        response (first character 'Y' or 'y')
C                                        and .FALSE. for 'N' or 'n'.  Other
C                                        inputs force a reprompt.
C
C        CR                  L      O    CR = .TRUE. if a null value (i.e.
C                                        Carriage Return only) was read,
C                                        else CR = .FALSE.
C        EOF                 L      O    EOF = .TRUE. if "End-of-File"
C                                        was detected (^Z under VMS; ^D
C                                        under Unix), else EOF = .FALSE.
C                                        
C
C     External devices:  See arguments SCREEN and KEYBRD.
C
C
C     External procedures:  UPCASE is used by the READC and READY options
C                           in place of the original in-line code - see Notes.
C
C
C     Environment:  DEC/OpenVMS Fortran 90
C                   SGI/IRIX    Fortran 90
C
C     Method:
C
C           The user is prompted for input, with the cursor left at the
C        end of the prompt if there is room for the expected response.  The
C        reply is read into a buffer as a character string and, except for
C        the READS and READC entries, left-justified with blanks squeezed
C        out and re-read internally in the requested format.  For READC, the
C        same repacking is done, then the result is converted to upper case.
C        No repacking or conversion is done in the case of READS.
C
C
C     Notes:
C
C        (1)  Conversions to upper case were originally done in-line by
C             adding the appropriate offset (= ICHAR ('A') - ICHAR ('a')).
C             However this has been shown to fail on some systems and has
C             been replaced by use of the UPCASE utility, which cannot fail.
C
C        (2)  The length of the string returned by READS and READC is set by
C             the calling program - excess characters read will be truncated
C             from the right.  No error flag is set for this condition.  The
C             maximum length for CHARS and STRING is MAXBUF.
C
C        (3)  The '$' carriage control character used for short prompts is
C             not compatible with Fortran 90 - this version uses ADVANCE='NO'
C             instead.
C
C
C     Author:  Robert Kennelly, NASA Ames/Sterling Software, Palo Alto, CA.
C
C
C     Development history:
C
C        27 Sep. 1983    RAK    Initial design and coding.
C         6 Oct. 1983    RAK    Extended to multiple data types.
C         2 Nov. 1983  RAK/LJC  Added READS for strings without modification.
C        16 Mar. 1984  RAK/LJC  Included "BN" specifier in formats to avoid
C                               unneccessary packing.
C        14 June 1984  DAS/RAK  Corrected header, TEST always defined, added
C                               CC, eliminated redundant trap on CASE, used
C                               shorter flags, e.g. 'R' instead of 'REAL'.
C         7 Sep. 1984    RAK    Leave STRING unchanged when CR entered,
C                               thus consistent with the other modes.
C        19 Oct. 1984    RAK    Initialize BUFFER prior to READ.
C        14 Dec. 1984    RAK    Make sure length of CHARS is not exceeded
C                               during repacking/case-conversion.
C        19 Dec. 1985  RAK/DBS  Addition of READY entry, and general 
C                               streamlining of code.
C        30 Dec. 1985  RAK/DAS  Edited header.
C        09 Sep. 1988    DAS    Sample usage added above; other cosmetics.
C        05 May  1989  DAS/RGL  Unix- and VMS-compatible now: ^Z/^D usage
C                               documented; STOP 'READER' message <=8 chars.;
C                               OFFSET not defined as a PARAMETER constant.
C                               (UNICOS displays '$' in column 1 at time of
C                               writing, but this may go away...)
C        01 Feb. 1990    DAS    '$' carriage control changed to (A,$) form
C                               (see Notes above; thanks to Scott Thomas).
C                               READS recommended over READC for file names
C                               on Unix systems.  Conversions to upper case
C                               now done via UPCASE - see Notes above.
C        23 July 1997     "     Replaced (A,$) with ADVANCE='NO' for Fortran 90.
C        24 July 1997     "     Eliminated ENTRY points as a work-around for an
C                               SGI compiler.  (See 6 ancillary routines below.)
C        19 May  1999     "     The CASE argument had not been described.
C        20 May  1999     "     Dexter Hermstad solved a control-D (EOF) problem
C                               encountered with IRIX 6 f90, making this version
C                               SGI-specific.
C        23 July 2004     "     Raised all the 80s to 132 for RDLIST purposes.
C        13 Feb  2006     "     IMPLICIT NONE was missing from the ancillary
C                               routines introduced to avoid ENTRY points.
C                               The missing declarations may overcome problems
C                               encountered with program CAVITY_MAP when it is
C                               compiled with the pgf90 compiler.
C-------------------------------------------------------------------------------


C     DECLARATIONS.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   INT4, KEYBRD, SCREEN
      REAL(KIND=4)                       ! RLC 20211025
     >   REAL4
      REAL(KIND=8)                     ! RLC 20211025
     >   REAL8
      CHARACTER
     >   CASE * 1, CHARS * (*), PROMPT * (*), STRING * (*)
      LOGICAL
     >   YES, CR, EOF

C     Local constants.  Note that the field lengths in the format strings
C     are hardcoded to MAXBUF.

      INTEGER
     >   MAXBUF, MAXLIN
      CHARACTER
     >   BLANK, FORMI * 10, FORMR * 12
      PARAMETER
     >  (BLANK  = ' ',
     >   FORMR  = '(BN, F132.0)',
     >   FORMI  = '(BN, I132)',
     >   MAXBUF = 132,
     >   MAXLIN = 132)

C     Local variables.

      LOGICAL
     >   ADVANCE
      INTEGER
     >   FILL, SCAN, STATUS, TEST
      CHARACTER
     >   BUFFER * (MAXBUF)

C     Procedures.

      EXTERNAL
     >   UPCASE


C     EXECUTION.
C     ----------

C     Set a typical maximum length for response to the non-character
C     entries.  This value and the length of the prompt string are used
C     to determine where the cursor is to be positioned after prompting.

      TEST = 16

      IF (CASE .EQ. 'S') THEN
         TEST = LEN (STRING)
      ELSE IF (CASE .EQ. 'C') THEN
         TEST = LEN (CHARS)
      ELSE IF (CASE .EQ. 'Y') THEN
         TEST = 5
      END IF

      ADVANCE = TEST + LEN (PROMPT) .GT. MAXLIN

C     Prompt user for input.
C     ----------------------

   30 CONTINUE
         BUFFER = BLANK

C        Issue prompt, then read up to MAXBUF characters from the keyboard.

         IF (ADVANCE) THEN
            WRITE (SCREEN, '(1X, A)') PROMPT
         ELSE
            WRITE (SCREEN, '(1X, A)', ADVANCE='NO') PROMPT
         END IF

         READ (KEYBRD, '(A)', IOSTAT = STATUS) BUFFER

C        Re-try on errors, exit on End-of-File, or continue.

         IF (STATUS .GT. 0) THEN
            WRITE (SCREEN, 1000)
            GO TO 30

         ELSE IF (STATUS .LT. 0) THEN ! SGI f90 won't allow further ^Ds
            EOF = .TRUE.              ! unless the keyboard is reopened
            CR = .FALSE.
            OPEN (KEYBRD, FILE='/dev/tty', STATUS='OLD') ! Dexter's fix
            GO TO 990

         ELSE
            EOF = .FALSE.
            CR = (BUFFER .EQ. BLANK)
         END IF


C        Process the data in the buffer, if any.
C        ---------------------------------------

         IF (.NOT. CR) THEN

            IF (CASE .EQ. 'I') THEN
               READ (BUFFER, FORMI, IOSTAT = STATUS) INT4

            ELSE IF (CASE .EQ. 'R') THEN
               READ (BUFFER, FORMR, IOSTAT = STATUS) REAL4

            ELSE IF (CASE .EQ. 'D') THEN
               READ (BUFFER, FORMR, IOSTAT = STATUS) REAL8

            ELSE IF (CASE .EQ. 'S') THEN

C              For literal strings, copy the buffered input directly to
C              the output variable.  If CR is true, we leave STRING alone
C              so that any default value set in the calling routine will
C              remain intact.  No error checking is required.

               STRING = BUFFER

            ELSE 

C              CASE is either 'C' or 'Y' - scan the input buffer, starting
C              from the left, to pack the data and count the non-blanks.

               FILL = 0
               DO SCAN = 1, MAXBUF
                  IF (BUFFER (SCAN : SCAN) .NE. BLANK) THEN
                     FILL = FILL + 1
                     BUFFER (FILL : FILL) = BUFFER (SCAN : SCAN)
                  END IF
               END DO

C              Convert to upper case.

               CALL UPCASE (BUFFER (1 : FILL))

               IF (CASE .EQ. 'C') THEN

                  CHARS = BUFFER (1 : FILL)

               ELSE

C                 CASE must be 'Y' - all we need is the first character.

                  IF (BUFFER (1 : 1) .EQ. 'Y') THEN
                     YES = .TRUE.
                  ELSE IF (BUFFER (1 : 1) .EQ. 'N') THEN
                     YES = .FALSE.
                  ELSE

C                    Buffer value is invalid - generate an ersatz input error
C                    to force a re-prompt.

                     STATUS = 999
                  END IF
               END IF
            END IF 

C           Re-try on errors, exit on End-of-File, else continue.
C           (Dropping through for 'C', 'S', 'Y' cases is OK and saves code.)

            EOF = (STATUS .LT. 0)
            IF (STATUS .GT. 0) THEN
               WRITE (SCREEN, 1000)
               GO TO 30

            END IF

      END IF


C     TERMINATION.
C     ------------

  990 RETURN

C     FORMATS.
C     --------

 1000 FORMAT (' Input error!  Please try again.')

      END SUBROUTINE READER


C     Ancillary subroutines used to avoid entry points in READER.
C     -----------------------------------------------------------
C
C     An SGI compiler also objects to omitted arguments (.., , ,),
C     so dummy arguments are passed as a work-around.

      SUBROUTINE READI (SCREEN, PROMPT, KEYBRD, INT4,   CR, EOF)

      IMPLICIT   NONE

C     Arguments:

      INTEGER    INT4, KEYBRD, SCREEN
      CHARACTER  PROMPT * (*)
      LOGICAL    CR, EOF

C     Local variables:

      REAL(KIND=4):: R4                            ! RLC 20211025
      REAL(KIND=8):: R8                                      ! RLC 20211025
      LOGICAL    Y

      CALL READER ('I', SCREEN, PROMPT, KEYBRD, CR, EOF,
     >             INT4, R4, R8, 'C', 'S', Y)
      END

      SUBROUTINE READR (SCREEN, PROMPT, KEYBRD, REAL4,  CR, EOF)

      IMPLICIT   NONE

C     Arguments:

      INTEGER    KEYBRD, SCREEN
      REAL(KIND=4):: REAL4                               ! RLC 20211025
      CHARACTER  PROMPT * (*)
      LOGICAL    CR, EOF

C     Local variables:

      INTEGER    I4
      REAL(KIND=8):: R8                                       ! RLC 20211025
      LOGICAL    Y

      CALL READER ('R', SCREEN, PROMPT, KEYBRD, CR, EOF,
     >             I4, REAL4, R8, 'C', 'S', Y)
      END

      SUBROUTINE READD (SCREEN, PROMPT, KEYBRD, REAL8,  CR, EOF)

      IMPLICIT   NONE

C     Arguments:

      INTEGER    KEYBRD, SCREEN
      REAL(KIND=8)::   REAL8                               ! RLC 20211025
      CHARACTER  PROMPT * (*)
      LOGICAL    CR, EOF

C     Local variables:

      INTEGER    I4
      REAL(KIND=4):: R4                                            ! RLC 20211025
      LOGICAL    Y

      CALL READER ('D', SCREEN, PROMPT, KEYBRD, CR, EOF,
     >             I4, R4, REAL8, 'C', 'S', Y)
      END

      SUBROUTINE READC (SCREEN, PROMPT, KEYBRD, CHARS,  CR, EOF)

      IMPLICIT   NONE

C     Arguments:

      INTEGER    KEYBRD, SCREEN
      CHARACTER  CHARS * (*), PROMPT * (*)
      LOGICAL    CR, EOF

C     Local variables:

      INTEGER    I4
      REAL(KIND=4):: R4                           ! RLC 20211025
      REAL(KIND=8):: R8                           ! RLC 20211025
      LOGICAL    Y

      CALL READER ('C', SCREEN, PROMPT, KEYBRD, CR, EOF,
     >             I4, R4, R8, CHARS, 'S', Y)
      END

      SUBROUTINE READS (SCREEN, PROMPT, KEYBRD, STRING, CR, EOF)

      IMPLICIT   NONE

C     Arguments:

      INTEGER    KEYBRD, SCREEN
      CHARACTER  PROMPT * (*), STRING * (*)
      LOGICAL    CR, EOF

C     Local variables:

      INTEGER    I4
      REAL(KIND=4):: R4                      ! RLC 20211025
      REAL(KIND=8):: R8                      ! RLC 20211025
      LOGICAL    Y

      CALL READER ('S', SCREEN, PROMPT, KEYBRD, CR, EOF,
     >             I4, R4, R8, 'C', STRING, Y)
      END

      SUBROUTINE READY (SCREEN, PROMPT, KEYBRD, YES,    CR, EOF)

      IMPLICIT   NONE

C     Arguments:

      INTEGER    KEYBRD, SCREEN
      CHARACTER  PROMPT * (*)
      LOGICAL    YES, CR, EOF

C     Logical variables:

      INTEGER    I4
      REAL(KIND=4):: R4                      ! RLC 20211025
      REAL(KIND=8):: R8                      ! RLC 20211025

      CALL READER ('Y', SCREEN, PROMPT, KEYBRD, CR, EOF,
     >             I4, R4, R8, 'C', 'S', YES)
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE ROTATE2D (N, X, Y, ANGLE, P, Q)
C
C  PURPOSE:
C     ROTATE2D rotates the given point(s) (X, Y) about the point (P, Q)
C     by the indicated angle.  Input values are overwritten.
C
C  ARGUMENTS:
C  ARG    TYPE  I/O/S   DIM   DESCRIPTION
C  N       I    I        -    Number of points to rotate
C  X,      R    I/O      N    Coordinates of point(s)
C   Y
C  ANGLE   R    I        -    Angle in degrees; positive is anticlockwise
C  P,      R    I        -    Coordinates of center of rotation
C   Q
C
C  ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C  HISTORY:  06/29/90  DAS  Initial implementation
C            08/11/94  DAS  Switched to COSD, SIND form.
C
C  AUTHOR: David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

C     Arguments:

      INTEGER
     >   N
      REAL
     >   X (N), Y (N), ANGLE, P, Q

C     Local variables:

      INTEGER
     >   I
      REAL
     >   CI, SI, XP, YP

C     Intrinsics:

      REAL       COSD, SIND

C     Execution:

      CI = COSD (ANGLE)
      SI = SIND (ANGLE)

      DO 100, I = 1, N
         XP = X (I) - P
         YP = Y (I) - Q
         X (I) = XP * CI - YP * SI + P
         Y (I) = XP * SI + YP * CI + Q
  100 CONTINUE

      RETURN
      END
C+---------------------------------------------------------------------
C
      SUBROUTINE RVERSE ( NX, X1, X2 )
C
C  PURPOSE:  RVERSE reverses the order of the elements of X1(*), re-
C            turning the result in X2(*). In-place reordering is OK.
C
C  ARGUMENTS:
C  ARG  DIM   TYPE I/O/S DESCRIPTION
C  NX   -       I    I   No. of elements in given array (odd or even).
C  X1   NX      R    I   Array whose ordering is to be reversed.
C  X2   NX      R    O   Reordered version of X1(*).  The calling pro-
C                        gram may pass the same array for X1 & X2.
C
C  ENVIRONMENT:  FORTRAN IV
C
C  AUTHOR: David Saunders, Informatics, Palo Alto, CA.  (07/30/82)
C
C----------------------------------------------------------------------
C
      DIMENSION  X1(NX), X2(NX)
C
      NXBY2 = (NX+1)/2
C
      DO 20 I = 1, NXBY2
         XI     = X1(I)
         II     = NX + 1 - I
         X2(I)  = X1(II)
         X2(II) = XI
 20   CONTINUE
C
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE SCAN2 (STRING, SEPS, FIRST, LAST, MARK)
C
C
C     Description and usage:
C
C           Looks for non-blank fields ("tokens") in a string, where the
C        fields are of arbitrary length and separated by any of a set of
C        user-specified separators (e.g., blanks or commas).  The position
C        of the end of the first token is also returned so that this
C        routine may be conveniently used within a loop to process an
C        entire line of text.
C
C           The procedure examines a substring, STRING (FIRST : LAST), which
C        may of course be the entire string (in which case just call SCAN2
C        with FIRST = 1 and LAST = LEN (STRING) ).  The indices returned
C        are relative to STRING itself, not the substring.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING    *         C    I      Text string containing data to be
C                                        scanned.
C        SEPS      *         C    I      String containing the separators.
C                                        Each character in SEPS counts as a
C                                        token delimiter.
C        FIRST               I    I/O    Index of first character of interest
C                                        in STRING.  FIRST >= 1 is assumed.
C                                        Output is index of the beginning of
C                                        the first non-separator, or 0 if no
C                                        token was found.
C        LAST                I    I/O    Index of last character of interest
C                                        in STRING.  LAST <= LEN (STRING) is
C                                        assumed.  Output is index of the end
C                                        of the last non-separator, or 0 if no
C                                        token was found.
C        MARK                I      O    Points to the end of the first token
C                                        found in the specified portion of
C                                        STRING.  Use FIRST = MARK + 2 for
C                                        further searches in the same string.
C                                        MARK is set to 0 if no token was found.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  FIRST >= 1 and LAST <= LEN (STRING) are assumed upon entry,
C             because these are the obvious bounds needed by the calling
C             program on the first scan of a given string - no need to
C             evaluate LEN (STRING) here with every call, while use of
C             MIN and MAX tends to disguise programming errors.
C
C        (3)  Entering with FIRST > LAST is an eventual normal consequence
C             of setting FIRST = MARK + 2 for further scans of the same string.
C             The zero DO loop iteration count feature of the language is
C             taken advantage of here - the result is simply to drop
C             out the bottom with FIRST = LAST = MARK = 0.
C
C        (4)  This version adjusts LAST using a BACKward search, which is
C             degenerate for all tokens after the first.  (Doing it in
C             the one forward loop means unnecessary repeated tokenizing
C             to find the end.)
C
C        (5)  Zeroing all of FIRST, LAST, and MARK (rather than just MARK)
C             when no token is found may, in retrospect, be undesirable in
C             that information is lost (esp. LAST), but there are too many
C             applications to change this now.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C         4 Mar 1986    RAK    Variation of SCANNR, which is hard-coded
C                              (for historical reasons and a modest speed
C                              advantage) for separators BLANK, TAB,
C                              COMMA, COLON, and EQUAL.
C
C         5 May  1988   DAS    Reverse search used to find LAST; MAX, MIN,
C                              and LEN also eliminated (see Notes).
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   FIRST, LAST, MARK
      CHARACTER
     >   SEPS * (*), STRING * (*)

C     Local variables.

      INTEGER
     >   HEAD, I, J, NSEPS, TAIL
      LOGICAL
     >   FOUND


C     Execution.
C     ----------

      NSEPS = LEN (SEPS)
      HEAD = FIRST
      TAIL = LAST

      FIRST = 0
      LAST = 0
      MARK = 0
      FOUND = .FALSE.

C     Look at each character in STRING (HEAD : TAIL) from left to right
C     until a token is found.  (Drop through if LAST > FIRST.)

      DO 30, I = HEAD, TAIL

C        Is the character a separator?

         DO 10, J = 1, NSEPS
            IF (STRING (I : I) .EQ. SEPS (J : J)) GO TO 20
   10    CONTINUE

C           Not a separator.  Check for the beginning of the FIRST token.

            IF (.NOT. FOUND) THEN
               FIRST = I
               FOUND = .TRUE.
            END IF
            GO TO 30

   20    CONTINUE

C           We found a separator.

            IF (FOUND) THEN

C              We just passed the "trailing edge" of the first token.

               MARK = I - 1
               GO TO 40
            END IF
   30 CONTINUE

C     We reached the last character while still seeking the first token.
C     Either this is part of the first token, or MARK and LAST are still zero.

      IF (FOUND) THEN
         MARK = TAIL
         LAST = TAIL
      END IF
      GO TO 99


   40 CONTINUE

C     We found the first token but haven't reached the end yet to adjust LAST.

      DO 60, I = TAIL, MARK + 1, -1

C        Keep searching as long as we have a separator.

         DO 50, J = 1, NSEPS
            IF (STRING (I : I) .EQ. SEPS (J : J)) GO TO 60
   50    CONTINUE

C           Not a separator, so we've found LAST.

            LAST = I
            GO TO 99
   60 CONTINUE

C     Dropped through, so there were no more tokens.

      LAST = MARK
     
C     Termination.
C     ------------

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE SCANNR (STRING, FIRST, LAST, MARK)
C
C
C     Description and usage:
C
C           Looks for significant fields ("tokens") in a string, where the
C        fields are of arbitrary length and separated by blanks, tabs, commas,
C        colons, or equal signs.  The position of the end of the first token
C        is also returned so that this routine may be conveniently used within
C        a loop to process an entire line of text.
C
C           The procedure examines a substring, STRING (FIRST : LAST), which
C        may of course be the entire string (in which case just call SCANNR
C        with FIRST = 1 and LAST = LEN (STRING) ).  The indices returned
C        are relative to STRING itself, not the substring.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING    *         C    I      Text string containing data to be
C                                        scanned.
C        FIRST               I    I/O    Index of first character of interest
C                                        in STRING.  FIRST >= 1 is assumed.
C                                        Output is index of the beginning of
C                                        the first token, or 0 if no token
C                                        was found.
C        LAST                I    I/O    Index of last character of interest
C                                        in STRING.  LAST <= LEN (STRING) is
C                                        assumed.  Output is index of the end
C                                        of the last non-separator, or 0 if no
C                                        token was found.
C        MARK                I      O    Points to the end of the first token
C                                        found in the specified portion of
C                                        STRING.  Use FIRST = MARK + 2 for
C                                        further searches in the same string.
C                                        MARK is set to 0 if no token was found.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.  Constant HT (Tab) is defined
C             in a non-standard way:  the CHAR function is not permitted
C             in a PARAMETER declaration (OK on VAX, though).  For Absoft
C             FORTRAN 77 on 68000 machines, use HT = 9.  In other cases, it
C             may be best to declare HT as a variable and assign HT = CHAR (9).
C
C        (2)  The pseudo-recursive structure of the original version has been
C             abandoned because the VAX compiler treated the SOLID statement
C             function as an in-line subroutine, with substantial penalty
C             in speed (factor of almost 3!).  The single-loop form used
C             later was almost as fast (especially for lines with only a few
C             tokens), and was more compact and easy to change if a different
C             set of delimiters was required.  However, ...
C
C        (3)  This version adjusts LAST using a BACKward search, which is
C             degenerate for all tokens after the first.  Repeated forward
C             scanning to locate LAST (in a calling program's loop over tokens)
C             amounts to an O(4N**2) operation count for N tokens in a string
C             versus the O(2N) that it should be.  (The 2 is in there because
C             the non-token fields are as numerous as the tokens and take a
C             similar time to scan.)  The price paid is some repeated code.
C
C        (4)  FIRST >= 1 and LAST <= LEN (STRING) are assumed upon entry,
C             because these are the obvious bounds needed by the calling
C             program on the first scan of a given string - no need to
C             evaluate LEN (STRING) here with every call, while use of
C             MIN and MAX tends to disguise programming errors.
C
C        (5)  Entering with FIRST > LAST is an eventual normal consequence
C             of setting FIRST = MARK + 2 for further scans of the same string.
C             The zero DO loop iteration count feature of the language is
C             taken advantage of here - the result is simply to drop
C             out the bottom with FIRST = LAST = MARK = 0.
C
C        (6)  Zeroing all of FIRST, LAST, and MARK (rather than just MARK)
C             when no token is found may, in retrospect, be undesirable in
C             that information is lost (esp. LAST), but there are too many
C             applications to change this now.
C
C        (7)  The variety of separators recognized limits the usefulness of
C             this routine somewhat.  The intent is to facilitate handling
C             such tokens as keywords or numerical values.  In other
C             applications, it might be necessary for ALL printing characters
C             to be significant.  Use SCAN2 (from the same author) if user-
C             supplied delimiters are appropriate.
C
C        (8)  Note that "null" characters are not treated here.  This should
C             not be a problem in that FORTRAN READs of short records, like
C             assignment of short strings, pads with blanks.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C        29 Dec. 1984    RAK    Initial design and coding, (very) loosely
C                               based on SCAN_STRING by Ralph Carmichael.
C        25 Feb. 1984    RAK    Added ':' and '=' to list of separators.
C        16 Apr. 1985    RAK    Defined SOLID in terms of variable DUMMY
C                               (previous re-use of STRING was ambiguous).
C         3 Mar. 1986    RAK    Restructured, without SOLID.  Checks for
C                               "state transitions" while executing a
C                               single forward DO loop.  Protect against
C                               funny input (FIRST > LAST).
C        13 May  1988    DAS    Introduced backward search for LAST;
C                               eliminated MIN, MAX and LEN.  (See Notes.)
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   FIRST, LAST, MARK
      CHARACTER
     >   STRING * (*)

C     Local constants.

      CHARACTER
     >   BLANK, EQUAL, COLON, COMMA, HT
      PARAMETER
     >   (BLANK = ' ',
     >    EQUAL = '=',
     >    COLON = ':',
     >    COMMA = ',',
     >    HT = CHAR (9))

C     Local variables.

      INTEGER
     >   HEAD, I, J, TAIL
      LOGICAL
     >   FOUND


C     Execution.
C     ----------

      HEAD = FIRST
      TAIL = LAST

      FIRST = 0
      LAST = 0
      MARK = 0
      FOUND = .FALSE.

C     Look at each character in STRING (HEAD : TAIL) from left to right
C     until a token is found.  (Drop through if LAST > FIRST.)

      DO 30, I = HEAD, TAIL

C        Is the character a separator?

         IF (STRING (I : I) .EQ. BLANK) GO TO 20
         IF (STRING (I : I) .EQ. HT   ) GO TO 20
         IF (STRING (I : I) .EQ. COMMA) GO TO 20
         IF (STRING (I : I) .EQ. EQUAL) GO TO 20
         IF (STRING (I : I) .EQ. COLON) GO TO 20

C           Not a separator.  Check for the beginning of the FIRST token.

            IF (.NOT. FOUND) THEN
               FIRST = I
               FOUND = .TRUE.
            END IF
            GO TO 30

   20    CONTINUE

C           We found a separator.

            IF (FOUND) THEN

C              We just passed the "trailing edge" of the first token.

               MARK = I - 1
               GO TO 40
            END IF
   30 CONTINUE

C     We reached the last character while still seeking the first token.
C     Either this is part of the first token, or MARK and LAST are still zero.

      IF (FOUND) THEN
         MARK = TAIL
         LAST = TAIL
      END IF
      GO TO 99


   40 CONTINUE

C     We found the first token but haven't reached the end yet to adjust LAST.

      DO 60, I = TAIL, MARK + 1, -1

C        Keep searching backwards as long as we have a separator.

         IF (STRING (I : I) .EQ. BLANK) GO TO 60
         IF (STRING (I : I) .EQ. HT   ) GO TO 60
         IF (STRING (I : I) .EQ. COMMA) GO TO 60
         IF (STRING (I : I) .EQ. EQUAL) GO TO 60
         IF (STRING (I : I) .EQ. COLON) GO TO 60

C           Not a separator, so we've found LAST.

            LAST = I
            GO TO 99
   60 CONTINUE

C     Dropped through, so there were no more tokens.

      LAST = MARK
     
C     Termination.
C     ------------

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE SELECT (PROMPT, NLINES, MENU, SUPRES, LUNCRT, LUNKBD,
     >                   SELNUM, SELNAM, QUIT)
C
C     Purpose:
C
C           SELECT modularizes the situation of choosing from a menu by
C        either NUMBER or NAME.  Its innovative (?) feature lies in its
C        parsing of the supplied menu text to identify valid responses,
C        thereby avoiding an additional dictionary as an argument. This
C        makes use of SELECT easy, at the expense of some complexity here.
C
C           If choice by NUMBER ALONE is adequate, it is probably better
C        to display the menu in-line and just prompt with, say, READER.
C        If device-dependent cursor control utilities are available and
C        acceptable from a portability point of view, a much different
C        approach could be taken.  This approach is as system-independent
C        as possible.
C
C           SELECT displays the menu unless the calling program suppresses
C        it, in which case the program user can still request the menu to
C        be shown if necessary.  A default selection must be defined by
C        the input values of SELNUM and SELNAM.  Valid keyboard entries
C        include all reasonable tokens or part-tokens from the menu text,
C        so long as they define a selection uniquely.  An ambiguous or
C        unidentified response leads to a diagnostic and a reprompt.
C
C           A numeric input is accepted right away if it is an in-range
C        item number, although the first tokens of each menu line are still
C        checked for possible duplication by a bad application (as well as
C        to identify the menu line number corresponding to the input).
C        (Descriptive menu text may contain numerical tokens which could
C        clash with the intended item number, so the normal search of all
C        tokens of every line is inappropriate for a numeric input.)
C
C           Two examples should clarify usage.
C
C           This version shows the user prompt and the default selection
C        in inverse video.  It also provides for having no default, via
C        a blank SELNAM on input, meaning both CR and CTRL Z mean quit
C        (handy for some applications).
C
C     Sample PROMPT:
C
C           'Select operating mode.'
C
C        The prompt seen on the screen will then be
C
C           Select operating mode.
C           ? = menu; CTRL Z = quit; <CR> = DISPLAY:
C
C     Corresponding sample menu text:
C
C           '0 = DISPLAY only (plot and/or tabulate - no changes)',
C           '1 = RECTIFY leading edge definition',
C           '2 = NORMALIZE or denormalize coordinates',
C           '3 = REDISTRIBUTE the abscissas',
C            :   :   :   :   :   :   :   :   :
C
C           (The calling program has this text in a DATA statement.)
C
C           Note that an item NUMBER is expected to be the FIRST token in
C        each menu line, and the corresponding NAME is taken as the SECOND
C        token, which should be unique to that item.  Tokens are separated
C        by any of the delimiters in SEPS below.  These include '.' so that
C        a menu could look like the following example, and so that decimal
C        points in the reply are ignored as probably unintentional.
C
C                Input MENU                       Output SELNUM, SELNAM
C
C           '   -1.  Return to previous menu.',              -1  'RETURN'
C           '    0.  Proceed.',                               0  'PROCEED'
C           '    1.  Display graphics parameters.',           1  'DISPLAY'
C                :   :   :   :   :   :   :   :                :   :   :
C
C        (Use of SELNAM in the calling program may not be appropriate with
C        this example because the mnemonics have not been well chosen.)
C        Other possible forms of item number are 1... and (1).
C
C           SYNONYMS may also be selected from the item descriptions as long
C        as they are unique across items.  An example of how this could be
C        handy is in item 2 above: 'denorm' could be a perfectly valid entry,
C        which would be interpreted as a request for operating mode 2.
C
C           Whatever valid entry (name or number) is identified, the returned
C        value of SELNAM will always be the corresponding SECOND token on the
C        line, fully expanded, in upper case, and left-justified. SELNAM should
C        therefore allow for the appropriate number of characters.
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        PROMPT      -        C     I    Prompt to appear first on the
C                                        screen (to which additional
C                                        standard text will be appended -
C                                        see example above).
C        NLINES      -        I     I    Number of lines of menu text.
C        MENU     NLINES      C     I    Menu, with each line in the form
C                                        shown by the above example - an
C                                        integer for token 1; a key name
C                                        for token 2; and possible synonym
C                                        names among further tokens.  Avoid
C                                        'HELP' for token 2 - it will be
C                                        interpreted as a request for the menu.
C        SUPRES      -        L     I    SUPRES=T means don't display the
C                                        menu (unless the user asks for it).
C        LUNCRT      -        I     I    Logical unit number for the screen.
C        LUNKBD      -        I     I    Logical unit number for the keyboard.
C        SELNUM      -        I    I/O   Input with default item number.
C                                        Output with item number selected from
C                                        the FIRST tokens on the menu lines.
C                                        SELNUM is NOT NECESSARILY within the
C                                        range [1:NLINES]. In the above example,
C                                        SELNUM = 0 is perfectly valid.
C        SELNAM      -        C    I/O   Input with default item name (which
C                                        should match the SECOND token of the
C                                        item referred to by SELNUM on input).
C                                        However a blank input (no default)
C                                        means both CR and CTRL Z will be
C                                        interpreted as meaning "QUIT."
C                                        Output with key name selected from
C                                        the SECOND tokens on the menu lines,
C                                        guaranteed upper case/left-justified.
C                                        May be preferred in calling program
C                                        over SELNUM, for mnemonic reasons.
C                                        Required length of SELNAM should not
C                                        exceed the parameter constant MXNAME
C                                        hard-coded below.
C        QUIT        -        L     O    QUIT = .TRUE. means Control-Z (EOF)
C                                        was entered and other outputs are
C                                        probably undefined.  See also
C                                        "SELNAM = blank" description.
C
C     Implementation notes:
C
C     1.    It would be convenient for the user and the programmer here not
C        to have to return a name as well as a number, because then the key
C        name would not have be token number TWO precisely, whereas if an
C        item number is entered and a corresponding name IS to be returned,
C        it must be in a known place, chosen as the second token.  The only
C        rationale is that the calling program may prefer to be mnemonic in
C        its use of the returned NAME rather than the returned NUMBER.  Here,
C        we provide for that choice, in the spirit of the lower-level LOOKUP
C        utility.
C
C     2.    Errors are handled in two main ways: they are either programmer
C        errors or user input errors.  If the menu text is found to break
C        the rules of the game, the application programmer has erred, and
C        execution halts here with a diagnostic.  If the selection from the
C        menu is too short to be unique, or is out of range or is otherwise
C        invalid, another selection is solicited as part of this routine's
C        purpose of ensuring meaningful inputs.
C
C           A third class of error consists of those that could conceivably
C        happen (such as two menu lines with the same leading integer, or
C        blank menu lines) but are more trouble than it's worth to handle.
C        Hooks are left in these cases, just in case someone has the urge
C        to make the routine more bullet-proof than it really needs to be.
C
C     3.    The SCAN2 version of SCANNR turned out to be the appropriate
C        utility here. TOKENS and LOOKUP could probably have been applied,
C        but the exact match for SELNUM, and the desired delimiters, would
C        cause difficulties.  Also, it was nice to avoid the local storage
C        that TOKENS would need to tokenize one menu line at a time.
C
C     Procedures:
C
C        NUMBER  Logical function for deciding if a string represents a number.
C        SCAN2   Finds first/last characters of next token; variable delimiters.
C        UPCASE  Converts text to upper case.
C
C     Environment:
C
C        VAX/VMS and IRIS 4D; FORTRAN 77 with some extensions:
C           >  IMPLICIT NONE
C           >  (A, $) carriage control
C           >  Escape sequences for inverse video
C
C     History:
C
C        09/24/86  DAS  Initial design and coding, for program PROFILE,
C                       using ideas from LOOKUP, etc., by R.A.Kennelly.
C        09/27/86  DAS  Switched from SCANNR to SCAN2; simplified defaults.
C        10/09/86  DAS  Added SUPRES argument for more flexibility, and
C                       inverse video at slight penalty to portability.
C        02/20/87  DAS  Provided for calling with SELNAM = "blank", as
C                       it is not always appropriate to have a default.
C                       In this case, both CR and CTRL Z mean "quit."
C        06/10/88  DAS  Use LAST = LEN (SELNAM), not MXNAME, since SCAN2
C                       was modified.
C        08/12/89  DAS  SCAN2 call on SELNAM uses SEPS now, not BLANK,
C                       to save scanning a default menu line externally.
C        02/07/90  DAS  Removed concatenation from I/O lists, and went
C                       to (A,$) in place of ('$',A) (both for IRIS 4D).
C                       Also moved escape sequences for inverse video
C                       out of parameter constants, where CHAR (27) and
C                       concatenation probably aren't portable.
C        03/19/97  DAS  The check for ambiguous inputs when a numeric
C                       selection is made can encounter unintended
C                       matches within the body of the menu text.
C                       Therefore, only FIRST tokens on each line should
C                       be checked for validity and ambiguity if the
C                       response is numeric.  (Is that clear?)
C
C     Author:
C
C        David Saunders, Sterling Software/NASA Ames, Mountain View, CA.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NLINES, LUNCRT, LUNKBD, SELNUM
      LOGICAL
     >   SUPRES, QUIT
      CHARACTER
     >   PROMPT * (*), MENU (NLINES) * (*), SELNAM * (*)

C     Local constants:

      INTEGER
     >   MXNAME
      CHARACTER
     >   BLANK, DFLTXT * 33, QUITXT * 28, SEPS * 7

      PARAMETER
     >  (BLANK   = ' ',
     >   DFLTXT  = ' ? = menu; CTRL Z = quit; <CR> = ',
     >   QUITXT  = ' ? = menu; CTRL Z or <CR> = ',
     >   MXNAME  = 16,
     >   SEPS    = ' =:.()/')

C     Local variables:

      INTEGER
     >   FIRST, I, I1, I2, LAST, LENMNU, LINE, MARK, N, NCHARS,
     >   NHITS, N1, N2, STATUS
      LOGICAL
     >   INTEGER, LISTIT, NODEFLT
      CHARACTER
     >   ESCAPE * 1, HELP * 4, INVERSE * 4, REPLY * (MXNAME), RESET * 4,
     >   TOKEN * (MXNAME)

C     Procedures:

      LOGICAL
     >   NUMBER
      EXTERNAL
     >   NUMBER, SCAN2, UPCASE

C     Storage:

      DATA
     >   HELP /'HELP'/


C     Execution.
C     ----------

C     Set up for inverse video (formerly done in VAX-dependent PARAMETERs).
C     2 formats and 3 I/O lists will need changing if this is removed.

      ESCAPE  = CHAR (27)
      INVERSE = ESCAPE // '[7m'
      RESET   = ESCAPE // '[0m'

C     Delimit the input default value for the selection name (may be blank).

      I1 = 1
      LAST = LEN (SELNAM)
      CALL SCAN2 (SELNAM, SEPS, I1, LAST, I2)

      NODEFLT = LAST .EQ. 0
      LISTIT = .NOT.SUPRES

 300  CONTINUE

C     Display menu only if calling program or program user asks for it.

      IF (LISTIT) THEN
         WRITE (LUNCRT, 1003) BLANK, MENU
      END IF

  310 CONTINUE

C     Issue supplied prompt plus standard default/quit/help info.

      WRITE (LUNCRT, 1005) INVERSE, PROMPT, RESET
 
      IF (NODEFLT) THEN
         WRITE (LUNCRT, 1006) QUITXT, INVERSE, 'quit', RESET, ': '
      ELSE
         WRITE (LUNCRT, 1006) DFLTXT, INVERSE, SELNAM (I1 : I2), RESET,
     >      ': '
      END IF

C     Accept the response, which should not exceed MXNAME characters.

      REPLY = BLANK
      READ (LUNKBD, 1001, IOSTAT=STATUS) REPLY

      IF (STATUS .GT. 0) THEN
         WRITE (LUNCRT, 1002) 'Input error -try again.'
         GO TO 310
      END IF

C     Check for CTRL Z or carriage return.

      IF (STATUS .LT. 0) THEN
         QUIT = .TRUE.
         GO TO 999
      END IF

      IF (REPLY .EQ. BLANK) THEN
         QUIT = NODEFLT
         GO TO 999
      END IF

      QUIT = .FALSE.

C     Response must be an item number, an item name, or a request to
C     display the menu, or it could be an invalid input.
C     Tokenize the response, and convert it to upper case, ignoring any
C     decimal point in what is presumably meant to be an item number.

      N1 = 1
      LAST = MXNAME
      CALL SCAN2 (REPLY, ' .', N1, LAST, N2)
      CALL UPCASE (REPLY (N1 : N2))
      N = N2 - N1 + 1

      IF (REPLY (N1 : N2) .EQ. '?' .OR.
     >    REPLY (N1 : N2) .EQ. HELP (1 : N)) THEN

C        (Some applications may find 'HELP' here clashes with an item name.
C        Just allowing for '?' here would be the way to go then.)

C        Display the menu, then reprompt.

         LISTIT = .TRUE.
         GO TO 300
      END IF

C     Try to match the reply with some token in the menu.  If the reply
C     is an item NUMBER, the match must be EXACT ('1' must not be claimed
C     to be a short form of '10', for instance).
C
C     The strategy is to isolate each token in the menu text and compare
C     its first N characters with the N significant characters of the
C     response string.  This has to be done in upper case, and since we
C     don't want to change the menu text in-place, we have to transfer
C     each token to a buffer of length MXNAME and UPCASE that.
C
C     Then ambiguities must be dealt with by searching the rest of the
C     menu for a second hit, in which case another response is solicited.

      INTEGER = NUMBER (REPLY (N1 : N2))
      LENMNU = LEN (MENU (1))
      NHITS = 0

C     For each line in the menu ...

      DO 600, I = 1, NLINES
         FIRST = 1
         LAST = LENMNU

C        For each token in the line (or just the FIRST if numeric) ...

  500    CONTINUE
            CALL SCAN2 (MENU (I), SEPS, FIRST, LAST, MARK)

            IF (MARK .GT. 0) THEN
C              Ignore tokens shorter than the reply.

               NCHARS = MARK - FIRST + 1
               IF (NCHARS .GE. N) THEN
                  IF (.NOT. INTEGER .OR.
     >               (INTEGER .AND. NCHARS .EQ. N)) THEN

                     TOKEN (1 : N) = MENU (I) (FIRST : MARK)
                     CALL UPCASE (TOKEN (1 : N))

                     IF (REPLY (N1 : N2) .EQ. TOKEN (1 : N)) THEN
                        NHITS = NHITS + 1
                        IF (NHITS .EQ. 1) THEN
C                          Ignore rest of this line, but go on to next line.
                           LINE = I
                           GO TO 600
                        ELSE
C                          Two hits:
                           WRITE (LUNCRT, 1002)
     >                        'Ambiguous selection - could be:',
     >                        MENU (LINE),
     >                        'or:',
     >                        MENU (I),
     >                        BLANK,
     >                        'Try again.'
                           GO TO 310
                        END IF
                     END IF
                  END IF
               END IF

               IF (.NOT. INTEGER) THEN
                  FIRST = MARK + 2
                  GO TO 500
               ELSE
                  ! Scan only FIRST tokens for numeric responses
               END IF

            END IF
  600 CONTINUE

      IF (NHITS .EQ. 0) THEN
         WRITE (LUNCRT, 1002) 'No such choice - try again.'
         GO TO 310
      END IF

C     Return valid choice in fully expanded/upper case form of tokens 1 and 2.

      FIRST = 1
      LAST = LENMNU
      CALL SCAN2 (MENU (LINE), SEPS, FIRST, LAST, MARK)
      READ (MENU (LINE) (FIRST : MARK), 1004) SELNUM

      FIRST = MARK + 2
      CALL SCAN2 (MENU (LINE), SEPS, FIRST, LAST, MARK)
      CALL UPCASE (MENU (LINE) (FIRST : MARK))
      SELNAM = MENU (LINE) (FIRST : MARK)


C     Termination.
C     ------------

  999 RETURN

C     Formats.
C     --------

 1001 FORMAT (A)
 1002 FORMAT (1X, A)
 1003 FORMAT (4X, A)
 1004 FORMAT (BN, I4)
 1005 FORMAT (/, 1X, 3A)
 1006 FORMAT (5A, $)

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE SOLVE(NDIM, N, A, B, IPVT)
C
      INTEGER NDIM, N, IPVT(N)
      REAL    A(NDIM,N),B(N)
C
C   SOLUTION OF LINEAR SYSTEM, A*X = B .
C   DO NOT USE IF DECOMP HAS DETECTED SINGULARITY.
C
C   INPUT..
C
C     NDIM = DECLARED ROW DIMENSION OF ARRAY CONTAINING A .
C
C     N = ORDER OF MATRIX.
C
C     A = TRIANGULARIZED MATRIX OBTAINED FROM DECOMP .
C
C     B = RIGHT HAND SIDE VECTOR.
C
C     IPVT = PIVOT VECTOR OBTAINED FROM DECOMP .
C
C   OUTPUT..
C
C     B = SOLUTION VECTOR, X .
C
C-----------------------------------------------------------------------
C
      INTEGER KB, KM1, NM1, KP1, I, K, M
      REAL    T
C
C     FORWARD ELIMINATION
C
      IF (N .EQ. 1) GO TO 50
      NM1 = N-1
      DO 20 K = 1, NM1
         KP1 = K+1
         M = IPVT(K)
         T = B(M)
         B(M) = B(K)
         B(K) = T
         DO 10 I = KP1, N
             B(I) = B(I) + A(I,K)*T
   10    CONTINUE
   20 CONTINUE
C
C     BACK SUBSTITUTION
C
      DO 40 KB = 1,NM1
         KM1 = N-KB
         K = KM1+1
         B(K) = B(K)/A(K,K)
         T = -B(K)
         DO 30 I = 1, KM1
             B(I) = B(I) + A(I,K)*T
   30    CONTINUE
   40 CONTINUE
   50 B(1) = B(1)/A(1,1)
      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION TABLE1 (NTABLE, XTABLE, YTABLE, INDEX, X, IER)
C
C ACRONYM: TABLE look-up, 1 dimension (linear inter-/extra-polation)
C          -----          -
C PURPOSE: TABLE1 returns a linearly-interpolated value from a given
C          1-D table for the given abscissa. It also returns the in-
C          dex of the table abscissa nearest to the given  abscissa,
C          on the low side, as this can speed subsequent look-ups in
C          some contexts.    The table's abscissas are assumed to be
C          monotonic increasing or monotonic decreasing.
C
C          This version permits extrapolation.
C
C          NOTE:  The subroutine form of TABLE1, called LINE1D,  may
C          be more appropriate if an array of target Xs is involved.
C
C METHOD:  Nothing more than linear interpolation  (2-point formula)
C          is attempted, because multi-point formulas demand special
C          treatment at the ends of the table, and because it is de-
C          sired to permit "safe" extrapolation.  (This is often due
C          to round-off or is otherwise intended, and is unlikely to
C          cause grief.)
C
C ARGUMENTS:
C    ARG    DIM   TYPE I/O/S DESCRIPTION
C  NTABLE    -      I    I   Length of the table; NTABLE > 1 assumed
C  XTABLE  NTABLE   R    I   Abscissas in table, assumed monotonic -
C                            increasing or decreasing
C  YTABLE  NTABLE   R    I   Ordinates in table
C  INDEX     -      I   I/O  Input as an estimate of the interval in
C                            the table which contains the target X -
C                            normally INDEX=1 on the first call in a
C                            sequence, and left as set by TABLE1 for
C                            subsequent calls in the sequence;
C                            output points to the element before X
C                            (or exactly at X) in the table;
C                            if extrapolation occurred to the LEFT
C                            of the table, output INDEX=1;
C                            if extrapolation occurred to the RIGHT,
C                            output INDEX=NTABLE-1 (not NTABLE),
C                            meaning the last two table values were
C                            used to do the extrapolation.
C    X       -      R    I   Abscissa at which interpolation is reqd.
C   IER      -      I    O   Error return code - see ERROR HANDLING.
C                            (IER=0 means no error detected.)
C  TABLE1    -      R    O   Interpolated or extrapolated value.
C                            TABLE1=0.0 if an error is detected (rather
C                            than being left undefined).
C EXTERNAL REFERENCES:
C    NAME     DESCRIPTION
C    INTERVAL Interpolation search for interval containing a point.
C
C ERROR HANDLING:
C          Once extrapolation was permitted, it was tempting to remove
C          the error return argument altogether.  However, the follow-
C          ing input checks may save programmers some grief:
C             IER = 0 if no error was detected;
C                 = 1 if NTABLE < 2;
C                 = 2 if INDEX < 1;
C                 = 3 if INDEX > NTABLE-1.
C          There is no check for monotonicity in XTABLE(*).
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   10/21/83   DAS    Initial design and code (sequential search; no
C                     extrapolation permitted).
C   06/06/86   DAS    Extrapolation now permitted.  Adapted so-called
C                     "interpolation search" from CSEVAL (by Robert
C                     Kennelly) to handle monotonic decreasing as well
C                     as increasing case.
C   08/08/86   DAS    TRIAL = MIN (NTABLE-2,...), not MIN (NTABLE-1,...)
C   10/20/87   RAK    Abstract the interval search (with revisions) as
C                     a separate module.
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   INDEX, NTABLE, IER
      REAL
     >   TABLE1, X, XTABLE (NTABLE), YTABLE (NTABLE)

C     Local variables:

      INTEGER
     >   LEFT
      REAL
     >   PM1

C     Execution:

      IER = 0
      TABLE1 = 0.E+0
      LEFT = INDEX

      IF (NTABLE .LT. 2) THEN
         IER = 1
         GO TO 99
      END IF

      IF (LEFT .LT. 1) THEN
         IER = 2
         GO TO 99
      END IF

      IF (LEFT .GT. NTABLE - 1) THEN
         IER = 3
         GO TO 99
      END IF

      PM1 = SIGN (1.E+0, XTABLE (LEFT+1) - XTABLE (LEFT))

C     Search for the best "left-hand" endpoint to use for interpolation.

      CALL INTERVAL (NTABLE, XTABLE, X, PM1, LEFT)

C     Linear interpolation or extrapolation.  Order does not matter.

      TABLE1 = YTABLE (LEFT)  + (YTABLE (LEFT+1) - YTABLE (LEFT)) *
     >   ((X - XTABLE (LEFT)) / (XTABLE (LEFT+1) - XTABLE (LEFT)))
      INDEX = LEFT

   99 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE TOKENS (STRING, NUMBER, LIST)
C
C
C     Description and usage:
C
C           An aid to parsing input data.  The individual "tokens" in a
C        character string are isolated, converted to uppercase, and stored
C        in an array.  Here, a token is a group of significant, contiguous
C        characters.  The following are NON-significant, and hence may
C        serve as separators:  blanks, horizontal tabs, commas, colons,
C        and equal signs.  See SCANNR for details.  (Use TOKEN2 and SCAN2
C        if you need a variable list of delimiters.)
C
C           Processing continues until the requested number of tokens have
C        been found or the end of the input string is reached.
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING              C    I      Input character string to be analyzed.
C        NUMBER              I    I/O    Number of tokens requested (input) and
C                                        found (output).
C        LIST    NUMBER      C      O    Array of tokens, changed to uppercase.
C
C
C     External references:
C
C        Name    Description
C        SCANNR  Finds positions of the first and last significant characters.
C        UPCASE  Converts a string to uppercase.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C        16 Jan. 1984    RAK    Initial design and coding.
C        16 Mar. 1984    RAK    Revised header to reflect full list of
C                               separators, repaired faulty WHILE clause
C                               in "10" loop.
C        18 Sep. 1984    RAK    Change elements of LIST to uppercase one
C                               at a time, leaving STRING unchanged.
C        12 Mar. 1986    RAK    Cross-referenced TOKEN2 variation.
C
C-----------------------------------------------------------------------


C     Variable declarations.
C     ----------------------

      IMPLICIT NONE

C     Parameters.

      CHARACTER
     >   BLANK
      PARAMETER
     >   (BLANK = ' ')

C     Variables.

      INTEGER
     >   COUNT, FIRST, I, LAST, MARK, NUMBER
      CHARACTER
     >   STRING * (*), LIST (NUMBER) * (*)

C     Procedures.

      EXTERNAL
     >   UPCASE, SCANNR


C     Executable statements.
C     ----------------------

C     WHILE there are tokens to find, loop UNTIL enough have been found.

      FIRST = 1
      LAST = LEN (STRING)

      COUNT = 0
   10 CONTINUE

C        Get delimiting indices of next token, if any.

         CALL SCANNR (STRING, FIRST, LAST, MARK)
         IF (LAST .GT. 0) THEN
            COUNT = COUNT + 1

C           Pass token to output string array, then change case.

            LIST (COUNT) = STRING (FIRST : MARK)
            CALL UPCASE (LIST (COUNT))
            FIRST = MARK + 2
            IF (COUNT .LT. NUMBER) GO TO 10

         END IF


C     Fill the rest of LIST with blanks and set NUMBER for output.

      DO 20 I = COUNT + 1, NUMBER
         LIST (I) = BLANK
   20 CONTINUE

      NUMBER = COUNT


C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE TRICPS ( N, D, L, R, S )
C
C ACRONYM: TRIdiagonal system (Cyclic, Positive definite, Symmetric)
C          ---                 -       -                  -
C
C PURPOSE: TRICPS solves one system  A s = r,  where A is a symmetric,
C          diagonally-dominant, irreducible matrix with positive diag-
C          onal elements, non-zero sub-/super-diagonals,  and non-zero
C          lower-left/upper-right elements, and zeros elsewhere. It is
C          suited to problems involving cyclic/periodic boundary  con-
C          ditions.
C
C METHOD:  The LDL' form of the Cholesky factorization of  A  is used.
C          L  in this case is unit lower bidiagonal with its last  row
C          also filled in:
C
C                 x x     x       1         x         1 x     x
C                 x x x           x 1         x         1 x   x
C                   x x x     =     x 1         x         1 x x
C                     x x x           x 1         x         1 x
C                 x     x x       x x x x 1         x         1
C
C          It is assumed that only one right-hand side is involved for
C          the given A,  so no attempt is made to separate the factor-
C          ization from the solution using the triangular factors.
C
C ARGUMENTS:
C   ARG  DIM TYPE I/O/S DESCRIPTION
C    N    -    I    I   Order of the system (>=3)
C    D    N    R   I/S  Main diagonal elements of A
C    L    N    R   I/S  Off-diagonals of A in L(1:N-1); A(N,1) in L(N)
C    R    N    R    I   Right-hand side elements (see next item)
C    S    N    R    O   Solution vector; may be the same argument as R
C                       in which case R is destroyed.  D and L are de-
C                       stroyed in either case.
C
C ERROR HANDLING:  None. A divide by zero can mean the matrix is sing-
C                  ular but is more likely to mean unwitting errors in
C                  the call to TRICPS.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   05/24/84   DAS    Initial design/coding, for a spline application.
C   05/22/86   RAK    Three divides by DI are worth replacing with one
C                     reciprocal and three multiplies in the main loop.
C
C AUTHOR: David Saunders, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER  N
      REAL     D(N), L(N), R(N), S(N)

C ... Local variables:

      INTEGER  I
      REAL     DI, DINV, LIM1

C ... Executable statements:

C ... The off-diagonals are modified in place.   The modified diagonals
C     do not have to be stored for the back-solve if the matrix factors
C     are treated as (LD)L' and not L(DL').  This opens the way for the
C     filled-in row of factor L to overwrite the original diagonals. It
C     is then not possible to include the I = N-1 case in the following
C     loop ( D(N-1) being the exceptional element ). Repeating the code
C     was considered preferable to an IF test inside the loop.

      DINV = 1.E+0 / D(1)
      D(1) = L(N) * DINV
      S(1) = R(1) * DINV
      D(N) = D(N) - L(N) * D(1)
      S(N) = R(N) - L(N) * S(1)

      DO 20 I = 2, N-2
         LIM1   = L(I-1)
         L(I-1) = LIM1 * DINV
         DI     = D(I) - LIM1 * L(I-1)
         DINV   = 1.E+0 / DI
         D(I)   = - LIM1 * D(I-1) * DINV
         S(I)   = ( R(I) - LIM1 * S(I-1) ) * DINV
         D(N)   = D(N) - ( DI * D(I) ) * D(I)
         S(N)   = S(N) - ( DI * D(I) ) * S(I)
   20 CONTINUE

      I      = N-1
      LIM1   = L(I-1)
      L(I-1) = LIM1 * DINV
      DI     = D(I) - LIM1 * L(I-1)
      D(I)   = ( L(I) - LIM1 * D(I-1) ) / DI
      S(I)   = ( R(I) - LIM1 * S(I-1) ) / DI
      D(N)   = D(N) - ( DI * D(I) ) * D(I)
      S(N)   = ( S(N) - ( DI * D(I) ) * S(I) ) / D(N)
      L(I)   = 0.E+0

      DO 30 I = N-1, 1, -1
         S(I) = S(I) - L(I) * S(I+1) - D(I)*S(N)
   30 CONTINUE

      RETURN
      END
C+---------------------------------------------------------------------
C
      SUBROUTINE TRDIAG (A, B, C, R, S, N)
C
C     TRDIAG solves a general N*N tridiagonal system which is
C     assumed to be suitably diagonally dominant...
C
C     A(*) is subdiagonal   (1st element unused here);
C     B(*) is main diagonal;
C     C(*) is superdiagonal (nth element unused here);
C     R(*) is the right-hand-side vector;
C     S(*) is returned with the solution.
C
C     B(*) is overwritten during the solution;
C     S(*) and R(*) may be the same locations, in which case
C     the RHS  R(*) is also overwritten...
C
C     REFERENCE:  Collected algorithms from CACM: Algorithm 24.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N
      REAL
     >   A (N), B (N), C (N), R (N), S (N)

C     Local variables:

      INTEGER
     >   J
      REAL
     >   W

C     Constants:

      REAL
     >   ONE
      PARAMETER
     >  (ONE = 1.E+0)

      W = ONE / B (1)
      S (1) = R (1) * W

      DO J = 1, N - 1
         B (J) = C (J) * W
         W =   ONE / (B (J + 1) - A (J + 1) * B (J))
         S (J + 1) = (R (J + 1) - A (J + 1) * S (J)) * W
      END DO

      DO J = N - 1, 1, -1
         S (J) = S (J) - B (J) * S (J + 1)
      END DO

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE UPCASE( STRING )
C
C PURPOSE:  UPCASE changes all lower case letters in the given
C           character string to upper case.
C
C METHOD:   Each character in STRING is treated in turn.  The intrinsic
C           function INDEX effectively allows a table lookup, with
C           the local strings LOW and UPP acting as two tables.
C           This method avoids the use of CHAR and ICHAR, which appear
C           to behave differently on ASCII and EBCDIC machines.
C
C ARGUMENTS
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C  STRING       *       C   I/O   Character string possibly containing
C                                 some lower-case letters on input;
C                                 strictly upper-case letters on output
C                                 with no change to any non-alphabetic
C                                 characters.
C
C EXTERNAL REFERENCES:
C  LEN    - Returns the declared length of a CHARACTER variable.
C  INDEX  - Returns the position of second string within first.
C
C ENVIRONMENT:  ANSI FORTRAN 77
C
C AUTHOR: Michael Saunders, Systems Optimization Lab., Stanford, 9/10/85
C         (rewrite of version by Hooper/Kennelly, Informatics, 1983)
C
C-----------------------------------------------------------------------

      CHARACTER      STRING * (*)
      INTEGER        I, J
      CHARACTER      C*1, LOW*26, UPP*26
      DATA           LOW /'abcdefghijklmnopqrstuvwxyz'/,
     $               UPP /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      DO 10 J = 1, LEN(STRING)
         C = STRING(J:J)
         IF (C .GE. 'a'  .AND.  C .LE. 'z') THEN
            I = INDEX( LOW, C )
            IF (I .GT. 0) STRING(J:J) = UPP(I:I)
         END IF
   10 CONTINUE

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE USESCALE (MODE, NDIM, NPTS, X, Y, Z, SCALE, SHIFT, IER)
C
C One-liner:  Normalizes or denormalizes 1, 2, or 3 dimensional data.
C
C Purpose:
C   USESCALE may be used to normalize or denormalize 1, 2 or 3D data.
C   Typical use is prior to or after application of numerical methods
C   which may be sensitive to data units.
C
C   USESCALE can be used with GETSCALE, which returns safeguarded scale
C   and shift factors for normalizing one or all of X, Y, Z to [0,1].
C
C Method:
C   Each of X[, Y[, Z]] is transformed in place as indicated below for X:
C
C   Normalize:            X <-- (X - SHIFT) / SCALE
C
C   Denormalize:          X <-- X * SCALE + SHIFT
C
C Arguments:
C    NAME    DIM    TYPE I/O/S DESCRIPTION
C   MODE      -     C*1    I   'N' to normalize or 'D' to denormalize
C   NDIM      -      I     I   Number of dimensions of the data
C   NPTS      -      I     I   Number of points to be transformed
C   X        NPTS    R    I/O  + X,Y,Z coordinates
C   Y         "      "     "   | of the input/output data
C   Z         "      "     "   + (see Notes)
C   SCALE    NDIM    R     I   Scale factor(s)
C   SHIFT     "      "     "   Shift factor(s)
C   IER       -      I     O   Error code = 0 for normal return
C                                         = -1 for bad MODE
C Notes:
C   For 1D and 2D cases, arguments Y and/or Z do not apply.  The calling
C   program should pass dummies, or pass X more than once.  The 1D case
C   is degenerate, but included for completeness.
C
C Environment:
C   VAX/VMS, FORTRAN 77 with:
C   IMPLICIT NONE
C   Names up to 8 characters
C
C Author: Michael Wong, Sterling Software, Palo Alto, CA
C
C History:
C   Oct. 88  D.A.Saunders  Initial design.
C   Nov. 88  M.D.Wong      Initial coding.
C   July 89  D.A.S.        Fixed glitches in header.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      CHARACTER  MODE*1
      INTEGER    NDIM, NPTS, IER
      REAL       X(NPTS), Y(NPTS), Z(NPTS), SCALE(NDIM), SHIFT(NDIM)

C     Local constants:

      REAL       ONE
      PARAMETER (ONE = 1.E+0)

C     Local variables:

      INTEGER    I
      REAL       TSCALE (3), TSHIFT (3)
      LOGICAL    NORMAL

C     Execution:

      NORMAL = MODE .EQ. 'N'

      IF (.NOT. NORMAL) THEN
         IER = -1
         IF (MODE .NE. 'D') GO TO 99
      END IF

      IER = 0

C     This short loop should save code vs. reuse of scalar temporararies...

      DO 10, I = 1, NDIM
         IF (NORMAL) THEN
            TSCALE (I) = ONE / SCALE (I)
            TSHIFT (I) = -SHIFT (I) * TSCALE (I)
         ELSE
            TSCALE (I) = SCALE (I)
            TSHIFT (I) = SHIFT (I)
         END IF
   10 CONTINUE

C     Transform the X coordinates common to all cases.

      DO 20, I = 1, NPTS
         X (I) = X (I) * TSCALE (1) + TSHIFT (1)
   20 CONTINUE

      IF (NDIM .GE. 2) THEN

C        Transform Y coordinates.

         DO 40, I = 1, NPTS
            Y (I) = Y (I) * TSCALE (2) + TSHIFT (2)
   40    CONTINUE

      END IF

      IF (NDIM .EQ. 3) THEN

C        Transform Z coordinates.

         DO 60 I = 1, NPTS
            Z (I) = Z (I) * TSCALE (3) + TSHIFT (3)
   60    CONTINUE

      END IF

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE WAGFIT (NX, X, Y, N, MAXWRK, WRK, C, RMSDEV, IER)
C
C  ACRONYM: WAGner series: FITting of first few terms by least squares
C           ---            ---
C  PURPOSE: WAGFIT  fits a linear combination of the  first  N  Wagner
C           functions in the least squares sense to the given dataset,
C           which is probably an airfoil surface or a  thickness  dis-
C           tribution or a meanline (camber) distribution. It may also
C           be a distribution of airfoil perturbations.   It must span
C           the interval [0.,1.], with Y(1) = Y(NX) = 0.
C
C           This version requires that the "ramp" function  needed  to
C           retain nonzero Y at X=1. be dealt with in the calling pro-
C           gram, for reasons of symmetry:  evaluation of the fit must
C           handle this nonzero case,  so  setting up the fit to allow
C           for it should be done in the same place.
C           
C           (WAGFIT originally included the ramp term  as  an  (N+1)th
C           coefficient to compute, but the preferred result is to re-
C           tain any original nonzero Y(NX) exactly.)
C
C  METHOD:  WAGFIT sets up the particular overdetermined system led to
C           by a Wagner series, one column at a time, then uses a gen-
C           eral linear least squares solver  to compute optimal coef-
C           ficients and a goodness-of-fit measure.
C
C  NOTES:   Wagner functions have found application in airfoil design,
C           as combinations of them can produce  airfoil-like  curves.
C           See subroutine BEVAL for details of the  Wagner functions,
C           and for evaluating the linear combination determined here.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S   DIM    DESCRIPTION
C    NX      I      I      -     Number of data points being fitted
C    X       R      I      NX    Abscissas of data points, in [0.,1.]
C    Y       R      I      NX    Corresponding ordinates: Y(1)=Y(NX)=0
C    N       I      I      -     Wagner functions defined by  1:N  are
C                                to be fitted as a linear combination,
C                                1 <= N <= NX. Values of N > 20 aren't
C                                permitted  (to  protect  over-zealous
C                                users from themselves).
C    MAXWRK  I      I      -     Maximum work-space provided.  Must be
C                                be at least NX*(N+2).
C    WRK     R      S    MAXWRK  Work-space for least squares solver
C    C       R      O      N     Computed coefs. C(I) applies to Wagner
C                                function I for I = 1:N.
C    RMSDEV  R      O      -     SQRT ( Optimal sum of squares / NX )
C    IER     I      O      -     0 means no errors;
C                                4 means no fit (N < 1 or > 20 );
C                                5 means no fit (N > NX);
C                                6 means no fit (MAXWRK too small).
C  EXTERNAL REFERENCES:
C    BEVAL      Evaluates a Wagner function for given N, X(*)
C    HDESOL     Linear least squares by Householder decomposition
C
C  ENVIRONMENT:  VAX/VMS FORTRAN
C
C  HISTORY:
C    03/05/85   DAS   Initial implementation, including "ramp" term.
C    07/17/86   DAS   Nonzero Y(NX) must be handled outside now.
C    04/05/92   DAS   Meaning of IER=5 above had not been updated.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER    N, NX, MAXWRK, IER
      REAL       C (N), RMSDEV, WRK (MAXWRK), X (NX), Y (NX)

C ... Local declarations:

      INTEGER    I, J, MAXN
      REAL       P (2)
      PARAMETER  (MAXN = 20)

C ... Externals:

      INTRINSIC  SQRT
      EXTERNAL   BEVAL, HDESOL

C ... Execution:

C ... Check for the more obvious erroneous calls:

      IF (N .LT. 1  .OR.  N .GT. MAXN) GO TO 940
      IF (N .GT. NX) GO TO 950
      IF (MAXWRK .LT. NX * (N + 2)) GO TO 960

C ... Set up the over-determined (or square) system by columns:

      P (2) = 1.E+0
      I = 1

      DO 200 J = 1, N
         P (1) = J

         CALL BEVAL ('WAGNER', 2, P, .FALSE., NX, X, WRK (I))

         I = I + NX
  200 CONTINUE

      J = I - 1

C ... Now for the right-hand-side:

      DO 300 I = 1, NX
         WRK (I+J) = Y (I)
  300 CONTINUE

      J = NX + J

C ... Solve the overdetermined system by orthogonal factorization.
C     The required coefs. are returned starting in WRK (J+1):

      CALL HDESOL (NX, NX, N+1, WRK (1), WRK (J+1), RMSDEV)

      DO 400 I = 1, N
         C (I) = WRK (J+I)
  400 CONTINUE

      RMSDEV = SQRT (RMSDEV / NX)

      IER = 0
      GO TO 999

C ... Error handling:

 940  CONTINUE
C ... Requested sum has too many terms to be a reasonable mathematical
C     model fitted by these techniques (i.e., N > MAXN), else N < 1:

      IER = 4
      GO TO 999

 950  CONTINUE
C ... N is too big relative to NX for the mathematical model to be fit
C     this way:

      IER = 5
      GO TO 999

 960  CONTINUE
C     Not enough work-space provided:

      IER = 6

 999  RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE XDERIVS (N, XP, XPP, YP, YPP, DYDX, D2YDX2, SCALE)
C
C  PURPOSE:
C
C        XDERIVS modularizes the calculation of Y vs. X derivatives
C     from parametric derivatives for X and Y vs. T.  The formulas are:
C
C           dY/dX    =  Y' / X'
C           d2Y/dX2  =  (Y" X'  -  Y' X") / X' ** 3
C
C     The case of dX/dT ~ 0. is safeguarded with a relative test.
C     (1.E-10 * SCALE is substituted for |X'| less than this value.)
C
C        See CURV2D for the (similar) curvature calculation.
C
C  ARGUMENTS:
C
C        Obvious from the description: N >= 1; all others except SCALE
C     are REAL (N).  SCALE should be the X data range, needed to make
C     the test against division by very small numbers relative.
C
C        Each 2nd derivative is calculated before each 1st derivative,
C     so the calling program may pass the same array for both DYDX and
C     D2YDX2 if only the first derivatives are required.
C
C  HISTORY:
C
C     10/25/91   DAS   Initial implementation, for program PROFILE.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    N
      REAL       XP (N), XPP (N), YP (N), YPP (N), DYDX (N), D2YDX2 (N),
     >           SCALE

C     Local constants:

      REAL       EPS
      PARAMETER (EPS = 1.E-10)

C     Local variables:

      INTEGER    I
      REAL       DXDT, EPSIL

C     Execution:

      EPSIL = EPS * SCALE
      DO 10, I = 1, N
         DXDT = XP (I)
         IF (ABS (DXDT) .LT. EPSIL) DXDT = SIGN (EPSIL, DXDT)
         D2YDX2 (I) = (YPP (I) - YP (I) * XPP (I) / DXDT) / DXDT ** 2
         DYDX (I) = YP (I) / DXDT    ! Avoid cubing DXDT to delay overflow
   10 CONTINUE

      RETURN
      END
C+---------------------------------------------------------------------
C
      SUBROUTINE XGRID (NPTS, MODE, XMIN, XMAX, X)
C
C  ONE-LINER:  Simple 1-D grid point distributions
C
C  PURPOSE:  XGRID generates NPTS abscissas in the range [XMIN, XMAX],
C            distributed according to MODE.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S   DIM     DESCRIPTION
C    NPTS     I     I      -      Desired number of abscissas.
C    MODE     I     I      -      MODE = 0 means equally-spaced points;
C                                      = 1 means sinusoidal bunching
C                                          towards XMIN;
C                                      = 2 means sinusoidal bunching
C                                          towards XMIN and XMAX;
C                                      = 3 means sinusoidal bunching
C                                          towards XMAX.
C  XMIN,XMAX  R     I      -      First and last values to include.
C    X        R     O     NPTS    Desired distribution of abscissas.
C
C  ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C  HISTORY:
C    04/30/82    PJT    Original design and coding.
C    01/20/84    DAS    Added MODE=2 option (for circular arc case).
C    02/03/86    DAS    Changed loops from 1:NPTS to 2:NPTS-1, in case
C                       X(1) and X(NPTS) are passed as XMIN, XMAX.
C    12/02/92    DAS    Added bunching at XMAX option.
C
C  AUTHOR:     Phillip J. Trosin, Informatics General, Palo Alto, CA.
C
C----------------------------------------------------------------------

      DIMENSION X (NPTS)

      PIBY2  = ASIN (1.0)
      DX     = XMAX - XMIN
      DTHETA = PIBY2 / (NPTS - 1)

      IF (MODE .EQ. 0) THEN   ! Uniform

         DX = DX / (NPTS - 1)

         DO 100 I = 2, NPTS - 1
            X (I) = XMIN + DX * (I - 1)
  100    CONTINUE

      ELSE IF (MODE .EQ. 1) THEN  ! Sinusoidal bunching towards XMIN

         DO 200 I = 2, NPTS - 1
            X (I) = XMIN + DX * (1.0 - COS ((I - 1) * DTHETA))
  200    CONTINUE

      ELSE IF (MODE .EQ. 2) THEN  ! Sinusoidal bunching at both ends:

         DTHETA = DTHETA + DTHETA
         DX = DX * 0.5

         DO 300 I = 2, NPTS - 1
            X (I) = XMIN + DX * (1.0 - COS ((I - 1) * DTHETA))
  300    CONTINUE

      ELSE IF (MODE. EQ. 3) THEN  ! Sinusoidal bunching towards XMAX

         DO 400 I = 2, NPTS - 1
            X (I) = XMIN + DX * SIN ((I - 1) * DTHETA)
  400    CONTINUE

      END IF

      X (1) = XMIN
      X (NPTS) = XMAX

      RETURN
      END
C+----------------------------------------------------------------------
C
      REAL FUNCTION ZEROIN ( AX, BX, F, TOL, IER )
      INTEGER       IER
      REAL          AX,BX,F,TOL
      EXTERNAL      F
C
C  PURPOSE:
C     ZEROIN estimates a zero of the function F(X) in the interval AX,BX.
C
C  INPUT:
C
C  AX     Left endpoint of initial interval
C  BX     Right endpoint of initial interval
C  F      FUNCTION subprogram which evaluates F(X) for any X in
C         the interval AX,BX
C  TOL    Desired length of the interval of uncertainty of the
C         final result (TOL >= zero)
C
C  OUTPUT:
C
C  ZEROIN Abcissa approximating a zero of F in the interval AX,BX
C  IER    0 means no error detected;
C         1 means F(AX) and F(BX) have the same sign;
C         2 means AX >= BX.
C
C     ZEROIN originally assumed that F(AX) and F(BX) have opposite signs
C  without a check.  This version introduces the IER argument.
C
C     ZEROIN returns a zero X in the given interval (AX,BX) to within
C  a tolerance 4 * MACHEPS * ABS(X) + TOL, where MACHEPS is the relative
C  machine precision.
C
C     This FUNCTION subprogram is a slightly modified translation of
C  the ALGOL 60 procedure ZERO given in Richard Brent, Algorithms for
C  Minimization Without Derivatives, Prentice-Hall, Inc. (1973).
C
C-----------------------------------------------------------------------
C
      REAL  A,B,C,D,E,EPS,FA,FB,FC,TOL1,XM,P,Q,R,S
      REAL  ABS,SIGN
C
C  COMPUTE EPS, THE RELATIVE MACHINE PRECISION
C
      EPS = 1.0E+0
   10 EPS = EPS/2.0E+0
      TOL1 = 1.0E+0 + EPS
      IF (TOL1 .GT. 1.0E+0) GO TO 10
C
C INITIALIZATION
C
      IER = 2
      A = AX
      B = BX
      IF ( A .GE. B ) GO TO 99
C
      FA = F(A)
      FB = F(B)
C
      IER = 1
      IF ( FA .EQ. 0.E+0 .AND. FB .NE. 0.E+0 ) GO TO 15
      IF ( FB .EQ. 0.E+0 .AND. FA .NE. 0.E+0 ) GO TO 15
      IF ( FA * FB .GE. 0.E+0 ) GO TO 99
C
   15 IER = 0
C
C BEGIN STEP
C
   20 C = A
      FC = FA
      D = B - A
      E = D
   30 IF (ABS(FC) .GE. ABS(FB)) GO TO 40
      A = B
      B = C
      C = A
      FA = FB
      FB = FC
      FC = FA
C
C CONVERGENCE TEST
C
   40 TOL1 = 2.0E+0*EPS*ABS(B) + 0.5E+0*TOL
      XM = 0.5E+0*(C - B)
      IF (ABS(XM) .LE. TOL1) GO TO 90
      IF (FB .EQ. 0.0E+0) GO TO 90
C
C IS BISECTION NECESSARY
C
      IF (ABS(E) .LT. TOL1) GO TO 70
      IF (ABS(FA) .LE. ABS(FB)) GO TO 70
C
C IS QUADRATIC INTERPOLATION POSSIBLE
C
      IF (A .NE. C) GO TO 50
C
C LINEAR INTERPOLATION
C
      S = FB/FA
      P = 2.0E+0*XM*S
      Q = 1.0E+0 - S
      GO TO 60
C
C INVERSE QUADRATIC INTERPOLATION
C
   50 Q = FA/FC
      R = FB/FC
      S = FB/FA
      P = S*(2.0E+0*XM*Q*(Q - R) - (B - A)*(R - 1.0E+0))
      Q = (Q - 1.0E+0)*(R - 1.0E+0)*(S - 1.0E+0)
C
C ADJUST SIGNS
C
   60 IF (P .GT. 0.0E+0) Q = -Q
      P = ABS(P)
C
C IS INTERPOLATION ACCEPTABLE
C
      IF ((2.0E+0*P) .GE. (3.0E+0*XM*Q - ABS(TOL1*Q))) GO TO 70
      IF (P .GE. ABS(0.5E+0*E*Q)) GO TO 70
      E = D
      D = P/Q
      GO TO 80
C
C BISECTION
C
   70 D = XM
      E = D
C
C COMPLETE STEP
C
   80 A = B
      FA = FB
      IF (ABS(D) .GT. TOL1) B = B + D
      IF (ABS(D) .LE. TOL1) B = B + SIGN(TOL1, XM)
      FB = F(B)
      IF ((FB*(FC/ABS(FC))) .GT. 0.0E+0) GO TO 20
      GO TO 30
C
C DONE
C
   90 ZEROIN = B
   99 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE ZERORC (AX, BX, X, F, TOL, NUMFUN, CALLER, LUNOUT,
     >                   HOLD, ISTAT)
C
C     Oneliner:
C
C     ZERO of a function (1-D; no derivatives; Reverse Communication)
C     ----                                     -       -
C
C     Description:
C
C        ZERORC estimates a zero of the function F(X) in the interval
C     (AX, BX).  Function derivatives are not required.
C
C        ZERORC is an adaptation of ZEROIN (see below), intended to
C     avoid the problems of passing as arguments the names of modules
C     with assumed calling sequences (such as ZEROIN's F(X)) at the
C     expense of having to return to the calling program for every
C     function evaluation.  Such reverse communication is usually
C     preferable to the use of common blocks forced on nontrivial
C     function routines by the likes of ZEROIN.
C
C
C     Arguments:
C
C     Name   Type   I/O/S   Description
C
C     AX      R     I       Left endpoint of initial interval.
C
C     BX      R     I       Right endpoint of initial interval.
C
C     X       R       O     X is output with the point at which the
C                           next function evaluation is required, until
C                           termination, when X is the best estimate
C                           of the location of the zero in (AX, BX).
C
C     F       R     I/O     After the initial call, F should be input
C                           with the function value corresponding to X.
C                           Upon termination, F is output with the
C                           function value corresponding to X, which
C                           should be essentially zero.
C
C     TOL     R     I       Desired length of the interval of uncertainty
C                           of the final result.  TOL >= 0.
C
C     NUMFUN  I     I/O     At the start of an iteration, NUMFUN should
C                           be input with the maximum number of function
C                           evaluations to be permitted.  On output, it
C                           contains the number of function evaluations so far.
C
C     CALLER  C*(*) I       Name of the calling program, printed if LUNOUT > 0.
C                           No more than 8 characters are significant.
C
C     LUNOUT  I     I       Logical unit number for output of the iteration
C                           history.  Use LUNOUT < 0 if you don't want the
C                           evaluations printed.  Any diagnostics go to
C                           the unit defined by ABS (LUNOUT).
C
C     HOLD (13)  R  I/O     Holds copies of the local real variables reused
C                           during the calls to ZERORC in an iteration.
C                           This permits finding a zero of a function which
C                           in turn involves finding a zero of another function.
C                           The contents should be manipulated only by ZERORC.
C
C     ISTAT   I     I/O     Status flag used as follows:
C                           ISTAT = +2 on input means this call initiates
C                                      a new iteration.
C                           ISTAT = +1 on return means the calling program
C                                      should evaluate the function at the
C                                      current value of X, and call ZERORC
C                                      again with this value in argument F
C                                      and ISTAT = +1 still.
C                           ISTAT =  0 on output means a zero has been
C                                      found to the specified tolerance.
C                           ISTAT = -1 on output means the limit on the
C                                      number of function evaluations
C                                      has been reached but the convergence
C                                      criteria have not been satisfied.
C                                      X, F are the best values found.
C                           ISTAT = -2 means AX >= BX.  This is checked only
C                                      at the start of an iteration.
C                           ISTAT = -3 means the input value of NUMFUN on
C                                      the first call is either non-positive
C                                      or absurdly large.
C                           ISTAT = -4 means the function has the same sign
C                                      at AX and BX.
C
C     Usage:
C
C        ISTAT = 2          ! Initialize the iteration
C        AX = ...
C        BX = ...
C        NUMFUN = ...       ! Max. no. of fn. evals. allowed
C                           <Define TOL and SUBNAME as well.>
C     10 CONTINUE
C
C           CALL ZERORC (AX, BX, X, F, TOL, NUMFUN, SUBNAME,
C       >                LUNOUT, HOLD, ISTAT)
C
C           IF (ISTAT .LT. -1) THEN      ! Fatal error
C
C              <Handle it>
C
C           ELSE IF (ISTAT .LT. 0) THEN  ! Iteration limit reached
C
C              <Unconverged, but X may still be usable>
C
C           ELSE IF (ISTAT .GT. 0) THEN  ! Evaluate the function
C
C              CALL fun (X, F)           ! Or whatever
C              GO TO 10
C
C           ! ELSE ISTAT = 0 - a zero has been found.
C
C           END IF
C
C 
C     Notes on the ZEROIN algorithm:
C
C       ZEROIN determines a zero X in the given interval (AX, BX) to within
C     a tolerance 4 * MACHEPS * ABS(X) + TOL, where MACHEPS is the relative
C     machine precision.
C
C       ZEROIN is discussed by Forsythe, Malcolm, and Moler in Computer
C     Methods for Mathematical Computations, Prentice-Hall, Inc. (1977).
C     It is a slightly modified translation of the ALGOL 60 procedure
C     ZERO given in Richard Brent, Algorithms for Minimization Without
C     Derivatives, Prentice-Hall, Inc. (1973).
C
C     History:
C
C     c. 1973  R. Brent            ALGOL 60 procedure ZERO.
C     c. 1977  M. Malcolm, et al.  ZEROIN version (FORTRAN 66).
C
C     04/11/92 D. Saunders         ZERORC version, to avoid the need for
C     Sterling Software/NASA Ames  common blocks for nontrivial functions.
C
C     10/29/10 D. Saunders         Matt MacLean reported that the machine
C              ERC, Inc./NASA ARC  epsilon iteration produced zero (!) with
C                                  the Sun/Oracle compiler switch for faster
C                                  but less accurate floating point arithmetic.
C                                  Therefore, use Fortran 90's intrinsic now.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   ISTAT, LUNOUT, NUMFUN
      REAL
     >   AX, BX, F, TOL, X, HOLD (13)
      CHARACTER
     >   CALLER * (*)

C     Local constants:

      REAL, PARAMETER ::
     >   HALF = 0.5E+0, ONE = 1.E+0, TWO = 2.E+0, THREE = 3.E+0,
     >   ZERO = 0.0E+0

C     Local variables:

      INTEGER
     >   LUNERR, MAXFUN
      REAL
     >   A, B, C, D, E, EPS, FA, FB, FC, TOL1, XM, P, Q, R, S
      LOGICAL
     >   BEGIN

C     System functions:

      INTRINSIC
     >   ABS, EPSILON, SIGN

      SAVE
     >   LUNERR, MAXFUN, EPS   ! Assume these are the same for any
                               ! iterations-within-iterations

C     Execution:

      IF (ISTAT .GT. 1) THEN   ! Initialize the iteration

         EPS = EPSILON (EPS)

C        Check input arguments:

         LUNERR = ABS (LUNOUT)
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1010) CALLER

         IF (AX .GE. BX) THEN
            ISTAT = -2
            WRITE (LUNERR, 1030) AX, BX, CALLER

         ELSE IF (NUMFUN .LE. 1 .OR. NUMFUN .GT. 200) THEN
            ISTAT = -3
            WRITE (LUNERR, 1040) NUMFUN, CALLER

         ELSE
            ISTAT = 1
            MAXFUN = NUMFUN
            NUMFUN = 0
            A = AX
            B = BX
            X = A
            HOLD (1) = A
            HOLD (2) = B

         END IF

         GO TO 99  ! Return early for the first function evaluation

      ELSE IF (NUMFUN .EQ. 0) THEN

C        Second call to ZERORC.  We have F(AX) but need F(BX).

         NUMFUN = 1
         FA = F
         HOLD (3) = FA
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020) NUMFUN, X, FA, CALLER

         X = BX
         GO TO 99

      ELSE IF (NUMFUN .EQ. 1) THEN
 
C        We now have F at both end points.
C        These function values must have opposite signs, or at least
C        one must be nonzero.

         FB = F
!!!      HOLD (4) = FB   ! Redundant
         FA = HOLD (3)

         IF (FA .EQ. ZERO .AND. FB .EQ. ZERO) THEN
            ISTAT = -4
         ELSE IF (FA * FB .GT. ZERO) THEN
            ISTAT = -4
         END IF

         IF (ISTAT .NE. 1) THEN
            WRITE (LUNERR, 1050) FA, FB, CALLER
            GO TO 99
         END IF

      END IF

      A  = HOLD (1)
      B  = HOLD (2)
      FA = HOLD (3)

      NUMFUN = NUMFUN + 1
      BEGIN  = NUMFUN .EQ. 2

      IF (.NOT. BEGIN) THEN   ! Restore local variables from previous call
!!!      FB = HOLD (4)        ! Redundant
         FC = HOLD (5)
         C  = HOLD (6)
         D  = HOLD (7)
         E  = HOLD (8)
         XM = HOLD (9)
         P  = HOLD (10)
         Q  = HOLD (11)
         R  = HOLD (12)
         S  = HOLD (13)

C        The next test had been at the end, after FB = F(B):

         FB = F
         BEGIN = FB * (FC /ABS (FC)) .GT. ZERO
      END IF

      IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020) NUMFUN, X, FB, CALLER

      IF (BEGIN) THEN  ! Begin step

         C  = A
         FC = FA
         D  = B - A
         E  = D

      END IF

      IF (ABS (FC) .LT. ABS (FB)) THEN
         A  = B
         B  = C
         C  = A
         FA = FB
         FB = FC
         FC = FA
      END IF

C     Convergence test:

      TOL1 = TWO * EPS * ABS (B) + HALF * TOL
      XM = HALF * (C - B)
      IF (ABS (XM) .LE. TOL1) GO TO 90
      IF (FB .EQ. ZERO) GO TO 90

C     Is bisection necessary?

      IF (ABS (E) .LT. TOL1) GO TO 70
      IF (ABS (FA) .LE. ABS (FB)) GO TO 70

C     Is quadratic interpolation possible?

      IF (A .EQ. C) THEN  ! No - do linear interpolation

         S = FB / FA
         P = TWO * XM * S
         Q = ONE - S

      ELSE                ! Inverse quadratic interpolation

         Q = FA / FC
         R = FB / FC
         S = FB / FA
         P = S * (TWO * XM * Q * (Q - R) - (B - A) * (R - ONE))
         Q = (Q - ONE) * (R - ONE) * (S - ONE)

      END IF

C     Adjust signs:

      IF (P .GT. ZERO) Q = -Q
      P = ABS (P)

C     Is interpolation acceptable?

      IF ((P + P) .GE. (THREE * XM * Q - ABS (TOL1 * Q))) GO TO 70
      IF (P .GE. ABS (HALF * E * Q)) GO TO 70

      E = D
      D = P / Q
      GO TO 80

C     Bisection:

   70 D = XM
      E = D

C     Complete step:

   80 A = B
      FA = FB
      IF (ABS (D) .GT. TOL1) B = B + D
      IF (ABS (D) .LE. TOL1) B = B + SIGN (TOL1, XM)

C*****FB = F(B)

C     Return for the next function value.

      X = B
      HOLD (1) = A
      HOLD (2) = B
      HOLD (3) = FA
      HOLD (4) = FB
      HOLD (5) = FC
      HOLD (6) = C
      HOLD (7) = D
      HOLD (8) = E
      HOLD (9) = XM
      HOLD (10) = P
      HOLD (11) = Q
      HOLD (12) = R
      HOLD (13) = S

      GO TO 99


C     Done:

   90 X = B
      F = FB
      ISTAT = 0


   99 RETURN


C     Formats:

 1010 FORMAT (/, ' ZERORC:  Iteration history for calls by ', A, //
     >        ' Iter', 15X, 'X', 24X, 'F(X)')
 1020 FORMAT (1X, I3, 1P, 2E14.6, 3X, A)
 1030 FORMAT (/, ' ZERORC:  AX >= BX is illegal: ', 1P, 2E15.6,
     >        '    Caller:  ', A)
 1040 FORMAT (/, ' ZERORC:  Suspicious input for (max.) NUMFUN: ', I10,
     >        '    Caller:  ', A)
 1050 FORMAT (/, ' ZERORC:  End-point function values have same sign: ',
     >        /, 1P, 2E15.6, '   Caller: ', A)

      END SUBROUTINE ZERORC
C+----------------------------------------------------------------------
C
      FUNCTION BESSEL (J, H, DEL)
C
C     One-liner: First derivative using central 3-point formula
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation using the central
C     3-point formula.  The data must be in the form of arrays containing
C     finite difference interval lengths and 2-point forward difference
C     derivatives.  BESSEL is intended to be used by PLSFIT for determin-
C     ing end conditions on an interval for (non-monotonic) interpolation
C     by piecewise cubics.  See the PLSFIT header for more details.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     BESSEL  R                 O    The function value is the adjusted
C                                    derivative.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), BESSEL

C     Local variables.

      REAL
     &   WEIGHT

C     Execution.
C     ----------

C     Estimate first derivative on left (J = 0) or right side (J = 1) of
C     an interval.

      WEIGHT = H (J) / (H (J) + H (J - 1))
      BESSEL = WEIGHT * DEL (J - 1) + (ONE - WEIGHT) * DEL (J)

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION BRODLIE (J, H, DEL)
C
C     One-liner: First derivative, adjusted for monotonicity
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        BRODLIE is intended to be used by PLSFIT for determining end
C     conditions on an interval for monotonic interpolation by piecewise
C     cubics. The data must be in the form of arrays containing finite
C     difference interval lengths and 2-point forward difference deriva-
C     tives. See the PLSFIT header for more details.
C
C        The method is due to Brodlie, Butland, Carlson, and Fritsch,
C     as referenced in the PLSFIT header.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C
C     BRODLIE R                 O    The function value is the adjusted
C                                    derivative.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C     04 Dec. 2002    DAS    SIGN work-around for Intel compiler bug.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE, THIRD
      PARAMETER
     &  (ZERO   = 0.0E+0,
     &   ONE    = 1.0E+0,
     &   THIRD  = ONE / 3.0E+0)

C     Arguments.

      INTEGER
     &   J
      REAL
     &   BRODLIE, H (-1:1), DEL (-1:1)

C     Local variables.

      REAL
     &   ALPHA, PRODUCT

C     Execution.
C     ----------

C     Compare the algebraic signs of the two DEL's.  Have to test that
C     at least one is positive to avoid a zero denominator (this fancy
C     test permits one term to be zero, but the answer below is zero
C     anyway in these cases).  The trick is to work around the SIGN
C     function, which returns positive even if its 2nd argument is zero.

C**** NO:  SIGN misbehaves on Intel systems when the 2nd argument is zero.

      PRODUCT = DEL (J - 1) * DEL (J)

      IF (PRODUCT == ZERO) THEN

         BRODLIE = ZERO

      ELSE IF (SIGN (ONE, -DEL (J - 1)) .NE. SIGN (ONE, DEL (J))) THEN

C        Form "weighted harmonic mean" of the two finite-difference
C        derivative approximations.  Note that we try to avoid overflow
C        by not multiplying them together directly.

         ALPHA   = THIRD * (ONE + H (J) / (H (J - 1) + H (J)))
         BRODLIE = PRODUCT / (ALPHA * DEL (J) +
     &      (ONE - ALPHA) * DEL (J - 1))
      ELSE

C        The signs differ, so make this point a local extremum.

         BRODLIE = ZERO
      END IF

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION BUTLAND (J, H, DEL)
C
C     One-liner: First derivative, non-central 3-point formula, adjusted
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation for PLSFIT over an
C     interval at a data boundary, using a modified forward or backward
C     3-point formula.  The data must be in the form of arrays containing
C     finite difference interval lengths and 2-point forward difference
C     derivatives, and the differencing direction is controlled by a flag.
C     See PLSFIT for more details, or THREEPT for the pure 3-pt. formula.
C
C        The "shape preserving adjustments" are from PCHIP, a monotone
C     piecewise cubic interpolation package by F. N. Fritsch.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C
C     BUTLAND R                 O    The function value is the adjusted
C                                    derivative.
C
C     Environment:  Fortran 90
C     ------------
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding, as THREEPT.
C     20 June 1991    DAS    Monotonic form renamed BUTLAND; THREEPT
C                            is now the pure 3-point formula.
C     04 Dec. 2002     "     SIGN work-arounds for Intel compiler bug.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), BUTLAND

C     Local constants.

      REAL, PARAMETER ::
     &   ZERO  = 0.0E+0,
     &   ONE   = 1.0E+0,
     &   THREE = 3.0E+0

C     Local variables.

      INTEGER
     &   STEP
      REAL
     &   DMAX, WEIGHT
      LOGICAL
     &   CONSTRAIN

C     Execution.
C     ----------

C     Estimate first derivative on a left-hand boundary using a 3-point
C     forward difference (STEP = +1), or with a backward difference for
C     the right-hand boundary (STEP = -1).

      STEP = 1 - J - J   ! J here is consistent with related modules.

C     In {H, DEL} form, the derivative looks like a weighted average.

C***  Avoid zero as the second argument of SIGN: Intel compiler misbehaves.

      IF (DEL (0) == ZERO) THEN

         BUTLAND = ZERO

      ELSE ! BUTLAND below cannot be zero

         WEIGHT  = -H (0) / (H (0) + H (STEP))
         BUTLAND = WEIGHT * DEL (STEP) + (ONE - WEIGHT) * DEL (0)

C        Shape-preserving adjustments.  Note that we try to avoid overflow
C        by not multiplying quantities directly.

         IF (SIGN (ONE, BUTLAND) /= SIGN (ONE, DEL (0))) THEN

C           Defer to the estimate closest to the boundary.

            BUTLAND = ZERO

C******  ELSE IF (SIGN (ONE, DEL (0)) .NE. SIGN (ONE, DEL (STEP))) THEN
         ELSE

            IF (DEL (STEP) == ZERO) THEN
               CONSTRAIN = DEL (0) < ZERO
            ELSE
               CONSTRAIN = SIGN (ONE, DEL (0)) /= SIGN (ONE, DEL (STEP))
            END IF

            IF (CONSTRAIN) THEN

C              If the monotonicity switches, may need to bound the estimate.

               DMAX = THREE * DEL (0)
               IF (ABS (BUTLAND) > ABS (DMAX)) BUTLAND = DMAX
            END IF

         END IF

      END IF

C     Termination.
C     ------------

      END
C+----------------------------------------------------------------------
C
      FUNCTION THREEPT (J, H, DEL)
C
C     One-liner: First derivative, non-central 3-point formula
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation for PLSFIT over an
C     interval at a data boundary, using a forward or backward 3-point
C     formula.  The data must be in the form of arrays containing finite
C     difference interval lengths and 2-point forward difference deriva-
C     tives, and the differencing direction is controlled by a flag. See
C     PLSFIT for more details.
C
C        See module BUTLAND for a version with "shape-preserving"
C     adjustments.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right. 
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     THREEPT R                 O    The function value is the derivative.
C
C     Environment:  VAX/VMS; FORTRAN 77
C     ------------
C
C     IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C     06 June 1991    DAS    Original THREEPT renamed BUTLAND; THREEPT
C                            now gives unmodified 1-sided 3-pt. results.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), THREEPT

C     Local constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Local variables.

      INTEGER
     &   STEP
      REAL
     &   WEIGHT

C     Execution.
C     ----------

C     Estimate first derivative on a left-hand boundary using a 3-point
C     forward difference (STEP = +1), or with a backward difference for
C     the right-hand boundary (STEP = -1).

      STEP = 1 - J - J   ! J here is consistent with related modules.

C     In {H, DEL} form, the derivative looks like a weighted average.

      WEIGHT  = -H (0) / (H (0) + H (STEP))
      THREEPT = WEIGHT * DEL (STEP) + (ONE - WEIGHT) * DEL (0)

C     Termination.
C     ------------

      RETURN
      END
