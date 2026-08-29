(* ::Package:: *)

(* ============================================================================
   GSE COMPANION NOTEBOOK -- all model figures in one place
   Brundage, Darity, Tavani: "Global Stratification Economics and
   International Extraction"
   Daniele Tavani, August 2026

   
   CONTENTS
   Part A  Theory model and Figures II-IV
           A1 primitives (retained share, extraction schedule, transition G)
           A2 fixed points and stability (bracketed bisection, no FindRoot)
           A3 validation block
           A4 exports: G_unique.dat, G_multiple.dat, G_family.dat
              (read by the pgfplots figures in the paper) + preview PDFs
   Part B  Calibration and Figures V-VII
           B1 data (Ricci 2019 x PWT 11.0) + import hook for annual WIOD data
           B2 estimation of the enforcement elasticity theta
           B3 exact multiplicity frontier (saddle-node characterization)
           B4 figure exports: fig_emp_schedule.pdf, fig_emp_diagnostic.pdf,
              fig_emp_2d.pdf
   Part C  Interactive explorer (Manipulate)

   HOW THE PAPER USES THESE OUTPUTS
   Figures II-IV are drawn natively by pgfplots and READ THE .dat FILES:
   upload G_unique.dat, G_multiple.dat, G_family.dat to Overleaf (Figures/).
   The preview PDFs of II-IV are for checking and slides, not for the paper.
   Figures V-VII are included as PDFs: upload the three fig_emp_*.pdf files.

   VALIDATION TARGETS (all verified against an independent Python
   implementation; a clean top-to-bottom run must reproduce them)
     Figure III equilibria:  Omega* = 1.7610 (G'=0.6687, stable)
                             Omega* = 3.2533 (G'=1.6825, unstable)
                             Omega* = 8.7174 (G'=0.4479, stable)
     Figure II equilibrium:  Omega* = 3.0261 (e*=0.2000, G'=0.4536)
     Saddle-node:            chi0crit = 1821.8
     Figure IV counts:       3 equilibria at chi0 = 202 and 700
                             1 equilibrium at chi0 = 1822 and 2400
     Estimation:             theta = 0.353, s.e. 0.116, R2 = 0.42,
                             CI = [0.126, 0.580]; grid range [0.334, 0.359]
     Exact frontier:         theta_min(0.90) = 0.543
   Runtime: Part A < 1 min; Part B ~ 2-4 min (frontier + 2-D contours).
   ==========================================================================*)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[] /. $Failed -> Directory[]];

(* ======================= PART A: THEORY MODEL ============================ *)

(* ---- A1. Primitives ----------------------------------------------------- *)

(* Parameters live in an Association so theory and calibration share code.  *)
theoryPar = <|"alpha" -> 0.5, "beta" -> 0.25, "sigma" -> 1.25,
   "kappa" -> 1.0, "lam" -> 0.5, "b" -> 0.5, "delta" -> 0.5,
   "chi0" -> 202., "theta" -> 3., "ebar" -> 0.5|>;

DD[e_, p_] := 1 - e + p["b"] e - p["delta"] e^2/2;      (* retained share  *)
mOf[p_] := 1 + p["lam"] (1 - p["b"]);
ldOf[p_] := p["lam"] p["delta"];

(* Extraction schedule, eq. (eq:extraction) with the cap ebar               *)
eSchedule[Om_?NumericQ, p_] := Module[{X = p["chi0"] Om^-p["theta"]},
   If[X <= ldOf[p] + mOf[p], p["ebar"],
    Min[mOf[p]/(X - ldOf[p]), p["ebar"]]]];

(* Transition function, eq. (eq: G)                                        *)
GG[Om_?NumericQ, p_] := Module[{e = eSchedule[Om, p], d, gam, inv},
   d = DD[e, p]; gam = p["beta"]/(1 - p["alpha"]);
   inv = 1/(1 - p["alpha"]);
   p["sigma"] (p["kappa"] Om^gam d^-inv + e/d)];

Gprime[Om_?NumericQ, p_] := (GG[Om + 10.^-6, p] - GG[Om - 10.^-6, p])/(2 10.^-6);

(* ---- A2. Fixed points and stability ------------------------------------ *)

(* Sign-scan on a geometric grid, then bisection on each bracket.
   No FindRoot anywhere in this file: every root is bracketed first.       *)
fixedPoints[p_, hi_: 60, n_: 20000] := Module[{Om, h, roots = {}, lo, hh, mid},
   Om = Exp@Subdivide[Log[10.^-3], Log[hi], n];
   h = (GG[#, p] - #) & /@ Om;
   Do[If[h[[i]] h[[i + 1]] < 0,
     {lo, hh} = {Om[[i]], Om[[i + 1]]};
     While[hh - lo > 10.^-9 (1. + hh),
      mid = (lo + hh)/2.;
      If[(GG[mid, p] - mid) h[[i]] > 0, lo = mid, hh = mid]];
     AppendTo[roots, (lo + hh)/2.]], {i, n}];
   roots];

equilibriumReport[p_, hi_: 60] := TableForm[
   {#, eSchedule[#, p], Gprime[#, p],
      If[Abs[Gprime[#, p]] < 1, "stable", "UNSTABLE"]} & /@
    fixedPoints[p, hi],
   TableHeadings -> {None, {"Omega*", "e*", "G'(Omega*)", ""}}];

(* Trap fixed point of the capped (chi0- and theta-free) branch, and the
   saddle-node value of eq. (eq: crit). Bracketed bisection: the capped
   branch a Om^gamma + c rises to a single peak at (a gamma)^(1/(1-gamma))
   and crosses the diagonal exactly once to the right of it.               *)
Omega3[p_] := Module[{gam, inv, a, c, lo, hi, mid},
   gam = p["beta"]/(1 - p["alpha"]); inv = 1/(1 - p["alpha"]);
   a = p["sigma"] p["kappa"] DD[p["ebar"], p]^-inv;
   c = p["sigma"] p["ebar"]/DD[p["ebar"], p];
   lo = Max[(a gam)^(1/(1 - gam)), 10.^-6];
   hi = Max[2 lo, 2.];
   While[a hi^gam + c - hi > 0, hi *= 2.];
   While[hi - lo > 10.^-10 (1. + hi),
    mid = (lo + hi)/2.;
    If[a mid^gam + c - mid > 0, lo = mid, hi = mid]];
   (lo + hi)/2.];

chi0crit[p_] := (ldOf[p] + mOf[p]/p["ebar"]) Omega3[p]^p["theta"];

(* ---- A3. Validation ----------------------------------------------------- *)

Print["=== Part A validation ==="];
Print["Figure III equilibria (expect 1.7610/3.2533/8.7174, \
G' 0.6687/1.6825/0.4479):"];
Print[equilibriumReport[theoryPar]];
Print["Figure II equilibrium (expect 3.0261, e*=0.2000, G'=0.4536):"];
Print[equilibriumReport[<|theoryPar, "chi0" -> 6.5, "theta" -> 0.|>, 20]];
Print["chi0crit (expect 1821.8): ", chi0crit[theoryPar]];
Print["Figure IV counts (expect 3,3,1,1): ",
  Table[Length@fixedPoints[<|theoryPar, "chi0" -> c|>, 120],
   {c, {202., 700., 1822., 2400.}}]];

(* ---- A4. Figure II-IV data exports and previews ------------------------- *)

figStyle = {Frame -> True, Axes -> False, FrameStyle -> GrayLevel[0.25],
   LabelStyle -> Directive[FontFamily -> "Times", 12],
   AspectRatio -> 0.72, ImageSize -> 480};

grid8  = Subdivide[10.^-3, 8., 220];
grid12 = Subdivide[10.^-3, 12., 220];
famChis = {202., 700., 1822., 2400.};

parII = <|theoryPar, "chi0" -> 6.5, "theta" -> 0.|>;

Export["G_unique.dat",
  Prepend[Transpose[{grid8, GG[#, parII] & /@ grid8}], {"Om", "G"}],
  "Table"];
Export["G_multiple.dat",
  Prepend[Transpose[{grid12, GG[#, theoryPar] & /@ grid12}], {"Om", "G"}],
  "Table"];
Export["G_family.dat",
  Prepend[Transpose[Join[{grid12},
     Table[GG[#, <|theoryPar, "chi0" -> c|>] & /@ grid12, {c, famChis}]]],
   {"Om", "G202", "G700", "G1822", "G2400"}], "Table"];
Print["Exported: G_unique.dat, G_multiple.dat, G_family.dat"];

(* Preview PDFs (checking/slides only; the paper draws these via pgfplots)  *)
eqII = fixedPoints[parII, 20]; eqIII = fixedPoints[theoryPar];


figThUnique = Show[
   Plot[{GG[Om, parII], Om}, {Om, 0.01, 8},
    PlotStyle -> {Directive[Blue, Thickness[0.004]],
      Directive[Black, Dashed, Thickness[0.002]]}],
   figStyle, PlotRange -> {{0, 8}, {0, 8}},
   FrameLabel -> {"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(t\)]\)", "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(t + 1\)]\)"},
   Epilog -> {PointSize[0.014], Point[{#, #}] & /@ eqII,
     Text[Style["\[CapitalOmega]*", 12], {First[eqII] - 0.45,
       First[eqII] + 0.35}]}]


Export["fig_theory_unique.pdf", figThUnique];



figThMultiple = Show[
   Plot[{GG[Om, theoryPar], Om}, {Om, 0.01, 12},
    PlotStyle -> {Directive[Blue, Thickness[0.004]],
      Directive[Black, Dashed, Thickness[0.002]]}],
   figStyle, PlotRange -> {{0, 12}, {0, 12}},
   FrameLabel -> {"\!\(\*SubscriptBox[
StyleBox[\"\[CapitalOmega]\",\nFontSlant->\"Italic\"], \(t\)]\)", "\!\(\*SubscriptBox[
StyleBox[\"\[CapitalOmega]\",\nFontSlant->\"Italic\"], \(t + 1\)]\)"},
   Epilog -> {PointSize[0.014],
     Point[{#, #}] & /@ {eqIII[[1]], eqIII[[3]]},
     {EdgeForm[Black], White,
      Rectangle[{eqIII[[2]] - 0.13, eqIII[[2]] - 0.13},
       {eqIII[[2]] + 0.13, eqIII[[2]] + 0.13}]}}]


Export["fig_theory_multiple.pdf", figThMultiple];



figThTrap = Show[
   Plot[Evaluate@Join[
      Table[GG[Om, <|theoryPar, "chi0" -> c|>], {c, famChis}], {Om}],
    {Om, 0.01, 12},
    PlotStyle -> {Directive[Blue, Thickness[0.004]],
      Directive[Purple, Thickness[0.004]],
      Directive[Red, Thickness[0.004]],
      Directive[Gray, Thickness[0.004]],
      Directive[Black, Dashed, Thickness[0.002]]},
    PlotLegends -> Placed[{"\[Chi]0=202", "\[Chi]0=700",
       "\[Chi]0=1822\[TildeTilde]\[Chi]0crit", "\[Chi]0=2400",
       "45\[Degree]"}, {0.18, 0.72}]],
   figStyle, PlotRange -> {{0, 12}, {0, 12}},
   FrameLabel -> {"\[CapitalOmega]t", "\[CapitalOmega]t+1"},
   Epilog -> {Red, Thickness[0.003], Circle[{8.7174, 8.7174}, 0.18]}]
Export["fig_theory_trapescape.pdf", figThTrap];


Print["Exported previews: fig_theory_unique.pdf, fig_theory_multiple.pdf, \
fig_theory_trapescape.pdf"]



(* ======================= PART B: CALIBRATION ============================= *)

(* ---- B1. Data ----------------------------------------------------------- *)

(* Ricci (2019, Table 2) extraction rates x PWT 11.0.
   Columns: year, region, e, OmegaPerWorker, OmegaAggregate, sigma          
calibData = {
   {1995, "East Europe",   0.168,  3.09, 23.85, 1.63},
   {1995, "Latin America", 0.062,  3.54, 12.74, 1.70},
   {1995, "China",         0.173, 26.58, 13.21, 0.76},
   {1995, "India",         0.209, 31.10, 27.85, 1.49},
   {1995, "Other Asia",    0.133,  7.44, 23.20, 1.20},
   {2000, "East Europe",   0.143,  2.96, 24.60, 1.35},
   {2000, "Latin America", 0.050,  3.38, 11.46, 1.41},
   {2000, "China",         0.099, 16.93,  8.41, 0.89},
   {2000, "India",         0.262, 21.31, 18.72, 1.33},
   {2000, "Other Asia",    0.152,  6.16, 18.38, 1.24},
   {2007, "East Europe",   0.075,  2.64, 22.85, 0.98},
   {2007, "Latin America", 0.064,  3.60, 10.82, 1.27},
   {2007, "China",         0.109,  8.96,  4.45, 0.71},
   {2007, "India",         0.170, 12.40, 10.18, 0.77},
   {2007, "Other Asia",    0.104,  4.38, 12.56, 0.98}};

(* Import hook for Brendan's annual estimates. Expected CSV layout:
   year, region, e, omega_pw, omega_agg[, sigma] with a header row. Keep
   the two WIOD vintages separate (the releases differ):
     calib9509 = importCalib["brendan_wiod2013_1995_2009.csv"];
     calib0014 = importCalib["brendan_wiod2016_2000_2014.csv"];
   then rerun B2 with the new dataset.                                     *)
importCalib[file_String] := Rest[Import[file, "CSV"]]; *)
(* ================= ANNUAL CALIBRATION DATA ================= *)

(* Import Brendan's CSVs and rearrange them into Daniele's internal format:
   {year, region, e, OmegaPerWorker, OmegaAggregate, sigma}
*)

importBrendan[file_String] := Module[
   {raw, header, rows, pos},

   raw = Import[file, "CSV"];
   header = First[raw];
   rows = Rest[raw];

   pos[name_] := First@FirstPosition[header, name];

   ({
       #[[pos["year"]]],
       #[[pos["region"]]],
       N@#[[pos["e"]]],
       N@#[[pos["Omega_pw"]]],
       N@#[[pos["Omega_agg"]]],
       N@#[[pos["sigma"]]]
       } &) /@ rows
   ];

(* Original WIOD release: 1995-2009 *)
calib9509 =
  importBrendan["brendan_wiod2013_1995_2009.csv"];

(* Subsequent WIOD release: 2000-2014 *)
calib0014 =
  importBrendan["brendan_wiod2016_2000_2014.csv"];

Print["1995-2009 N = ", Length[calib9509]];
Print["1995-2009 regions = ",
  DeleteDuplicates[calib9509[[All, 2]]]];

Print["2000-2014 N = ", Length[calib0014]];
Print["2000-2014 regions = ",
  DeleteDuplicates[calib0014[[All, 2]]]];



(* ---- B2. Estimation ----------------------------------------------------- *)

(* Estimating equation from the interior schedule (eq:extraction):
     ln( m/e + lambda delta ) = ln chi0 - theta ln Omega
   so -theta is the OLS slope and H0: theta = 0 is the null of no
   power-based enforcement.                                                *)
estimateTheta[data_, lam_: 0.5, b_: 0.5, delta_: 0.5, col_: 4] :=
  Module[{m = 1 + lam (1 - b), e, Om, fit, th, se},
   e = data[[All, 3]]; Om = data[[All, col]];
   fit = LinearModelFit[Transpose[{Log[Om], Log[m/e + lam delta]}], t, t];
   th = -fit["BestFitParameters"][[2]];
   se = fit["ParameterErrors"][[2]];
   <|"theta" -> th, "se" -> se, "R2" -> fit["RSquared"],
     "chi0hat" -> Exp[fit["BestFitParameters"][[1]]],
     "CI" -> {th - 1.96 se, th + 1.96 se}|>];

(*basePW  = estimateTheta[calibData];
baseAgg = estimateTheta[calibData, 0.5, 0.5, 0.5, 5];
gridRange = MinMax@Flatten@Table[
    estimateTheta[calibData, lam, b, d]["theta"],
    {lam, {0.2, 0.5, 0.8}}, {b, {0.3, 0.5, 0.9}}, {d, {0.3, 0.5, 0.9}}];

Print["=== Part B validation ==="];
Print["per-worker (expect theta 0.353, se 0.116, R2 0.42): ", basePW];
Print["aggregate  (expect theta 0.353, se 0.242, R2 0.14): ", baseAgg];
Print["grid range (expect {0.334, 0.359}): ", gridRange];*)
(* Aggregate Omega is column 5 in Daniele's internal data structure *)

base9509 =
  estimateTheta[calib9509, 0.5, 0.5, 0.5, 5];

base0014 =
  estimateTheta[calib0014, 0.5, 0.5, 0.5, 5];

Print["=== ANNUAL AGGREGATE CALIBRATION ==="];

Print["WIOD 2013, 1995-2009: ", base9509];

Print["WIOD 2016, 2000-2014: ", base0014];



(* ---- B3. Exact multiplicity frontier ------------------------------------ *)

(* theta_min(alpha+beta) = smallest theta such that SOME chi0 delivers three
   equilibria. No chi0 scan is needed: under weak gains (b<1), G is
   pointwise decreasing in chi0 on the interior branch (Proposition 1),
   while the capped branch -- hence the trap fixed point Omega3 -- does not
   depend on chi0 or theta. The deepest admissible dip of G below the
   45-degree line therefore occurs at the saddle-node value chi0crit of
   eq. (eq: crit), and theta_min is the root of that dip in theta.
   Sanity anchor: with theoryPar this machinery returns Omega3 = 8.7174
   and chi0crit = 1821.8 (Part A).                                         *)

empPar = <|"alpha" -> 0.578, "beta" -> 0.322, "sigma" -> 0.705,
   "kappa" -> 1.175, "lam" -> 0.5, "b" -> 0.5, "delta" -> 0.5,
   "chi0" -> 1., "theta" -> 1., "ebar" -> 0.2618|>;

dip[p_, theta_, n_: 30000] := Module[{q, O3, Om, gam, inv, X, ee, d, G},
   q = <|p, "theta" -> theta|>;
   O3 = Omega3[q];
   q["chi0"] = chi0crit[q];
   gam = q["beta"]/(1 - q["alpha"]); inv = 1/(1 - q["alpha"]);
   Om = Exp@Subdivide[Log[10.^-4], Log[0.99 O3], n]; (* exclude pinned end *)
   X = q["chi0"] Om^-theta;
   ee = Map[If[# <= ldOf[q] + mOf[q], q["ebar"],
       Min[mOf[q]/(# - ldOf[q]), q["ebar"]]] &, X];
   d = DD[ee, q];
   G = q["sigma"] (q["kappa"] Om^gam d^-inv + ee/d);
   Min[G - Om]];

thetaMin[ab_, p_: empPar] := Module[{q, lo = 0.02, hi = 10.0, mid},
   q = <|p, "beta" -> ab - p["alpha"]|>;
   If[dip[q, hi] > 0, Return[Missing["NoMultiplicity"]]];
   While[hi - lo > 0.002, mid = (lo + hi)/2.;
    If[dip[q, mid] < 0, hi = mid, lo = mid]];
   hi];

(* Verified reference values (empirical parameterization):
     alpha+beta : 0.80   0.85   0.90   0.92   0.94   0.96
     theta_min  : 0.940  0.752  0.543  0.449  0.346  0.231
   Frontier enters the CI (upper bound 0.580) at alpha+beta ~ 0.885 and
   crosses the point estimate (0.353) at alpha+beta ~ 0.939.               *)
Print["frontier check at alpha+beta = 0.90 (expect ~0.543): ",
  thetaMin[0.90]]



(* ---- B4. Figures V-VII -------------------------------------------------- *)

(* CONFIRM before circulating: the shaded band in Figure V is described in
   the caption as the Hickel et al. (2021) drain estimate. Set the interval
   to whatever the July figure used; [0.20, 0.2618] is a placeholder. Note
   ebar = 0.2618 equals the sample-maximum extraction rate (India 2000),
   which is a different object from the Hickel band.                       *)
hickelBand = {0.07, 0.09};
paperFigStyle = {
   Frame -> {{True, False}, {True, False}},
   Axes -> False,
   FrameStyle -> Directive[Black, Thickness[0.0015]],
   LabelStyle -> Directive[
     FontFamily -> "Latin Modern Roman",
     FontSize -> 11,
     Black
     ],
   BaseStyle -> Directive[
     FontFamily -> "Latin Modern Roman",
     FontSize -> 10
     ],
   Background -> White
   };

regionColors = <|
   "China" -> RGBColor[31/255, 119/255, 180/255],
   "Latin America" -> RGBColor[255/255, 127/255, 14/255],
   "India" -> RGBColor[44/255, 160/255, 44/255],
   "East Europe" -> RGBColor[214/255, 39/255, 40/255],
   "Other Asia" -> RGBColor[148/255, 103/255, 189/255]
   |>;

regionLabelOffsets = <|
   "China" -> {0.025, 0.004},
   "East Europe" -> {0.025, 0.010},
   "India" -> {0.025, 0.004},
   "Latin America" -> {0.025, -0.006},
   "Other Asia" -> {0.025, 0.004}
   |>;

(* -- Figure V: data vs fitted extraction schedule (log-x by coordinate
      transform; ScalingFunctions does not survive Show over mixed
      graphics) ------------------------------------------------------------ *)
(*figSchedule = Module[{m, ld, th, c0, lx = Log10, ticks, regLines, dataPts},
   m = 1 + 0.5 (1 - 0.5); ld = 0.5*0.5;
   th = basePW["theta"]; c0 = basePW["chi0hat"];
   ticks = {lx[#], #} & /@ {2, 3, 5, 10, 20, 30, 40};
   regLines = Table[{lx[#[[4]]], #[[3]]} & /@
      SortBy[Select[calibData, #[[2]] == r &], First],
     {r, DeleteDuplicates[calibData[[All, 2]]]}];
   dataPts = Table[{lx[#[[4]]], #[[3]]} & /@
      Select[calibData, #[[1]] == yr &], {yr, {1995, 2000, 2007}}];
   Show[
    Plot[m/(c0 (10^u)^-th - ld), {u, lx[1.6], lx[42]},
     PlotStyle -> Directive[RGBColor[0.69, 0.19, 0.19], Thickness[0.005]],
     Prolog -> {RGBColor[0.91, 0.91, 0.95],
       Rectangle[{lx[1.6], hickelBand[[1]]}, {lx[42], hickelBand[[2]]}]}],
    ListLinePlot[regLines,
     PlotStyle -> Directive[GrayLevel[0.78], Thickness[0.0018]]],
    ListPlot[dataPts,
     PlotMarkers -> {{"\[EmptyCircle]", 16}, {"\[EmptySquare]", 14},
       {"\[EmptyUpTriangle]", 15}},
     PlotStyle -> RGBColor[0.12, 0.31, 0.47],
     PlotLegends -> Placed[{"1995", "2000", "2007"}, {0.14, 0.78}]],
    figStyle,
    FrameTicks -> {{Automatic, None}, {ticks, None}},
    FrameLabel -> {"\[CapitalOmega] (capital per worker, core/periphery)",
      "extraction rate e"},
    PlotRange -> {{Log10[1.6], Log10[42]}, {0, 0.30}},
    Epilog -> {GrayLevel[0.45],
      Text[Style["Hickel et al. band", 9], {Log10[28], 0.245}]}]]*)
      
      makeSchedule[data_, base_, markYears_, tag_] :=
 Module[
  {m, ld, th, c0, lx = Log10, regions, regLines,
   pointPlots, yearMarkers, firstYear, regionLabels,
   logTicks, yTicks, fitFun, loFun, hiFun, legendBox, fig},

  m = 1 + 0.5 (1 - 0.5);
  ld = 0.5*0.5;

  th = base["theta"];
  c0 = base["chi0hat"];

  regions = DeleteDuplicates[data[[All, 2]]];

  yearMarkers = <|
    markYears[[1]] -> "\[FilledCircle]",
    markYears[[2]] -> "\[FilledSquare]",
    markYears[[3]] -> "\[FilledUpTriangle]"
    |>;

  (* Log-axis styling like the paper: label 10, keep smaller ticks unlabeled *)
  logTicks = Join[
    {{lx[10], Superscript[10, 1]}},
    ({lx[#], ""} &) /@ {2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40}
    ];

  yTicks = Table[{y, NumberForm[y, {3, 2}]}, {y, 0, 0.30, 0.05}];

  (* Colored within-region annual paths *)
  regLines =
   Table[
    ListLinePlot[
     ({lx[#[[5]]], #[[3]]} &) /@
      SortBy[Select[data, #[[2]] == r &], First],
     PlotStyle ->
      Directive[
       regionColors[r],
       Opacity[0.35],
       Thickness[0.0012]
       ]
     ],
    {r, regions}
    ];

  (* Filled year markers; color identifies region *)
  pointPlots =
   Flatten[
    Table[
     With[
      {pts =
        ({lx[#[[5]]], #[[3]]} &) /@
         Select[data, #[[1]] == yr && #[[2]] == r &]},
      If[
       Length[pts] == 0,
       Nothing,
       ListPlot[
        pts,
        PlotMarkers -> {yearMarkers[yr], 10},
        PlotStyle -> regionColors[r]
        ]
       ]
      ],
     {yr, markYears},
     {r, regions}
     ]
    ];

  (* Label each region at first-year observation *)
  firstYear = Min[data[[All, 1]]];

  regionLabels =
   Graphics[
    Table[
     With[
      {row =
        First@Select[
          data,
          #[[1]] == firstYear && #[[2]] == r &
          ]},
      Text[
       Style[r, 9, regionColors[r]],
       {lx[row[[5]]], row[[3]]} + regionLabelOffsets[r],
       {-1, 0}
       ]
      ],
     {r, regions}
     ]
    ];

  (* Fitted schedule and CI schedules *)
  fitFun[om_] :=
   Min[m/(c0 om^-th - ld), empPar["ebar"]];

  loFun[om_] :=
   Min[m/(c0 om^-base["CI"][[1]] - ld), empPar["ebar"]];

  hiFun[om_] :=
   Min[m/(c0 om^-base["CI"][[2]] - ld), empPar["ebar"]];

  (* Manual legend to mimic paper *)
  legendBox =
   Inset[
    Grid[
     {
      {
       Style["\[FilledCircle]", 10, GrayLevel[0.35]],
       ToString[markYears[[1]]],
       Spacer[14],
       Style["\[FilledUpTriangle]", 10, GrayLevel[0.35]],
       ToString[markYears[[3]]]
       },
      {
       Style["\[FilledSquare]", 10, GrayLevel[0.35]],
       ToString[markYears[[2]]],
       Spacer[14],
       Graphics[
        {
         GrayLevel[0.20],
         Thickness[0.08],
         Line[{{0, 0}, {1, 0}}]
         },
        ImageSize -> {28, 8}
        ],
       Row[{
         "fitted, ",
         "\!\(\*OverscriptBox[\(\[Theta]\), \(^\)]\)=",
         NumberForm[th, {3, 2}]
         }]
       }
      },
     Alignment -> Left,
     Spacings -> {0.45, 0.25}
     ],
    Scaled[{0.18, 0.91}]
    ];

  fig =
   Show[

    (* Hickel band *)
    Graphics[
     {
      GrayLevel[0.92],
      Rectangle[
       {lx[1.6], hickelBand[[1]]},
       {lx[42], hickelBand[[2]]}
       ]
      }
     ],

    (* Dashed 95% CI *)
    Plot[
     {
      loFun[10^u],
      hiFun[10^u]
      },
     {u, lx[1.6], lx[42]},
     PlotStyle -> {
       Directive[GrayLevel[0.55], Dashed, Thickness[0.002]],
       Directive[GrayLevel[0.55], Dashed, Thickness[0.002]]
       }
     ],

    (* Main fitted schedule *)
    Plot[
     fitFun[10^u],
     {u, lx[1.6], lx[42]},
     PlotStyle ->
      Directive[
       GrayLevel[0.20],
       Thickness[0.004]
       ]
     ],

    Sequence @@ regLines,
    Sequence @@ pointPlots,
    regionLabels,

    paperFigStyle,

    AspectRatio -> 0.55,
    ImageSize -> 520,

    FrameTicks -> {
      {yTicks, None},
      {logTicks, None}
      },

    FrameLabel -> {
      "Relative power \[CapitalOmega]  (core / region, aggregate capital)",
      "Extraction rate  e"
      },

    PlotRange -> {
      {lx[1.6], lx[42]},
      {0, 0.30}
      },

    Epilog -> {
      legendBox,

      Text[
       Style["Hickel et al. band", 8, GrayLevel[0.45]],
       {lx[3.0], Mean[hickelBand]}
       ],

      Text[
       Style["95% CI", 8, GrayLevel[0.45]],
       {lx[29], loFun[29] + 0.007}
       ]
      }
    ];

  Export[
   "fig_emp_schedule_" <> tag <> ".pdf",
   fig
   ];

  fig
  ];


Export["fig_emp_schedule.pdf", figSchedule];



figSchedule9509 =
 makeSchedule[
  calib9509,
  base9509,
  {1995, 2002, 2009},
  "1995_2009"
  ];

figSchedule0014 =
 makeSchedule[
  calib0014,
  base0014,
  {2000, 2007, 2014},
  "2000_2014"
  ];


(* -- Figure VI: frontier vs estimate -------------------------------------- 
abGrid = Range[0.80, 0.965, 0.005];
frontierCurve = {#, thetaMin[#]} & /@ abGrid;

figDiagnostic = Show[
   ListLinePlot[frontierCurve,
    PlotStyle -> Directive[Black, Thickness[0.005]],
    Filling -> Top, FillingStyle -> RGBColor[0.95, 0.89, 0.89]],
   figStyle,
   FrameLabel -> {"returns to scale \[Alpha]+\[Beta]",
     "enforcement elasticity \[Theta]"},
   PlotRange -> {{0.80, 0.965}, {0, 1.6}},
   Epilog -> {
     {RGBColor[0.85, 0.90, 0.96], Opacity[0.55],
      Rectangle[{0.80, basePW["CI"][[1]]}, {0.965, basePW["CI"][[2]]}]},
     {RGBColor[0.12, 0.31, 0.47], Thickness[0.0035],
      Line[{{0.80, basePW["theta"]}, {0.965, basePW["theta"]}}]},
     Text[Style["multiple equilibria", 11, Italic], {0.86, 1.3}],
     Text[Style["95% CI", 9, GrayLevel[0.35]], {0.812, 0.50}],
     Text[Style[
       "\!\(\*OverscriptBox[\(\[Theta]\), \(^\)]\)=0.353", 10,
       RGBColor[0.12, 0.31, 0.47]], {0.948, basePW["theta"] + 0.06}]}]*)
       
       abGrid = Range[0.80, 0.965, 0.005];
frontierCurve = {#, thetaMin[#]} & /@ abGrid;

makeDiagnostic[base_, tag_] :=
 Module[{xTicks, yTicks, fig},

  xTicks =
   Table[{x, NumberForm[x, {3, 2}]}, {x, 0.80, 0.96, 0.02}];

  yTicks =
   Table[{y, NumberForm[y, {2, 1}]}, {y, 0, 1.6, 0.2}];

  fig =
   Show[

    (* Multiplicity region *)
    ListLinePlot[
     frontierCurve,
     PlotStyle ->
      Directive[Black, Thickness[0.0035]],
     Filling -> Top,
     FillingStyle -> RGBColor[0.95, 0.89, 0.89],
     PlotRange -> {{0.80, 0.965}, {0, 1.6}}
     ],

    (* Confidence interval *)
    Graphics[
     {
      RGBColor[0.85, 0.90, 0.96],
      Opacity[0.65],
      Rectangle[
       {0.80, base["CI"][[1]]},
       {0.965, base["CI"][[2]]}
       ]
      }
     ],

    (* theta point estimate *)
    Graphics[
     {
      RGBColor[0.12, 0.31, 0.47],
      Thickness[0.003],
      Line[
       {
        {0.80, base["theta"]},
        {0.965, base["theta"]}
        }
       ]
      }
     ],

    paperFigStyle,

    AspectRatio -> 0.67,
    ImageSize -> 500,

    FrameTicks -> {
      {yTicks, None},
      {xTicks, None}
      },

    FrameLabel -> {
      "returns to scale \[Alpha]+\[Beta]",
      "enforcement elasticity \[Theta]"
      },

    PlotRange -> {
      {0.80, 0.965},
      {0, 1.6}
      },

    Epilog -> {

      Text[
       Style["multiple equilibria", 10, Italic],
       {0.865, 1.27}
       ],

    Text[
 Style["95% CI", 8, GrayLevel[0.35]],
 {
  0.812,
  base["theta"] +
   0.65 (base["CI"][[2]] - base["theta"])
  }
 ],

      Text[
       Style[
        "\!\(\*OverscriptBox[\(\[Theta]\), \(^\)]\) = " <>
         ToString[NumberForm[base["theta"], {4, 3}]],
        9,
        RGBColor[0.12, 0.31, 0.47]
        ],
       {0.948, base["theta"] + 0.035}
       ],

      Inset[
       LineLegend[
        {Black},
        {"multiplicity frontier \!\(\*SubscriptBox[\(\[Theta]\), \(min\)]\)"},
        LegendMarkerSize -> 22,
        LabelStyle -> Directive[
          FontFamily -> "Latin Modern Roman",
          FontSize -> 8
          ]
        ],
       Scaled[{0.78, 0.94}]
       ]
      }
    ];

  Export[
   "fig_emp_diagnostic_" <> tag <> ".pdf",
   fig
   ];

  fig
  ];


Export["fig_emp_diagnostic.pdf", figDiagnostic];


figDiagnostic9509 =
 makeDiagnostic[base9509, "1995_2009"];

figDiagnostic0014 =
 makeDiagnostic[base0014, "2000_2014"];


(* -- Figure VII: (b, alpha+beta) plane ------------------------------------ *)
contourAB[b_, target_] := Module[{lo = 0.80, hi = 0.985, mid, f},
   f[ab_] := With[{t = thetaMin[ab, <|empPar, "b" -> b|>]},
     If[MissingQ[t], Infinity, t]];
   If[f[lo] < target, Return[lo]];
   If[f[hi] > target, Return[Missing["None"]]];
   While[hi - lo > 0.003, mid = (lo + hi)/2.;
    If[f[mid] > target, lo = mid, hi = mid]];
   (lo + hi)/2.];

(*bGrid = Range[0.05, 0.95, 0.05];
cHat = {#, contourAB[#, basePW["theta"]]} & /@ bGrid;
cCI  = {#, contourAB[#, basePW["CI"][[2]]]} & /@ bGrid;

figTwoD = Show[
   ListLinePlot[{DeleteMissing[cHat, 1, 2], DeleteMissing[cCI, 1, 2]},
    PlotStyle -> {Directive[Black, Thickness[0.005]],
      Directive[Black, Dashed, Thickness[0.003]]},
    Filling -> {1 -> Top}, FillingStyle -> RGBColor[0.95, 0.89, 0.89],
    PlotLegends -> Placed[
      {"frontier at \!\(\*OverscriptBox[\(\[Theta]\), \(^\)]\)=0.353",
       "frontier at CI upper bound"}, {0.72, 0.16}]],
   figStyle,
   FrameLabel -> {"integration gains b", "returns to scale \[Alpha]+\[Beta]"},
   PlotRange -> {{0.05, 0.95}, {0.80, 0.99}},
   Epilog -> {Text[Style["multiple equilibria", 10, Italic], {0.22, 0.955}]}]*)
  makeTwoD[base_, tag_] :=
 Module[{bGrid, cHat, cCI, xTicks, yTicks, fig},

  bGrid = Range[0.05, 0.95, 0.05];

  cHat =
   {#, contourAB[#, base["theta"]]} & /@ bGrid;

  cCI =
   {#, contourAB[#, base["CI"][[2]]]} & /@ bGrid;

  xTicks =
   Table[
    {x, NumberForm[x, {2, 1}]},
    {x, 0.1, 0.9, 0.1}
    ];

  yTicks =
   Table[
    {y, NumberForm[y, {4, 3}]},
    {y, 0.80, 0.975, 0.025}
    ];

  fig =
   Show[

    ListLinePlot[
     {
      DeleteMissing[cHat, 1, 2],
      DeleteMissing[cCI, 1, 2]
      },

     PlotStyle -> {
       Directive[Black, Thickness[0.0035]],
       Directive[Black, Dashed, Thickness[0.0025]]
       },

     Filling -> {1 -> Top},

     FillingStyle ->
      RGBColor[0.95, 0.89, 0.89],

     PlotRange -> {
       {0.05, 0.95},
       {0.80, 0.99}
       },

     PlotLegends ->
      Placed[
       {
        "frontier at \!\(\*OverscriptBox[\(\[Theta]\), \(^\)]\) = " <>
         ToString[NumberForm[base["theta"], {4, 3}]],

        "frontier at CI upper bound"
        },
       {0.72, 0.16}
       ]
     ],

    paperFigStyle,

    AspectRatio -> 0.67,
    ImageSize -> 480,

    FrameTicks -> {
      {yTicks, None},
      {xTicks, None}
      },

    FrameLabel -> {
      "integration gains b",
      "returns to scale \[Alpha]+\[Beta]"
      },

    PlotRange -> {
      {0.05, 0.95},
      {0.80, 0.99}
      },

    Epilog -> {
      Text[
       Style[
        "multiple equilibria",
        9,
        Italic
        ],
       {0.22, 0.955}
       ]
      }
    ];

  Export[
   "fig_emp_2d_" <> tag <> ".pdf",
   fig
   ];

  fig
  ];


Export["fig_emp_2d.pdf", figTwoD];

Print["Exported: fig_emp_schedule.pdf, fig_emp_diagnostic.pdf, \
fig_emp_2d.pdf"];


figTwoD9509 =
 makeTwoD[base9509, "1995_2009"];

figTwoD0014 =
 makeTwoD[base0014, "2000_2014"];


(* ======================= PART C: EXPLORER ================================ *)

(* Interactive phase diagram. Coarser grids for responsiveness; use Part A
   functions for publication numbers.                                       *)
Manipulate[Module[{p, eqs},
   p = <|"alpha" -> a, "beta" -> bb, "sigma" -> s, "kappa" -> k,
     "lam" -> l, "b" -> bpar, "delta" -> d, "chi0" -> c0, "theta" -> th,
     "ebar" -> eb|>;
   eqs = Quiet@fixedPoints[p, 40, 1500];
   Plot[{GG[Om, p], Om}, {Om, 0.01, 15},
    PlotStyle -> {Blue, {Black, Dashed}}, PlotRange -> {0, 15},
    AxesLabel -> {"\[CapitalOmega]t", "\[CapitalOmega]t+1"},
    Epilog -> {PointSize[Large], Point[{#, #}] & /@ eqs},
    PlotLabel -> Row[{"equilibria: ", Round[eqs, .01]}]]],
  {{a, .5, "\[Alpha]"}, .2, .7}, {{bb, .25, "\[Beta]"}, .1, .45},
  {{s, 1.25, "\[Sigma]=sC/sP"}, .5, 2}, {{k, 1, "\[Kappa]"}, .5, 2},
  {{l, .5, "\[Lambda]"}, 0, 1}, {{bpar, .5, "b"}, 0, 2},
  {{d, .5, "\[Delta]"}, .1, 1}, {{c0, 202., "\[Chi]0"}, 5, 2500},
  {{th, 3., "\[Theta]"}, 0, 5}, {{eb, .5, "ebar"}, .1, .9}]

