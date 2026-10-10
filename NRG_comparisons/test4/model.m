(* ASQ_variational: physical L/R basis, Gamma = GammaL + GammaR.
   Snapshot this file and modeloperators.m in the calculation directory.
   WN and WS are energies, not rho W. bcsgap is a signed pairing coefficient. *)

def2ch[1];

If[ExpandFileName[$InputFileName] =!= FileNameJoin[{Directory[], "model.m"}],
  MyError["ASQ requires the calculation-local model.m snapshot."];
];
If[BAND =!= "flat" || DISCRETIZATION =!= "Z" || SYMTYPE =!= "SPU1" ||
   POLARIZED || POL2x2 || RUNGS || paramdefaultbool["substeps", False] || Ninit =!= 0,
  MyError["ASQ supports only flat/Z/SPU1, Ninit=0, no polarization, rungs or substeps."];
];
If[!finiteRealNumberQ[lambda] || !TrueQ[lambda >= 4],
  MyError["ASQ requires Lambda >= 4."];
];
If[VARIANT =!= "" || Sort[loptions] =!= {"COMPLEX", "SCHUR"},
  MyError["ASQ requires options=COMPLEX SCHUR and no variant."];
];
If[Sort[First /@ listdata["extra"]] =!= Sort[{"Bz", "WN", "WS", "nu", "phi"}],
  MyError["ASQ extra fields must be exactly nu, Bz, WN, WS, phi."];
];
If[!(And @@ (finiteRealNumberQ[paramnum[#, "extra"]] & /@ {"nu", "Bz", "WN", "WS", "phi"})),
  MyError["ASQ parameters must be finite and real."];
];
If[!positiveFiniteRealNumberQ[bandrescale] ||
   !finiteRealNumberQ[realU] || !TrueQ[realU >= 0] ||
   !finiteRealNumberQ[realGamma] || !TrueQ[realGamma >= 0] ||
   !paramexists["bcsgap"] || !TrueQ[paramnum["bcsgap"] < 0] ||
   !finiteRealNumberQ[paramnum["bcsgap"]] ||
   (Or @@ (paramexists[#] & /@ {"bcsgap1", "bcsgap2", "gap", "hook_bcs", "shift0"})),
  MyError["ASQ requires U,Gamma >= 0, D > 0 and one negative bcsgap=-Delta/D for both leads."];
];

(* Keep these symbolic: initial.m appends the extra-field replacement rules
   AFTER loading this model. OwnValues would destroy the phase derivative. *)
snegrealconstants[nu, Bz, WN, WS, phi];
H1 = U/2 pow[number[d[]] - nu, 2] + Bz spinz[d[]];
Hdotlead = (gammaPol/Sqrt[2]) (
  hopphi[f[0], d[], -phi/4] + hopphi[f[1], d[], phi/4]);
Hbridge = genhop[(WN + I WS) Exp[-I phi/2], f[0], f[1], UP] +
          genhop[(WN - I WS) Exp[-I phi/2], f[0], f[1], DO];
Hc = Hdotlead + Hbridge;
H = H0 + H1 + Hc;
asqCurrentSymbolic = D[Hdotlead + Hbridge, phi];

MyPrint["ASQ_variational physical L/R model: ", ExpandFileName[$InputFileName]];
MyPrint["Gamma is TOTAL; gammaPol=sqrt(2 D Gamma/Pi), each lead gammaPol/sqrt(2)."];
MyPrint["Dot energy is exactly U/2 (n-nu)^2 + Bz Sz; WN/WS are dimensional."];
MyPrint["Both pairing coefficients are -Delta/D; physical Delta is positive."];
MyPrint["asqJ=dH/dphi; asqJ/Delta gives I in units 2 e Delta/hbar."];
