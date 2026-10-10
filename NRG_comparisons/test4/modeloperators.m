(* Return an operator table, not Null. No ambient current-operator dependency. *)
If[calcopq["asqJ"],
  If[ExpandFileName[$InputFileName] =!= FileNameJoin[{Directory[], "modeloperators.m"}],
    MyError["ASQ requires the calculation-local modeloperators.m snapshot."];
  ];
  MyPrint["ASQ scalar operator asqJ from ", ExpandFileName[$InputFileName]];
  Module[{asqNumericCurrent = asqCurrentSymbolic /. params},
    Block[{MPVCFAST = False},
      mtSingletOp["asqJ", asqNumericCurrent]
    ]
  ],
  {}
]
