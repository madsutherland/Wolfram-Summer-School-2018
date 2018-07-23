(* ::Package:: *)

(* Wolfram Language Package *)


(* $Id: ForceFieldMolecule.m,v 1.3 2017/11/05 13:10:56 Bob Exp $ *)

(* :Title: ForceFieldMolecule.m *)

(* :Context: ForceFieldMolecule` *)

(* :Author: Robert B. Nachbar *)

(* :Summary:
   Functions for generating molecular models with 3D coordinates, reporting geometries, 
     visualization.
 *)

(* :Copyright: (C) 2017 by Robert B. Nachbar *)

(* :Package Version: $Revision: 1.3 $ *)

(* :Mathematica Version: 11.2 *)

(* :History:
*)

(* :Keywords: template, skeleton, package *)

(* :Sources:
*)

(* :Warnings:
   <description of global effects, incompatibilities>
*)

(* :Limitations:
   <special cases not handled, known problems>
*)

(* :Discussion:
   <description of algorithm, information for experts>
*)

(* :Requirements:
*)

(* :Examples:
   <sample input that demonstrates the features of this package>
*)


(* set up the package context, including public imports *)

BeginPackage["ForceFieldMolecule`"]
ClearAll[Evaluate[Context[]<>"*"]]
(* Exported symbols added here with SymbolName::usage *)

Molecule::usage = "Molecule is the Head of a WL molecule."; 
Bond::usage = "Bond is the Head of a WL Molecule bond."; 

FromMolecule::usage = "FromMolecule[mol, newProperties] returns a force field molecule \
Association given a WL Molecule mol. newProperties is a sequence of Rules."; 

ToMolecule::usage = "ToMolecule[mol] retruns a WL Molecule given a force field molecule \
Association mol."; 

\[ScriptCapitalR]::usage = "\[ScriptCapitalR] represents a symbolic bond length."
\[ScriptCapitalA]::usage = "\[ScriptCapitalA] represents a symbolic bond angle."

MakeRingMolecule::usage = 
  "MakeRingMolecule[q\[Phi], params] returns the atoms, coordinates, bonds, \
angles, torsions, nonbonds, and cross terms for an n-membered \
aliphatic cyclic molecule given the out-of-plane distortion \
coordinates q\[Phi], which has the form \
{\!\(\*SubscriptBox[\(\[Phi]\), \(0\)]\),\!\(\*SubscriptBox[\(q\), \
\(1\)]\),\!\(\*SubscriptBox[\(\[Phi]\), \
\(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\
\(\[Phi]\), \(2\)]\),\[Ellipsis]}, where n is the length of q\[Phi]. \
Reference bond lengths and angles are taken from the list params. 
MakeRingMolecule[q\[Phi], {}] will leave the molecular structure \
and energy parameters in symbolic form.";

MakeZmatrixMolecule::usage = 
  StringForm[
    "MakeZmatrixMolecule[Z, params] returns the atoms, coordinates, bonds, \
angles, torsions, nonbonds, and cross terms for a hydrocarbon \
molecule given the input Z-matrix. Z has the form ``, where \
\[ScriptI], \[ScriptJ], \[ScriptK], and \[ScriptM] are unique indices \
numbering the atoms, a is an atomic symbol (String), b is a bond \
order (\"\[ScriptS]\", \"\[ScriptD]\", or \"\[ScriptT]\"), \
r is a bond distance, \[Theta] is an angle in \
degrees, \[Phi] is a signed torsion angle in degrees, and i, j, and k \
are indicies of atoms sequentially defining the distance, angle, and \
torsion, respectively. \
Reference bond lengths and angles are taken from the list params. 
MakeZmatrixMolecule[Z, {}] will leave the \
molecular structure and energy parameters in symbolic form.", 
    {
       {Subscript["a", "\[ScriptI]"], "\[ScriptI]", 0, 0, 0, 0, 0, 0, 0},
       {Subscript["a", "\[ScriptJ]"], "\[ScriptJ]", 
       Subscript["i", "\[ScriptJ]"], 0, 0, 
       Subscript["b", "\[ScriptJ]", Subscript["i", "\[ScriptJ]"]], 
       Subscript["r", "\[ScriptJ]", Subscript["i", "\[ScriptJ]"]], 0, 0},
       {Subscript["a", "\[ScriptK]"], "\[ScriptK]", 
       Subscript["i", "\[ScriptK]"], Subscript["j", "\[ScriptK]"], 0, 
       Subscript["b", "\[ScriptK]", Subscript["i", "\[ScriptK]"]], 
       Subscript["r", "\[ScriptK]", Subscript["i", "\[ScriptK]"]], 
       Subscript["\[Theta]", "\[ScriptK]", Subscript["i", "\[ScriptK]"], 
        Subscript["j", "\[ScriptK]"]], 0},
       {Subscript["a", "\[ScriptL]"], "\[ScriptL]", 
       Subscript["i", "\[ScriptL]"], Subscript["j", "\[ScriptL]"], 
       Subscript["k", "\[ScriptL]"], 
       Subscript["b", "\[ScriptL]", Subscript["i", "\[ScriptL]"]], 
       Subscript["r", "\[ScriptL]", Subscript["i", "\[ScriptL]"]], 
       Subscript["\[Theta]", "\[ScriptL]", Subscript["i", "\[ScriptL]"], 
        Subscript["j", "\[ScriptL]"]], 
       Subscript["\[Phi]", "\[ScriptL]", Subscript["i", "\[ScriptL]"], 
        Subscript["j", "\[ScriptL]"], Subscript["k", "\[ScriptL]"]]},
       {Subscript["a", "\[ScriptM]"], "\[ScriptM]", 
       Subscript["i", "\[ScriptM]"], Subscript["j", "\[ScriptM]"], 
       Subscript["k", "\[ScriptM]"], 
       Subscript["b", "\[ScriptM]", Subscript["i", "\[ScriptM]"]], 
       Subscript["r", "\[ScriptM]", Subscript["i", "\[ScriptM]"]], 
       Subscript["\[Theta]", "\[ScriptM]", Subscript["i", "\[ScriptM]"], 
        Subscript["j", "\[ScriptM]"]], 
       Subscript["\[Phi]", "\[ScriptM]", Subscript["i", "\[ScriptM]"], 
        Subscript["j", "\[ScriptM]"], Subscript["k", "\[ScriptM]"]]},
       {"\[VerticalEllipsis]", "\[VerticalEllipsis]", "\[VerticalEllipsis]", 
       	"\[VerticalEllipsis]", "\[VerticalEllipsis]", "\[VerticalEllipsis]", 
       	"\[VerticalEllipsis]", "\[VerticalEllipsis]", "\[VerticalEllipsis]"}
      } // MatrixForm] // ToString[#, StandardForm] &;

AddHydrogens::usage = "AddHydrogens[mol, params] returns a new molecule with \
vacant valences filled with hydrogen atoms."; 

AddConstraints::usage = "AddConstraints[mol, newConstraints] replaces the constraints \
in the molecule mol with newConstraints, which has the form \
{{type, atoms, value, forceConstant}...} where type is the type of \
constraint (\"distance\", \"angle\", \"torsion\", \"dihedral\"), atoms is a list of atoms defining \
the distance, angle, torsion, or dihedral constraint, and value is the desired value. 
AddConstraints[mol, {}] effectively removes all constraints"; 

ShowMolecule::usage = "ShowMolecule[mol] returns a Graphics3D rendering of the molecule mol."; 

DisplayMolecule::usage = "ShowMolecule[mol] returns a Graphics3D rendering of the molecule mol \
in a convenient Manipulate."; 

StructureReport::usage = 
"StructureReport[mol] returns a report of the structure of the \
molecule mol including the cartesian coordinates, bond \
lengths, angles, torsions, and short nonbonded distance.";

Distance::usage := "Distance[mol, i, j] returns the Euclidean distance \
between atoms i and j in molecule mol."

Angle::usage := "Angle[mol, i, j, k] returns the plane angle (in degrees) \
subtended by atoms i-j-k in molecule mol."

Torsion::usage := "Torsion[mol, i, j, k, l] returns the signed torsion angle (in degrees) \
subtended by atoms i-j-k-l in molecule mol. The planes i-j-k and j-k-l define the angle."

Dihedral::usage := "Dihedral[mol, i, j, k, l] returns the signed out-of-plane angle (in degrees) \
subtended by atoms j-i(-k)-l molecule mol. The planes j-i-k and j-i-l define the angle."

PointGroupSymmetry::usage = "PointGroupSymmetry[mol] returns the point group symmetry of \
molecule mol in the Sch\[ODoubleDot]nflies notation.";

Schoenflies::usage = "Schoenflies[mol]."

AddHydrogens::bndord = "Bond orders `` must be a list of \"\[ScriptS]\", \"\[ScriptD]\", or \"\[ScriptT]\"."; 

AddHydrogens::badhyb = "Hybridization for atoms `` could not be computed."; 

ShowMolecule::nncoord = "Nonnumeric value encountered in AtomCoordinates."; 

Begin["`Private`"] (* Begin Private Context *) 
ClearAll[Evaluate[Context[]<>"*"]]

norm[x_] := x/Sqrt[x.x]

normC = Compile[{{a, _Real, 1}}, 
	a/Sqrt[a.a]
	]

cross[{a1_, a2_, a3_}, {b1_, b2_, b3_}] := 
	{-a3 b2 + a2 b3, 
	a3 b1 - a1 b3, 
	-a2 b1 + a1 b2}

crossC = Compile[{{a, _Real, 1}, {b, _Real, 1}}, 
	{-a[[3 ]] b[[2]] + a[[2 ]] b[[3]], 
	a[[3]] b[[1]] - a[[1]] b[[3]], 
	-a[[2]] b[[1]] + a[[1]] b[[2]]}
	]


FromMolecule[mol_Molecule, newProperties___Rule] :=
	Module[{atoms, bonds, properties, coords, newMol},
		{atoms, bonds} = List @@ mol[[{1, 2}]];
		properties = If[Length[mol] == 3, Association@mol[[3]], {}];
  
		bonds = Replace[bonds, 
			{
			Bond[a : {_, _}] :> {a, "\[ScriptS]"}, 
			Bond[a : {_, _}, t_] :> {a, Switch[t, 
				"Single", "\[ScriptS]", 
				"Double", "\[ScriptD]", 
				"Triple", "\[ScriptT]", 
				_, Missing[t]]}
			}, {1}];
  
		newMol = <|"Atoms" -> atoms, "Bonds" -> bonds|>;
  
		coords = properties["AtomPositions"];
		If[MatchQ[coords, _?MatrixQ], 
			newMol["AtomCoordinates"] = coords/100. (* pm -> Angstrom *)
			];
  
		Join[newMol, <|newProperties|>]
		]

ToMolecule[mol_Association] := 
	Module[{atoms, bonds, coords, properties}, 
		atoms = mol["Atoms"]; 
		bonds = mol["Bonds"]; 
		
		bonds = Replace[bonds, 
			{
			{{i_, j_}, "\[ScriptS]"} :> Bond[{i, j}], 
			{{i_, j_}, "\[ScriptD]"} :> Bond[{i, j}, "Double"], 
			{{i_, j_}, "\[ScriptT]"} :> Bond[{i, j}, "Triple"]
			}, {1}];
		
		coords = mol["AtomCoordinates"];
		If[MatchQ[coords, _?MatrixQ], 
			coords *= 100. (* Angstrom -> pm *)
			];

		properties = DeleteCases[Normal @ mol, ("Atoms" | "Bonds" | "AtomCoordinates") -> _]; 
		properties = Prepend[properties, "AtomPositions" -> coords]; 

		Molecule[atoms, bonds, properties]
		]

qPhiRules[qPhi_] := 
	With[{n = Length[qPhi]}, 
		Thread[Table[Piecewise[{{Subscript[\[Phi], (i - 1)/2], OddQ[i]}, 
			{Subscript[q, i/2], True}}], {i, n}] -> qPhi]]

ringCarbons[qPhi_] := 
	Module[{n = Length[qPhi], r}, 
		r = (1/2)*Subscript[\[ScriptCapitalR], "CC"]*Csc[Pi/n]; 
		Piecewise[{
			{Table[{r*Cos[((2*Pi)*(j - 1))/n], r*Sin[((2*Pi)*(j - 1))/n], 
				Sqrt[2/n]*Sum[Subscript[q, m]*Cos[Subscript[\[Phi], m] + (2*Pi*m*(j - 1))/n], 
				{m, 2, (n - 1)/2}]}, {j, n}], OddQ[n]}, 
			{Table[{r*Cos[((2*Pi)*(j - 1))/n], r*Sin[((2*Pi)*(j - 1))/n], 
				Sqrt[2/n]*Sum[Subscript[q, m]*Cos[Subscript[\[Phi], m] + (2*Pi*m*(j - 1))/n], 
				{m, 2, n/2 - 1}] + (Subscript[q, n/2]*(-1)^(j - 1))/Sqrt[n]}, {j, n}], True}
			}] /. qPhiRules[qPhi]
		]

ringHydrogens[c_] := 
 Module[{Cijk = Partition[c, 3, 1, {2, 2}], Vji, Vjk, norm, 
   axis, alpha = 1/2 Subscript[\[ScriptCapitalA], "HCH"]}, 
  Vji = Normalize /@ (#[[1]] - #[[2]] &) /@ Cijk; 
  Vjk = Normalize /@ (#[[3]] - #[[2]] &) /@ Cijk; 
  norm = Normalize /@ MapThread[cross, {Vji, Vjk}]; 
  axis = Normalize /@ (-(Vji + Vjk)); 
  MapThread[{#1 + 
      Subscript[\[ScriptCapitalR], "CH"] Normalize[
        Cos[alpha] #2 + Sin[alpha] #3], #1 + 
      Subscript[\[ScriptCapitalR], "CH"] Normalize[
        Cos[alpha] #2 - Sin[alpha] #3]} &, {c, axis, norm}]]

(* TODO: "AtomTypes" *)
Options[MakeRingMolecule] = {
	"BondOrders" -> Automatic
	};

MakeRingMolecule[qPhi_, params_, opts:OptionsPattern[]] := 
	Module[{n = Length[qPhi], c, h, bc, bh, coords, atoms, bonds, bondOrders}, 
		c = ringCarbons[qPhi] //. params; 
		atoms = ConstantArray["C", Length @ c]; 
		bc = Partition[Range[n], 2, 1, {1, 1}] // Thread[{#, "\[ScriptS]"}] &;

		bondOrders = OptionValue["BondOrders"]; 
		bc = Switch[bondOrders, 
			Automatic, 
				bc, 
			{("\[ScriptS]" | "\[ScriptD]" | "\[ScriptT]")..}, 
(* TODO: check that bondOrders has correct length *)
				Thread[{First /@ bc, bondOrders}], 
			{({_, _} -> ("\[ScriptS]" | "\[ScriptD]" | "\[ScriptT]"))...}, 
(* TODO: check that atom indices are in range *)
				Fold[With[{ij=First@#2}, Replace[#1, {a:(ij | Reverse@ij), _} :> {a, Last@#2}, {1}]]&, 
					bc, bondOrders], 
			_, 
				Message[MakeRingMolecule:erropts, bondOrders, "BondOrders"];
				Return[$Failed]
			];

(*
		h = Join @@ ringHydrogens[c] //. forceFieldParameters; 
		bh = Join @@ Thread /@ Thread@{Range[n], Partition[Range[n + 1, n + 2 n], 2, 2]} // 
			Thread[{#, "\[ScriptS]"}] &; atoms = Join["C" & /@ c, "H" & /@ h]; 

		coords = Join[c, h]; 
		bonds = Join[bc, bh]; 
*)
		coords = c; 
		bonds = bc; 

		AddHydrogens[<| "Atoms" -> atoms, "Bonds" -> bonds, "AtomCoordinates" -> coords |>, params]
		]


ZmatrixAtoms[Zin_, params_:{}] := 
 	Module[{A = Zin[[All, 1]], B = Zin[[All, 2 ;; 6]], 
   Z = N[Zin[[All, -3 ;; All]]], 
   			bt, z, a, b, c, d, t, r, theta, phi, bc, n, M}, 
  		(* need to check B and Z for proper structure *)
  		
  Z[[All, {2, 3}]] = (Pi/180.)*Z[[All, {2, 3}]]; 
  		z = Z; 
  		If[Length[B] >= 1, 
   			{d, c, b, a, t} = B[[1]]; 
   			{r, theta, phi} = z[[1]]; 
   			Z[[d]] = {0., 0., 0.}
   			]; 
  		If[Length[B] >= 2, 
   			{d, c, b, a, t} = B[[2]]; 
   			{r, theta, phi} = z[[2]]; 
   			Z[[d]] = Z[[c]] + {-r, 0., 0.}
   			]; 
  		If[Length[B] >= 3, 
   			{d, c, b, a, t} = B[[3]]; 
   			{r, theta, phi} = z[[3]]; 
   			bc = Z[[c]] - Z[[b]] // norm; 
   			n = cross[Z[[b]] - {0., 1., 0.}, bc] // norm;  
   			M = {bc, cross[n, bc], n} // Transpose; 
   			Z[[d]] = Z[[c]] + M.{-r Cos[theta], r Sin[theta], 0}
   			]; 
  		Do[
   			{d, c, b, a, t} = B[[i]]; 
   			Which[
    (*
    				A[[i]] \[Equal] Null, 
    					, 
    *)
    				b > 0 && a > 0, 
    					{r, theta, phi} = z[[i]]; 
    					bc = Z[[c]] - Z[[b]] // norm; 
    					n = cross[Z[[b]] - Z[[a]], bc] // norm; 
    					M = {bc, cross[n, bc], n} // Transpose; 
    					Z[[d]] = 
     Z[[c]] + 
      M.{-r Cos[theta], r Cos[phi] Sin[theta], 
        r Sin[phi] Sin[theta]}, 
    				True, 
    					Z[[i]] = Null
    				], 
   			{i, 4, Length[B]}
   			]; 
 
		{DeleteCases[A, Null], 
		Thread[{B[[All, {1, 2}]], 
		B[[All, 5]](* /. {1 -> "\[ScriptS]", 2 -> "\[ScriptD]", 3 -> "\[ScriptT]"}*)}] // 
			DeleteCases[#, {{_, 0}, _}] &, 
   		DeleteCases[Z, Null]} //. params
  		]

Options[MakeZmatrixMolecule] = {
	"RingBonds" -> {}
	};

MakeZmatrixMolecule[Z_, params_:{}, OptionsPattern[]] := 
	Module[{atoms, coords, bonds, ringBonds}, 
		ringBonds = OptionValue["RingBonds"]; 
		ringBonds = Switch[ringBonds, 
			{({_, _} -> ("\[ScriptS]" | "\[ScriptD]" | "\[ScriptT]"))...}, 
				ringBonds, 
			_, 
				Message[MakeZmatrixMolecule::optvg, "RingBonds", ringBonds, 
					"a list of rules of the form {_, _} -> (\"\[ScriptS]\" | \"\[ScriptD]\" | \"\[ScriptT]\")."]; 
				Return[$Failed]
			];

		{atoms, bonds, coords} = (*Echo @ *)ZmatrixAtoms[Z, params]; 
		AddHydrogens[<| "Atoms" -> atoms, "Bonds" -> Join[bonds, List @@@ ringBonds], "AtomCoordinates" -> coords |>, params]
		]


hypridizationRules = {
	{"P", {"\[ScriptD]", "\[ScriptS]", "\[ScriptS]", "\[ScriptS]"}} | 
		{"S", {"\[ScriptD]", "\[ScriptD]", "\[ScriptS]", "\[ScriptS]"}} -> "sp3", 
	{"C" | "N" | "O" | "F" | "Si" | "P" | "S" | "Cl" | "Br", {"\[ScriptS]" ..}} -> "sp3", 
	{"C" | "N" | "O" | "Si" | "P" | "S", {"\[ScriptD]", "\[ScriptS]" ...}} -> "sp2", 
	{"C" | "N" | "Si" | "P", {"\[ScriptS]"..., "\[ScriptT]"}} -> "sp", 
	{"C" | "Si", {"\[ScriptD]", "\[ScriptD]"}} -> "sp", 
	{"H", {"\[ScriptS]"}} -> "s"
	};

debugEcho[expr_, args___] := If[TrueQ[Global`$debug], Echo[expr, args], expr]

expr : AddHydrogens[mol_Association, params_:{}] /; Length[mol["AtomCoordinates"]] >= 3 := 
	Module[{atomsC, bondsC, coordsC, nC, atomsH, bondsH, coordsH, nH, 
		bo, hybrid, adj, k, Vji, Vjk, normal, axis, 
		alpha = N[(1/2)*ArcCos[-1/3]], beta = N[(1/2)*ArcCos[-1/2]], 
		p}, 

		atomsC = mol["Atoms"];
		bondsC = mol["Bonds"]; 
		coordsC = mol["AtomCoordinates"]; 
		nC = Length[coordsC]; 

  		bo = Join @@ Thread /@ bondsC // SortBy[#, First] & // SplitBy[#, First] & // 
  			Map[Last, #, {2}] & // Map[Sort, #] &(* // Echo*);
  		If[!MatchQ[bo, {{("\[ScriptS]" | "\[ScriptD]" | "\[ScriptT]")...}..}], 
  			Message[AddHydrogens::bndord, bo];
  			Return[HoldForm[expr]]
  			];
		hybrid = debugEcho[#, "hybrid"]& @ Replace[debugEcho @ Thread[{atomsC, bo}], hypridizationRules, {1}]; 
		If[!MatchQ[hybrid, {("s" | "sp" | "sp2" | "sp3")..}], 
			p = First /@ Position[hybrid, Except["s" | "sp" | "sp2" | "sp3"], {1}, Heads -> False]; 
			Message[AddHydrogens::badhyb, Thread[p -> hybrid[[p]]]]; 
			Return[$Failed]
			];
		adj = (*Echo @ *)(Join[bondsC[[All, 1]], Reverse /@ bondsC[[All, 1]]] // 
			Sort // SplitBy[#1, First] & // Map[{#1[[1, 1]], #1[[All, 2]]} &, #1] &); 
			(* need to make sure Union[Join @@ bondsC[[All,{1, 2}]]] == Range[nC] *)
		nH = nC; 

  		{atomsH, bondsH, coordsH} = Join @@ MapThread[Function[{ij, h}, 
      			Switch[debugEcho[#, "h"]& @ h, 
       				"sp3",
       					Switch[debugEcho[#, "length[ij]"]& @ Length[ij[[2]]], 
        						1, 
									k = Cases[adj, {First[ij[[2]]], jList_} :> 
										First[DeleteCases[jList, ij[[1]]]]] // First; 
        							Vji = (coordsC[[ij[[1]]]] - coordsC[[ij[[2, 1]]]]) // norm; 
        							Vjk = (coordsC[[k]] - coordsC[[ij[[2, 1]]]]) // norm; 
        							normal = cross[Vji, Vjk] // norm; 
        							axis = -(Vji + Vjk) // norm; 
        							Switch[debugEcho[#, "atomsC[ij]"]& @ atomsC[[ij[[1]]]], 
         								"C" | "Si", 
         									{{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] - Subscript[\[ScriptCapitalR], "CH"] Vjk}, 
          									{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] + Subscript[\[ScriptCapitalR], "CH"] * 
													norm[(-Cos[alpha]) axis + Sin[alpha] normal]}, 
          									{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] + Subscript[\[ScriptCapitalR], "CH"] * 
													norm[(-Cos[alpha]) axis - Sin[alpha] normal]}}, 
         								"N" | "P", 
         									{{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] - Subscript[\[ScriptCapitalR], "CH"] Vjk}, 
          									{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] + Subscript[\[ScriptCapitalR], "CH"] * 
													norm[(-Cos[alpha]) axis + Sin[alpha] normal]}}, 
         								"O" | "S", 
         									{{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] - Subscript[\[ScriptCapitalR], "CH"] Vjk}}, 
         								"H" | "F" | "Cl" | "Br" | _, 
         									{}
         								], 
        						2, 
        							Vji = norm /@ (coordsC[[ij[[1]]]] - # &) /@ coordsC[[ij[[2]]]]; 
        							normal = cross @@ Vji // norm; 
        							axis = Total[Vji] // norm; 
        							Switch[debugEcho[#, "atomsC[ij]"]& @ atomsC[[ij[[1]]]], 
         								"C" | "Si", 
         									{{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] + Subscript[\[ScriptCapitalR], "CH"] * 
													norm[Cos[alpha] axis + Sin[alpha] normal]}, 
          									{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] + Subscript[\[ScriptCapitalR], "CH"] * 
													norm[Cos[alpha] axis - Sin[alpha] normal]}}, 
         								"N" | "P", 
         									{{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] + Subscript[\[ScriptCapitalR], "CH"] * 
													norm[Cos[alpha] axis + Sin[alpha] normal]}}, 
         								_, 
         									{}
         								], 
        						3, 
									Vji = norm /@ (coordsC[[ij[[1]]]] - # &) /@ coordsC[[ij[[2]]]] // Total // norm; 
        							Switch[debugEcho[#, "atomsC[ij]"]& @ atomsC[[ij[[1]]]], 
         								"C" | "Si", 
         									{{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] + Subscript[\[ScriptCapitalR], "CH"] Vji}}, 
         								_, 
         									{}
         								], 
        						_, 
        							{}
        						],
       				"sp2", 
       					Switch[debugEcho[#, "length[ij]"]& @ Length[ij[[2]]], 
        						1, 
									k = Cases[adj, {First[ij[[2]]], jList_} :> 
										First[DeleteCases[jList, ij[[1]]]]] // First; 
        							Vji = (coordsC[[ij[[1]]]] - coordsC[[ij[[2, 1]]]]) // norm; 
        							Vjk = (coordsC[[k]] - coordsC[[ij[[2, 1]]]]) // norm; 
        							normal = cross[Vji, Vjk] // norm; 
        							axis = -(Vji + Vjk) // norm; 
        							Switch[debugEcho[#, "atomsC[ij]"]& @ atomsC[[ij[[1]]]], 
         								"C" | "Si", 
         									{{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] - Subscript[\[ScriptCapitalR], "CH"] Vjk}, 
          									{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] + Subscript[\[ScriptCapitalR], "CH"] * 
													norm[(-Cos[beta]) axis + Sin[0] normal]}}, 
         								"N" | "P", 
         									{{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] - Subscript[\[ScriptCapitalR], "CH"] Vjk}}, 
         								_, 
         									{}
         								], 
        						2, 
        							Vji = norm /@ (coordsC[[ij[[1]]]] - # &) /@ coordsC[[ij[[2]]]]; 
        							normal = cross @@ Vji // norm; 
        							axis = Total[Vji] // norm; 
        							Switch[debugEcho[#, "atomsC[ij]"]& @ atomsC[[ij[[1]]]], 
         								"C" | "Si", 
         									{{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] + Subscript[\[ScriptCapitalR], "CH"] * (*
													norm[Cos[beta] axis + Sin[beta] normal]*)axis}}, 
         								_, 
         									{}
         						], 
        						_, 
        							{}
        						], 
       				"sp", 
       					Switch[debugEcho[#, "length[ij]"]& @ Length[ij[[2]]], 
        						1, 
									k = Cases[adj, {First[ij[[2]]], jList_} :> 
										First[DeleteCases[jList, ij[[1]]]]] // First; 
        							Vji = (coordsC[[ij[[1]]]] - coordsC[[ij[[2, 1]]]]) // norm; 
        							Vjk = (coordsC[[k]] - coordsC[[ij[[2, 1]]]]) // norm; 
        							normal = cross[Vji, Vjk] // norm; 
        							axis = -(Vji + Vjk) // norm; 
        							Switch[debugEcho[#, "atomsC[ij]"]& @ atomsC[[ij[[1]]]], 
         								"C" | "Si", 
         									{{"H", {{ij[[1]], ++nH}, "\[ScriptS]"}, 
												coordsC[[ij[[1]]]] - Subscript[\[ScriptCapitalR], "CH"] Vjk}}, 
         								_, 
         									{}
         								], 
        						_, 
        							{}
        						], 
       				_ (* "s" *), 
       					{}
       				]
      			], {adj, hybrid}] //. params // Switch[#, {}, {{}, {}, {}}, _, Transpose[#]]&;

		<|"Atoms" -> Join[atomsC, atomsH], "Bonds" -> Join[bondsC, bondsH], "AtomCoordinates" -> Join[coordsC, coordsH]|>
		]

AddConstraints::badtyp = "Constraint type `` in `` is not recognized."; 
AddConstraints::badatm = "Constrained atom `` in `` is not recognized."; 
AddConstraints::badref = "Constrained `` `` in `` is not in allowed range ``."; 
AddConstraints::badfrc = "Constrained force constant `` in `` is not in allowed range ``."; 

(* need more checking for proper syntax *)
AddConstraints[mol_Association, c:{_String, {__}, _?NumberQ, _?Positive}] := 
	AddConstraints[mol, {c}]
AddConstraints[mol_Association, newConstraints:{{_String, {__}, _?NumberQ, _?Positive} ...}] := 
	Module[{na, constraints, ok = True}, 
		na = Length @ mol["Atoms"];
  		Scan[Function[constraint, 
  			Scan[If[#<1 || #>na, 
  				Message[AddConstraints::badatm, #, constraint]; 
  				ok = False]&, constraint[[2]]]; 
  			Switch[constraint[[1]], 
				"distance", If[constraint[[3]] <= 0, 
  					Message[AddConstraints::badref, constraint[[1]], constraint[[3]], constraint, 
  						Subscript["c", 0] > 0]; 
  					ok = False], 
				"angle", If[constraint[[3]] <= 0 || constraint[[3]] > 180, 
  					Message[AddConstraints::badref, constraint[[1]], constraint[[3]], constraint, 
  						0 < Subscript["c", 0] <= 180]; 
  					ok = False], 
				"torsion", If[constraint[[3]] <= -360 || constraint[[3]] > 360, 
  					Message[AddConstraints::badref, constraint[[1]], constraint[[3]], constraint, 
  						-360 < Subscript["c", 0] <= 360]; 
  					ok = False], 
				"dihedral", If[constraint[[3]] <= -180 || constraint[[3]] > 180, 
  					Message[AddConstraints::badref, constraint[[1]], constraint[[3]], constraint, 
  						-180 < Subscript["c", 0] <= 180]; 
  					ok = False], 
				_, Message[AddConstraints::badtyp, constraint[[1]], constraint]; ok = False
				]; 
			If[constraint[[4]] < 0, 
  				Message[AddConstraints::badfrc, constraint[[4]], constraint, 
  					0 <= Subscript["K", "c"]]; 
				ok = False
				];
			], newConstraints]; 

  		constraints = If[ok, newConstraints, {}]; 

  		Append[mol, "Constraints" -> constraints]
  		]

Options[ShowMolecule] = {
	"ShowHydrogens" -> Automatic, 
	"ShowLabels" -> False, 
	"ShowAxes" -> False, 
	"Width" -> 6
	};

ShowMolecule[mol_Association, OptionsPattern[{ShowMolecule, Graphics3D}]] := 
	Module[{atoms, bonds, coords, showHydrogens, showLabels, showLabelsH, showAxes, 
			width, centroid, axisLength}, 
		atoms = mol["Atoms"];
		bonds = mol["Bonds"];
		coords = mol["AtomCoordinates"]; 

		showHydrogens = OptionValue["ShowHydrogens"]; 
		showHydrogens = Switch[showHydrogens, 
			True | False, showHydrogens, 
			Automatic, Count[atoms, Except["H"]] <= 3, 
			_, 
				Message[ShowMolecule::opttfa, "ShowHydrogens", showHydrogens]; 
				Return[$Failed]
			];
		
		showLabels = If[OptionValue["ShowLabels"], 
			{#1, {White, Inset[#2, First@#1 + {0, 0, 0.50}]}}&, 
			#1&, 
				Message[ShowMolecule::opttfa, "ShowLabels", showLabels]; 
				Return[$Failed]
			];
		showLabelsH = If[OptionValue["ShowLabels"], 
			{#1, {White, Inset[#2, First@#1 + {0, 0, 0.35}]}}&, 
			#1&
			];
		width = OptionValue["Width"];

		showAxes = OptionValue["ShowAxes"]; 
		showAxes = Switch[showAxes, 
			True | False, showAxes, 
			Automatic, Count[atoms, Except["H"]] <= 3, 
			_, 
				Message[ShowMolecule::opttfa, "ShowAxes", showAxes]; 
				Return[$Failed]
			];

		If[!MatchQ[coords, {{__?NumberQ}..}], 
			Message[ShowMolecule::nncoord];
			Return[$Failed]
			];

		centroid = Mean @ coords;
		axisLength = Max[-0.75 Subtract @@@ MinMax /@ Transpose @ coords]; 

		Graphics3D[{
			{{ColorData["Atoms", atoms[[#]]], showLabels[Sphere[coords[[#]], 0.35], #]} & /@ 
				Select[Range[Length[atoms]], ! MatchQ[atoms[[#]], "H"] &]}, 
			{ColorData["Atoms", "Ag"], Cylinder[coords[[#]], 0.10] & /@ 
				Select[bonds[[All, 1]], ! MatchQ[atoms[[#]], {_, "H"} | {"H", _}] &]}, 
			If[showHydrogens, 
				{{ColorData["Atoms", "H"], showLabelsH[Sphere[coords[[#]], 0.20], #] & /@ 
					Select[Range[Length[coords]], MatchQ[atoms[[#]], "H"] &]}, 
				{ColorData["Atoms"]["Ag"], Cylinder[coords[[#]], 0.10] & /@ 
					Select[bonds[[All, 1]], MatchQ[atoms[[#]], {_, "H"} | {"H", _}] &]}
				}, 
			  (* else *)
				{}
				], 
			If[showAxes, 
				{White, Thickness[Large], 
					Line[centroid  #& /@ axisLength{{-1, 0, 0}, {1, 0, 0}}], Text["X", axisLength{1, 0, 0}, {-1, -1}], 
					Line[centroid  #& /@ axisLength{{0, -1, 0}, {0, 1, 0}}], Text["Y", axisLength{0, 1, 0}, {-1, -1}], 
					Line[centroid  #& /@ axisLength{{0, 0, -1}, {0, 0, 1}}], Text["Z", axisLength{0, 0, 1}, {-1, -1}]
					}, 
			  (* esle *)
				{}
				]
			}, 
			Lighting -> "Neutral", 
			Background -> ColorData["WebSafe"]["#0066FF"], 
			ImageSize -> width 72 {1, 3/4}, 
			Boxed -> False(*,
			SphericalRegion\[Rule]True,PlotRange\[Rule]Automatic,
			BoxRatios\[Rule]{1,1,1}*)
			]
		]

DisplayMolecule[mol_Association] := DisplayMolecule[{mol}]
DisplayMolecule[mols:{__Association}] :=
	DynamicModule[{molNames, molList, mol, hydrogens = True, labels = False, axes = False, width = 6},
		molNames = MapIndexed[Switch[#1["Name"], _Missing, First@#2, _, #1["Name"]]&, mols]; 
		molList = AssociationThread[molNames -> mols]; 
		
		Manipulate[
			ShowMolecule[molList[mol], "ShowHydrogens" -> hydrogens, 
				 "ShowLabels" -> labels, "ShowAxes" -> axes, "Width" -> width],
   
			Grid[{
				{PopupMenu[Dynamic[mol], molNames], ""},
				{Button[RadioButton[Dynamic[hydrogens], True], hydrogens = ! hydrogens, 
					Appearance -> Frameless], "Hydrogens"},
				{Button[RadioButton[Dynamic[labels], True], labels = ! labels, 
					Appearance -> Frameless], "Labels"},
				{Button[RadioButton[Dynamic[axes], True], axes = ! axes, 
					Appearance -> Frameless], "Axes"},
				{PopupMenu[Dynamic[width], {4, 6, 8, 10}], "Width"}
				},
				Alignment -> {{Right, Left}}
				],
   
			TrackedSymbols :> {mol, hydrogens, labels, axes, width},
			ContentSize -> Dynamic[24 + 72 width],
			SaveDefinitions -> True
			]
		]

makeStructureLists[atoms_, bondsIn_] := 
	Module[{n = Length[atoms], adj, idx, adjTemp, bonds, 
		angles = {}, torsions = {}, oops = {}, nonbonds = {}},

		idx = Range[n];

(*Print["bondsIn = ", bondsIn];*)
		bonds = Sort /@ First /@ bondsIn // Sort; 

		adjTemp = {} & /@ idx;
		Scan[(AppendTo[adjTemp[[#[[1]]]], #[[2]]]; 
			AppendTo[adjTemp[[#[[2]]]], #[[1]]]) &, bonds];
		adjTemp = Sort /@ adjTemp;
		ClearAll[adj];
		MapThread[(adj[#1] = #2) &, {idx, adjTemp}];
(*
Print["adj = ",TableForm[adj/@idx//PadRight[#, Automatic, ""]&,
  TableHeadings -> {idx,Automatic}]];
*)
  
		angles = Join @@ Function[{j}, Insert[#, j, 2] & /@ Sort /@ Subsets[adj[j], {2}]] /@ idx // 
			SortBy[#, #[[{2, 1, 3}]]&]&;

		torsions = Join @@ Function[{j, k}, Join @@ Outer[{#1, j, k, #2} &, DeleteCases[adj[j], k], 
			DeleteCases[adj[k], j]]] @@@ bonds // 
			SortBy[#, #[[{2, 3, 1, 4}]]&]&;
  

(*TODO: use same atom order as in energy function *)
		oops = Function[{j, k}, Module[{il}, il = Sort @ DeleteCases[adj[j], k]; 
				If[Length[il] == 2, {j, k, First@il, Last@il}, Unevaluated@Sequence[]]]] @@@ 
			Join @@ ({#, Reverse /@ #}& @Cases[bondsIn, {jk_, "\[ScriptD]"} :> jk]) // 
			SortBy[#, #[[{1, 2, 3, 4}]]&]&;
  
		nonbonds = Subsets[idx, {2}] // DeleteCases[#, {i_, j_} /; MemberQ[bonds, {i, j} | {j, i}]] & // 
			DeleteCases[#, {i_, j_} /; MemberQ[angles, {i, _, j} | {j, _, i}]] & // Sort;

		{bonds, angles, torsions, oops, nonbonds}
		]

Distance[mol_Association, i_, j_] := 
	With[{X=mol["AtomCoordinates"]}, 
		distance[X[[i]], X[[j]]]
		]

distance[a_, b_] := With[{d = b - a}, Sqrt[d.d]]

Angle[mol_Association, i_, j_, k_] := 
	With[{X=mol["AtomCoordinates"]}, 
		angle[X[[i]], X[[j]], X[[k]]] / Degree
		]

angle[a_, b_, c_] := With[{ba = a - b // norm, bc = c - b // norm}, ArcCos[ba.bc]]

Torsion[mol_Association, i_, j_, k_, l_] := 
	With[{X=mol["AtomCoordinates"]}, 
		torsion[X[[i]], X[[j]], X[[k]], X[[l]]] / Degree
		]

torsion[a_, b_, c_, d_] := Module[{ba = a - b // norm, bc = c - b // norm, cd = d - c // norm, mu, nu}, 
	mu = cross[ba, bc]; 
	nu = cross[-bc, cd]; 
	ArcTan[mu.nu, cross[mu, nu].bc]
	]

Dihedral[mol_Association, i_, j_, k_, l_] := 
	With[{X=mol["AtomCoordinates"]}, 
		dihedral[X[[i]], X[[j]], X[[k]], X[[l]]] / Degree
		]

(*dihedral[a_, b_, c_, d_] := Module[{ab = b - a // norm, ac = c - a // norm, ad = d - a // norm, mu, nu}, 
	mu = cross[ac, ab] // norm; 
	nu = cross[ab, ad] // norm; 
	ArcCos[mu.nu]
	]*)
dihedral[a_, b_, c_, d_] := With[{phi=torsion[c, a, b, d]}, Sign[phi] Pi - phi]

Options[StructureReport] = {
	"ShowHydrogens" -> Automatic
	}; 

StructureReport[mol_Association, OptionsPattern[]] := 
	Module[{atoms, coords, bonds, angles, torsions, oops, nonbonds, 
			constraints, cartesian, internal, rows, includeHydrogens}, 

		includeHydrogens = OptionValue["ShowHydrogens"]; 
		atoms = mol["Atoms"]; 
		bonds = mol["Bonds"]; 
		coords = mol["AtomCoordinates"]; 
		constraints = mol["Constraints"] /. _Missing -> {}; 

		{bonds, angles, torsions, oops, nonbonds} = 
			makeStructureLists[atoms, bonds]; 

		includeHydrogens = Switch[includeHydrogens, 
			True | False, includeHydrogens, 
			Automatic, Count[atoms, Except["H"]] <= 3, 
			_, 
				Message[StructureReport::opttfa, "ShowHydrogens", includeHydrogens]; 
				Return[$Failed]
			];
		
		rows = If[includeHydrogens, True & /@ atoms, # =!= "H" & /@ atoms]; 
		cartesian = Grid[{
			{"Coordinates"}, 
			Grid[Prepend[MapThread[Join, {{Pick[Range[Length[atoms]], rows], 
				Pick[atoms, rows]} // Thread, Chop[#]}], {"", "Atom", "X", "Y","Z"}], 
				Spacings -> {2, 0.5}, 
				Dividers -> {{False, True, False}, {False, True, False}}, 
				Alignment -> {{Right, Center, {Decimal}}, Automatic}, 
				ItemSize -> All] & /@ 
				{Pick[coords, rows]}
			}, 
			Spacings -> {4}
			]; 

		internal = Grid[{
			{"Bonds", "Angles", "Torsions", "Out-of-Planes", "Nonbonds < 3.0 \[Angstrom]", "Constraints"}, 
			{Grid[Append[#, Chop[distance @@ coords[[#]]]] & /@ 
				Select[bonds, includeHydrogens || FreeQ[atoms[[#]], "H"] &] // 
				Prepend[#, {"i", "j", "r (\[Angstrom])"}] &, 
				Alignment ->{{{Right}, Decimal}, Automatic, {{1, 3} -> " "}}, 
				Dividers -> {None, {False, True, False}}, 
				ItemSize -> All], 
			Grid[Append[#, Chop[(angle @@ coords[[#]]) 180/Pi]] & /@ 
				Select[angles, includeHydrogens || FreeQ[atoms[[#]], "H"] &] // 
				Prepend[#, {"i", "j", "k", "\[Theta] (\[Degree])"}] &, 
				Alignment ->{{{Right}, Decimal}, Automatic, {{1, 4} -> " "}}, 
				Dividers -> {None, {False, True, False}}, 
				ItemSize -> All], 
			Grid[Append[#, Chop[(torsion @@ coords[[#]]) 180/Pi]] & /@ 
				Select[torsions, includeHydrogens || FreeQ[atoms[[#]], "H"] &] // 
				Prepend[#, {"i", "j", "k", "l", "\[Phi] (\[Degree])"}] &, 
				Alignment ->{{{Right}, Decimal}, Automatic, {{1, 5} -> " "}}, 
				Dividers -> {None, {False, True, False}}, 
				ItemSize -> All], 
			Grid[Append[#, Chop[(dihedral @@ coords[[#]]) 180/Pi]] & /@ 
				Select[oops, includeHydrogens || FreeQ[atoms[[#]], "H"] &] // 
				Prepend[#, {"i", "j", "k", "l", "\[Chi] (\[Degree])"}] &, 
				Alignment ->{{{Right}, Decimal}, Automatic, {{1, 5} -> " "}}, 
				Dividers -> {None, {False, True, False}}, 
				ItemSize -> All], 
			Grid[Append[#, Chop[distance @@ coords[[#]]]] & /@ 
				Select[nonbonds, includeHydrogens || FreeQ[atoms[[#]], "H"] &] // 
				DeleteCases[#, {_, _, r_} /; r >= 3.0] & // 
				Prepend[#, {"i", "j", "r (\[Angstrom])"}] &, 
				Alignment ->{{{Right}, Decimal}, Automatic, {{1, 3} -> " "}}, 
				Dividers -> {None, {False, True, False}}], 
			Grid[(Module[{func = #[[1]] /. {"distance" -> distance, "angle" -> angle, 
					"torsion" -> torsion, "dihedral" -> dihedral}, 
					s, units}, 
				{s, units} = Switch[func, 
					distance, {1, "\[Angstrom]"}, 
					angle | torsion | dihedral, {1/Degree, "\[Degree]"}, 
					_, {1, ""}]; 
				Join[PadRight[#[[2]], 4, ""], {#[[1]], (*StringForm["````", s func @@ coords[[#[[2]]]], units]*)s func @@ coords[[#[[2]]]], 
					(*StringForm["````", #[[3]], units]*)#[[3]], #[[4]]}]] & /@ constraints) // 
				Prepend[#, {"i", "j", "k", "l", "type", "q", Subscript["q", 0], Subscript["K", "q"]}] &, 
				Alignment ->{{{Right}, Center, Decimal}, Automatic}, 
				Dividers -> {None, {False, True, False}}, 
				ItemSize -> All]
			}}, 
			Spacings -> {2}, 
			Alignment -> {Center, Top}
			]; 

		Column[{cartesian, internal}, Spacings -> 3]
		]

$nominalMasses = {
	"H" -> "H1", "He" -> "He4", 
	"B" -> "B11", "C" -> "C12", "N" -> "N14", "O" -> "O16", "F" -> "F19", "Ne" -> "Ne20", 
	"Al" -> "Al27", "Si" -> "Si28", "P" -> "P31", "S" -> "S32", "Cl" -> "Cl35", "Ar" -> "Ar40", 
	"Br" -> "Br79"
	}; 

(* TODO: check all uses of RotationMatrix and ReflectionMatrix for Transpose *)
		
alignAxes[mass_?VectorQ, coords_?MatrixQ, atomClasses:{{__Integer}..}, toler_?Positive] /; 
		Length[coords] == Length[mass] := 
	Module[{xyz, n, centerOfMass, inertiaTensor, d, Ia, Ib, Ic, axes, 
			thirdAxis, firstAxis, newCoords, theta, class, atoms, atoms2, point, point2, points, 
			mid, truncatedClass, rotationalTop}, 
		xyz = DiagonalMatrix[ConstantArray[1, {3}]]; 
		n = Length @ coords; 

		centerOfMass = (Total @ MapThread[Times, {coords, mass}]) / Total[mass](* // Echo[#,"com"]&*);
		newCoords = (#1 - centerOfMass) & /@ coords(* // Echo*);

		inertiaTensor = Transpose[newCoords].DiagonalMatrix[mass].newCoords(* // Echo[#, "sum mxyz", MatrixForm]&*); 
		d = Diagonal[inertiaTensor](* // Echo*); 
		inertiaTensor = DiagonalMatrix[d] - inertiaTensor + DiagonalMatrix[ConstantArray[Total@d, Length@d] - d](* // Echo[#, "I", MatrixForm]&*); 
		{{Ia, Ib, Ic}, axes} = Eigensystem[inertiaTensor](* // Echo*); 
		If[Cross[axes[[1]], axes[[2]]].axes[[3]] < 0, 
			axes[[3]] = -axes[[3]]
			];
debugEcho[{Ia, Ib, Ic}, "moments"];

		rotationalTop = Which[

			(2 Abs[Ia-Ib]/(Ia+Ib)) <= toler && (2 Abs[Ib-Ic]/(Ib+Ic)) <= toler, 
				(* these symmetries all have sets of 3 mutually orthogonal 2-fold rotations *)
				(* T, Td, and Th have 1 set; O and Oh have 4 sets; I and Ih have 5 sets *)

				class = (*Echo @ *)FirstCase[atomClasses, {_, _, __}]; 
				(* TODO: handle Missing *)
				{atoms, point} = Catch@(Do[
					mid = Mean[{newCoords[[class[[1]]]], newCoords[[class[[i]]]]}];
					If[Total@mid^2 > toler^2 && 
						superimposableQ[atomClasses, newCoords, (newCoords.RotationMatrix[Pi, mid]), toler], 
debugEcho[class[[{1, i}]], "found first 2-fold"];
						Throw[{class[[{1, i}]], mid }]
						],
					{i, 2, Length[class]}
					];
				(*need fatal Throw here*)
				{Missing[], Missing[]}
				); 

				truncatedClass = (*Echo @ *)DeleteCases[class, First @ atoms]; 
				{atoms2, point2} = Catch@(Do[
					mid = Mean[{newCoords[[truncatedClass[[1]]]], newCoords[[truncatedClass[[i]]]]}];
					If[Total@mid^2 > toler^2 &&
						(point.mid)^2 < toler^2 &&
						superimposableQ[atomClasses, newCoords, (newCoords.RotationMatrix[Pi, mid]), toler], 
debugEcho[truncatedClass[[{1, i}]], "found second 2-fold"];
						Throw[{truncatedClass[[{1, i}]], mid }]
						],
					{i, 2, Length[truncatedClass]}
					];
				(*need fatal Throw here*)
				{Missing[], Missing[]}
				); 

				(* align first 2-fold with X axis *)
				newCoords = Join[newCoords, {point, point2}]; 
				newCoords = newCoords.Transpose@Quiet[RotationMatrix[{point, xyz[[1]]}], RotationMatrix::spln]; 

				(* align second 2-fold with Y axis *)
				newCoords = newCoords.Transpose@RotationMatrix[-ArcTan @@ newCoords[[n + 2, {2, 3}]], xyz[[1]]]; 
				newCoords = Drop[newCoords, -2];
				"Spherical",

			Ic <= toler && 2 Abs[Ia-Ib]/(Ia+Ib) <= toler, 
				newCoords = Join[newCoords, axes]; 
				(* align 0-momentum axis with Z axis *)
				newCoords = With[{r = Quiet[RotationMatrix[{xyz[[3]], newCoords[[n + 3]]}], RotationMatrix::spln]}, 
					newCoords.r]; 
				newCoords = Drop[newCoords, -3];
				"Linear",

			(thirdAxis = 2 Abs[Ia-Ib]/(Ia+Ib) <= toler) || (firstAxis = 2 Abs[Ib-Ic]/(Ib+Ic) <= toler), 
(*
				newCoords = Join[newCoords, If[thirdAxis, axes, RotateLeft@axes]]; 
				(* align nondegenerate axis with Z axis *)
				newCoords = With[{r = Quiet[RotationMatrix[{xyz[[3]], newCoords[[n + 3]]}], RotationMatrix::spln]}, 
					newCoords.r]; 
				(* align one of the degenerate axes with X axis *)
				newCoords = With[{r = RotationMatrix[ArcTan @@ newCoords[[n + 2, {1, 2}]], xyz[[3]]]}, 
					newCoords.r]; 
				newCoords = Drop[newCoords, -3];
*)
				(* align nondegenerate axis with Z axis *)
				newCoords = With[{r = Quiet[RotationMatrix[{xyz[[3]], If[thirdAxis, axes[[3]], axes[[1]]]}], RotationMatrix::spln]}, 
					newCoords.r]; 
debugEcho[newCoords, "aligned with Z axis", MatrixForm];
				(* align orthogonal 2-fold axis or embedding mirror plane with X axis *)
				point = Catch @ (
					Scan[
						Function[class, 
				(* find an atom in an off-axis, special position *)
debugEcho[class[[1]], "testing"];
							If[newCoords[[class[[1]], 1]]^2 + newCoords[[class[[1]], 2]]^2 > toler^2, 
debugEcho[class[[1]], "off axis"];
								theta = ArcTan[newCoords[[class[[1]], 1]], newCoords[[class[[1]], 2]]]; 
								If[newCoords[[class[[1]], 3]]^2 < toler^2, 
debugEcho["testing 2-fold"];
									If[superimposableQ[atomClasses, newCoords, (newCoords.C2XY[theta]), toler], 
debugEcho[class[[1]], "on 2-fold axis"];
										Throw[newCoords[[class[[1]]]]]
										], 
								  (* else *)
debugEcho["testing mirror"];
									If[superimposableQ[atomClasses, newCoords, (newCoords.sigmaXY[theta]), toler], 
debugEcho[class[[1]], "in mirror plane"];
										Throw[Append[newCoords[[class[[1]], {1, 2}]], 0.]]
										]
									]
								]; 
							If[Total[newCoords[[class, 3]]^2] > Length[class] toler^2,  
								If[MatchQ[points = SplitBy[SortBy[newCoords[[class]], Last], Sign[Last[#]]&], {{__}, {__}}], 
									(* if atoms present at +Z and -Z, look for 2-fold axis in XY plane *)
									Scan[Function[point2, 
											mid = Mean[{points[[1, 1]], point2}]; 
											If[Total@mid^2 > toler^2, 
debugEcho["testing 2-fold"];
												If[superimposableQ[atomClasses, newCoords, 
													(newCoords.C2XY[ArcTan[mid[[1]], mid[[2]]]]), toler], 
debugEcho[mid, "found 2-fold axis"];
													Throw[mid]
													]
												]
											], 
										points[[2]]
										], 
								  (* else *)
									(* if atoms present at +Z or -Z, look for mirror plane perpendicular to XY plane *)
									Scan[Function[point2, 
											mid = Mean[{points[[1, 1]], point2}]; 
											If[Total@mid^2 > toler^2, 
debugEcho["testing mirror"];
												If[superimposableQ[atomClasses, newCoords, 
													(newCoords.sigmaXY[ArcTan[mid[[1]], mid[[2]]]]), toler], 
debugEcho[mid, "found mirror plane"];
													Throw[ReplacePart[mid, 3 -> 0.]]
													]
												]
											], 
										points[[1, 2;;]]
										]
									]
								]
							], 
						Cases[atomClasses, {_, __}]
						]; 
debugEcho["no atom for alignment with X axis"];
					Missing[]
					);
				If[MatchQ[point, _?VectorQ], 
				(* align the point with the X axis *)
debugEcho[point, "aligning with X axis"];
					newCoords = With[{r = Transpose @ Quiet[RotationMatrix[debugEcho@{point, xyz[[1]]}], RotationMatrix::spln]}, 
					newCoords.r]
					]; 
				"Symmetric",

			True, 
				newCoords = Join[newCoords, axes]; 
				(* align largest moment axis with Z axis *)
				newCoords = With[{r = Quiet[Transpose @ RotationMatrix[{newCoords[[n + 1]], xyz[[3]]}], RotationMatrix::spln]}, 
					newCoords.r]; 
				(* align second largest moment axis with X axis *)
				newCoords = With[{r = Transpose @ RotationMatrix[-ArcTan @@ newCoords[[n + 2, {1, 2}]], xyz[[3]]]}, 
					newCoords.r]; 
				newCoords = Drop[newCoords, -3];
				"Asymmetric"
			]; 

		inertiaTensor = Transpose[newCoords].DiagonalMatrix[mass].newCoords; 
		d = Diagonal[inertiaTensor]; 
		inertiaTensor = DiagonalMatrix[d] - inertiaTensor + DiagonalMatrix[ConstantArray[Total@d, Length@d] - d]; 
		{{Ia, Ib, Ic}, axes} = Eigensystem[inertiaTensor](* // Echo*); 

		{newCoords, rotationalTop}
		]

classify[p_, toler_] := 
	With[{u = Gather[p, (2 Abs[#1 - #2] / (#1 + #2) <= toler&)]}, 
		(First @ FirstPosition[u, #1] &) /@ p
		]

orbits[A_, toler_] := 
	Module[{i = 1, deg, oldClass, newClass}, 
		deg = Plus @@ (A^i) // debugEcho[#, "degree"]&; 
		oldClass = classify[deg, toler] // debugEcho[#, "class"]&; 
		While[
			deg = Plus @@ A^(++i) // debugEcho[#, "degree"]&; 
			newClass = classify[deg, toler] // debugEcho[#, "class"]&; 
			newClass != oldClass, 

			oldClass = newClass
			]; 
		(Flatten[Position[newClass, #1]] &) /@ Union[newClass]
		]

(* might want special function rotation-inversion matrice, can use inplace of reflection and inversion *)

RotationReflectionMatrix[theta_, axis_] := ReflectionMatrix[axis].RotationMatrix[theta, axis] // Transpose

(* n-fold rotation about Z axis *)
CnZ[n_] := RotationMatrix[2 Pi / n, {0, 0, 1}] // Transpose

(* 2-fold rotation in X,Y plane, at angle theta from X axis *)
C2XY[theta_] := RotationMatrix[Pi, RotationMatrix[theta, {0, 0, 1}].{1, 0, 0}] // Transpose

(* reflection normal to Z axis *)
sigmaZ[] := ReflectionMatrix[{0, 0, 1}] // Transpose

(* reflection normal to axis in X,Y plane, at angle theta from X axis *)
sigmaXY[theta_] := ReflectionMatrix[RotationMatrix[theta, {0, 0, 1}].{0, 1, 0}] // Transpose

(* n-fold rotation reflection about Z axis *)
SnZ[n_] := RotationMatrix[2 Pi / n, {0, 0, 1}].ReflectionMatrix[{0, 0, 1}] // Transpose

(* 4-fold rotations for octahedral symmetry *)
C4X = Transpose @ RotationMatrix[2 Pi / 4, {1, 0, 0}];
C4Y = Transpose @ RotationMatrix[2 Pi / 4, {0, 1, 0}];
C4Z = Transpose @ RotationMatrix[2 Pi / 4, {0, 0, 1}];

(* 5-fold rotations for icosahedral symmetry *)
With[{tau = (1 + Sqrt[5])/2}, 
	C5a = (1/2)*{{1, tau, 1/tau}, {-tau, 1/tau, 1}, {1/tau, -1, tau}};
	C5b = (1/2)*{{1, -(1/tau), tau}, {-(1/tau), tau, 1}, {-tau, -1, 1/tau}}; 
]

highestOrderRotation[coords_, atomClasses_, toler_] := 
	Module[{divisors, n}, 
		divisors = Reverse @ Divisors @ Max[Length /@ atomClasses]; 
		n = Catch @ Scan[
			If[superimposableQ[atomClasses, coords, coords.CnZ[#], toler], Throw[#]]&, 
			divisors]; 
		n
		]

superimposableQ[classes_, c1_, c2_, toler_] := 
	FreeQ[Function[class, FirstCase[class, a_ /; EuclideanDistance[c2[[a]], #] <= toler]& /@ 
			c1[[class]]] /@ classes, 
		_Missing]

Options[Schoenflies] = {
	Tolerance -> 10.^(-4)
	}; 

(* TODO: return coordinates in standard orientaion *)

Schoenflies[mol_Association, opts:OptionsPattern[]] := 
	Module[{toler, atoms, coords, mass, top, distanceMatrix, atomClasses, n, symmetry}, 
		toler = OptionValue[Tolerance]; 

		atoms = mol["Atoms"]; 
		
		If[Length @ atoms == 1, 
			Return[Subscript["K", "h"]]
			]; 
		
		coords = mol["AtomCoordinates"]; 
debugEcho[coords, "before", MatrixForm];
		mass = QuantityMagnitude @ IsotopeData[#, "AtomicMass"]& /@ 
			(atoms /. $nominalMasses) // debugEcho[#, "mass"]&; 

(*
		distanceMatrix = DistanceMatrix[debugEcho[#, "mass-weighted coords", MatrixForm]& @ MapThread[Times, {coords, mass}], 
			DistanceFunction -> EuclideanDistance]; 
*)
(*
		distanceMatrix = DistanceMatrix[coords, DistanceFunction -> EuclideanDistance] + 
			DiagonalMatrix[mass]; 
*)
		distanceMatrix = DistanceMatrix[coords, DistanceFunction -> EuclideanDistance] + 
			DiagonalMatrix[mass];
debugEcho[distanceMatrix, "distance matrix", TableForm[#, TableHeadings->Automatic]&]; 
		atomClasses = debugEcho[#, "atomClasses"]& @ orbits[distanceMatrix, toler];

		{coords, top} = alignAxes[mass, coords, atomClasses, toler]; 
debugEcho[coords, "after", MatrixForm];
debugEcho[top, "top"];

		symmetry = Switch[top, 
			"Spherical", (* Ih, I, Oh, O, Th, Td, T *)
				(* make final realignment for octahedral symmetry *)
				coords = If[superimposableQ[atomClasses, coords, coords.C4X, toler], 
					If[superimposableQ[atomClasses, coords, coords.C4Y, toler], 
						coords (* already aligned *), 
					  (* else *)
						coords.Transpose@RotationMatrix[Pi / 4, {1, 0, 0}]
						],
					If[superimposableQ[atomClasses, coords, coords.C4Y, toler], 
						coords.Transpose@RotationMatrix[Pi / 4, {0, 1, 0}], 
					  (* else *)
						If[superimposableQ[atomClasses, coords, coords.C4Z, toler], 
							coords.Transpose@RotationMatrix[Pi / 4, {0, 0, 1}], 
						  (* else *)
							coords (* not octahedral *)
							]
						]
					];
				If[superimposableQ[atomClasses, coords, coords.Transpose @ RotationMatrix[Pi / 2, {0, 0, 1}], toler], 
					If[superimposableQ[atomClasses, coords, coords.DiagonalMatrix[{-1,-1,-1}], toler], 
						Subscript["O", "h"], 
					  (* else *)
						"O"
						], 
				  (* else *)
					If[superimposableQ[atomClasses, coords, coords.C5a, toler] || 
							superimposableQ[atomClasses, coords, coords.C5b, toler], 
						If[superimposableQ[atomClasses, coords, coords.DiagonalMatrix[{-1,-1,-1}], toler], 
							Subscript["I", "h"], 
						  (* else *)
							"I"
							], 
					  (* else *)
						If[superimposableQ[atomClasses, coords, coords.DiagonalMatrix[{-1,-1,-1}], toler], 
							Subscript["T", "h"], 
						  (* else *)
							If[superimposableQ[atomClasses, coords, coords.sigmaXY[Pi / 4], toler], 
								Subscript["T", "d"], 
							  (* else *)
								"T"
								]
							]
						]
					], 
			"Linear", (* Dooh, Coov *)
				(*Echo[top -> atomClasses]; *)
				If[debugEcho[#, "atomClasses singletons"]& @ MatchQ[atomClasses, {{_}..}], 
					Subscript["C", Row[{Infinity, "v"}]], 
				  (* else *)
					Subscript["D", Row[{Infinity, "h"}]]
					], 
			"Symmetric", (* Dnh, Dnd, Dn, Cnh, Cnv, Cn, Sn, n>=3; D2d *)
				(*top -> atomClasses*)
				n = debugEcho[#, "n"]& @ highestOrderRotation[coords, atomClasses, toler];
				
				(*If[qqq, 
					xxx, 
				  (* else *)
					yyy
					]*)
				If[debugEcho[#, StringForm["SnZ[``]", 2 n]]& @ superimposableQ[atomClasses, coords, coords.SnZ[2 n], toler], 
					If[debugEcho[#, StringForm["sigmaXY[``]", 0]]& @ superimposableQ[atomClasses, coords, coords.sigmaXY[0], toler] || 
							debugEcho[#, StringForm["sigmaXY[``]", Pi/(2 n)]]& @ superimposableQ[atomClasses, coords, coords.sigmaXY[Pi/(2 n)], toler], 
						Subscript["D", Row[{n, "d"}]], 
					  (* else *)
						Subscript["S", 2 n]
						], 
				  (* else *)
					If[debugEcho[#, StringForm["C2XY[``]", 0]]& @ superimposableQ[atomClasses, coords, coords.C2XY[0], toler], 
						If[debugEcho[#, "sigmaZ[]"]& @ superimposableQ[atomClasses, coords, coords.sigmaZ[], toler], 
							Subscript["D", Row[{n, "h"}]], 
						  (* else *)
							Subscript["D", n]
							], 
					  (* else *)
						If[debugEcho[#, "sigmaXY[0]"]& @ superimposableQ[atomClasses, coords, coords.sigmaXY[0], toler], 
							Subscript["C", Row[{n, "v"}]], 
						  (* else *)
							If[debugEcho[#, "sigmaZ[]"]& @ superimposableQ[atomClasses, coords, coords.sigmaZ[], toler], 
								Subscript["C", Row[{n, "h"}]], 
							  (* else *)
								Subscript["C", n]
								]
							]
						]
					],
			_(* "Asymmetric" *),  (* D2h, D2d, D2, C2h, C2v, C2, Cs, Ci (S2), C1 *)
				(*Echo[top -> atomClasses]; *)
				If[debugEcho[#, "CnZ[2]"]& @ superimposableQ[atomClasses, coords, coords.CnZ[2], toler], 
					If[debugEcho[#, "C2XY[0]"]& @ superimposableQ[atomClasses, coords, coords.C2XY[0], toler], 
						If[debugEcho[#, "sigmaZ[]"]& @ superimposableQ[atomClasses, coords, coords.sigmaZ[], toler], 
							Subscript["D", Row[{2, "h"}]], 
						  (* else *)
							Subscript["D", 2]
							], 
					  (* else *)
						If[debugEcho[#, "sigmaZ[]"]& @ superimposableQ[atomClasses, coords, coords.sigmaZ[], toler], 
							Subscript["C", Row[{2, "h"}]], 
						  (* else *)
							If[debugEcho[#, "sigmaXY[0]"]& @ superimposableQ[atomClasses, coords, coords.sigmaXY[0], toler], 
								Subscript["C", Row[{2, "v"}]], 
							  (* else *)
							  	Subscript["C", 2]
								] 
							]
						], 
				  (* else *)
					If[debugEcho[#, "C2XY[0]"]& @ superimposableQ[atomClasses, coords, coords.C2XY[0], toler], 
						coords = coords.RotationMatrix[Pi/2, {0, 1, 0}]; 
						debugEcho["rotate to Z"];
						If[debugEcho[#, "sigmaZ[]"]& @ superimposableQ[atomClasses, coords, coords.sigmaZ[], toler], 
							Subscript["C", Row[{2, "h"}]], 
						  (* else *)
							If[debugEcho[#, "sigmaXY[0]"]& @ superimposableQ[atomClasses, coords, coords.sigmaXY[0], toler], 
								Subscript["C", Row[{2, "v"}]], 
							  (* else *)
							  	Subscript["C", 2]
								] 
							], 
					  (* else *)
						If[debugEcho[#, "C2XY[Pi/2]"]& @ superimposableQ[atomClasses, coords, coords.C2XY[Pi/2], toler], 
							coords = coords.RotationMatrix[Pi/2, {1, 0, 0}]; 
							debugEcho["rotate to Z"];
							If[debugEcho[#, "sigmaZ[]"]& @ superimposableQ[atomClasses, coords, coords.sigmaZ[], toler], 
								Subscript["C", Row[{2, "h"}]], 
							  (* else *)
								If[debugEcho[#, "sigmaXY[0]"]& @ superimposableQ[atomClasses, coords, coords.sigmaXY[0], toler], 
									Subscript["C", Row[{2, "v"}]], 
								  (* else *)
								  	Subscript["C", 2]
									] 
								], 
						  (* else *)
							If[debugEcho[#, "sigmaZ[]"]& @ superimposableQ[atomClasses, coords, coords.sigmaZ[], toler], 
								Subscript["C", "s"], 
							  (* else *)
								If[debugEcho[#, "sigmaXY[0]"]& @ superimposableQ[atomClasses, coords, coords.sigmaXY[0], toler], 
									coords = coords.RotationMatrix[Pi/2, {0, 1, 0}]; 
									debugEcho["rotate to Z"];
									Subscript["C", "s"], 
								  (* else *)
									If[debugEcho[#, "sigmaXY[Pi/2]"]& @ superimposableQ[atomClasses, coords, coords.sigmaXY[Pi/2], toler], 
										coords = coords.RotationMatrix[Pi/2, {1, 0, 0}]; 
										debugEcho["rotate to Z"];
										Subscript["C", "s"], 
									  (* else *)
										If[debugEcho[#, "SnZ[2]"]& @ superimposableQ[atomClasses, coords, coords.SnZ[2], toler], 
											Subscript["C", "i"], 
										  (* else *)
											Subscript["C", 1]
											]
										]
									]
								]
							]
						]
					]
			]; 
debugEcho[coords, "final", MatrixForm];
		symmetry
		]

PointGroupSymmetry[mol_Association] := Missing["NotImplemented"]

End[] (* End Private Context *)

EndPackage[]
