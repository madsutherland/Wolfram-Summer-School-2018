(* ::Package:: *)

(* :Title: SphericalMultivariatePlot.m -- a package template *)

(* :Context: SphericalMultivariatePlot` *)

(* :Author: Roman E. Maeder *)

(* :Summary:
   The SphericalMultivariatePlot package is a syntactically correct framework for package
   development.
 *)

(* :Copyright: (C) <year> by <name or institution> *)

(* :Package Version: 2.0 *)

(* :Mathematica Version: 3.0 *)

(* :History:
   2.0 for Programming in Mathematica, 3rd ed.
   1.1 for Programming in Mathematica, 2nd ed.
   1.0 for Programming in Mathematica, 1st ed.
*)

(* :Keywords: template, SphericalMultivariatePlot, package *)

(* :Sources:
   Roman E. Maeder. Programming in Mathematica, 3rd ed. Addison-Wesley, 1996.
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
   ProgrammingInMathematica/Package1.m
   ProgrammingInMathematica/Package2.m
   ProgrammingInMathematica/Package3.m
*)

(* :Examples:
   <sample input that demonstrates the features of this package>
*)


(* set up the package context, including public imports *)

BeginPackage["SphericalMultivariatePlot`"(*, 
	{"ProgrammingInMathematica`Package1`", "ProgrammingInMathematica`Package2`"}*)]
ClearAll[Evaluate[Context[]<>"*"]]

(* usage messages for the exported functions and the context itself *)

SphericalMultivariatePlot::usage = "SphericalMultivariatePlot.m is a package that does nothing."

fxSphericalToCartesian::usage = "fxSphericaltoCartesian[expression] \
returns expression converted from spherical to Cartesian coordinates. \
\[Theta] is measured from positive z axis and defined from 0 to \
\[Pi], and it is assumed {x,y,z} \[Element] Reals.";

(* error messages for the exported objects *)

fxSphericalToCartesian::nortp = "Input does not contain r, \[Theta], or \[Phi].";
  
Begin["`Private`"]    (* begin the private context (implementation part) *)
ClearAll[Evaluate[Context[]<>"*"]]

(*Needs["ProgrammingInMathematica`Package3`"]*)    (* read in any hidden imports *)

(* unprotect any system functions for which definitions will be made *)

(*protected = Unprotect[ Sin, Cos ]*)

(* definition of auxiliary functions and local (static) variables *)

(*Aux[f_] := Do[something]*)

(*staticvar = 0*)
(*use for example, for a physical constant you refer to in the package*)

(* definition of the exported functions *)

fxSphericalToCartesian[expr_] := (
	If[And @@ (FreeQ[expr, #] &) /@ {Global`r,Global`\[Theta],Global`\[Phi]}, 
		Message[fxSphericalToCartesian::nortp]
		];
	Simplify[expr //. <|Global`r -> Sqrt[Global`x^2 + Global`y^2 + Global`z^2], 
		Global`\[Theta] -> ArcCos[Global`z/Global`r], 
		Global`\[Phi] -> ArcTan[Global`x, Global`y]|>, {Global`x, Global`y, Global`z} \[Element] Reals
		]
	)
    
(* definitions for system functions *)

(*Sin/: Sin[x_]^2 := 1 - Cos[x]^2*)

(*Protect[ Evaluate[protected] ]*)     (* restore protection of system symbols *)

End[ ]         (* end the private context *)

(*Protect[ Function1, Function2 ]*)    (* protect exported symbols - do not do this until you're done debugging *)

EndPackage[ ]  (* end the package context *)



