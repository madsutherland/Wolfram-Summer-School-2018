(* ::Package:: *)

(* :Title: DataImportandNormalization.m -- a package template *)

(* :Context: DataImportandNormalization.m` *)

(* :Author: Madeleine Sutherland *)

(* :Summary:
   This package contains commonly-needed functions for getting data imported from Excel into
   an acceptable format for analysis.
 *)

(* :Copyright: (C) <2018> by <MIT> *)

(* :Mathematica Version: 11.3 *)

(* :Keywords: data, normalization, package *)

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

BeginPackage["`DataImportandNormalization`"] 
	(*{"ProgrammingInMathematica`Package1`", "ProgrammingInMathematica`Package2`"}*)
ClearAll[Evaluate[Context[]<>"*"]]

(* usage messages for the exported functions and the context itself *)

DataImportandNormalization::usage = "DataImportandNormalization.m is a package that manages data imported 
from Excel."

deleteZeroRows::usage = "deleteZeroRows[matrix] removes rows of zeros and rows of nulls"
deleteZeroColumns::usage = "deleteZeroColumns[matrix] removes columns of zeros and columns of nulls"
ratioNormalize::usage = "ratioNormalize[list,position] takes the ratio of each element of {list} to
the specified element"
logNormalize::usage = "logNormalize[n,c] takes the natural log of n+c. Can be applied to matrices where c
is a positive number of magnitude greater than the lowest negative number in the data set"
Begin["`Private`"]    (* begin the private context (implementation part) *)
ClearAll[Evaluate[Context[]<>"*"]]

deleteZeroRows[m_]:=DeleteCases[m,{0..}|{Null..},1]
deleteZeroColumns[m_]:=Block[{Global`s},
Global`s=DeleteCases[Transpose[m],{0..}|{Null..},1];
Transpose[Global`s]]
ratioNormalize[list_,pos_:43]:=list/list[[pos]]
logNormalize[n_,c_]:=Log[n+c]
End[ ]         (* end the private context *)

EndPackage[ ]  (* end the package context *)





















