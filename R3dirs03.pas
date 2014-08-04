(* Programid      : RnMod3d                                                                *)
(* Versionid      : Version 0.8 (Sep. 15, 1997 - July 18, 2000)                            *)
(* Description    : Radon transport model - 3D time-dependent - finite volume              *)
(* Programmer1    : Claus E. Andersen, Risoe National Laboratory                           *)
(* Programmer2    : DK-4000 Roskilde, Denmark. claus.andersen@risoe.dk                     *)
(* Copyright      : Risoe National Laboratory, Denmark                                     *)
(* Documentation  : User's Guide to RnMod3d, Risoe-R-1201(EN)                              *)
(* Files          : R3Defi03.pas = Global definitions                                      *)
(*                  R3Main03.pas = Main program                                            *)
(*                  R3Delp03.pas = Special Delphi 3.0 / Borland Pascal 7.0 code            *)
(*                  R3Dirs03.pas = Compiler directives                                     *)
(*                  R3Writ03.pas = Mainly procedures for writing text                      *)



{$DEFINE Delphi}
{$apptype console}
(**)

{$DEFINE imax100}
{$DEFINE jmax100}
{$DEFINE kmax200}

{$DEFINE do_not_check_maxtime_every_iteration}

{$N+} (* 80287 numerical processor       *)
{$R-} (* Range checking                  *)
{$S+} (* Stack-overflow checking         *)
{$F+}

