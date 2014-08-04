unit R3Delp03;
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


interface
procedure getdate(var year,month,day,dayofweek:word);
procedure gettime(var hour,minute,second,sec100:word);

implementation

uses sysutils;

procedure getdate(var year,month,day,dayofweek:word);
var td1:TdateTime;
begin
td1:=Date;
DecodeDate(td1,year,month,day);
dayofweek:=0; (* Not implemented *)
end;

procedure gettime(var hour,minute,second,sec100:word);
var tt1:TdateTime; milsec:word;
begin
tt1:=time;
DecodeTime(tt1,hour,minute,second,milsec);
sec100:=milsec div 10;
end;

begin
end.
