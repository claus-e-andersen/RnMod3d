unit R3Main03;
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

{$I R3dirs03}

interface

{$IFDEF Delphi}
uses sysutils,R3Delp03,R3Writ03,R3Defi03; (* Delphi *)
{$ELSE}
uses dos,R3Writ03,R3Defi03; (* Borland Pascal *)
{$ENDIF}


procedure change_con(i:itype;j:jtype;k:ktype;
                     dir:dirtype;var conOld:nodecontype; conNew:nodecontype);
procedure change_node(i:itype;j:jtype;k:ktype;
                      nodetypNew:nodetyptype;
                      wconNew,econNew,
                      sconNew,nconNew,
                      bconNew,tconNew:nodecontype);
procedure set_node(i:itype;j:jtype;k:ktype;nodetypNew:nodetyptype);
procedure error_node(idst:string;i:itype;j:jtype;k:ktype;message:string);
procedure update_flxval(flx:flxtype;dir:dirtype;i:itype;j:jtype;k:ktype;sign:signtype);
procedure initialize_vars;

procedure create_GP(var GP:gridtype);
procedure dispose_GP(var GP:gridtype);
procedure initialize_GP(var GP:gridtype);
procedure copy_GP(var GP1,GP2:gridtype);
procedure swap_GP(var GP1,GP2:gridtype);

procedure dispose_Fieldbuffer(cBUFname:use_fieldbuffer_type);
procedure create_Fieldbuffer(var cBUF:fieldtype);
procedure save_Field_in_buffer(var cBUF:fieldtype);
procedure restore_Field_from_buffer(var cBUF:fieldtype);

procedure dispose_qBUF(var qBUF:flowfieldtype);
procedure create_qBUF(var qBUF:flowfieldtype);
procedure save_FlowField_in_qBUF(var qBUF:flowfieldtype);
procedure restore_FlowField_from_qBUF(var qBUF:flowfieldtype);

procedure initialize_nodes_and_connectors;
procedure convert_deadnodes_to_NOPs;
procedure check_node_connections;
procedure set_materials;
procedure set_coefficients;
procedure check_coefficients;
function  mat_string(m:mattype):string;
procedure wr_count_nodes(var OM:text);
procedure find_field;
procedure find_better_field_Gauss_Seidel;
procedure find_better_field_Thomas;
procedure reset_fluxes;
procedure wr_residuals(var OM:text);
procedure wr_fluxes(var OM:text);
procedure wr_obses(var OM:text);
procedure wr_all_nodedata(var OM:text);
procedure wr_all_coefficients(var OM:text);
procedure export_field_proc(filename:string);
procedure run_model;
procedure close_model;
procedure update_plotfile(plt:pltfiletype;wdir:dirtype);
procedure solve_iline(j0:jtype;k0:ktype);
procedure solve_jline(i0:itype;k0:ktype);
procedure solve_kline(i0:itype;j0:jtype);
procedure set_cvsize(i:itype;j:jtype;k:ktype;
                     var ArW,ArE,ArS,ArN,ArB,ArT,dV:datatype);
function  FixId(dir:dirtype;h:htype;var wFixfound:FixType):boolean;
function  nodetyp_string(x:nodetyptype):string;
function  nodecon_string(x:nodecontype):string;
procedure find_residuals;
procedure calc_fluxes(iter_period:longint);
procedure calc_obses(iter_period:longint);
function  war_string(w:warningtype):string;

procedure extract_time_and_daystrings(t:timetype; var tst,dayst,dnost:string);
procedure get_t(var t:timetype);
function  time_difference(t1,t2:timetype):datatype;
procedure wr_t(var OM:text; t:timetype);
function  node_flux(dir:dirtype;i:itype;j:jtype;k:ktype):datatype;

implementation
function  converged:boolean; forward;
procedure set_boundary_conditions(bc_proc:aijkproctype); forward;
procedure boundary_conditions_2D(i:itype;j:jtype;k:ktype); forward;

procedure extract_time_and_daystrings(t:timetype; var tst,dayst,dnost:string);
begin
{$IFDEF Delphi}
tst:=TimeToStr(t);
dayst:=FormatDateTime('dd/mm/yyyy',t);
dnost:='';
{$ELSE}
tst:=' ';
dayst:=' ';
dnost:=' ';
{$ENDIF}
end;

{$IFDEF Delphi}
procedure get_t(var t:timetype);
begin
t:=now;
end;
{$ELSE}
procedure get_t(var t:timetype);
begin
t:=0;
end;
{$ENDIF}

function dt(t1,t2:timetype):datatype;
const ID = 'dt';
begin
dt:=(t2-t1)*24*3600; (* seconds *)
end;

function time_difference(t1,t2:timetype):datatype;
(* Time is a renamed version of dt for external use *)
(* dt is too common a name for the purpose          *)
(* June 19, 2000                                    *)
begin
time_difference:=dt(t1,t2);
end;


{$IFDEF Delphi}
procedure wr_t(var OM:text; t:timetype);
var tst,dst,dnost:string;
begin
extract_time_and_daystrings(t,tst,dst,dnost);
writeln(OM,'* Time = ',dnost,' ',dst,' ',tst);
end;
{$ELSE}
function LeadingZero(w :Word): String;
var s:String;
begin
Str(w:0,s);
if Length(s) = 1 then s := '0' + s;
LeadingZero := s;
end;
function get_datetimestring:string;
var Hour, Minute, Second, Sec100: Word;
    Year, Month, Day, DayOfWeek: Word;
    tst,yearst,monthst,dayst:string;
begin
GetTime(Hour, Minute, Second, Sec100);
GetDate(Year, Month, Day, DayOfWeek);
tst:=LeadingZero(hour)+':'+LeadingZero(minute)+':'+LeadingZero(second);
str(year,yearst);
str(month,monthst);
str(day,dayst);
get_datetimestring:=dayst+'-'+monthst+'-'+yearst+' '+tst;
end;
procedure wr_t(var OM:text; t:timetype);
var tst,dst,dnost:string;
begin
dst:='there';
dnost:='';
writeln(OM,'* Time = ',get_datetimestring);
end;
{$ENDIF}



function nodecon_string(x:nodecontype):string;
var st:string;
begin
case x of
  nill:      st:='nill   ';
  std:       st:='std    ';
  noflow:    st:='noflow ';
  ConX:      st:='unchanged ';
else
  st:='Unknown !!';
end; (* case *)
nodecon_string:=st;
end;

function nodetyp_string(x:nodetyptype):string;
var st:string;
begin
case x of
  NOP:        st:='NOP    ';
  free:       st:='free   ';
  fixed1:     st:='fixed1 ';
  fixed2:     st:='fixed2 ';
  fixed3:     st:='fixed3 ';
  fixed4:     st:='fixed4 ';
  fixed5:     st:='fixed5 ';
  NodX: st:='unchanged ';
else
  st:='Unknown !!';
end; (* case *)
nodetyp_string:=st;
end;

function mat_string(m:mattype):string;
var st:string;
begin
case m of
  mat_undefined : st:='mat? ';
  mat1 : st:='mat1 ';
  mat2 : st:='mat2 ';
  mat3 : st:='mat3 ';
  mat4 : st:='mat4 ';
  mat5 : st:='mat5 ';
  mat6 : st:='mat6 ';
  mat7 : st:='mat7 ';
  mat8 : st:='mat8 ';
  mat9 : st:='mat9 ';
  mat10: st:='mat10';
  mat11: st:='mat11';
  mat12: st:='mat12';
  mat13: st:='mat13';
  mat14: st:='mat14';
  mat15: st:='mat15';
else
  st:='Unknown !!';
end; (* case *)
mat_string:=st;
end;

function geometry_string(g:geometrytype):string;
var st:string;
begin
case g of
  cartesian3D   : st:='cartesian3D   ';
  cartesian2D   : st:='cartesian2D   ';
  cylindrical2D : st:='cylindrical2D ';
  else
  st:='Unknown !!';
end; (* case *)
geometry_string:=st;
end;


function war_string(w:warningtype):string;
var st:string;
begin
case w of
  war_none:          st:='war_none';
  war_first:         st:='war_first';
  war_interpolation: st:='war_interpolation';
  war_other:         st:='war_other';
  war_fileimport:    st:='war_fileimport';
  war_convergence:   st:='war_convergence';
  war_residual:      st:='war_residual';
 else
  st:='war_Unknown !!';
end; (* case *)
war_string:=st;
end;

procedure change_con(i:itype;j:jtype;k:ktype;
                     dir:dirtype;var conOld:nodecontype; conNew:nodecontype);
const ID='change_con';
var err:integer;
    conResult:nodecontype;
begin
err:=0;
conresult:=nill; (* Arbitary initialization just to make compiler happy *)
case conOld of
  Nill:
    begin
      case conNew of
        Nill: err:=err+10000;
        NoFlow,ConX: conResult:=conOld; (* Do nothing *)
        Std: err:=err+1;
      else
        err:=err+10;
      end; (* case *)
    end; (* nill *)
  NoFlow:
    begin (* noFlow *)
      case conNew of
        Nill: err:=err+10000;
        NoFlow,ConX: conResult:=conOld; (* Do nothing *)
        Std: begin (* reestablish connection to/from the other node *)
               conResult:=conNew;
               case dir of
                 west:    GP[i-1]^[j]^[k].econ:=std;
                 east:    GP[i+1]^[j]^[k].wcon:=std;
                 south:   GP[i]^[j-1]^[k].ncon:=std;
                 north:   GP[i]^[j+1]^[k].scon:=std;
                 bottom:  GP[i]^[j]^[k-1].tcon:=std;
                 top:     GP[i]^[j]^[k+1].bcon:=std;
               else
                 err:=err+100;
               end; (* dir case *)
             end; (* std *)
        else
          err:=err+10;
        end; (* case conNew *)
    end; (* noflow *)
  Std:
    begin (* std *)
      case conNew of
        nill: err:=err+10000;
        Std,ConX: conResult:=conOld; (* Do nothing *)
        noFlow: begin (* break connection to/from the other node *)
                  conResult:=conNew;
                  case dir of
                    west:    GP[i-1]^[j]^[k].econ:=noflow;
                    east:    GP[i+1]^[j]^[k].wcon:=noflow;
                    south:   GP[i]^[j-1]^[k].ncon:=noflow;
                    north:   GP[i]^[j+1]^[k].scon:=noflow;
                    bottom:  GP[i]^[j]^[k-1].tcon:=noflow;
                    top:     GP[i]^[j]^[k+1].bcon:=noflow;
                  else
                   err:=err+100;
                  end; (* dir case *)
             end; (* std *)
        else
          err:=err+1;
        end; (* case conNew *)
    end; (* std *)
else (* incl. ConX *)
   err:=5000;
end; (* main case *)

if err<>0 then
  begin
    writeln('* Old connection    : ',nodecon_string(conOld));
    writeln('* Wanted connection : ',nodecon_string(conNew));
    writeln('* Direction         : ',dirname(dir));
    writeln('* Error code        : ',err:10);
    writeln('*     1 = Cannot change Nill');
    writeln('*    10 = Unknown New connection');
    writeln('*   100 = Illegal direction (use west, east, south, north, top, bottom)');
    writeln('*  1000 = Cannot change Std to Nill');
    writeln('*  5000 = Illegal old connection');
    writeln('* 10000 = Nill is an Illegal new connection');
    error_node(ID,i,j,k,'Some problem--see above.');
  end
else
  conOld:=conResult;
end;

procedure change_node(i:itype;j:jtype;k:ktype;
                      nodetypNew:nodetyptype;
                      wconNew,econNew,
                      sconNew,nconNew,
                      bconNew,tconNew:nodecontype);
const ID='change_node';
begin
with(GP[i]^[j]^[k]) do
  begin
    if (nodetypNew<>NodX) then nodetyp:=nodetypNew;
    if (nodetyp=NOP) then GP[i]^[j]^[k].valid_fieldvalue:=false
    else GP[i]^[j]^[k].valid_fieldvalue:=true;

    if (nodetypNew=NOP) and
       ( (wconNew<>noFlow) or (econNew<>noFlow) or
         (sconNew<>noFlow) or (nconNew<>noFlow) or
         (bconNew<>noFlow) or (tconNew<>noFlow) ) then
      begin
        writeln('* w-con : ',nodecon_string(wconNew));
        writeln('* e-con : ',nodecon_string(econNew));
        writeln('* s-con : ',nodecon_string(sconNew));
        writeln('* n-con : ',nodecon_string(nconNew));
        writeln('* b-con : ',nodecon_string(bconNew));
        writeln('* t-con : ',nodecon_string(tconNew));
        error_node(ID,i,j,k,'When you try to set nodetyp=NOP,'+
                            ' all connections should be set to noFlow');
      end;
    change_con(i,j,k,west,  wcon,wconNew);
    change_con(i,j,k,east,  econ,econNew);
    change_con(i,j,k,south, scon,sconNew);
    change_con(i,j,k,north, ncon,nconNew);
    change_con(i,j,k,bottom,bcon,bconNew);
    change_con(i,j,k,top,   tcon,tconNew);
  end;
end; (* change_node *)

procedure set_node(i:itype;j:jtype;k:ktype;nodetypNew:nodetyptype);
begin
if (nodetypNew=NOP) then
  change_node(i,j,k,NOP,noFlow,noFlow,noFlow,noFlow,noFlow,noFlow)
else
  change_node(i,j,k,nodetypNew,ConX,ConX,ConX,ConX,ConX,ConX);
end;


procedure set_nodes_region_ijk(iA,iB:itype;ireg:regsettype;
                               jA,jB:jtype;jreg:regsettype;
                               kA,kB:ktype;kreg:regsettype;
                               nodetyp0:nodetyptype;
                               wcon0,econ0,
                               scon0,ncon0,
                               bcon0,tcon0:nodecontype);
const ID='set_nodes_region_ijk';
var i:itype;
    j:jtype;
    k:ktype;
begin
if not (iA<=iB) then error_std(ID,'iA>iB');
if not (jA<=jB) then error_std(ID,'jA>jB');
if not (kA<=kB) then error_std(ID,'kA>kB');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      if in_region_ijk(i,iA,iB,ireg,j,jA,jB,jreg,k,kA,kB,kreg) then
        change_node(i,j,k,nodetyp0,wcon0,econ0,scon0,ncon0,bcon0,tcon0);
end; (* set_nodes_region_ijk *)

procedure set_nodes_region(xFixA,xFixB:Fixtype;xreg:regsettype;
                           yFixA,yFixB:Fixtype;yreg:regsettype;
                           zFixA,zFixB:Fixtype;zreg:regsettype;
                           nodetyp0:nodetyptype;
                           wcon0,econ0,
                           scon0,ncon0,
                           bcon0,tcon0:nodecontype);
const ID='set_nodes_region';
var   iA,iB:itype;
      jA,jB:jtype;
      kA,kB:ktype;
begin
check_fix_region_param(ID,xFixA,xFixB,yFixA,yFixB,zFixA,zFixB);
iA:=wFixVal[xFixA].h; iB:=wFixVal[xFixB].h;
jA:=wFixVal[yFixA].h; jB:=wFixVal[yFixB].h;
kA:=wFixVal[zFixA].h; kB:=wFixVal[zFixB].h; (* ERROR: kB was missing april 7, 1998 *)
set_nodes_region_ijk(iA,iB,xreg,jA,jB,yreg,kA,kB,zreg,
                     nodetyp0,
                     wcon0,econ0,scon0,ncon0,bcon0,tcon0);
end; (* set_nodes_region_fix *)

procedure error_node(idst:string;i:itype;j:jtype;k:ktype;message:string);
begin
writeln('Error_node message from ',idst);
writeln(message);
writeln;
writeln('Nodedata : ');
wr_nodedata(output,i,j,k,true);
writeln('Press enter');
close(LOG);
if press_enter_wanted then readln;
halt;
end;

procedure check_node_connections;
const id='check_node_connections';
var i:itype;
    j:jtype;
    k:ktype;
    err:integer;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      with GP[i]^[j]^[k] do
        begin
          if (econ<>std) and (wcon<>std) and
             (scon<>std) and (ncon<>std) and
             (bcon<>std) and (tcon<>std) and
             (nodetyp<>NOP) and (not (nodetyp in FixedBCs)) then
             error_node(id,i,j,k,'This node is completely disconnected.'+cr+
             'Why has it not been set to NOP ?');
         err:=0;
         if (wcon=std) and (GP[i-1]^[j]^[k].econ<>std) then err:=err+1;
         if (econ=std) and (GP[i+1]^[j]^[k].wcon<>std) then err:=err+1;
         if (scon=std) and (GP[i]^[j-1]^[k].ncon<>std) then err:=err+1;
         if (ncon=std) and (GP[i]^[j+1]^[k].scon<>std) then err:=err+1;
         if (bcon=std) and (GP[i]^[j]^[k-1].tcon<>std) then err:=err+1;
         if (tcon=std) and (GP[i]^[j]^[k+1].bcon<>std) then err:=err+1;

         if (wcon=noFlow) and (GP[i-1]^[j]^[k].econ<>noFlow) then err:=err+1;
         if (econ=noFlow) and (GP[i+1]^[j]^[k].wcon<>noFlow) then err:=err+1;
         if (scon=noFlow) and (GP[i]^[j-1]^[k].ncon<>noFlow) then err:=err+1;
         if (ncon=noFlow) and (GP[i]^[j+1]^[k].scon<>noFlow) then err:=err+1;
         if (bcon=noFlow) and (GP[i]^[j]^[k-1].tcon<>noFlow) then err:=err+1;
         if (tcon=noFlow) and (GP[i]^[j]^[k+1].bcon<>noFlow) then err:=err+1;

         if err<>0 then error_node(ID,i,j,k,'Problem std--non-std!!')

        end;

end;

procedure wr_count_nodes(var OM:text);
const id='wr_count_nodes';
var i:itype;
    j,jusemin,jusemax:jtype;
    k:ktype;
    typ:nodetyptype;
    totsum:longint;
    sum:array[nodetyptype] of longint;
begin
wr_line(OM);
writeln(OM,ID);
if (geometry=cartesian3D) then
  begin
    jusemin:=1;
    jusemax:=jmax;
  end
else
  begin
    jusemin:=2;
    jusemax:=2;
  end;
totsum:=0;
for typ:=Nop to Nodx do
  sum[typ]:=0;
for i:=1 to imax do
  for j:=jusemin to jusemax do
    for k:=1 to kmax do
       begin
         inc(totsum);
         inc(sum[GP[i]^[j]^[k].nodetyp]);
       end;
writeln(OM,'* Type and number of nodes incl. boundary conditions : ');
for typ:=Nop to Nodx do
  if typ in FixedBCs then
    writeln(OM,'*',nodetyp_string(typ):15,' ',sum[typ]:6,' value = ',cBC[typ]:20)
  else
    writeln(OM,'*',nodetyp_string(typ):15,' ',sum[typ]:6);
writeln(OM,'*','Total':15,' ',totsum:6);
end;


function FlxName(Flx:FlxType):string;
const ID='FlxName';
var st:string;
begin
case Flx of
  flx1 : st:='Flx1  ';
  flx2 : st:='Flx2  ';
  flx3 : st:='Flx3  ';
  flx4 : st:='Flx4  ';
  flx5 : st:='Flx5  ';
(*flx6 : st:='Flx6  ';
  flx7 : st:='Flx7  ';
  flx8 : st:='Flx8  ';
  flx9 : st:='Flx9  ';
  flx10: st:='Flx10 ';*)
else
  st:='Unknown!';
end;
FlxName:=st;
end; (* FlxName *)

function ObsName(Obs:ObsType):string;
const ID='ObsName';
var st:string;
begin
case Obs of
  obs1 : st:='Obs1  ';
  obs2 : st:='Obs2  ';
  obs3 : st:='Obs3  ';
  obs4 : st:='Obs4  ';
  obs5 : st:='Obs5  ';
else
  st:='Unknown!';
end;
ObsName:=st;
end; (* ObsName *)


function FixId(dir:dirtype;h:htype;var wFixfound:FixType):boolean;
const ID='FixId';
var wFix,FixStart,FixStop:FixType;
found:boolean;
begin
FixStart:=xFix1; (* Arbitary initialization just to keep compiler happy! *)
FixStop :=xFix1; (* Arbitary initialization just to keep compiler happy! *)

wFixfound:=wFixlast;
case dir of
  xdir: begin FixStart:=xFix1; FixStop:=pred(yFix1) end;
  ydir: begin FixStart:=yFix1; FixStop:=pred(zFix1) end;
  zdir: begin FixStart:=zFix1; FixStop:=pred(wFixLast) end;
else
  error_std(ID,'Unknown dir');
end;
found:=false;
for wFix:=FixStart to FixStop do
  if (wFixVal[wFix].defined) and (wFixVal[wFix].h=h) then
    begin
      found:=true;
      wFixfound:=wFix;
    end;
FixID:=found;
end;


procedure set_cvsize(i:itype;j:jtype;k:ktype;
                     var ArW,ArE,ArS,ArN,ArB,ArT,dV:datatype);
const ID='set_cvsize';
begin
case geometry of
  cartesian3D:
    begin
      ArW:=(y[j+1]-y[j])*(z[k+1]-z[k]);
      ArE:=ArW;
      ArS:=(x[i+1]-x[i])*(z[k+1]-z[k]);
      ArN:=ArS;
      ArB:=(x[i+1]-x[i])*(y[j+1]-y[j]);
      ArT:=ArB;
      dV:=(x[i+1]-x[i])*(y[j+1]-y[j])*(z[k+1]-z[k]);
    end;
  cartesian2D:
    begin
      ArW:=Ly*(z[k+1]-z[k]);
      ArE:=ArW;
      ArS:=(x[i+1]-x[i])*(z[k+1]-z[k]);
      ArN:=ArS;
      ArB:=(x[i+1]-x[i])*Ly;
      ArT:=ArB;
      dV:=(x[i+1]-x[i])*Ly*(z[k+1]-z[k]);
    end;
  cylindrical2D:
    begin
      ArW:=2*pi*x[i]*(z[k+1]-z[k]);
      ArE:=2*pi*x[i+1]*(z[k+1]-z[k]);
      ArS:=1;
      ArN:=1;
      ArB:=pi*(sqr(x[i+1])-sqr(x[i]));
      ArT:=ArB;
      if j=2 then dV:=pi*(sqr(x[i+1])-sqr(x[i]))*(z[k+1]-z[k])
      else dV:=0;
    end;
else
  error_std(ID,'Unknown geometry');
end;
end; (* set_cvsize *)


procedure convert_deadnodes_to_NOPs;
const ID='convert_deadnodes_to_NOPs';
var sum,ArW,ArE,ArS,ArN,ArB,ArT,dV:datatype;
    i:itype;j:jtype;k:ktype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      with GP[i]^[j]^[k] do
        begin
          set_cvsize(i,j,k,ArW,ArE,ArS,ArN,ArB,ArT,dV);
          if (ArW=0) and (wcon=std) then
            change_con(i,j,k,west,wcon,noflow);
          if (ArE=0) and (econ=std) then
            change_con(i,j,k,east,econ,noflow);
          if (ArS=0) and (scon=std) then
            change_con(i,j,k,south,scon,noflow);
          if (ArN=0) and (ncon=std) then
            change_con(i,j,k,north,ncon,noflow);
          if (ArB=0) and (bcon=std) then
            change_con(i,j,k,bottom,bcon,noflow);
          if (ArT=0) and (tcon=std) then
            change_con(i,j,k,top,tcon,noflow);

          sum:=0;
          if (wcon=std) then sum:=sum+1;
          if (econ=std) then sum:=sum+1;
          if (scon=std) then sum:=sum+1;
          if (ncon=std) then sum:=sum+1;
          if (bcon=std) then sum:=sum+1;
          if (tcon=std) then sum:=sum+1;

          if (sum=0) and (nodetyp=free) then (* Fixed nodes may well be disconnected ! *)
            begin
              change_node(i,j,k,NOP,noFlow,noFlow,noFlow,noFlow,noFlow,noFlow);
              if wr_details then writeln(LOG,ID,' Warning. Node at (i j k)=(',i:5,' ',j:5,' ',k:5,') was set to NOP');
            end;

        end;
end;

function Apower(P:datatype):datatype;
const ID='Apower';
      zero=1E-196;
var A:datatype;
begin
case scheme of
  powerlaw:
    begin
      if (1-0.1*abs(P)>zero) then A:=max(0,exp(5*ln(1-0.1*abs(P))))
      else A:=0;
    end;
  central: A:=1-0.5*abs(P);
  upwind:  A:=1;
  hybrid:  A:=max(0,1-0.5*abs(P));
  exact :  begin
             if abs(exp(abs(P))-1)>zero then
               A:=abs(P)/(exp(abs(P))-1)
             else
               A:=1/(1+0.5*abs(P));
           end;
 else
   begin
     A:=0;
     error_std(ID,'Unknown scheme type!');
   end;
 end;
Apower:=A;
end;

function Deff(D1,D2,f:datatype):datatype;
(* See Patankar p. 45 or Risoe-R-623(EN) p. 15 *)
begin
if (D1=0) or (D2=0) then error_std('Deff','D1=0 or D2=0');
if ((1.0-f)/D1+f/D2)=0 then error_std('Deff','((1-f)/D1+f/D2)=0');
Deff:=1.0/((1.0-f)/D1+f/D2);
end;


function a_coeff(Ar,q:datatype;dir:dirtype;i:itype;j:jtype;k:ktype):datatype;
const ID='a_coeff';
      PeMax=100;
var aa,ds,f,D1,D2,Y,QQ,Pe:datatype;
begin
Pe:=0; (* Aribtary initialization just to keep compiler happy    ! *)
ds:=0; (* Aribtary initialization just to keep compiler happy    ! *)
f:=0;  (* Aribtary initialization just to keep compiler happy    ! *)
D1:=0; (* Aribtary initialization just to keep compiler happy    ! *)
D2:=0; (* Aribtary initialization just to keep compiler happy    ! *)
QQ:=0; (* Aribtary initialization just to keep compiler happy    ! *)

case dir of
  west: begin
          ds:=0.5*(dx[i-1]+dx[i]);
          if ds=0 then error_node(ID,i,j,k,'ds=0 '+DirName(dir));
          f:=dx[i]/ds/2;
          D1:=D_def(xdir,i-1,j,k);
          D2:=D_def(xdir,i,j,k);
          QQ:=q;
        end;
  east: begin
          ds:=0.5*(dx[i]+dx[i+1]);
          if ds=0 then error_node(ID,i,j,k,'ds=0 '+DirName(dir));
          f:=dx[i+1]/ds/2;
          D1:=D_def(xdir,i,j,k);
          D2:=D_def(xdir,i+1,j,k);
          QQ:=-q;
        end;
  south: begin
          ds:=0.5*(dy[j-1]+dy[j]);
          if ds=0 then error_node(ID,i,j,k,'ds=0 '+DirName(dir));
          f:=dy[j]/ds/2;
          D1:=D_def(ydir,i,j-1,k);
          D2:=D_def(ydir,i,j,k);
          QQ:=q;
        end;
  north: begin
          ds:=0.5*(dy[j]+dy[j+1]);
          if ds=0 then error_node(ID,i,j,k,'ds=0 '+DirName(dir));
          f:=dy[j+1]/ds/2;
          D1:=D_def(ydir,i,j,k);
          D2:=D_def(ydir,i,j+1,k);
          QQ:=-q;
        end;
  bottom: begin
          ds:=0.5*(dz[k-1]+dz[k]);
          if ds=0 then error_node(ID,i,j,k,'ds=0 '+DirName(dir));
          f:=dz[k]/ds/2;
          D1:=D_def(zdir,i,j,k-1);
          D2:=D_def(zdir,i,j,k);
          QQ:=q;
        end;
  top: begin
          ds:=0.5*(dz[k]+dz[k+1]);
          if ds=0 then error_node(ID,i,j,k,'ds=0 '+DirName(dir));
          f:=dz[k+1]/ds/2;
          D1:=D_def(zdir,i,j,k);
          D2:=D_def(zdir,i,j,k+1);
          QQ:=-q;
        end;
else
  error_std(ID,'Illegal direction '+DirName(dir));
end;

Y:=Deff(D1,D2,f)*Ar/ds;

if Ar<0 then
  error_node(ID,i,j,k,'Ar<0. Direction : '+DirName(dir));
if Y<0 then
  error_node(ID,i,j,k,'Y<0. Direction : '+DirName(dir));
if (Y=0) and (q>0) then Pe:=PeMax;
if (Y=0) and (q<0) then Pe:=-PeMax;
if (Y=0) and (q=0) then
  error_std(ID,'Y=0 and q=0, Cannot calculate Peclet number !');
if (Y<>0) then Pe:=q/Y;
if Pe>PeMax then Pe:=PeMax;
if Pe<-Pemax then Pe:=-PeMax;
aa:=Y*Apower(Pe)+max(0,QQ);
a_coeff:=aa;
end;

procedure set_coefficients;
const ID='set_coefficients';
var ArW,ArE,ArS,ArN,ArB,ArT,dV:datatype;
    i:itype;j:jtype;k:ktype;
    ap_new_dt:datatype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      begin
        set_cvsize(i,j,k,ArW,ArE,ArS,ArN,ArB,ArT,dV);
        with GP[i]^[j]^[k] do
          begin

            if nodetyp=NOP then
              begin
                aW:=0; aE:=0; aSS:=0; aN:=0; aB:=0; aT:=0;
                aP_old_dt:=0;
                b:=dumReal;
                aP:=1;
                valid_fieldvalue:=false;
              end;

            if nodetyp in fixedBCs then
              begin
                aW:=0; aE:=0; aSS:=0; aN:=0; aB:=0; aT:=0;
                aP_old_dt:=0;
                b:=cBC[nodetyp];
                aP:=1;
                valid_fieldvalue:=true;
              end;

           if nodetyp=free then
             begin

               case Wcon of
                 std: aW:=a_coeff(ArW,qw,west,i,j,k);
                 nill,noflow: aW:=0;
               else
                 error_node(ID,i,j,k,'Unknown node connector (W)');
               end; (* Wcon *)

               case Econ of
                 std: aE:=a_coeff(ArE,qE,east,i,j,k);
                 nill,noflow: aE:=0;
               else
                 error_node(ID,i,j,k,'Unknown node connector (E)');
               end; (* Econ *)

               case Scon of
                 std: aSS:=a_coeff(ArS,qS,south,i,j,k);
                 nill,noflow: aSS:=0;
               else
                 error_node(ID,i,j,k,'Unknown node connector (S)');
               end; (* Scon *)

               case Ncon of
                 std: aN:=a_coeff(ArN,qN,north,i,j,k);
                 nill,noflow: aN:=0;
               else
                 error_node(ID,i,j,k,'Unknown node connector (N)');
               end; (* Ncon *)

               case Bcon of
                 std: aB:=a_coeff(ArB,qB,bottom,i,j,k);
                 nill,noflow: aB:=0;
               else
                 error_node(ID,i,j,k,'Unknown node connector (B)');
               end; (* Bcon *)

               case Tcon of
                 std: aT:=a_coeff(ArT,qT,top,i,j,k);
                 nill,noflow: aT:=0;
               else
                 error_node(ID,i,j,k,'Unknown node connector (T)');
               end; (* Tcon *)


               (* The "density" beta may now change from one time step to the
                  next. aP_old should correspond to the time when the last
                  field was calculated. If we ignore transport, generation
                  and decay we have: beta(1)*ca(1) = beta(0)*ca(0). In terms
                  of coefficients this means, that ap_new*ca = ap_old*ca_old) *)

               case solution of
                 unsteady: begin (* Revised June 6, 1999 *)
                             if dtim<=0 then error_std(ID,'dtim<=0 !!');
                             aP_new_dt:=beta_def(i,j,k)*dV;
                             if not unsteady_has_been_started then aP_old_dt:=aP_new_dt;
                           end;
                 steady: begin (* Revised November 1, 1999 *)
                           (* Now a steady-state field can be used as a start of am unstaedy calc. *)
                           ap_new_dt:=beta_def(i,j,k)*dV;;
                           aP_old_dt:=0
                         end;
               else
                 ap_new_dt:=0; (* Just to keep compiler happy *)
                 error_std(ID,'This solution type is invalid.')
               end;

               if solution=unsteady then
                 begin
                   b:=e_def(i,j,k)*dV*G_def(i,j,k)+aP_old_dt/dtim*GP[i]^[j]^[k].c;
                   aP:=aW+aE+aSS+aN+aB+aT+aP_new_dt/dtim+beta_def(i,j,k)*dV*lambda_def(i,j,k);
                 end
               else
                 begin
                   b:=e_def(i,j,k)*dV*G_def(i,j,k);
                   aP:=aW+aE+aSS+aN+aB+aT+beta_def(i,j,k)*dV*lambda_def(i,j,k);
                 end;


               aP_old_dt:=aP_new_dt; (* To be used in next timne step *)

               valid_fieldvalue:=true;

            end; (* nodetyp=free *)

          if not (nodetyp in [NOP,free]+fixedBCs) then
            error_std(ID,'Illegal or undefined nodetyp.');

          end;
      end;

if (solution=unsteady) then
  unsteady_has_been_started:=true
else
  unsteady_has_been_started:=false;
end; (* set_coefficients *)

procedure revise_BC_running_coefficients;
const ID='revise_BC_running_coefficients';
var ArW,ArE,ArS,ArN,ArB,ArT,dV:datatype;
    i:itype;j:jtype;k:ktype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      begin
        set_cvsize(i,j,k,ArW,ArE,ArS,ArN,ArB,ArT,dV);
        with GP[i]^[j]^[k] do
          begin

            if nodetyp in fixedBCs then
              begin
                aW:=0; aE:=0; aSS:=0; aN:=0; aB:=0; aT:=0;
                aP_old_dt:=0;
                b:=cBC[nodetyp];
                aP:=1;
                valid_fieldvalue:=true;
              end;

          end;

      end;
end; (* revise_BC_running_coefficients *)

procedure check_coefficients;
const ID='check_coefficients';
var   i:itype;j:jtype;k:ktype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      if GP[i]^[j]^[k].ap=0 then
        error_node(ID,i,j,k,'ap=0. Try another scheme (e.g. upwind)!');
end;

procedure initialize_vars;
const ID='initialize_vars';
var nodetyp:nodetyptype;
    i:itype;j:jtype;k:ktype;
    wFix:FixType;
    obs:obstype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
for i:=1 to imaxTot do
  begin
    x[i]:=dumReal;
    dx[i]:=dumReal;
    dcdxMax[i]:=0;
  end;
for j:=1 to jmaxTot do
  begin
    y[j]:=dumReal;
    dy[j]:=dumReal;
    dcdyMax[j]:=0;
  end;
for k:=1 to kmaxTot do
  begin
    z[k]:=dumReal;
    dz[k]:=dumReal;
    dcdzMax[k]:=0;
  end;
for wFix:=xFix1 to wFixlast do
  begin
    wFixVal[wFix].w:=DumReal;
    wFixVal[wFix].h:=0;
    wFixVal[wFix].defined:=false;
  end;

imax:=1;
jmax:=1;
kmax:=1;

for nodetyp:=free to NodX do
  cBC[nodetyp]:=0;

for obs:=obs1 to obs_last do
  begin
    Obsval[obs]:=0.0;
    Obsval_old[obs]:=0.0;
    Obsval_change[obs]:=0.0;
  end;


residual_b_sum:=0;
residual_sum:=0;
residual_max:=0;
i_residual_max:=1;
j_residual_max:=1;
k_residual_max:=1;
end; (* initialize_vars; *)

procedure initialize_vars_session;
const ID='initialize_vars_session';
var war:warningtype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
for war:=war_first to war_last do
  warning_table[war]:=0;
end; (* initialize_vars_session *)

procedure initialize_GP(var GP:gridtype);
const ID='initialize_GP';
var   i:itype;j:jtype;k:ktype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      with GP[i]^[j]^[k] do
        begin
          c:=0.0;
          ap:=0;
          aP_old_dt:=0;
          b:=0;
          aw:=0; ae:=0; ass:=0; an:=0; ab:=0; at:=0;
          qw:=0; qe:=0; qs:=0; qn:=0; qb:=0; qt:=0;
          nodetyp:=free;
          econ:=std; wcon:=std;
          scon:=std; ncon:=std;
          bcon:=std; tcon:=std;
          valid_fieldvalue:=false;
          mat:=mat_undefined;
        end;
end;

procedure copy_GP(var GP1,GP2:gridtype);
const ID='copy_GP';
var   i:itype;j:jtype;k:ktype;
begin
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      begin
        if (GP1[i]=nil) then error_std(ID,'GP1 does not exist!');
        if (GP1[i]^[j]=nil) then error_std(ID,'GP1 does not exist!');
      (*GP2[i]^[j]^[k]:=GP1[i]^[j]^[k];*)
      end;
end;

procedure swap_GP(var GP1,GP2:gridtype);
const ID='swap_GP';
var GP0:gridtype;
begin
GP0:=GP1;
GP1:=GP2;
GP2:=GP0;
end;

procedure reset_flowfield;
const ID='reset_flowfield';
var   i:itype;j:jtype;k:ktype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      with GP[i]^[j]^[k] do
        begin
          qw:=0; qe:=0; qs:=0; qn:=0; qb:=0; qt:=0;
        end;
end;


procedure initialize_nodes_and_connectors;
const ID='initialize_nodes_and_connectors';
var   i:itype;j:jtype;k:ktype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      with GP[i]^[j]^[k] do
        begin
          nodetyp:=free; (* The following may be overwritten below *)
          wcon:=std;
          econ:=std;
          scon:=std;
          ncon:=std;
          bcon:=std;
          tcon:=std;

          (* x-sides *)
          if in_plane_ijk([inside,eqAB],i,1,1,j,1,jmax,k,1,kmax) then
            wcon:=nill;
          if in_plane_ijk([inside,eqAB],i,imax,imax,j,1,jmax,k,1,kmax) then
            econ:=nill;

          (* y-sides *)
          if in_plane_ijk([inside,eqAB],i,1,imax,j,1,1,k,1,kmax) then
            scon:=nill;
          if in_plane_ijk([inside,eqAB],i,1,imax,j,jmax,jmax,k,1,kmax) then
            ncon:=nill;

          (* z-sides*)
          if in_plane_ijk([inside,eqAB],i,1,imax,j,1,jmax,k,1,1) then
            bcon:=nill;
          if in_plane_ijk([inside,eqAB],i,1,imax,j,1,jmax,k,kmax,kmax) then
            tcon:=nill;

      end
end; (* initialize_grid *)

procedure set_materials;
var i:itype;
    j:jtype;
    k:ktype;
begin
if (@materials_def=nil) then error_std('set_materials','materials_def = nil is illegal.'
                                       +cr+'A function with materials must be defined.');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      GP[i]^[j]^[k].mat:=materials_def(i,j,k);
end;

procedure create_GP(var GP:gridtype);
const ID='create_GP';
var   i:itype;j:jtype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');

if (imax=1) and (jmax=1) and (kmax=1) then
  begin
    writeln(' imax = ',imax);
    writeln(' jmax = ',jmax);
    writeln(' kmax = ',kmax);
    error_std(ID,'You are creating a GP before run_model has been run the first time!'+cr+
    'Hence imax etc. are not yet defined.'+cr+
    'Please move create_GP to a place after the first run_model.');
  end;

if wr_details then
  begin
    writeln('sizeof(GP[i]^)         = ',sizeof(GP[1]^));
    writeln('sizeof(GP[i]^[j]^)     = ',sizeof(GP[1]^[1]^));
    writeln('sizeof(GP[i]^[j]^[k])  = ',sizeof(GP[1]^[1]^[1]));
  end;

{$IFDEF Delphi}
for i:=1 to imax do
  new(GP[i]);
for i:=1 to imax do
  for j:=1 to jmax do
    new(GP[i]^[j])
{$ELSE}
for i:=1 to imax do
  begin
    if MaxAvail > sizeof(GP[i]^) then
      new(GP[i])
    else
      begin
        writeln('MaxAvail = ',MaxAvail);
        writeln('i=',i);
        error_std(ID,'Too little dynamic space.');
      end
  end;
for i:=1 to imax do
  for j:=1 to jmax do
    begin
      if MaxAvail > sizeof(GP[i]^[j]^) then
        new(GP[i]^[j])
      else
        begin
          writeln('MaxAvail = ',MaxAvail);
          writeln('i=',i,' j=',j);
          error_std(ID,'Too little dynamic space.');
        end;
    end;
{$ENDIF}
end; (* create_GP *)

procedure dispose_GP(var GP:gridtype);
const ID='dispose_GP';
var   i:itype;j:jtype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
for i:=1 to imax do
  for j:=1 to jmax do
        dispose(GP[i]^[j]);
for i:=1 to imax do dispose(GP[i]);
if wr_details then
  begin
    wr_line(output);
    writeln(output,ID,' : Setting: imax:=1, jmax:=1 and kmax:=1');
    wr_line(LOG);
    writeln(LOG,ID,' : Setting: imax:=1, jmax:=1 and kmax:=1');
  end;
imax:=1;
jmax:=1;
kmax:=1;
end; (* dispose_GP *)

procedure create_Fieldbuffer(var cBUF:fieldtype);
const ID='create_Fieldbuffer';
var   i:itype;j:jtype;k:ktype;
begin
if wr_main_procedure_id then writeln(ID,'...');
cBUF.imax:=imax;
cBUF.jmax:=jmax;
cBUF.kmax:=kmax;
for i:=1 to imax do
  new(cBUF.GP[i]);
for i:=1 to imax do
  for j:=1 to jmax do
    new(cBUF.GP[i]^[j]);
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      begin
        cBUF.GP[i]^[j]^[k].c:=0;
        cBUF.GP[i]^[j]^[k].aP_old_dt:=0;
      end;
end; (* create_fieldbuffer *)

procedure save_Field_in_buffer(var cBUF:fieldtype);
const ID='save_Field_in_buffer';
var   i:itype;j:jtype;k:ktype;
begin
if wr_details then writeln(ID);
if LOG_file_is_open then writeln(LOG,ID);
if (imax<>cBUF.imax) or (jmax<>cBUF.jmax) or (kmax<>cBUF.kmax) then
  error_std(ID,'Current field does not match the fieldbuffer. '+cr+
  'Call dispose_Fieldbuffer before this run.');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      begin
        cBUF.GP[i]^[j]^[k].c:=GP[i]^[j]^[k].c;
        cBUF.GP[i]^[j]^[k].aP_old_dt:=GP[i]^[j]^[k].aP_old_dt;
      end;
end; (* save_Field_in_buffer *)

procedure restore_Field_from_buffer(var cBUF:fieldtype);
const ID='restore_Field_from_buffer';
var   i:itype;j:jtype;k:ktype;
begin
if wr_details then writeln(ID);
if LOG_file_is_open then writeln(LOG,ID);
if (imax<>cBUF.imax) or (jmax<>cBUF.jmax) or (kmax<>cBUF.kmax) then
  error_std(ID,'Current field does not match the fieldbuffer. Make a new buffer.');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      begin
        GP[i]^[j]^[k].c:=cBUF.GP[i]^[j]^[k].c;
        GP[i]^[j]^[k].aP_old_dt:=cBUF.GP[i]^[j]^[k].aP_old_dt;
      end;
end; (* restore_Field_from_buffer *)

procedure dispose_fieldbuffer(cBUFname:use_fieldbuffer_type);
const ID='dispose_fieldbuffer';
var   i:itype;j:jtype;
var   cBUF:fieldtype;
begin
case cBUFname of
  cBUF1: cBUF:=cBUF1v;
  cBUF2: cBUF:=cBUF2v;
else
  error_std(ID,'Unknown cBUF');
end;
if wr_main_procedure_id then writeln(ID,'...');
for i:=1 to cBUF.imax do
  for j:=1 to cBUF.jmax do
    if (cBUF.GP[i]^[j]<>nil) then dispose(cBUF.GP[i]^[j]);
for i:=1 to cBUF.imax do
  if (cBuf.GP[i]<>nil) then dispose(cBuf.GP[i]);
cBUF.imax:=1;
cBUF.jmax:=1;
cBUF.kmax:=1;
case cBUFname of
  cBUF1: cBUF1_has_been_created:=false;
  cBUF2: cBUF2_has_been_created:=false;
else
  error_std(ID,'Unknown cBUF');
end;
end; (* dispose_fieldbuffer *)


procedure dispose_qBUF(var qBUF:flowfieldtype);
const ID='dispose_qBUF';
var   i:itype;j:jtype;
begin
if wr_main_procedure_id then writeln(ID,'...');
for i:=1 to qBUF.imax do
  for j:=1 to qBUF.jmax do
    if (qBUF.q[i]^[j]<>nil) then dispose(qBUF.q[i]^[j]);
for i:=1 to qBUF.imax do
  if (qBuf.q[i]<>nil) then dispose(qBuf.q[i]);
qBUF.imax:=1;
qBUF.jmax:=1;
qBUF.kmax:=1;
end; (* dispose_Flowfieldbuffer *)

procedure create_qBUF(var qBUF:flowfieldtype);
const ID='create_qBUF';
var   i:itype;j:jtype;
begin
if wr_main_procedure_id then writeln(ID,'...');
QBUF.imax:=imax;
QBUF.jmax:=jmax;
QBUF.kmax:=kmax;
for i:=1 to imax do
  new(QBUF.q[i]);
for i:=1 to imax do
  for j:=1 to jmax do
    new(QBUF.q[i]^[j]);
qBuf_has_been_created:=true;
end; (* create_qBUF *)

procedure wr_field;
var i:itype;j:jtype;k:ktype;
begin
j:=2;
for k:=kmax downto 1 do
  begin
    for i:=1 to imax do write(GP[i]^[j]^[k].c:5:0,' ');
    writeln;
  end;
end;

procedure wr_field_nodes(kref:ktype);
var i:itype;j:jtype;
begin
for j:=jmax downto 1 do
  begin
    for i:=1 to imax do
      if in_cube([eqAB],i,xFix2,xFix3,j,yFix2,yFix3,2,zFix1,zFix2)
        then
           write('True  ')
        else
           write('False ');
        (*write(nodetyp_string(GP[i]^[j]^[kref].nodetyp):5,' ');*)
    writeln;
  end;
end;

procedure Thomas(Neq:htype);
(* Neq is the number of equations. Hence use only: h=1,2,3..Neq *)
const ID='Thomas';
var alpha: datatype;
    h:htype;
begin
for h:=1 to Neq-1 do
  begin
    if BB[h]=0 then error_std(ID,'BB[h]=0 (1)');
    alpha := AA[h+1]/BB[h];
    BB[h+1]:= BB[h+1]-alpha*CC[h];
    DD[h+1]:= DD[h+1]-alpha*DD[h];
  end;
for h:=Neq downto 2 do
 begin
   if BB[h]=0 then
     begin
     writeln('h  = ',h);
     writeln('AA = ',AA[h]:13);
     writeln('BB = ',BB[h]:13);
     writeln('CC = ',CC[h]:13);
     writeln('DD = ',DD[h]:13);
     error_std(ID,'BB[h]=0 (2)');
     end;
   alpha:=CC[h-1]/BB[h];
   DD[h-1]:=(DD[h-1]-alpha*DD[h]);
 end;
for h:=1 to Neq do
  DD[h]:=DD[h]/BB[h];
end;

procedure solve_iline(j0:jtype;k0:ktype);
const ID='solve_iline';
var i:itype;
    cS,cN,cB,cT,c_old,c_new: datatype;
begin
if wr_all_procedure_id then writeln(ID,' j= ',j0,' k= ',k0);
if wr_all_procedure_id then writeln(LOG,ID,' j= ',j0,' k= ',k0);
for i:=1 to imax do
  with GP[i]^[j0]^[k0] do
    begin
      if j0=1    then cS:=0 else cS:=GP[i]^[j0-1]^[k0].c;
      if j0=jmax then cN:=0 else cN:=GP[i]^[j0+1]^[k0].c;
      if k0=1    then cB:=0 else cB:=GP[i]^[j0]^[k0-1].c;
      if k0=kmax then cT:=0 else cT:=GP[i]^[j0]^[k0+1].c;
      AA[i]:=-aW;
      BB[i]:=aP;
      CC[i]:=-aE;
      DD[i]:=b+aSS*cS+aN*cN+aB*cB+aT*cT;
    end;
Thomas(imax);
for i:=1 to imax do
  begin
    c_old     :=GP[i]^[j0]^[k0].c;
    c_new     :=DD[i];
    if GP[i]^[j0]^[k0].nodetyp<>nop then
      GP[i]^[j0]^[k0].c:=c_old+relax_factor*(c_new-c_old)
    else  (* Do not relax NOP's, May 30, 1999 *)
      GP[i]^[j0]^[k0].c:=c_new;
  end;
end;

procedure solve_jline(i0:itype;k0:ktype);
const ID='solve_jline';
var j:jtype;
    cW,cE,cB,cT,c_old,c_new: datatype;
begin
if wr_all_procedure_id then writeln(ID,' i= ',i0,' k= ',k0);
if wr_all_procedure_id then writeln(LOG,ID,' i= ',i0,' k= ',k0);
for j:=1 to jmax do
  with GP[i0]^[j]^[k0] do
    begin
      if i0=1    then cW:=0 else cW:=GP[i0-1]^[j]^[k0].c;
      if i0=imax then cE:=0 else cE:=GP[i0+1]^[j]^[k0].c;
      if k0=1    then cB:=0 else cB:=GP[i0]^[j]^[k0-1].c;
      if k0=kmax then cT:=0 else cT:=GP[i0]^[j]^[k0+1].c;
      AA[j]:=-aSS;
      BB[j]:=aP;
      CC[j]:=-aN;
      DD[j]:=b+aW*cW+aE*cE+aB*cB+aT*cT;
    end;
Thomas(jmax);
for j:=1 to jmax do
  begin
    c_old     :=GP[i0]^[j]^[k0].c;
    c_new     :=DD[j];
    if GP[i0]^[j]^[k0].nodetyp<>nop then
       GP[i0]^[j]^[k0].c:=c_old+relax_factor*(c_new-c_old)
    else (* Do not relax NOP's, May 30, 1999 *)
       GP[i0]^[j]^[k0].c:=c_new;
  end;
end;

procedure solve_kline(i0:itype;j0:jtype);
const ID='solve_kline';
var k:ktype;
    cW,cE,cS,cN,c_old,c_new: datatype;
begin
if wr_all_procedure_id then writeln(ID,' i= ',i0,' j= ',j0);
if wr_all_procedure_id then writeln(LOG,ID,' i= ',i0,' j= ',j0);

for k:=1 to kmax do
  with GP[i0]^[j0]^[k] do
    begin
      if i0=1    then cW:=0 else cW:=GP[i0-1]^[j0]^[k].c;
      if i0=imax then cE:=0 else cE:=GP[i0+1]^[j0]^[k].c;
      if j0=1    then cS:=0 else cS:=GP[i0]^[j0-1]^[k].c;
      if j0=jmax then cN:=0 else cN:=GP[i0]^[j0+1]^[k].c;
      AA[k]:=-aB;
      BB[k]:=aP;
      CC[k]:=-aT;
      DD[k]:=b+aW*cW+aE*cE+aSS*cS+aN*cN;
    end;
Thomas(kmax);
for k:=1 to kmax do
  begin
    c_old     :=GP[i0]^[j0]^[k].c;
    c_new     :=DD[k];
    if GP[i0]^[j0]^[k].nodetyp<>nop then
      GP[i0]^[j0]^[k].c:=c_old+relax_factor*(c_new-c_old)
    else (* Do not relax NOP's, May 30, 1999 *)
      GP[i0]^[j0]^[k].c:=c_new;
  end;
end;

procedure Find_better_field_Thomas;
const ID='find_better_field_Thomas';
var i:itype;j:jtype;k:ktype;
begin
for j:=1 to jmax do
  for k:=1 to kmax do
    solve_iline(j,k);

if (geometry=cartesian3D) then
  for i:=1 to imax do
    for k:=1 to kmax do
      solve_jline(i,k);

for i:=1 to imax do
  for j:=1 to jmax do
    solve_kline(i,j);
end; (* find_better_field_Thomas *)


procedure Find_better_field_Gauss_Seidel;
const ID='find_better_field_Gauss_Seidel';
var cW,cE,cS,cN,cB,cT,cold:datatype;
    jUseMin,jUseMax:jtype;
    i:itype;j:jtype;k:ktype;
begin
if (geometry=cartesian3D) then
  begin
    jUseMin:=1;
    jUseMax:=jmax
  end
  else
    begin
      jUseMin:=2;
      jUseMax:=2
    end;

for i:=1 to imax do
  for j:=jUseMin to jUseMax do
    for k:=1 to kmax do
      with GP[i]^[j]^[k] do
        begin
          if nodetyp=NOP then c:=DumReal;

          if nodetyp in fixedBCs then c:=cBC[nodetyp];

          if nodetyp=free then
            begin
              cold:=c;
              if Wcon=std then cW:=GP[i-1]^[j]^[k].c else cW:=0;
              if Econ=std then cE:=GP[i+1]^[j]^[k].c else cE:=0;
              if Scon=std then cS:=GP[i]^[j-1]^[k].c else cS:=0;
              if Ncon=std then cN:=GP[i]^[j+1]^[k].c else cN:=0;
              if Bcon=std then cB:=GP[i]^[j]^[k-1].c else cB:=0;
              if Tcon=std then cT:=GP[i]^[j]^[k+1].c else cT:=0;
              if ap>0 then
                c:=(aW*cw+aE*cE+aSS*cS+aN*cN+aB*cB+aT*cT+b)/aP
              else
                begin
                  writeln('* ap = ',ap);
                  error_node(ID,i,j,k,'ap<=0');
                end;
              c:=cold+relax_factor*(c-cold)
            end;

          if not (nodetyp in [NOP,free]+fixedBCs) then
            error_std(ID,'Illegal or undefined nodetyp.');
        end;

end; (* find_better_field_Gauss_Seidel *)

procedure find_field;
const ID='find_field';
var flxcalcno,obscalcno:longint;
    ddt:datatype;
    conv:boolean;
    fresh_convergence_evaluation:boolean;
    iter_BC:longint;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
if (BC_running) and (@BC_running_update_of_cBCs_def=nil) then
  error_std(ID,'No BC_running_update_of_cBCs_def defined!');
get_t(t1);
t2:=t1;
ddt:=0;
calc_fluxes(conv_evaluation_period);
flxcalcno:=0;
calc_obses(conv_evaluation_period);
obscalcno:=0;
iter:=0;
iter_BC:=0;
fresh_convergence_evaluation:=false;
repeat

    if ((iter mod conv_evaluation_period)=0) then
      begin (* time to evaluate convergence... *)
        find_residuals;
        if wr_iteration_line_screen then
          begin
            wr_line(output);
            writeln('Iteration = ',iter:7,' (',max_iterations,')',
                    ' Time = ',ddt/60:6:2,' min (',max_time/60:6:2,')',
                    ' Residual = ',residual_sum:11);
          end;
        if wr_iteration_line_log then
          begin
            wr_line(LOG);
             writeln(LOG,'Iteration = ',iter:7,' (',max_iterations,')',
                        ' Time = ',ddt/60:6:2,' min (',max_time/60:6:2,')',
                        ' Residual = ',residual_sum:11);
          end;
        if wr_residual_during_calc_screen then wr_residuals(output);
        if wr_residual_during_calc_log then wr_residuals(LOG);

        if (@flux_def<>nil) then
          begin
            inc(flxcalcno);
            calc_fluxes(conv_evaluation_period);
            if wr_flux_during_calc_log then wr_fluxes(LOG);
            if wr_flux_during_calc_screen then wr_fluxes(output);
          end;
       if (@probe_def<>nil) then
        begin
          inc(obscalcno);
          calc_obses(conv_evaluation_period);
          if wr_probes_during_calc_log then wr_obses(LOG);
          if wr_probes_during_calc_screen then wr_obses(output);
        end;

        close(LOG);
        append(LOG);
        get_t(t2);
        ddt:=dt(t1,t2);
        fresh_convergence_evaluation:=true;
      end; (* ... time to evaluate convergence *)

     if (BC_running and
         (iter_BC>BC_running_min_iterations) and
         (residual_sum<BC_running_max_residual_sum_before_new_BC) and
         fresh_convergence_evaluation) then

       begin (* time to re-evaluate boundary conditions *)
         iter_BC:=0;
         if wr_BC_running_messages_screen then
           begin
             wr_line(output);
             writeln('Re-evaluating boundary conditions ... (iter=',iter:10,')');
          end;
        if wr_BC_running_messages_log then
           begin
             wr_line(LOG);
             writeln(LOG,'Re-evaluating boundary conditions ... (iter=',iter:10,')');
          end;
         (* We evaluate fluxes and probes because these may be used as part of *)
         (* the running boundary condition.                                    *)
         if (@flux_def<>nil) then calc_fluxes(1);
         if (@probe_def<>nil) then calc_obses(1);
         if (@BC_running_update_of_cBCs_def<>nil) then BC_running_update_of_cBCs_def;
         revise_BC_running_coefficients;
         if wr_BC_running_messages_screen then wr_count_nodes(output); (* bondary conditions may have changed *)
         if wr_BC_running_messages_log then wr_count_nodes(LOG);
         find_residuals;
         fresh_convergence_evaluation:=false;
       end;

    if (@user_procedure_each_iter_def<>nil) then
      user_procedure_each_iter_def;

    if (@solver_def<>nil) then
      solver_def; (* find_better_field *)

    iter_BC:=iter_BC+1;
    iter:=iter+1;
    {$IFDEF check_maxtime_every_iteration}
      get_t(t2);
      ddt:=dt(t1,t2);
    {$ENDIF}
    if ((@flux_def=nil) or (flux_convset=[]) or (flxcalcno>=min_flxcalcno)) and
       ((@probe_def=nil) or (probe_convset=[]) or (obscalcno>=min_obscalcno)) and
       (converged) then
         conv:=true else conv:=false;
    if (BC_running) and
       (@BC_running_convergence_def<>nil) and
       (not BC_running_convergence_def) then conv:=false;

    if conv then (* check also the residual *)
      begin
        if (residual_sum>max_residual_sum) then conv:=false;
      end;

until (iter>min_iterations) and
      ((iter>max_iterations) or
       (ddt>max_time) or
       (conv and fresh_convergence_evaluation));
get_t(t2);
ddt:=dt(t1,t2);
if wr_final_results_screen then
  begin
    wr_line(output);
    writeln('**** RUN ID    : ',runid);
    writeln('**** RUN TITLE : ',runtitle);
  end;
if wr_final_results_log then
  begin
    wr_line(LOG);
    writeln(LOG,'**** RUN ID    : ',runid);
    writeln(LOG,'**** RUN TITLE : ',runtitle);
  end;
if conv then
  begin
    if wr_final_results_log or wr_iteration_line_log then writeln(LOG,'***** Converged');
    if wr_final_results_screen or wr_iteration_line_screen then writeln('***** Converged');
  end;
if iter>max_iterations then
  begin
    warning_std(ID,'***** Max number of iterations was reached',war_other);
  end;
if ddt>max_time then
  begin
    warning_std(ID,'***** Max time reached',war_other);
  end;
if not conv then
  begin
    warning_std(ID,'***** This run did not converge !',war_convergence);
  end;
if wr_final_results_screen or wr_iteration_line_screen then
  begin
    wr_line(output);
    writeln('Iteration = ',iter:7,' (',max_iterations,')',
            ' Time = ',ddt/60:6:2,' min (',max_time/60:6:2,')',
            ' Residual = ',residual_sum:11);
  end;
if wr_final_results_log or wr_iteration_line_log then
  begin
    wr_line(LOG);
    writeln(LOG,'Iteration = ',iter:7,' (',max_iterations,')',
               ' Time = ',ddt/60:6:2,' min (',max_time/60:6:2,')',
               ' Residual = ',residual_sum:11);
  end;
end; (* Find_field *)

procedure Find_residuals;
const ID='find_residuals';
var cW,cE,cS,cN,cB,cT:datatype;
    jUseMin,jUseMax:jtype;
    i:itype;j:jtype;k:ktype;
    residual:datatype;
    residual_sum_old,
    residual_max_old:datatype;
begin
residual_sum_old:=residual_sum;
residual_max_old:=residual_max;
residual_b_sum:=0;
residual_max:=0;
residual_sum:=0;
i_residual_max:=1;
j_residual_max:=1;
k_residual_max:=1;
if (geometry=cartesian3D) then
  begin
    jUseMin:=1;
    jUseMax:=jmax
  end
  else
    begin
      jUseMin:=2;
      jUseMax:=2
    end;

for i:=1 to imax do
  for j:=jUseMin to jUseMax do
    for k:=1 to kmax do
      with GP[i]^[j]^[k] do
        begin

          if nodetyp=free then
            begin
              if Wcon=std then cW:=GP[i-1]^[j]^[k].c else cW:=0;
              if Econ=std then cE:=GP[i+1]^[j]^[k].c else cE:=0;
              if Scon=std then cS:=GP[i]^[j-1]^[k].c else cS:=0;
              if Ncon=std then cN:=GP[i]^[j+1]^[k].c else cN:=0;
              if Bcon=std then cB:=GP[i]^[j]^[k-1].c else cB:=0;
              if Tcon=std then cT:=GP[i]^[j]^[k+1].c else cT:=0;
              if ap>0 then
                residual:=abs((aW*cw+aE*cE+aSS*cS+aN*cN+aB*cB+aT*cT+b)-ap*c)
              else
                begin
                  residual:=dumreal;
                  writeln('* ap = ',ap);
                  error_node(ID,i,j,k,'ap<=0');
                end;
              residual_sum:=residual_sum+residual;
              residual_b_sum:=residual_b_sum+abs(b);
              if residual>residual_max then
                begin
                  residual_max:=residual;
                  i_residual_max:=i;
                  j_residual_max:=j;
                  k_residual_max:=k;
                end;
            end;
        end;

residual_sum_change:=0;
residual_max_change:=0;
if residual_sum_old>0 then
  residual_sum_change:=(residual_sum-residual_sum_old)/residual_sum_old;
if residual_max_old>0 then
  residual_max_change:=(residual_max-residual_max_old)/residual_max_old;
end; (* find_residuals *)


function node_flux(dir:dirtype;i:itype;j:jtype;k:ktype):datatype;
const ID='node_flux';
var ArW,ArE,ArS,ArN,ArB,ArT,dV:datatype;
    res:datatype;
    con:nodecontype;
    found:boolean;
begin
con:=nill; (* Arbitary initialization just to keep compiler happy! *)
res:=0;    (* Arbitary initialization just to keep compiler happy! *)
set_cvsize(i,j,k,ArW,ArE,ArS,ArN,ArB,ArT,dV);
with GP[i]^[j]^[k] do
  begin
    case dir of
      west  : con:=wcon;
      east  : con:=econ;
      south : con:=scon;
      north : con:=ncon;
      bottom: con:=bcon;
      top   : con:=tcon;
   else
     error_std(ID,'Illegal direction. '+Dirname(dir));
   end; (* dir case *)

   if not (nodetyp in [NOP,free,fixed1..NodX]-[NodX]) then
     error_std(ID,'Illegal type of node '+nodetyp_string(nodetyp));

   found:=false;

   if (nodetyp = NOP ) or (con in [noflow,nill]) then
     begin
       res:=0;
       found:=true;
     end;

   if (not found) and (nodetyp=free) then
      begin
        found:=true;
        case dir of
          west:   res:=qw*c+aw*(GP[i-1]^[j]^[k].c-c);
          east:   res:=qe*c+ae*(c-GP[i+1]^[j]^[k].c);

          south:  res:=qs*c+ass*(GP[i]^[j-1]^[k].c-c);
          north:  res:=qn*c+an*(c-GP[i]^[j+1]^[k].c);

          bottom: res:=qb*c+ab*(GP[i]^[j]^[k-1].c-c);
          top:    res:=qt*c+at*(c-GP[i]^[j]^[k+1].c);
        else
          error_std(ID,'Illegal direction. '+Dirname(dir));
        end; (* case *)
     end; (* Free *)

   if (not found) and (nodetyp in [fixed1..NodX]-[NodX]) then
      begin
        (* found:=true; Of no use ! *)

        if (con<>std) then
          error_std(ID,'Connection must be std!'+nodecon_string(con));

        case dir of
          west: with GP[i-1]^[j]^[k] do
                  res:=qe*c+ae*(c-GP[i]^[j]^[k].c);
          east: with GP[i+1]^[j]^[k] do
                  res:=qw*c+aw*(GP[i]^[j]^[k].c-c);

          south: with GP[i]^[j-1]^[k] do
                   res:=qn*c+an*(c-GP[i]^[j]^[k].c);
          north: with GP[i]^[j+1]^[k] do
                   res:=qs*c+ass*(GP[i]^[j]^[k].c-c);

          bottom: with GP[i]^[j]^[k-1] do
                    res:=qt*c+at*(c-GP[i]^[j]^[k].c);
          top: with GP[i]^[j]^[k+1] do
                 res:=qb*c+ab*(GP[i]^[j]^[k].c-c);
        else
          error_std(ID,'Illegal direction. '+Dirname(dir));
        end; (* case *)
     end; (* Free *)

end; (* with *)
node_flux:=res;
end; (* node_flux *)

procedure reset_fluxes;
const ID='reset_fluxes';
var flx:flxtype;
begin
for flx:=flx1 to flx_last do
  begin
    flxval[flx].j:=0.0;
    flxval[flx].q:=0.0;
  end;
end; (* reset_fluxes *)

procedure reset_obses;
const ID='reset_obses';
var obs:obstype;
begin
for obs:=obs1 to obs_last do
  obsval[obs]:=0.0;
end; (* reset_obses *)

procedure update_flxval(flx:flxtype;
                        dir:dirtype;
                        i:itype;j:jtype;k:ktype;
                        sign:signtype);
const ID='update_flxval';
var s,q:datatype;
begin
q:=0; (* Arbitary initialization just to keep compiler happy! *)
case dir of
  west  : q:=GP[i]^[j]^[k].qw;
  east  : q:=GP[i]^[j]^[k].qe;
  south : q:=GP[i]^[j]^[k].qs;
  north : q:=GP[i]^[j]^[k].qn;
  bottom: q:=GP[i]^[j]^[k].qb;
  top   : q:=GP[i]^[j]^[k].qt
else
  error_std(ID,'Illegal dir '+DirName(dir));
end; (* case *)
if sign=plus then s:=1.0 else s:=-1.0;
flxval[flx].j:=flxval[flx].j+s*node_flux(dir,i,j,k);
flxval[flx].q:=flxval[flx].q+s*q;
end;

procedure calc_fluxes(iter_period:longint);
const ID='calc_fluxes';
var i:itype;
    j:jtype;
    k:ktype;
    flx:flxtype;
begin
if wr_details then writeln(ID,'...');
if wr_details then writeln(LOG,ID,'...');
if iter_period<=0 then error_std(ID,'iter_period<=0');
reset_fluxes;
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      if (@flux_def<>nil) then flux_def(i,j,k);
for flx:=flx1 to flx_last do
  begin
    if (iter_period*FlxVal_old[flx].j)<>0 then
      flxval_change[flx].j:=
        (FlxVal[flx].j-FlxVal_old[flx].j)/FlxVal_old[flx].j/iter_period
    else
      flxval_change[flx].j:=0.0;
    flxval_old[flx]:=flxval[flx];
    flxval_change[flx].q:=0.0;
  end;
end; (* calc_flux *)

procedure calc_obses(iter_period:longint);
const ID='calc_obses';
var obs:obstype;
begin
if wr_details then writeln(ID,'...');
if wr_details then writeln(LOG,ID,'...');
estimate_fieldvalues_at_lost_NOPs;
reset_obses;
if (@probe_def<>nil) then probe_def;

for obs:=obs1 to obs_last do
  begin
    if (iter_period*ObsVal_old[obs])<>0 then
      Obsval_change[obs]:=
        (ObsVal[Obs]-ObsVal_old[Obs])/ObsVal_old[Obs]/iter_period
    else
      Obsval_change[obs]:=0.0;
    Obsval_old[obs]:=obsval[obs];
  end;
end; (* calc_obs *)

function converged:boolean;
var flx:flxtype;
    flxres:boolean;
    obs:obstype;
    obsres:boolean;
begin
flxres:=true;
if (@flux_def<>nil) then
  for flx:=flx1 to pred(flx_last) do
    if (flx in flux_convset) and (abs(flxval_change[flx].j)>max_change) then
      flxres:=false;
obsres:=true;
if (@probe_def<>nil) then
  for obs:=obs1 to pred(obs_last) do
    if (obs in probe_convset) and (abs(obsval_change[obs])>max_change) then
      obsres:=false;
if flxres and obsres then
  converged:=true
else
  converged:=false;
end;

procedure wr_residuals(var OM:text);
const ID='wr_residuals';
begin
writeln(OM,'* Abs. sum of bs          = ',residual_b_sum:14);
writeln(OM,'* Abs. sum of residuals   = ',residual_sum:14,' (change = ',residual_sum_change:14,')');
writeln(OM,'* Max residual            = ',residual_max:14,' (change = ',residual_max_change:14,')');
if (geometry=cartesian3D) then
  begin
    writeln(OM,'* Max residual at (i,j,k) = (',i_residual_max:4,',',j_residual_max:4,',',k_residual_max:4,')');
    writeln(OM,'* Max residual at (x,y,z) = (',xnod(i_residual_max):12,',',
                                               ynod(j_residual_max):12,',',
                                               znod(k_residual_max):12,')');
  end
else
  begin
    writeln(OM,'* Max residual at (i,k) = (',i_residual_max:4,',',k_residual_max:4,')');
    writeln(OM,'* Max residual at (x,z) = (',xnod(i_residual_max):12,',',
                                             znod(k_residual_max):12,')');
  end

end;


procedure wr_fluxes(var OM:text);
const ID='wr_fluxes';
var flx:flxtype;
begin
for flx:=flx1 to pred(flx_last) do
    writeln(OM,FlxName(flx),
               ': J = ',FlxVal[flx].j:16,
               ' ( change = ',flxval_change[flx].j:16,' )'+
               ' Q = ',FlxVal[flx].q:16);
end;

procedure wr_obses(var OM:text);
const ID='wr_obses';
var obs:obstype;
begin
for obs:=obs1 to pred(obs_last) do
    writeln(OM,ObsName(obs),
               ': c = ',ObsVal[obs]:16,
               ' ( change = ',Obsval_change[obs]:16,' )');
end;

procedure wr_all_nodedata(var OM:text);
const ID='wr_all_nodedata';
var i:itype;
    j:jtype;
    k:ktype;
begin
wr_line(OM);
writeln(OM,ID);
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      if (i=1) and (j=1) and (k=1) then
        wr_nodedata(OM,i,j,k,true)
      else
        wr_nodedata(OM,i,j,k,false);
end;

procedure wr_all_coefficients(var OM:text);
const ID='wr_all_coefficients';
var i:itype;
    j:jtype;
    k:ktype;
begin
wr_line(OM);
writeln(OM,ID);
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      with GP[i]^[j]^[k] do
        (* if (Nodetyp<>NOP) then *)
           writeln(OM,i:3,' ',j:3,' ',k:3,' ',mat_string(mat),
                                       ' ap=',aP:12 ,' b=',b:12,
                                       ' aw=',aW:12 ,' ae=',aE:12,
                                       ' as=',aSS:12,' an=',aN:12,
                                       ' ab=',aB:12 ,' at=',aT:12);
end;

procedure export_field_proc(filename:string);
const ID='export_field';
var FLD:text;
    i:itype;
    j:jtype;
    k:ktype;
begin
writeln(ID,' (',filename,') ... ');
writeln(LOG,ID,' (',filename,') ... ');
assign(FLD,filename);
rewrite(FLD);
writeln(FLD,imax:10,' ',jmax:10,' ',kmax:10);
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      writeln(FLD,GP[i]^[j]^[k].c:27);
close(FLD);
end;

procedure save_FlowField_in_qBUF(var qBUF:flowfieldtype);
const ID='save_FlowField_in_qBUF';
var i:itype;
    j:jtype;
    k:ktype;
begin
if wr_details then writeln(ID);
if LOG_file_is_open then writeln(LOG,ID);
if (imax<>qBUF.imax) or (jmax<>qBUF.jmax) or (kmax<>qBUF.kmax) then
  begin
    wr_line(output);
    write('Creating new buffer for qBUF (flowfield) w. ');
    writeln('imax = ',imax,' jmax = ',jmax,' kmax = ',kmax);
    wr_line(output);
    wr_line(LOG);
    write(LOG,'Creating new buffer for qBUF (flowfield) w. ');
    writeln(LOG,'imax = ',imax,' jmax = ',jmax,' kmax = ',kmax);
    wr_line(LOG);
    dispose_qBuf(qBUF); (* dispose the old one *)
    create_qBUF(qBUF);
  end;
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      with GP[i]^[j]^[k] do
        begin
          qBUF.q[i]^[j]^[k].qW:=node_flux(west,i,j,k);
          qBUF.q[i]^[j]^[k].qE:=node_flux(east,i,j,k);
          qBUF.q[i]^[j]^[k].qS:=node_flux(south,i,j,k);
          qBUF.q[i]^[j]^[k].qN:=node_flux(north,i,j,k);
          qBUF.q[i]^[j]^[k].qB:=node_flux(bottom,i,j,k);
          qBUF.q[i]^[j]^[k].qT:=node_flux(top,i,j,k);
        end;
end;

procedure restore_FlowField_from_qBUF(var qBUF:flowfieldtype);
const ID='restore_FlowField_from_qBUF';
var i:itype; j:jtype; k:ktype;
begin
if wr_details then writeln(ID);
if LOG_file_is_open then writeln(LOG,ID);
if not qBUF_has_been_created then
  error_std(ID,'No flowfield has yet been saved in a buffer in this session.'+cr+
  'Import therefore failed.'+cr+
  'Import of flowfields from session to session can only be done with a file.');
if (imax<>qBUF.imax) or (jmax<>qBUF.jmax) or (kmax<>qBUF.kmax) then
  error_std(ID,'Flowfield in the buffer does not match the current grid!.');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      with GP[i]^[j]^[k] do
        begin
          qW:=flowfactor*qBUF.q[i]^[j]^[k].qW;
          qE:=flowfactor*qBUF.q[i]^[j]^[k].qE;
          qS:=flowfactor*qBUF.q[i]^[j]^[k].qS;
          qN:=flowfactor*qBUF.q[i]^[j]^[k].qN;
          qB:=flowfactor*qBUF.q[i]^[j]^[k].qB;
          qT:=flowfactor*qBUF.q[i]^[j]^[k].qT;
        end;
end; (* restore_FlowField_from_qBUF *)

procedure export_flowfield_proc(filename:string);
const ID='export_flowfield';
var FLD:text;
    i:itype;
    j:jtype;
    k:ktype;
begin
writeln(ID,' (',filename,') ... ');
assign(FLD,filename);
rewrite(FLD);
writeln(FLD,imax:10,' ',jmax:10,' ',kmax:10);
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      with GP[i]^[j]^[k] do (* export j-fluxes as q-fluxes ! *)
      writeln(FLD,node_flux(west,i,j,k):17  ,' ',node_flux(east,i,j,k):17,' ',
                  node_flux(south,i,j,k):17 ,' ',node_flux(north,i,j,k):17,' ',
                  node_flux(bottom,i,j,k):17,' ',node_flux(top,i,j,k):17);
close(FLD);
end;


procedure import_flowfield_proc(filename:string);
const ID='import_flowfield';
var FLD:text;
    i:itype;
    j:jtype;
    k:ktype;
    imax_rd,jmax_rd,kmax_rd:integer;
begin
write(ID,' (',filename,') ...');
writeln(LOG,ID,' (',filename,')');
if not file_exists(filename) then
  begin
    writeln;
    warning_std(ID,'No such file!. Import ignored. '+filename,war_fileimport)
  end
else
  begin
    assign(FLD,filename);
    reset(FLD);
    readln(FLD,imax_rd,jmax_rd,kmax_rd);
    if (imax_rd<>imax) or (jmax_rd<>jmax) or (kmax_rd<>kmax) then
      begin
        writeln;
        writeln('File : imax  = ',imax_rd:10,' jmax = ',jmax_rd:10,' kmax = ',kmax_rd:10);
        writeln('Grid : imax  = ',imax:10,' jmax = ',jmax:10,' kmax = ',kmax:10);
        writeln('WARNING : This file does not match the current grid!');
        writeln('Import ignored');
        warning_std(ID,'This file does not match the current grid! '+filename,war_fileimport);
      end
    else
      begin
        for i:=1 to imax do
          for j:=1 to jmax do
            for k:=1 to kmax do
              with GP[i]^[j]^[k] do
                begin
                  readln(FLD,qw,qe,qs,qn,qb,qt);
                  qw:=flowfactor*qw;
                  qe:=flowfactor*qe;
                  qs:=flowfactor*qs;
                  qn:=flowfactor*qn;
                  qb:=flowfactor*qb;
                  qt:=flowfactor*qt;
                end;
        writeln(' ... Succeeded.')
      end;
    close(FLD);
  end;
end;


procedure update_plotfile(plt:pltfiletype;wdir:dirtype);
begin
pltfileset:=pltfileset+[plt];
pltfiledirs[plt]:=wdir;
end;

procedure wr_plotfiles_proc;
const ID='wr_plofiles_proc';
var FLD:text;
    i:itype;
    j:jtype;
    k:ktype;
    a:integer;
    ast,bst:string[2];
    plt:pltfiletype;
    plotfile_name:string;
    mt:longint;
    mtst:string;
    firstline:array[pltfiletype] of boolean;
begin
a:=1;
for plt:=plt1 to pred(plt_last) do
  firstline[plt]:=false;
for plt:=plt1 to pred(plt_last) do
  begin
    ast:='xx';
    str(a,ast);
    bst:='00';
    if (ord(ast[2])<ord('0')) or (ord(ast[2])>ord('9')) then
      bst[2]:=ast[1]
    else
      bst:=ast;
    plotfile_name:='f'+runid+'_'+bst+'.dat';
    writeln(ID,' (',plotfile_name,') ... ');
    assign(FLD,plotfile_name);
    rewrite(FLD);

    for i:=1 to imax do
      for j:=1 to jmax do
        for k:=1 to kmax do
           begin
             pltfileset:=[];
             plotfiles_def(i,j,k);
             if plt in pltfileset then
               begin
                 if (@materials_def<>nil) then
                   begin
                     mt:=ord(materials_def(i,j,k))+1;
                     mtst:=mat_string(materials_def(i,j,k))
                   end
                 else
                   begin
                     mt:=0;
                     mtst:='undef'
                   end;
                 if (pltfiledirs[plt]=xdir) and not(firstline[plt]) then
                     begin (* Header line with labels *)
                       writeln(FLD,'j':3,' ',
                                   'k':3,' ',
                                   'y':15,' ',
                                   'z':15,' ',
                                   'c':17,' ',
                                   'nodeno':6,' ',
                                   'matno' :6,' ',
                                   'node   ':7,' ',
                                   'mat  ':5);
                       firstline[plt]:=true;
                     end;
                 if (pltfiledirs[plt]=ydir) and not(firstline[plt]) then
                     begin (* Header line with labels *)
                       writeln(FLD,'i':3,' ',
                                   'k':3,' ',
                                   'x':15,' ',
                                   'z':15,' ',
                                   'c':17,' ',
                                   'nodeno':6,' ',
                                   'matno' :6,' ',
                                   'node   ':7,' ',
                                   'mat  ':5);
                       firstline[plt]:=true;
                     end;
                 if (pltfiledirs[plt]=zdir) and not(firstline[plt]) then
                     begin (* Header line with labels *)
                       writeln(FLD,'i':3,' ',
                                   'j':3,' ',
                                   'x':15,' ',
                                   'y':15,' ',
                                   'c':17,' ',
                                   'nodeno':6,' ',
                                   'matno' :6,' ',
                                   'node   ':7,' ',
                                   'mat  ':5);
                       firstline[plt]:=true;
                     end;
                 if pltfiledirs[plt]=xdir then
                     writeln(FLD,j:3,' ',
                                 k:3,' ',
                                 y[j]+0.5*dy[j]:15,' ',
                                 z[k]+0.5*dz[k]:15,' ',
                                 GP[i]^[j]^[k].c:17,' ',
                                 ord(GP[i]^[j]^[k].nodetyp):6,' ',
                                 mt:6,' ',
                                 nodetyp_string(GP[i]^[j]^[k].nodetyp):7,' ',
                                 mtst:5);
                 if pltfiledirs[plt]=ydir then
                   writeln(FLD,i:3,' ',
                               k:3,' ',
                               x[i]+0.5*dx[i]:15,' ',
                               z[k]+0.5*dz[k]:15,' ',
                               GP[i]^[j]^[k].c:17,' ',
                               ord(GP[i]^[j]^[k].nodetyp):6,' ',
                               mt:6,' ',
                               nodetyp_string(GP[i]^[j]^[k].nodetyp):7,' ',
                               mtst:5);
                 if pltfiledirs[plt]=zdir then
                   writeln(FLD,i:3,' ',
                               j:3,' ',
                               x[i]+0.5*dx[i]:15,' ',
                               y[j]+0.5*dy[j]:15,' ',
                               GP[i]^[j]^[k].c:17,' ',
                               ord(GP[i]^[j]^[k].nodetyp):6,' ',
                               mt:6,' ',
                               nodetyp_string(GP[i]^[j]^[k].nodetyp):7,' ',
                               mtst:5);
               end
           end;
   close(FLD);
   inc(a);
  end;
end;

procedure import_field_proc(filename:string);
const ID='import_field';
var FLD:text;
    i:itype;
    j:jtype;
    k:ktype;
    imax_rd,jmax_rd,kmax_rd:integer;
begin
write(ID,' (',filename,') ...');
writeln(LOG,ID,' (',filename,')');
if not file_exists(filename) then
  begin
    writeln;
    warning_std(ID,'No such file!. Import ignored. '+filename,war_fileimport)
  end
else
  begin
    assign(FLD,filename);
    reset(FLD);
    readln(FLD,imax_rd,jmax_rd,kmax_rd);
    if (imax_rd<>imax) or (jmax_rd<>jmax) or (kmax_rd<>kmax) then
      begin
        writeln;
        writeln('File : imax  = ',imax_rd:10,' jmax = ',jmax_rd:10,' kmax = ',kmax_rd:10);
        writeln('Grid : imax  = ',imax:10,' jmax = ',jmax:10,' kmax = ',kmax:10);
        writeln('WARNING : This file does not match the current grid!');
        writeln('Import ignored');
        warning_std(ID,'This file does not match the current grid! '+filename,war_fileimport);
      end
    else
      begin
        for i:=1 to imax do
          for j:=1 to jmax do
            for k:=1 to kmax do
              readln(FLD,GP[i]^[j]^[k].c);
        writeln(' ... Succeeded.')
      end;
    close(FLD);
  end;
end;

procedure set_initialfield;
const ID='set_initialfield';
var i:itype;
    j:jtype;
    k:ktype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
       if (@initialfield_def<>nil) and
          (GP[i]^[j]^[k].nodetyp=free) then (* No impact on fixed nodes *)
         begin
           GP[i]^[j]^[k].c:=initialfield_def(i,j,k);
         end
end;

procedure boundary_conditions_2D(i:itype;j:jtype;k:ktype);
var xFix,zFix,xFixMax,zFixmax:Fixtype;
begin
xFix:=xFix1;
repeat
  xFixMax:=xFix;
  xFix:=succ(xFix)
until (xFix=yFix1) or (not wFixVal[xFix].defined);
zFix:=zFix1;
repeat
  zFixMax:=zFix;
  zFix:=succ(zFix)
until (zFix=wFixlast) or (not wFixVal[zFix].defined);

if in_plane([inside,eqAB],
            i,xFix1,xFixMax,
            j,yFix1,yFix1,
            k,zFix1,zFixMax) then
  change_node(i,j,k,NOP,noFlow,noFlow,noFlow,noFlow,noFlow,noFlow);
if in_plane([inside,eqAB],
            i,xFix1,xFixMax,
            j,yFix2,yFix2,
            k,zFix1,zFixMax) then
  change_node(i,j,k,NOP,noFlow,noFlow,noFlow,noFlow,noFlow,noFlow);
end;

procedure set_boundary_conditions(bc_proc:aijkproctype);
var i:itype;j:jtype; k:ktype;
begin
if (@bc_proc<>nil) then
  for i:=1 to imax do
    for j:=1 to jmax do
      for k:=1 to kmax do bc_proc(i,j,k);
end;


procedure evaluate_grid;
const ID='evaluate_grid';
var i:itype;j:jtype; k:ktype;
begin
if wr_details then writeln(ID,'...');
if wr_details then writeln(LOG,ID,'...');
dcdwmaxall:=0;
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      with GP[i]^[j]^[k] do
        begin
          if (nodetyp=free) and (econ=std) then
            dcdxMax[i]:=max(dcdxMax[i],abs(c-GP[i+1]^[j]^[k].c));
          if (nodetyp=free) and (ncon=std) then
            dcdyMax[j]:=max(dcdyMax[j],abs(c-GP[i]^[j+1]^[k].c));
          if (nodetyp=free) and (tcon=std) then
            dcdzMax[k]:=max(dcdzMax[k],abs(c-GP[i]^[j]^[k+1].c));
         dcdwMaxAll:=max(dcdwMaxAll,dcdxMax[i]);
         dcdwMaxAll:=max(dcdwMaxAll,dcdyMax[j]);
         dcdwMaxAll:=max(dcdwMaxAll,dcdzMax[k]);
        end;
end;


procedure new_grid;
const ID='new_grid (run_model)';
begin
if (@grid_def=nil) then
  error_std(ID,'Grid_def = nil is illegal.'+cr+'A grid must be defined.');
warning_std(ID,'Redefining the grid!',war_other);
wr_memory_status(LOG,'');
wr_memory_status(output,'');
dispose_GP(GP);
number_of_nodes:=0;
initialize_vars;
wr_memory_status(LOG,'');
wr_memory_status(output,'');
if geometry<>cartesian3D then
  begin (* 2D-system *)
    writeln('Geometry = ',geometry_string(geometry),' (preset y-axis)');
    set_FixVal(yFix1,0.00);
    set_FixVal(yFix2,Ly);
    set_axis_single(yFix1,yFix2,1,FocusA,1);
  end;

grid_def;
grid_def_last:=grid_def;

if (geometry=cylindrical2D) and
   (wFixVal[xFix1].defined) and
   (wFixVal[xFix1].w<0)
then
  error_std(ID,'xFix1 (r) < 0 is not possible for cylindrical coordinates!');

check_axes;
number_of_nodes:=imax;
number_of_nodes:=number_of_nodes*jmax;
number_of_nodes:=number_of_nodes*kmax;

create_GP(GP);
initialize_GP(GP);
initialize_nodes_and_connectors;
wr_memory_status(LOG,'');
wr_memory_status(output,'');
wr_line(output);

end; (* New grid *)

procedure first_runproc;
const ID='first_runproc (run_model)';
begin
    warnings_were_issued:=false;
    max_issued_warning_priority:=war_none;
    log_name:='f'+runid+'LOG.dat';
    res_name:='f'+runid+'RES.dat';
    assign(LOG,log_name);
    assign(RES,res_name);
    rewrite(LOG);
    rewrite(RES);
    LOG_file_is_open:= true;
    RES_file_is_open:= true;
    wr_header(output);
    wr_header(LOG);
    get_t(t1);
    wr_t(LOG,t1);
    writeln('**** LOG File  : ',log_name);
    writeln('**** RES File  : ',res_name);
    writeln('**** RUN ID    : ',runid);
    writeln('**** RUN TITLE : ',runtitle);
    writeln(LOG,'**** LOG File  : ',log_name);
    writeln(LOG,'**** RES File  : ',res_name);
    writeln(LOG,'**** RUN ID    : ',runid);
    writeln(LOG,'**** RUN TITLE : ',runtitle);
    initialize_vars;
    initialize_vars_session;

    wr_memory_status(LOG,'');

    geometry_original:=geometry;

    if geometry<>cartesian3D then
      begin (* 2D-systen *)
        wr_line(output);
        writeln('Geometry = ',geometry_string(geometry),' (preset y-axis)');
        wr_line(LOg);
        writeln(LOG,'Geometry = ',geometry_string(geometry),' (preset y-axis)');
        set_FixVal(yFix1,0.00);
        set_FixVal(yFix2,Ly);
        set_axis_single(yFix1,yFix2,1,FocusA,1);
      end;

    if (@grid_def=nil) then
      error_std(ID,'Grid_def = nil is illegal.'+cr+'A grid must be defined.');

    grid_def;

    grid_def_last:=grid_def;

    if (geometry=cylindrical2D) and
       (wFixVal[xFix1].defined) and
       (wFixVal[xFix1].w<0) then
       error_std(ID,'xFix1 (r) < 0 is not possible for cylindrical coordinates!');

    check_axes;
    number_of_nodes:=imax;
    number_of_nodes:=number_of_nodes*jmax;
    number_of_nodes:=number_of_nodes*kmax;
    create_GP(GP);
    initialize_GP(GP);
    initialize_nodes_and_connectors;
    qBUF_has_been_created:=false;
    wr_memory_status(LOG,'');
    wr_memory_status(output,'');
    wr_line(output);
    close(LOG);
    append(LOG);
    close(RES);
    append(RES);
    first_run:=false;
end;

procedure run_model;
const ID='run_model';
begin
if wr_main_procedure_id then writeln(ID,'...');
if LOG_file_is_open and wr_main_procedure_id then writeln(LOG,ID,'...');

(* Observe, that if the specified buffer has not yet been created, then
   nothing is transfered to GP. *)
case use_fieldbuffer of
  cBUF1: begin
           if wr_details then writeln(ID,' Now recalling an old state from cBUF1');
           if LOG_file_is_open and wr_details then writeln(LOG,ID,' Now recalling an old state from cBUF1');
           if cBUF1_has_been_created then
             restore_Field_from_buffer(cBUF1v); (* transfer cBUF1 to GP *)
         end;
  cBUF2: begin
           if wr_details then writeln(ID,' Now recalling an old state from cBUF1');
           if LOG_file_is_open and wr_details then writeln(LOG,ID,' Now recalling an old state from cBUF1');
           if cBUF2_has_been_created then
             restore_Field_from_buffer(cBUF2v); (* transfer cBUF2 to GP *)
         end;
  no_cBUF: (*do nothing *);
  else
    error_std(ID,'Unknown cBUF!')
end; (* case *)

if ((@grid_def_last<>nil) and (@grid_def_last<>@grid_def)) or
   ((force_new_grid_in_every_run) and (not first_run))
then new_grid;


if first_run then first_runproc;

if (not first_run) and (geometry<>geometry_original) then
  error_std(Id,'The geometry cannot be redefined within a session. '+cr+
               'Do not change geometry!');

if import_field_name='' then
  import_field_name:='f'+runid+'_00.dat';
if export_field_name='' then
  export_field_name:='f'+runid+'_00.dat';
if temp_field_name='' then
  temp_field_name:='f'+runid+'tmp.dat';
if flowfield_name='' then
  flowfield_name:='f'+runid+'flw.dat';

case flowfield of
  import           : import_flowfield_proc(flowfield_name);
  import_from_qbuf : restore_flowfield_from_qbuf(qbuf);
else
  reset_flowfield;
end;

initialize_nodes_and_connectors; (* Always start from here. June 12, 1999 *)
set_boundary_conditions(boundary_conditions_def);

if geometry<>cartesian3D then
  set_boundary_conditions(boundary_conditions_2D);
convert_deadnodes_to_NOPs;
check_node_connections;

(* Observe: These fields go into the coefficients *)
if (solution=steady) and
   ((import_initialfield) or (@initialfield_def<>nil)) then
  warning_std(ID,'Do not try to import an initial field with solution=steady.',
                 war_other);
if (import_initialfield) and (@initialfield_def<>nil) then
  warning_std(ID,'You have requested that a field is imported from file '+cr+
                 'AND that a field is set by initialfield_def !!',
                 war_other);

if unsteady_has_been_started and ((import_initialfield) or (@initialfield_def<>nil)) then
  error_std(ID,'You have tried to import an initial field after the first (unsteady) run_model.'+cr+
                 'You should set initialfield_def:=nil or'+cr+
                 'import_initialfield:=false');

if import_initialfield then import_field_proc(import_field_name);
if (@initialfield_def<>nil) then set_initialfield;

if (@D_def=nil)      then error_std(ID,'D_def = nil is illegal.'+cr+'A function with diffusivity must be defined.');
if (@e_def=nil)      then error_std(ID,'e_def = nil is illegal.'+cr+'A function with porosity must be defined.');
if (@beta_def=nil)   then error_std(ID,'beta_def = nil is illegal.'
                                        +cr+'A function with partition-corrected-porosity must be defined.');
if (@G_def=nil)      then error_std(ID,'G_def = nil is illegal.'+cr+'A function with generation rate must be defined.');
if (@lambda_def=nil) then error_std(ID,'lambda_def = nil is illegal.'+cr+'A function with decay constant must be defined.');

if (@boundary_conditions_def=nil) then error_std(ID,'boundary_conditions_def = nil is illegal.'
                                                    +cr+'A procedure with boundary conditions must be defined.');

set_materials;
set_coefficients;
check_coefficients;
if wr_node_numbers then wr_count_nodes(LOG);
if wr_nodes then wr_all_nodedata(LOG);
if wr_node_sizes then wr_cvsize(LOG);
if wr_coefficients then wr_all_coefficients(LOG);
if import_finalfield_guess then import_field_proc(import_field_name);
if (flux_convset=[]) and (probe_convset=[]) then
  warning_std(ID,'Both flx_convset and obs_convset are empty. '+cr+
                 'Check of convergence, therefore will be invalid.',
                 war_other);
close(LOG);
append(LOG);
close(RES);
append(RES);

find_field;

if wr_final_results_log then wr_residuals(LOG);
if wr_final_results_screen then wr_residuals(output);

if (@flux_def<>nil) then
        begin
          calc_fluxes(1);
          if wr_final_results_log then wr_fluxes(LOG);
          if wr_final_results_screen then wr_fluxes(output);
        end;
if (@probe_def<>nil) then
        begin
          calc_obses(1);
          if wr_final_results_log then wr_obses(LOG);
          if wr_final_results_screen then wr_obses(output);
        end;

if (residual_sum>max_residual_sum) or
   ((residual_sum*residual_sum_warning_limit>residual_b_sum) and
    (residual_b_sum<>0)) then
   warning_std(ID,'The sum of abs. residuals seems to be too large compared to the source term (b)'+cr+
                  'Solution: Do more iterations or decrease the diffusivity.',
                  war_residual);

if wr_materials_volumes then
  wr_material_volumes_etc(LOG);
if wr_BC_running_messages_log then
  wr_count_nodes(LOG); (* bondary conditions may have changed *)

evaluate_grid;
if wr_axes then wr_axes_proc(LOG);
if export_field then
  export_field_proc(export_field_name);
if (@plotfiles_def<>nil) then
  wr_plotfiles_proc;
case flowfield of
  export         : export_flowfield_proc(flowfield_name);
  export_to_qbuf : save_flowfield_in_qbuf(qbuf);
end;

case use_fieldbuffer of
  cBUF1: begin
           if not cBUF1_has_been_created then create_Fieldbuffer(cBUF1v);
           cBUF1_has_been_created:=true;
           save_Field_in_Buffer(cBUF1v);
         end;
  cBUF2: begin
           if not cBUF2_has_been_created then create_Fieldbuffer(cBUF2v);
           cBUF2_has_been_created:=true;
           save_Field_in_Buffer(cBUF2v);
         end;
  no_cBUF: (* do nothing *);
  else
    error_std(ID,'Unknown cBUF!')
end; (* case *)


close(LOG);
append(LOG);
close(RES);
append(RES);
end; (* run_model *)

procedure close_model;
begin
if warnings_were_issued then
  begin
    wr_line(output);
    writeln('OBSERVE : Warnings were issued during this session. ');
    wr_warning_table(output);
    wr_line(output);
  end;
if warnings_were_issued then
  begin
    wr_line(LOG);
    writeln(LOG,'OBSERVE : Warnings were issued during this session. ');
    wr_warning_table(LOG);
    wr_line(LOG);
  end;
get_t(t1);
wr_line(LOG);
wr_t(LOG,t1);
wr_line(output);
writeln('**** LOG File  : ',log_name);
writeln('**** RES File  : ',res_name);
writeln('**** RUN ID    : ',runid);
writeln('**** RUN TITLE : ',runtitle);
wr_line(output);
wr_t(output,t1);
writeln('ByeBye (',programid,')');
if press_enter_wanted then readln;
close(LOG);
close(RES);
end;


begin (* main *)
imax                      := 1;
jmax                      := 1;
kmax                      := 1;
first_run                 := true;
unsteady_has_been_started := false;
grid_def_last             := nil;
pltfileset                := [];
temp_field_name           := '';
LOG_file_is_open          := false;
RES_file_is_open          := false;
cBUF1_has_been_created    := false;
cBUF2_has_been_created    := false;
set_control_variables_to_defaults;
end. (* main *)
