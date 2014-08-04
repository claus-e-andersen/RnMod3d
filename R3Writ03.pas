unit R3Writ03;
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
uses sysutils,R3Defi03; (* Delphi *)
{$ELSE}
uses dos,R3defi03; (* Borland Pascal *)
{$ENDIF}

procedure error_std(idst:string;message:string);
procedure warning_std(idst:string;message:string;war:warningtype);

procedure wr_memory_status(var OM:text;st:string);
procedure wr_material_volumes_etc(var OM:text);
procedure wr_axis(var OM:text;dir:dirtype);
procedure wr_axes_proc(var OM:text);
procedure wr_nodedata(var OM:text; i:itype;j:jtype;k:ktype; header:boolean);
procedure wr_cvsize(var OM:text);
procedure get_fieldvalue(xp,dx,yp,dy,zp,dz:datatype; var c,dc:datatype; var valid:boolean);
procedure get_fieldvalue2D(xp,dx,zp,dz:datatype; var c,dc:datatype; var valid:boolean);
function fieldvalue(xp,yp,zp:datatype; var valid:boolean):datatype;
function fieldvalue2D(xp,zp:datatype; var valid:boolean):datatype;
function xnod(i:itype):datatype;
function ynod(j:jtype):datatype;
function znod(k:ktype):datatype;
procedure find_i_intpol(xp:datatype; var iA1,iB1,iA2,iB2:itype);
procedure find_j_intpol(yp:datatype; var jA1,jB1,jA2,jB2:jtype);
procedure find_k_intpol(zp:datatype; var kA1,kB1,kA2,kB2:ktype);
function find_i_face(xp:datatype):itype;
function find_j_face(yp:datatype):jtype;
function find_k_face(zp:datatype):ktype;
function find_i_cv(xp:datatype):itype;
function find_j_cv(yp:datatype):jtype;
function find_k_cv(zp:datatype):ktype;
procedure get_avgfield(x1,x2,y1,y2,z1,z2,ddd:datatype; var c,dc:datatype);
procedure get_avgfield2D(x1,x2,z1,z2,ddd:datatype; var c,dc:datatype);
procedure get_activity(mset: mattypeset; var act_tot,vol:datatype);
procedure wr_gridfiles;

Implementation

uses R3Main03;

procedure error_std(idst:string;message:string);
begin
writeln('Error_std message from : ',idst);
writeln(message);
writeln('Press enter');
if LOG_file_is_open then close(LOG);
if RES_file_is_open then close(RES);
if press_enter_wanted then readln;
halt;
end;

procedure warning_std(idst:string;
                      message:string;
                      war:warningtype);
begin
inc(warning_table[war]);
warnings_were_issued:=true;
if ord(war)>ord(max_issued_warning_priority) then
  max_issued_warning_priority:=war;
if ord(war)>=ord(warning_priority_screen) then
  begin
    wr_line(output);
    writeln('Warning_std message from : ',idst,' (warning = ',war_string(war),')');
    writeln(message);
    wr_line(output);
  end;
if (LOG_file_is_open) and (ord(war)>=ord(warning_priority_log)) then
  begin
    wr_line(LOG);
    writeln(LOG,'Warning_std message from : ',idst,' (warning = ',war_string(war),')');
    writeln(LOG,message);
    wr_line(LOG);
  end;
end;


procedure wr_memory_status(var OM:text;st:string);
const ID='wr_memory_status';
begin
wr_line(OM);
writeln(OM,ID);
{$IFNDEF Delphi}
writeln(OM,'* MaxAvail = ',MaxAvail/1000/1000:12:4,' Mb '+st,' Nodes    = ',number_of_nodes:12);
{$ENDIF}
writeln(OM,'  imax     = ',imax:6   ,' jmax   = ',jmax:6   ,' kmax    = ',kmax:6);
writeln(OM,'  imaxTot  = ',imaxTot:6,' jmaxTot= ',jmaxTot:6,' kmaxTot = ',kmaxTot:6);

{$IFDEF Delphi}
(* These functions can be used to tell something about the memory status under Delphi *)
(*
writeln(OM,'TotalAddrSpace   = ',GetHeapStatus.TotalAddrSpace);
writeln(OM,'TotalUncommitted = ',GetHeapStatus.TotalUncommitted);
writeln(OM,'Totalcommitted   = ',GetHeapStatus.Totalcommitted);
writeln(OM,'Total allocated  = ',GetHeapStatus.TotalAllocated);
writeln(OM,'Total free       = ',GetHeapStatus.TotalFree);
writeln(OM,'Free small       = ',GetHeapStatus.FreeSmall);
writeln(OM,'Free big         = ',GetHeapStatus.FreeBig);
writeln(OM,'Unused           = ',GetHeapStatus.Unused);
writeln(OM,'Overhead         = ',GetHeapStatus.overhead);
writeln(OM,'HeapErrorCode    = ',GetHeapStatus.HeapErrorCode);
*)
{$ENDIF}
end;


procedure get_activity(mset: mattypeset; var act_tot,vol:datatype);
var conc,dum,dV:datatype;
    i:itype;
    j,jusemin,jusemax:jtype;
    k:ktype;
begin
act_tot:=0;
vol:=0;
if geometry=cartesian3D then
  begin
    jusemin:=1;
    jusemax:=jmax;
  end
else
  begin
    jusemin:=2;
    jusemax:=2;
  end;
for i:=1 to imax do
  for j:=jusemin to jusemax do
    for k:=1 to kmax do
      if (materials_def(i,j,k) in mset) then
          begin
            set_cvsize(i,j,k,dum,dum,dum,dum,dum,dum,dV);
            vol:=vol+dV;
            if (dV<>0) and (GP[i]^[j]^[k].valid_fieldvalue) then
              begin
                conc:=GP[i]^[j]^[k].c;
                act_tot:=act_tot+beta_def(i,j,k)*conc*dV;
                (* Observe:                                                   *)
                (* conc is the radon concentration in the air-filled pores    *)
                (* whereas act_tot is the total activity regardless of phase. *)
              end;
          end;
end; (* get_activity *)

procedure wr_material_volumes_etc(var OM:text);
const ID='wr_material_volumes_etc';
var i:itype;
    j,jusemin,jusemax:jtype;
    k:ktype;
    ArW,ArE,ArS,ArN,ArB,ArT,dV,conc:datatype;
    m:mattype;
    first:array[mattype] of boolean;
    num,N_invalid:array[mattype] of longint;
    act,avgconc,sumconc,minconc,maxconc,vol:array[mattype] of datatype;
    act_tot,avgconc_tot,sumconc_tot,vol_tot:datatype;
    i_minconc,i_maxconc:array[mattype] of itype;
    j_minconc,j_maxconc:array[mattype] of jtype;
    k_minconc,k_maxconc:array[mattype] of ktype;
    matset: set of mattype;
begin
wr_line(OM);
writeln(OM,ID,' (volume-averaged field values)');
if geometry=cartesian3D then
  begin
    jusemin:=1;
    jusemax:=jmax;
  end
else
  begin
    jusemin:=2;
    jusemax:=2;
  end;
vol_tot:=0;
act_tot:=0;
sumconc_tot:=0;
for m:=mat_undefined to mat_last do
  begin
    first[m]:=true;
    num[m]:=0;
    avgconc[m]:=0.0;
    sumconc[m]:=0.0;
    minconc[m]:=0.0;
    maxconc[m]:=0.0;
    i_minconc[m]:=1;
    j_minconc[m]:=1;
    k_minconc[m]:=1;
    i_maxconc[m]:=1;
    j_maxconc[m]:=1;
    k_maxconc[m]:=1;
    vol[m]:=0;
    act[m]:=0;
    N_invalid[m]:=0;
  end;
for i:=1 to imax do
  for j:=jusemin to jusemax do
    for k:=1 to kmax do
      begin
        set_cvsize(i,j,k,ArW,ArE,ArS,ArN,ArB,ArT,dV);
        m:=GP[i]^[j]^[k].mat;
        vol[m]:=vol[m]+dV;
        vol_tot:=vol_tot+dV;

        (* if (dV<>0) and (not GP[i]^[j]^[k].valid_fieldvalue) then *)
        if (dV=0) or (not GP[i]^[j]^[k].valid_fieldvalue) then (* Changed May 30, 1999 *)
          begin
            inc(N_invalid[m]); (* Changed June 25, 1998 from: error_node(ID,i,j,k,'dV<>0 and fieldvalue=false !'); *)
          end
        else
          begin
            conc:=GP[i]^[j]^[k].c;

            sumconc[m]:=sumconc[m]+conc*dV;
            sumconc_tot:=sumconc_tot+conc*dV;
            act[m]:=act[m]+beta_def(i,j,k)*conc*dV;
            act_tot:=act_tot+beta_def(i,j,k)*conc*dV;
            if (dV>0) then
              begin
                inc(num[m]);
                if (first[m]) or (conc<minconc[m]) then
                  begin
                    minconc[m]:=conc;
                    i_minconc[m]:=i;
                    j_minconc[m]:=j;
                    k_minconc[m]:=k;
                  end;
                if (first[m]) or (conc>maxconc[m]) then
                  begin
                    maxconc[m]:=conc;
                    i_maxconc[m]:=i;
                    j_maxconc[m]:=j;
                    k_maxconc[m]:=k;
                  end;
                if first[m] then first[m]:=false;
            end;
          end; (* valid point *)
      end;
if vol[mat_last]<>0 then
  error_std(ID,'Some volumes have been assigned to mat_last.');
for m:=mat_undefined to pred(mat_last) do
  if (vol[m]<>0) then avgconc[m]:=sumconc[m]/vol[m];
if (vol_tot<>0) then
  avgconc_tot:=sumconc_tot/vol_tot;

matset:=[];
for m:=mat_undefined to pred(mat_last) do
  if (num[m]+N_invalid[m]<>0) then matset:=matset+[m];


writeln(OM,'mat ':6,' ',
           'Avg(conc)':18,' ',
           'Activity':18,' ',
           'Volume':18,' ',
           'N':10,' ',
           'N_invalid':10);

for m:=mat_undefined to pred(mat_last) do
  if m in matset then
    writeln(OM,mat_string(m):6,' ',
               avgconc[m]:18,' ',
               act[m]:18,' ',
               vol[m]:18,' ',
               num[m]:10,' ',
               N_invalid[m]:10);

if geometry=cartesian3D then
  begin
    writeln(OM,'mat ':6,' ',
               'Min(conc)':18,' ',
               'i':3,' ',
               'j':3,' ',
               'k':3,' ',
               'x':12,' ',
               'y':12,' ',
               'z':12);
    for m:=mat_undefined to pred(mat_last) do
      if m in matset then
        writeln(OM,mat_string(m):6,' ',
                   minconc[m]:18,' ',
                   i_minconc[m]:3,' ',
                   j_minconc[m]:3,' ',
                   k_minconc[m]:3,' ',
                   xnod(i_minconc[m]):12,' ',
                   ynod(j_minconc[m]):12,' ',
                   znod(k_minconc[m]):12);

    writeln(OM,'mat ':6,' ',
               'Max(conc)':18,' ',
               'i':3,' ',
               'j':3,' ',
               'k':3,' ',
               'x':12,' ',
               'y':12,' ',
               'z':12);
    for m:=mat_undefined to pred(mat_last) do
      if m in matset then
        writeln(OM,mat_string(m):6,' ',
                   maxconc[m]:18,' ',
                   i_maxconc[m]:3,' ',
                   j_maxconc[m]:3,' ',
                   k_maxconc[m]:3,' ',
                   xnod(i_maxconc[m]):12,' ',
                   ynod(j_maxconc[m]):12,' ',
                   znod(k_maxconc[m]):12);

  end
else
  begin (* 2D only *)
    writeln(OM,'mat ':6,' ',
               'Min(conc)':18,' ',
               'i':3,' ',
               'k':3,' ',
               'x':12,' ',
               'z':12);
    for m:=mat_undefined to pred(mat_last) do
      if m in matset then
        writeln(OM,mat_string(m):6,' ',
                   minconc[m]:18,' ',
                   i_minconc[m]:3,' ',
                   k_minconc[m]:3,' ',
                   xnod(i_minconc[m]):12,' ',
                   znod(k_minconc[m]):12);


    writeln(OM,'mat ':6,' ',
               'Max(conc)':18,' ',
               'i':3,' ',
               'k':3,' ',
               'x':12,' ',
               'z':12);
    for m:=mat_undefined to pred(mat_last) do
      if m in matset then
        writeln(OM,mat_string(m):6,' ',
                   maxconc[m]:18,' ',
                   i_maxconc[m]:3,' ',
                   k_maxconc[m]:3,' ',
                   xnod(i_maxconc[m]):12,' ',
                   znod(k_maxconc[m]):12);

  end;
writeln(OM,'Total geometric volume      = ',vol_tot:20);
writeln(OM,'Total activity              = ',act_tot:20);
writeln(OM,'Overall mean concentration  = ',avgconc_tot:20);
end;

procedure wr_axis(var OM:text;dir:dirtype);
const ID='wr_axis';
var h,hmax:htype;
    w,dw,dcdwMax:axistype;
    fixst,st,sth:string;
    wFix:FixType;
begin
hmax:=0; (* Arbitary initialization just to keep compiler happy! *)
case dir of
  xdir: begin hmax:=imax; w:=x; dw:=dx; sth:='i'; st:='x'; dcdwMax:=dcdxMax end;
  ydir: begin hmax:=jmax; w:=y; dw:=dy; sth:='j'; st:='y'; dcdwMax:=dcdyMax end;
  zdir: begin hmax:=kmax; w:=z; dw:=dz; sth:='k'; st:='z'; dcdwMax:=dcdzMax end;
else
  error_std(ID,'Illegal dir.');
end;

writeln(OM,'axis',' ',
           sth:3,' ',
           st+'['+sth+']':11,' ',
           st+'['+sth+'+1]':11,' ',
           'd'+st+'['+sth+']':12,' ',
           'dcd'+st:12,' ',
           'dcd'+st+'norm':10,' ',
           'Fixpts':7);

for h:=1 to hmax do
  begin
    if FixID(dir,h,wFix) then
      FixSt:=FixName(wFix)
    else
      FixSt:='-';
    write(OM,st:4,' ',
             h:3,' ',
             w[h]:11:5,' ',
             w[h+1]:11:5,' ',
             dw[h]:12:7,' ',
             dcdwMax[h]:12,' ');
    if dcdwMaxAll<>0 then
      write(OM,dcdwMax[h]/dcdwMaxAll:10:8,' ')
    else
      write(OM,0.0:10:8);
    writeln(OM,FixSt:7);
  end;
end; (* wr_axis *)

procedure wr_axes_proc(var OM:text);
const ID='wr_axes_proc';
begin
wr_line(OM);
writeln(OM,ID);
wr_axis(OM,xdir);
writeln(OM);
wr_axis(OM,ydir);
writeln(OM);
wr_axis(OM,zdir);
writeln(OM);
end;

procedure wr_nodedata(var OM:text; i:itype;j:jtype;k:ktype; header:boolean);
const ID='wr_nodedata';
begin
wr_line(OM);
writeln(OM,ID);
if header then
  begin
    wr_line(OM);
    writeln(OM,'i':3,' ','j':3,' ','k':3,' ','c':12,' ',
               'nodetyp',' ',
               'west   ','east   ',
               'south  ','north  ',
               'bottom ','top    ');
    wr_line(OM);
  end;
with GP[i]^[j]^[k] do
writeln(OM,i:3,' ',j:3,' ',k:3,' ',c:12:3,' ',
nodetyp_string(nodetyp),' ',
nodecon_string(wcon),nodecon_string(econ),
nodecon_string(scon),nodecon_string(ncon),
nodecon_string(bcon),nodecon_string(tcon));
end;

procedure wr_cvsize(var OM:text);
const ID='wr_cvsize';
var ArW,ArE,ArS,ArN,ArB,ArT,dV:datatype;
    i:itype;j:jtype;k:ktype;
begin
wr_line(OM);
writeln(OM,ID);
for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      begin
        set_cvsize(i,j,k,ArW,ArE,ArS,ArN,ArB,ArT,dV);
        writeln(OM,i:3,' ',j:3,' ',k:3,' ArW=',ArW:12,' ArE=',ArE:12,
                                       ' ArS=',ArS:12,' ArN=',ArN:12,
                                       ' ArB=',ArB:12,' ArT=',ArT:12,
                                       ' dV=',dV:12);
      end;
end;

function xnod(i:itype):datatype;
begin
xnod:=x[i]+0.5*dx[i];
end;

function ynod(j:jtype):datatype;
begin
ynod:=y[j]+0.5*dy[j];
end;

function znod(k:ktype):datatype;
begin
znod:=z[k]+0.5*dz[k];
end;

procedure find_i_intpol(xp:datatype; var iA1,iB1,iA2,iB2:itype);
const ID='find_i_intpol';
var found:boolean;
    i:itype;
begin
found:=false;
for i:=1 to imax do
  if (not found) and (abs(xp-xnod(i))<=0) then (* First check if xp is exactly on top of a node ... *)
    begin
      found:=true;
      if i=1 then
        begin (* Unique *)
          iA1:=1; iB1:=2;
          iA2:=1; iB2:=2;
        end;
      if (2<=i) and (i<=imax-1) then
        begin (* Two solutions *)
          iA1:=i-1; iB1:=i;
          iA2:=i;   iB2:=i+1;
        end;
      if (i=imax) then
        begin (* Unique *)
          iA1:=imax-1; iB1:=imax;
          iA2:=imax-1; iB2:=imax;
        end;
    end;
if not found then (* ... then look between nodes *)
  for i:=1 to imax-1 do
    if (xnod(i)<xp) and (xp<xnod(i+1)) then
      begin
        found:=true;
        iA1:=i; iB1:=i+1;
        iA2:=i; iB2:=i+1;
      end;
if not found then
  begin
    writeln('x = ',xp);
    error_std(ID,'Could not find x on this axis !!')
  end;
end; (* find_i_intpol *)

procedure find_j_intpol(yp:datatype; var jA1,jB1,jA2,jB2:jtype);
const ID='find_j_intpol';
var found:boolean;
    j:jtype;
begin
found:=false;
for j:=1 to jmax do
  if (not found) and (abs(yp-ynod(j))<=0) then (* First check if yp is exactly on top of a node ... *)
    begin
      found:=true;
      if j=1 then
        begin (* Unique *)
          jA1:=1; jB1:=2;
          jA2:=1; jB2:=2;
        end;
      if (2<=j) and (j<=jmax-1) then
        begin (* Two solutions *)
          jA1:=j-1; jB1:=j;
          jA2:=j;   jB2:=j+1;
        end;
      if (j=jmax) then
        begin (* Unique *)
          jA1:=jmax-1; jB1:=jmax;
          jA2:=jmax-1; jB2:=jmax;
        end;
    end;
if not found then (* ... then look between nodes *)
  for j:=1 to jmax-1 do
    if (ynod(j)<yp) and (yp<ynod(j+1)) then
      begin
        found:=true;
        jA1:=j; jB1:=j+1;
        jA2:=j; jB2:=j+1;
      end;
if not found then
  begin
    writeln('y = ',yp);
    error_std(ID,'Could not find y on this axis !!')
  end;
end; (* find_j_intpol *)

procedure find_k_intpol(zp:datatype; var kA1,kB1,kA2,kB2:ktype);
const ID='find_k_intpol';
var found:boolean;
    k:ktype;
begin
found:=false;
for k:=1 to kmax do
  if (not found) and (abs(zp-znod(k))<=0) then (* First check if zp is exactly on top of a node ... *)
    begin
      found:=true;
      if k=1 then
        begin (* Unique *)
          kA1:=1; kB1:=2;
          kA2:=1; kB2:=2;
        end;
      if (2<=k) and (k<=kmax-1) then
        begin (* Two solutions *)
          kA1:=k-1; kB1:=k;
          kA2:=k;   kB2:=k+1;
        end;
      if (k=kmax) then
        begin (* Unique *)
          kA1:=kmax-1; kB1:=kmax;
          kA2:=kmax-1; kB2:=kmax;
        end;
    end;
if not found then (* ... then look between nodes *)
  for k:=1 to kmax-1 do
    if (znod(k)<zp) and (zp<znod(k+1)) then
      begin
        found:=true;
        kA1:=k; kB1:=k+1;
        kA2:=k; kB2:=k+1;
      end;
if not found then
  begin
    writeln('z = ',zp);
    error_std(ID,'Could not find z on this axis !!')
  end;
end; (* find_k_intpol *)

function find_i_cv(xp:datatype):itype;
const ID='find_i_cv';
(* Gives the number of the control volume where xp belongs *)
var ip:itype;
begin
if (xp<x[1]) or (xp>x[imax]) then
  begin
    writeln('xp = ',xp:14,' Valid range : ',x[1]:14,' to ',x[imax]:14);
    error_std(ID,'Out of range');
  end;
ip:=imax;
repeat
  ip:=ip-1
until ((x[ip]<=xp) and (xp<=x[ip+1]));
find_i_cv:=ip;
end;


function find_j_cv(yp:datatype):jtype;
const ID='find_j_cv';
(* Gives the number of the control volume where yp belongs *)
var jp:jtype;
begin
if (yp<y[1]) or (yp>y[jmax]) then
  begin
    writeln('yp = ',yp:14,' Valid range : ',y[1]:14,' to ',y[jmax]:14);
    error_std(ID,'Out of range');
  end;
jp:=jmax;
repeat
  jp:=jp-1
until ((y[jp]<=yp) and (yp<=y[jp+1]));
find_j_cv:=jp;
end;


function find_k_cv(zp:datatype):ktype;
const ID='find_k_cv';
(* Gives the number of the control volume where xp belongs *)
var kp:ktype;
begin
if (zp<z[1]) or (zp>z[kmax]) then
  begin
    writeln('zp = ',zp:14,' Valid range : ',z[1]:14,' to ',z[kmax]:14);
    error_std(ID,'Out of range');
  end;
kp:=kmax;
repeat
  kp:=kp-1
until ((z[kp]<=zp) and (zp<=z[kp+1]));
find_k_cv:=kp;
end;

function find_i_face(xp:datatype):itype;
const ID='find_i_face';
(* Observe if xp is between the x-coordinate of node no. i1 and i2
this function returns i=i1 although xp may be within control volume
i2. *)
var ip:itype;
begin
if (xp<xnod(1)) or (xp>xnod(imax)) then
  begin
    writeln('xp = ',xp:14,' Valid range : ',xnod(1):14,' to ',xnod(imax):14);
    error_std(ID,'Out of range');
  end;
ip:=imax;
repeat
  ip:=ip-1
until ((xnod(ip)<=xp) and (xp<=xnod(ip+1)));
find_i_face:=ip;
end;

function find_j_face(yp:datatype):jtype;
const ID='find_j_face';
(* Observe if yp is between the y-coordinate of node no. j1 and j2
this function returns j=j1 although yp may be within control volume
j2. *)
var jp:jtype;
begin
if (yp<ynod(1)) or (yp>ynod(jmax)) then
  begin
    writeln('yp = ',yp:14,' Valid range : ',ynod(1):14,' to ',ynod(jmax):14);
    error_std(ID,'Out of range');
  end;
jp:=jmax;
repeat
  jp:=jp-1
until ((ynod(jp)<=yp) and (yp<=ynod(jp+1)));
find_j_face:=jp;
end;

function find_k_face(zp:datatype):ktype;
const ID='find_k_face';
(* Observe if zp is between the z-coordinate of node no. k1 and k2
this function returns k=k1 although zp may be within control volume
k2. *)
var kp:ktype;
begin
if (zp<znod(1)) or (zp>znod(kmax)) then
  begin
    writeln('zp = ',zp:14,' Valid range : ',znod(1):14,' to ',znod(kmax):14);
    error_std(ID,'Out of range');
  end;
kp:=kmax;
repeat
  kp:=kp-1
until ((znod(kp)<=zp) and (zp<=znod(kp+1)));
find_k_face:=kp;
end;

function fieldvalue_old(xp,yp,zp:datatype):datatype;
const ID='fieldvalue_old';
var xA,xB,yA,yB,zA,zB,
    f1,f2,f3,f4,f5,f6,f7,f8,
    fp,
    t1,t2,t3:datatype;
    ip:itype; jp:jtype; kp:ktype;
begin
ip:=find_i_face(xp);
jp:=find_j_face(yp);
kp:=find_k_face(zp);

xA:=xnod(ip);
xB:=xnod(ip+1);
yA:=ynod(jp);
yB:=ynod(jp+1);
zA:=znod(kp);
zB:=znod(kp+1);

if (xB-xA)<>0 then t1:=(xp-xA)/(xB-xA) else t1:=0;
if (yB-yA)<>0 then t2:=(yp-yA)/(yB-yA) else t2:=0;
if (zB-zA)<>0 then t3:=(zp-zA)/(zB-zA) else t3:=0;

if ((GP[ip]^  [jp]^  [kp]  .nodetyp in [NOP,NodX]) or
    (GP[ip+1]^[jp]^  [kp]  .nodetyp in [NOP,NodX]) or
    (GP[ip]^  [jp+1]^[kp]  .nodetyp in [NOP,NodX]) or
    (GP[ip+1]^[jp+1]^[kp]  .nodetyp in [NOP,NodX]) or
    (GP[ip]^  [jp]^  [kp+1].nodetyp in [NOP,NodX]) or
    (GP[ip+1]^[jp]^  [kp+1].nodetyp in [NOP,NodX]) or
    (GP[ip]^  [jp+1]^[kp+1].nodetyp in [NOP,NodX]) or
    (GP[ip+1]^[jp+1]^[kp+1].nodetyp in [NOP,NodX])) then
      begin
        wr_nodedata(output,ip,  jp,  kp,true);
        wr_nodedata(output,ip+1,jp,  kp,false);
        wr_nodedata(output,ip,  jp+1,kp,false);
        wr_nodedata(output,ip+1,jp+1,kp,false);
        wr_nodedata(output,ip,  jp,  kp+1,false);
        wr_nodedata(output,ip+1,jp,  kp+1,false);
        wr_nodedata(output,ip,  jp+1,kp+1,false);
        wr_nodedata(output,ip+1,jp+1,kp+1,false);
        writeln('xp = ',xp:12);
        writeln('yp = ',yp:12);
        writeln('zp = ',zp:12);
        error_std(ID,'One or more of the coner points are undefined');
      end;
f1:=GP[ip]^  [jp]^  [kp].c;
f2:=GP[ip+1]^[jp]^  [kp].c;
f3:=GP[ip]^  [jp+1]^[kp].c;
f4:=GP[ip+1]^[jp+1]^[kp].c;
f5:=GP[ip]^  [jp]^  [kp+1].c;
f6:=GP[ip+1]^[jp]^  [kp+1].c;
f7:=GP[ip]^  [jp+1]^[kp+1].c;
f8:=GP[ip+1]^[jp+1]^[kp+1].c;

fp:=(1-t1)*(1-t2)*(1-t3)*f1+
    (t1)*  (1-t2)*(1-t3)*f2+
    (1-t1)*(t2)*  (1-t3)*f3+
    (t1)*  (t2)*  (1-t3)*f4+
    (1-t1)*(1-t2)*(t3)*  f5+
    (t1)*  (1-t2)*(t3)*  f6+
    (1-t1)*(t2)*  (t3)*  f7+
    (t1)*  (t2)*  (t3)*  f8;

fieldvalue_old:=fp;
end; (* fieldvalue *)


procedure interpolate(xp,yp,zp:datatype;
                      i1,i2,j1,j2,k1,k2:integer;
                      var fp:datatype;
                      var valid:boolean);
var xA,xB,yA,yB,zA,zB,
    f1,f2,f3,f4,f5,f6,f7,f8,
    t1,t2,t3:datatype;
begin
valid:=true;
fp:=dumreal;
if (i1>=i2) or (i1<1) or (i2>imax) then valid:=false;
if (j1>=j2) or (j1<1) or (j2>jmax) then valid:=false;
if (k1>=k2) or (k1<1) or (k2>kmax) then valid:=false;

if valid then
  begin (* Do not do these tests, if you are outside valid ranges ! *)
    if (GP[i1]^[j1]^[k1].valid_fieldvalue = false) or
       (GP[i2]^[j1]^[k1].valid_fieldvalue = false) or
       (GP[i1]^[j1]^[k2].valid_fieldvalue = false) or
       (GP[i2]^[j1]^[k2].valid_fieldvalue = false) or
       (GP[i1]^[j2]^[k1].valid_fieldvalue = false) or
       (GP[i2]^[j2]^[k1].valid_fieldvalue = false) or
       (GP[i1]^[j2]^[k2].valid_fieldvalue = false) or
       (GP[i2]^[j2]^[k2].valid_fieldvalue = false)
    then
       valid:=false;
  end;

if valid then
  begin (* Still valid after all tests ! *)
    xA:=xnod(i1);
    xB:=xnod(i2);
    yA:=ynod(j1);
    yB:=ynod(j2);
    zA:=znod(k1);
    zB:=znod(k2);

    if (xB-xA)<>0 then t1:=(xp-xA)/(xB-xA) else t1:=0;
    if (yB-yA)<>0 then t2:=(yp-yA)/(yB-yA) else t2:=0;
    if (zB-zA)<>0 then t3:=(zp-zA)/(zB-zA) else t3:=0;

    f1:=GP[i1]^[j1]^[k1].c;
    f2:=GP[i2]^[j1]^[k1].c;
    f3:=GP[i1]^[j2]^[k1].c;
    f4:=GP[i2]^[j2]^[k1].c;
    f5:=GP[i1]^[j1]^[k2].c;
    f6:=GP[i2]^[j1]^[k2].c;
    f7:=GP[i1]^[j2]^[k2].c;
    f8:=GP[i2]^[j2]^[k2].c;

    fp:=(1-t1)*(1-t2)*(1-t3)*f1+
        (t1)*  (1-t2)*(1-t3)*f2+
        (1-t1)*(t2)*  (1-t3)*f3+
        (t1)*  (t2)*  (1-t3)*f4+
        (1-t1)*(1-t2)*(t3)*  f5+
        (t1)*  (1-t2)*(t3)*  f6+
        (1-t1)*(t2)*  (t3)*  f7+
        (t1)*  (t2)*  (t3)*  f8;

  end;

end; (* interpolate *)

function fieldvalue(xp,yp,zp:datatype; var valid:boolean):datatype;
const ID='fieldvalue';
var iA1,iB1,iA2,iB2,icv:itype;
    jA1,jB1,jA2,jB2,jcv:jtype;
    kA1,kB1,kA2,kB2,kcv:ktype;
    found:boolean;
    res:datatype;
    dum,dV:datatype;
begin
icv:=find_i_cv(xp);
jcv:=find_j_cv(yp);
kcv:=find_k_cv(zp);
set_cvsize(icv,jcv,kcv,dum,dum,dum,dum,dum,dum,dV);
if (dV<>0) and (GP[icv]^[jcv]^[kcv].nodetyp=NOP) then
  begin
    valid:=false;
    if ord(warning_priority_screen)<=ord(war_other) then writeln('x = ',xp:14:7,'y = ',yp:14:7,'  z = ',zp:14:7);
    if ord(warning_priority_log)<=ord(war_other) then writeln(LOG,'x = ',xp:14:7,'y = ',yp:14:7,'  z = ',zp:14:7);
    warning_std(ID,'This point is within a NOP node of non-zero volume',
                war_other);
  end;
find_i_intpol(xp,iA1,iB1,iA2,iB2); (* NB: iB1=iA1+1 and iB2:=iA2+1 *)
find_j_intpol(yp,jA1,jB1,jA2,jB2);
find_k_intpol(zp,kA1,kB1,kA2,kB2);
found:=false;
(* if not GP[icv]^[jcv]^[kcv].valid_fieldvalue then
  begin
    valid:=false;
    found:=true;
    res:=dumreal;
  end; !!!!!!!!!!!!!!!! claus *)
(* Try four different squares. Normally these squares are all equal ! *)
if not found then interpolate(xp,yp,zp,iA1,iB1,jA1,jB1,kA1,kB1,res,found);
if not found then interpolate(xp,yp,zp,iA2,iB2,jA1,jB1,kA1,kB1,res,found);
if not found then interpolate(xp,yp,zp,iA1,iB1,jA2,jB2,kA1,kB1,res,found);
if not found then interpolate(xp,yp,zp,iA2,iB2,jA2,jB2,kA1,kB1,res,found);

if not found then interpolate(xp,yp,zp,iA1,iB1,jA1,jB1,kA2,kB2,res,found);
if not found then interpolate(xp,yp,zp,iA2,iB2,jA1,jB1,kA2,kB2,res,found);
if not found then interpolate(xp,yp,zp,iA1,iB1,jA2,jB2,kA2,kB2,res,found);
if not found then interpolate(xp,yp,zp,iA2,iB2,jA2,jB2,kA2,kB2,res,found);
if not found then
  begin
    found:=true;
    res:=GP[icv]^[jcv]^[kcv].c;
    if ord(warning_priority_screen)<=ord(war_interpolation) then
      writeln('x = ',xp:14:7,'y = ',yp:14:7,'  z = ',zp:14:7,' c = ',res:14);
    if ord(warning_priority_log)<=ord(war_interpolation) then
      writeln(LOG,'x = ',xp:14:7,'y = ',yp:14:7,'  z = ',zp:14:7,' c = ',res:14);
    warning_std(ID,'Could not find anything by bi-linear interpolation. '+cr+
                   'The field-value of the control volume is used.',
                   war_interpolation);
  end;
valid:=found;
fieldvalue:=res;
end; (* fieldvalue *)

procedure interpolate2D(xp,zp:datatype;
                        i1,i2,k1,k2:integer;
                        var fp:datatype;
                        var valid:boolean);
var xA,xB,zA,zB,
    f1,f2,f3,f4,
    t1,t3:datatype;
begin
valid:=true;
fp:=dumreal;
if (i1>=i2) or (i1<1) or (i2>imax) then valid:=false;
if (k1>=k2) or (k1<1) or (k2>kmax) then valid:=false;

if valid then
  begin (* Do not do these tests, if you are outside valid ranges ! *)
    if (GP[i1]^[2]^[k1].valid_fieldvalue = false) or
       (GP[i2]^[2]^[k1].valid_fieldvalue = false) or
       (GP[i1]^[2]^[k2].valid_fieldvalue = false) or
       (GP[i2]^[2]^[k2].valid_fieldvalue = false) then
       valid:=false;
  end;

if valid then
  begin (* Still valid after all tests ! *)
    xA:=xnod(i1);
    xB:=xnod(i2);
    zA:=znod(k1);
    zB:=znod(k2);

    if (xB-xA)<>0 then t1:=(xp-xA)/(xB-xA) else t1:=0;
    if (zB-zA)<>0 then t3:=(zp-zA)/(zB-zA) else t3:=0;

    f1:=GP[i1]^[2]^[k1].c;
    f2:=GP[i2]^[2]^[k1].c;
    f3:=GP[i1]^[2]^[k2].c;
    f4:=GP[i2]^[2]^[k2].c;

    fp:=(1-t1)*(1-t3)*f1+
        (t1)*  (1-t3)*f2+
        (1-t1)*(t3)*  f3+
        (t1)*  (t3)*  f4;
  end;
end; (* interpolate2D *)


function fieldvalue2D(xp,zp:datatype; var valid:boolean):datatype;
const ID='fieldvalue2D';
var iA1,iB1,iA2,iB2,icv:itype;
    kA1,kB1,kA2,kB2,kcv:ktype;
    found:boolean;
    res:datatype;
    dum,dV:datatype;
begin
icv:=find_i_cv(xp);
kcv:=find_k_cv(zp);
set_cvsize(icv,2,kcv,dum,dum,dum,dum,dum,dum,dV);
if (dV<>0) and (GP[icv]^[2]^[kcv].nodetyp=NOP) then
  begin
    valid:=false;
    if ord(warning_priority_screen)<=ord(war_interpolation) then writeln('x = ',xp:14:7,'  z = ',zp:14:7);
    if ord(warning_priority_log)<=ord(war_interpolation) then writeln(LOG,'x = ',xp:14:7,'  z = ',zp:14:7);
    warning_std(ID,'This point is within a NOP node of non-zero volume',
    war_interpolation);
  end;
find_i_intpol(xp,iA1,iB1,iA2,iB2); (* NB: iB1=iA1+1 and iB2:=iA2+1 *)
find_k_intpol(zp,kA1,kB1,kA2,kB2);
found:=false;
if not GP[icv]^[2]^[kcv].valid_fieldvalue then
  begin
    valid:=false;
    found:=true;
    res:=dumreal;
  end;
(* Try four different squares. Normally these squares are all equal ! *)
if not found then interpolate2D(xp,zp,iA1,iB1,kA1,kB1,res,found);
if not found then interpolate2D(xp,zp,iA2,iB2,kA1,kB1,res,found);
if not found then interpolate2D(xp,zp,iA1,iB1,kA2,kB2,res,found);
if not found then interpolate2D(xp,zp,iA2,iB2,kA2,kB2,res,found);
if not found then
  begin
    found:=true;
    res:=GP[icv]^[2]^[kcv].c;
    if ord(warning_priority_screen)<=ord(war_interpolation) then writeln('x = ',xp:14:7,'  z = ',zp:14:7);
    if ord(warning_priority_log)<=ord(war_interpolation) then writeln(LOG,'x = ',xp:14:7,'  z = ',zp:14:7);
    warning_std(ID,'Could not find anything by bi-linear interpolation. '+cr+
                   'The field-value of the control volume is used.',
                   war_interpolation);
  end;
valid:=found;
fieldvalue2D:=res;
end; (* fieldvalue2D *)

procedure get_fieldvalue(xp,dx,yp,dy,zp,dz:datatype;
                         var c,dc:datatype;
                         var valid:boolean);
const ID='get_fieldvalue';
      q=10; (* Resolution factor for numerical derivate calculations *)
var ddx,ddy,ddz,
    dc_ddx,dc_ddy,dc_ddz,
    x_minus,x_plus,
    y_minus,y_plus,
    z_minus,z_plus:datatype;
    ip:itype;
    jp:jtype;
    kp:ktype;
(* The field value is found by bi-linear interpolation as described
   in Numerical Recipes in Pascal by Press et al. pp. 108. The uncertainty
   of the field value is estimated from derivatives df_dx calculated
   numerically. Analytical expressions are avoided because of the
   problems of control volumes of zero thickness (located at fix points). *)
begin
if geometry=cylindrical2D then
  error_std(ID,'Use special 2D-cylindrical procedure !!');

ddx:=dx/q;
ddy:=dy/q;
ddz:=dz/q;

if ((dx<0) or (dy<0) or (dz<0)) then
  begin
  writeln('xp = ',xp:12,' dx = ',dx:12);
  writeln('yp = ',yp:12,' dy = ',dy:12);
  writeln('zp = ',zp:12,' dz = ',dz:12);
  error_std(ID,'dx, dy or dz <0')
  end;

if (xp<x[1]) or (xp>x[imax]) or
   (yp<y[1]) or (yp>y[jmax]) or
   (zp<z[1]) or (zp>z[kmax]) then
   begin
     writeln('xp = ',xp:14,' Valid range: ',x[1]:14,' to ',x[imax]:14);
     writeln('yp = ',yp:14,' Valid range: ',y[1]:14,' to ',y[jmax]:14);
     writeln('zp = ',zp:14,' Valid range: ',z[1]:14,' to ',z[kmax]:14);
     error_std(ID,'xp, yp or zp out or range !');
   end;

ip:=find_i_face(xp);
jp:=find_j_face(yp);
kp:=find_k_face(zp);

if ip>1 then x_minus:=max(xnod(ip-1),xp-0.5*ddx)
  else x_minus:=max(xnod(1),xp-0.5*ddx);
if ip<imax then x_plus:=min(xnod(ip+1),xp+0.5*ddx)
  else x_plus:=min(xnod(imax),xp+0.5*ddx);
if (x_minus=x_plus) then dc_ddx:=0
else
  dc_ddx:=(fieldvalue(x_plus,yp,zp,valid)-fieldvalue(x_minus,yp,zp,valid))/(x_plus-x_minus);

if jp>1 then y_minus:=max(ynod(jp-1),yp-0.5*ddy)
  else y_minus:=max(ynod(1),yp-0.5*ddy);
if jp<jmax then y_plus:=min(ynod(jp+1),yp+0.5*ddy)
  else y_plus:=min(ynod(jmax),yp+0.5*ddy);
if (y_minus=y_plus) then
  dc_ddy:=0
else
  dc_ddy:=(fieldvalue(xp,y_plus,zp,valid)-fieldvalue(xp,y_minus,zp,valid))/(y_plus-y_minus);

if kp>1 then z_minus:=max(znod(kp-1),zp-0.5*ddz)
  else z_minus:=max(znod(1),zp-0.5*ddz);
if kp<kmax then z_plus:=min(znod(kp+1),zp+0.5*ddz)
  else z_plus:=min(znod(kmax),zp+0.5*ddz);
if (z_minus=z_plus) then
  dc_ddz:=0
else
  dc_ddz:=(fieldvalue(xp,yp,z_plus,valid)-fieldvalue(xp,yp,z_minus,valid))/(z_plus-z_minus);

c:=fieldvalue(xp,yp,zp,valid);
dc:=sqrt(sqr(dc_ddx*dx)+sqr(dc_ddy*dy)+sqr(dc_ddz*dz));
end; (* get_fieldvalue *)

procedure get_fieldvalue2D(xp,dx,zp,dz:datatype; var c,dc:datatype; var valid:boolean);
const ID='get_fieldvalue2D';
      q=10;
var ddx,ddz,
    dc_ddx,dc_ddz,
    x_minus,x_plus,
    z_minus,z_plus:datatype;
    ip:itype;
    kp:ktype;
begin
valid:=true;
if (geometry<>cartesian2D) and (geometry<>cylindrical2D) then
  error_std(ID,'This procedure is only ok for the 2D geometries !!');

ddx:=dx/q;
ddz:=dz/q;

if ((dx<0) or (dz<0)) then
  begin
    writeln('xp = ',xp:12,' dx = ',dx:12);
    writeln('zp = ',zp:12,' dz = ',dz:12);
    error_std(ID,'dx or dz <0')
  end;

if (xp<x[1]) or (xp>x[imax]) or
   (zp<z[1]) or (zp>z[kmax]) then
   begin
     writeln('xp = ',xp:14,' Valid range: ',x[1]:14,' to ',x[imax]:14);
     writeln('zp = ',zp:14,' Valid range: ',z[1]:14,' to ',z[kmax]:14);
     error_std(ID,'xp or zp out or range !');
   end;

ip:=find_i_face(xp);
kp:=find_k_face(zp);

if ip>1 then x_minus:=max(xnod(ip-1),xp-0.5*ddx)
  else x_minus:=max(xnod(1),xp-0.5*ddx);
if ip<imax then x_plus:=min(xnod(ip+1),xp+0.5*ddx)
  else x_plus:=min(xnod(imax),xp+0.5*ddx);
if (x_minus=x_plus) then dc_ddx:=0
else
  dc_ddx:=(fieldvalue2D(x_plus,zp,valid)-fieldvalue2D(x_minus,zp,valid))/(x_plus-x_minus);

if kp>1 then z_minus:=max(znod(kp-1),zp-0.5*ddz)
  else z_minus:=max(znod(1),zp-0.5*ddz);
if kp<kmax then z_plus:=min(znod(kp+1),zp+0.5*ddz)
  else z_plus:=min(znod(kmax),zp+0.5*ddz);
if (z_minus=z_plus) then
  dc_ddz:=0
else
  dc_ddz:=(fieldvalue2D(xp,z_plus,valid)-fieldvalue2D(xp,z_minus,valid))/(z_plus-z_minus);

c:=fieldvalue2D(xp,zp,valid);
dc:=sqrt(sqr(dc_ddx*dx)+sqr(dc_ddz*dz));
end; (* get_fieldvalue2D *)

procedure get_avgfield(x1,x2,y1,y2,z1,z2,ddd:datatype; var c,dc:datatype);
const ID='get_avgfield';
var xx,yy,zz,val,sum,sum2:datatype;
    Nsum:longint;
    err:integer;
    valid:boolean;
begin
if wr_details then
  begin
    writeln(ID);
    writeln('x1 = ',x1:14,' x2 = ',x2:14);
    writeln('y1 = ',y1:14,' y2 = ',y2:14);
    writeln('z1 = ',z1:14,' z2 = ',z2:14);
  end;
err:=0;
if (x1>x2) or (x1<x[1]) or (x2>x[imax]) then err:=1;
if (y1>y2) or (y1<y[1]) or (y2>y[jmax]) then err:=10;
if (z1>z2) or (z1<z[1]) or (z2>z[kmax]) then err:=100;
if err<>0 then
  begin
    writeln('errcode   1: x-problem');
    writeln('errcode  10: y-problem');
    writeln('errcode 100: z-problem');
    writeln('error = ',err,' (see above for explanation)');
    error_std(ID,'(x1>x2) or (x1<x[1]) or (x2>x[imax]) or similar'+
                 'problem for y or z.');
  end;

if (ddd<=0) then
  error_std(ID,'ddd<=0 !');

Nsum:=0;
sum:=0;
sum2:=0;
xx:=x1;
repeat
  yy:=y1;
  repeat
    zz:=z1;
    repeat
      val:=Fieldvalue(xx,yy,zz,valid);
      if not valid then
        begin
          writeln('Problem at x = ',xx:14:7,' y = ',yy:14:7,' z = ',zz:14:7);
          error_std(ID,'NodeTyp=NOP and volume<>0. '+
                       'Do not try to average over this point');
        end
      else
        begin (* valid point *)
          inc(Nsum);
          sum:=sum+val;
          sum2:=sum2+sqr(val);
        end;
      zz:=zz+ddd;
    until zz>z2;
    yy:=yy+ddd;
  until yy>y2;
  xx:=xx+ddd;
until xx>x2;
if Nsum>0 then c:=sum/Nsum else c:=0;
if Nsum>1 then dc:=sqrt((sum2-Nsum*sqr(c))/(Nsum-1)) else dc:=0;
end;  (* get_avgfield *)

procedure get_avgfield2D(x1,x2,z1,z2,ddd:datatype; var c,dc:datatype);
const ID='get_avgfield2D';
var xx,zz,val,dum,dV,sum,sum2:datatype;
    Nsum:longint;
    err:integer;
    icv:itype;
    jcv:jtype;
    kcv:ktype;
    valid:boolean;
begin
if wr_details then
  begin
    writeln(ID);
    writeln('x1 = ',x1:14,' x2 = ',x2:14);
    writeln('z1 = ',z1:14,' z2 = ',z2:14);
  end;
err:=0;
if (x1>x2) or (x1<x[1]) or (x2>x[imax]) then err:=1;
if (z1>z2) or (z1<z[1]) or (z2>z[kmax]) then err:=100;
if err<>0 then
  begin
    writeln('errcode   1: x-problem');
    writeln('errcode 100: z-problem');
    writeln('error = ',err,' (see above for explanation)');
    error_std(ID,'(x1>x2) or (x1<x[1]) or (x2>x[imax]) or similar'+
                 'problem for z.');
  end;
Nsum:=0;
sum:=0;
sum2:=0;
xx:=x1;
repeat
  zz:=z1;
  repeat
    icv:=find_i_cv(xx);
    jcv:=2;
    kcv:=find_k_cv(zz);
    set_cvsize(icv,jcv,kcv,dum,dum,dum,dum,dum,dum,dV);
    if (GP[icv]^[jcv]^[kcv].nodetyp=NOP) and (dV<>0) then
      begin
        writeln('Problem at x = ',xx:14:7,' z = ',zz:14:7);
        error_std(ID,'NodeTyp=NOP and volume<>0. '+
                     'Do not try to average over this point');
      end;
    val:=Fieldvalue2D(xx,zz,valid);
    if not valid then
      begin
        writeln('Problem at x = ',xx:14:7,' z = ',zz:14:7);
        error_std(ID,'NodeTyp=NOP and volume<>0. '+
                     'Do not try to average over this point');
      end;
    inc(Nsum);
    sum:=sum+val;
    sum2:=sum2+sqr(val);
    zz:=zz+ddd;
  until zz>z2;
xx:=xx+ddd;
until xx>x2;
if Nsum>0 then c:=sum/Nsum else c:=0;
if Nsum>1 then dc:=sqrt((sum2-Nsum*sqr(c))/(Nsum-1)) else dc:=0;
end;  (* get_avgfield2D *)


procedure wr_gridfiles;
(* These procedures write grid files needed to visualize the computations *)
(* For example, the data may be read and plotted with the S-plus script   *)
(* called RnMod3d_readgrid.                                               *)
(* This procedure is not decriober in the Users Guide                     *)
(* What to do: Simply include the statement wr_gridfiles in the job file  *)
(* and four extra files will be created.                                  *)
(* Date: September 1, 1999                                                *)
(* Revised: September 1, 1999                                             *)
const ID = 'wr_gridfiles';

procedure wr_gridfile1;
(* Overall values *)
var GRD:text;
    filename:string;
begin
filename:='f'+runid+'_'+'gr1.dat';
assign(GRD,filename);
write(' '+filename+' ');
rewrite(GRD);
writeln(GRD,'imax',', ',imax:5);
writeln(GRD,'jmax',', ',jmax:5);
writeln(GRD,'kmax',', ',kmax:5);
writeln(GRD,'geometry',', ',ord(geometry):5);
close(GRD);
end;

procedure wr_gridfile2;
(* Grid lines *)
var GRD:text;
    filename:string;
    i:itype; j:jtype; k:ktype;
begin
filename:='f'+runid+'_'+'gr2.dat';
assign(GRD,filename);
write(' '+filename+' ');
rewrite(GRD);
writeln(GRD,'dir':3,', ','h':5,', ','w':15);
for i:=1 to imax do
  writeln(GRD,'i':1,',   ',i:5,', ',x[i]:15);
for j:=1 to jmax do
  writeln(GRD,'j':1,',   ',j:5,', ',y[j]:15);
for k:=1 to kmax do
  writeln(GRD,'k':1,',   ',k:5,', ',z[k]:15);
close(GRD);
end;

procedure wr_gridfile3;
(* Fix lines *)
var GRD:text;
    filename:string;
    i:itype; j:jtype; k:ktype;
    wFix:FixType;
begin
filename:='f'+runid+'_'+'gr3.dat';
assign(GRD,filename);
write(' '+filename+' ');
rewrite(GRD);
writeln(GRD,'dir':3,', ','h':5,', ','w':19);
for i:=1 to imax do
  if FixID(xdir,i,wFix) then
    writeln(GRD,'i':1,',   ',i:5,', ',x[i]:19);
for j:=1 to jmax do
  if FixID(ydir,j,wFix) then
    writeln(GRD,'j':1,',   ',j:5,', ',y[j]:19);
for k:=1 to kmax do
  if FixID(zdir,k,wFix) then
    writeln(GRD,'k':1,',   ',k:5,', ',z[k]:19);
close(GRD);
end;

procedure wr_gridfile4;
(* Control volumes *)
var GRD:text;
    filename:string;
    i:itype; j:jtype; k:ktype;
    no:longint;
begin
filename:='f'+runid+'_'+'gr4.dat';
assign(GRD,filename);
write(' '+filename+' ');
rewrite(GRD);
writeln(GRD,'no':8,', ',
            'i':5,', ',
            'j':5,', ',
            'k':5,', ',
            'xmin':19,', ',
            'xmax':19,', ',
            'ymin':19,', ',
            'ymax':19,', ',
            'zmin':19,', ',
            'zmax':19,', ',
            'xnod':19,', ',
            'ynod':19,', ',
            'znod':19,', ',
            'nodetyp':19,', ',
            'mat':19,', ',
            'c':19);
for k:=1 to kmax do
  for j:=1 to jmax do
    for i:=1 to imax do
      begin
        no:=i+(j-1)*imax+(k-1)*imax*jmax;
        writeln(GRD,no:8,', ',
                    i:5,', ',
                    j:5,', ',
                    k:5,', ',
                    x[i]:19,', ',
                    x[i+1]:19,', ',
                    y[j]:19,', ',
                    y[j+1]:19,', ',
                    z[k]:19,', ',
                    z[k+1]:19,', ',
                    xnod(i):19,', ',
                    ynod(j):19,', ',
                    znod(k):19,', ',
                    nodetyp_string(GP[i]^[j]^[k].nodetyp),', ',
                    ord(materials_def(i,j,k)):19,', ',
                    GP[i]^[j]^[k].c:19)
      end;
close(GRD);
end;
begin (* wr_gridfiles *)
write(ID+':');
wr_gridfile1;
wr_gridfile2;
wr_gridfile3;
wr_gridfile4;
writeln;
end;

begin
end.