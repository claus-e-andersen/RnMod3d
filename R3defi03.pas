unit R3Defi03;
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

Interface
const programid      = 'RnMod3d';
      versionid      = 'Version 0.8 (Sep. 15, 1997 - July 18, 2000)';
      description    = 'Radon and soil gas transport model';
      programmer1    = 'Claus E. Andersen, Risoe National Laboratory, ';
      programmer2    = 'DK-4000 Roskilde, Denmark. claus.andersen@risoe.dk';
      documentation  = 'User''s Guide to RnMod3d, Risoe-R-1201(EN)';
      copyright      = 'Copyright, Risoe National Laboratory, Denmark';

      {$IfDEF imax3}   imaxTot=  3; {$ENDIF} (* Max number of nodes on x-axis *)
      {$IfDEF imax10}  imaxTot= 10; {$ENDIF}
      {$IfDEF imax50}  imaxTot= 50; {$ENDIF}
      {$IfDEF imax100} imaxTot=100; {$ENDIF}
      {$IfDEF imax150} imaxTot=150; {$ENDIF}
      {$IfDEF imax200} imaxTot=200; {$ENDIF}
      {$IfDEF imax250} imaxTot=250; {$ENDIF}
      {$IfDEF imax300} imaxTot=300; {$ENDIF}
      {$IfDEF imax350} imaxTot=350; {$ENDIF}
      {$IfDEF imax400} imaxTot=400; {$ENDIF}
      {$IfDEF imax450} imaxTot=450; {$ENDIF}
      {$IfDEF imax500} imaxTot=500; {$ENDIF}

      {$IfDEF jmax3}   jmaxTot= 3;  {$ENDIF} (* Max number of nodes on y-axis *)
      {$IfDEF jmax10}  jmaxTot= 10; {$ENDIF}
      {$IfDEF jmax50}  jmaxTot= 50; {$ENDIF}
      {$IfDEF jmax100} jmaxTot=100; {$ENDIF}
      {$IfDEF jmax150} jmaxTot=150; {$ENDIF}
      {$IfDEF jmax200} jmaxTot=200; {$ENDIF}
      {$IfDEF jmax250} jmaxTot=250; {$ENDIF}
      {$IfDEF jmax300} jmaxTot=300; {$ENDIF}
      {$IfDEF jmax350} jmaxTot=350; {$ENDIF}
      {$IfDEF jmax400} jmaxTot=400; {$ENDIF}
      {$IfDEF jmax450} jmaxTot=450; {$ENDIF}
      {$IfDEF jmax500} jmaxTot=500; {$ENDIF}

      {$IfDEF kmax3}   kmaxTot=  3; {$ENDIF} (* Max number of nodes on z-axis *)
      {$IfDEF kmax10}  kmaxTot= 10; {$ENDIF}
      {$IfDEF kmax50}  kmaxTot= 50; {$ENDIF}
      {$IfDEF kmax100} kmaxTot=100; {$ENDIF}
      {$IfDEF kmax150} kmaxTot=150; {$ENDIF}
      {$IfDEF kmax200} kmaxTot=200; {$ENDIF}
      {$IfDEF kmax250} kmaxTot=250; {$ENDIF}
      {$IfDEF kmax300} kmaxTot=300; {$ENDIF}
      {$IfDEF kmax350} kmaxTot=350; {$ENDIF}
      {$IfDEF kmax400} kmaxTot=400; {$ENDIF}
      {$IfDEF kmax450} kmaxTot=450; {$ENDIF}
      {$IfDEF kmax500} kmaxTot=500; {$ENDIF}

      hmaxTot        = imaxTot+jmaxTot+kmaxTot;

      cr             = chr(10)+chr(13);  (* Carriage return character *)
      dumReal        = -9.99e1;          (* Dummy real    *)
      dumInt         = -99;              (* Dummy integer *)

      min_flxcalcno  = 3; (* First calculations of flux changes are invalid *)
      min_obscalcno  = 3; (* First calculations of obs changes are invalid *)
      residual_sum_warning_limit = 100;  (* Limit used by war_residual *)
      Lambda_Rn222   = 2.09838e-6;       (* Decay constant for Rn-222 *)

type datatype = extended; (* Main data type for all floats *)

     itype    = 1..imaxTot;
     jtype    = 1..jmaxTot;
     ktype    = 1..kmaxTot;

     htype    = 0..hmaxTot; (* Generic index coordinate (i,j,k) *)

     axistype = array[htype] of datatype;

     (* Fix points *)
     fixtype =  (xfix1,xfix2,xfix3,xfix4,xfix5,xfix6,xfix7,
                 yfix1,yfix2,yfix3,yfix4,yfix5,yfix6,yfix7,
                 zfix1,zfix2,zfix3,zfix4,zfix5,zfix6,zfix7,wfixLast);

     (* Directions *)
     dirtype=(east,west,north,south,top,bottom,xdir,ydir,zdir,nodir);

     (* Nodes *)
     nodetyptype=(nop,free,fixed1,fixed2,fixed3,fixed4,fixed5,NodX);

     (* Connectors *)
     nodecontype=(nill,std,noflow,ConX);

     (* Materials *)
     mattype=(mat_undefined,
              mat1,mat2,mat3,mat4,mat5,
              mat6,mat7,mat8,mat9,mat10,
              mat11,mat12,mat13,mat14,mat15,
              mat_last);
     mattypeset=set of mattype;

     (* Flux probes *)
     Flxtype=(Flx1,Flx2,Flx3,Flx4,Flx5,Flx_last);

     (* Field value probes *)
     Obstype=(Obs1,Obs2,Obs3,Obs4,Obs5,Obs_last);

     (* Plot files *)
     pltfiletype=(plt1,plt2,plt3,plt_last);

     (* Geometry *)
     geometrytype=(cartesian3D,cartesian2D,cylindrical2D);

     (* Schemes *)
     schemetype=(powerlaw,central,upwind,hybrid,exact);

     (* In-functions *)
     regtype=(inside,outside,eqA,eqB,eqAB);
     regsettype=set of regtype;

     (* Type needed in the set_axis procedures *)
     funcidtype=(FocusA,FocusB);

     (* Type needed for flux calculations *)
     signtype=(plus,minus);

     (* Actions related to files *)
     filehandletype=(none,import,export,import_from_qbuf,export_to_qbuf);

     (* Type of solutions *)
     solutiontype=(steady,unsteady);

     (* Buffers *)
     use_fieldbuffer_type=(cBUF1,cBUF2,no_cBUF);

     (* Tag for each run *)
     runidtype=string[4];
     runtitletype=string;

     (* Warnings *)
     warningtype=(war_none,
                  war_first,
                  war_interpolation,
                  war_other,
                  war_fileimport,
                  war_convergence,
                  war_residual,
                  war_last);

     (* Data type used in GP for each single control volume *)
     nodedatatype=record
                    c:datatype;
                    aE,aW,aN,aSS,aT,aB,b,ap,ap_old_dt:datatype;
                    qE,qW,qN,qS,qT,qB:datatype;
                    nodetyp:nodetyptype;
                    Wcon,Econ,Scon,Ncon,Bcon,Tcon:nodecontype;
                    mat:mattype;
                    valid_fieldvalue:boolean;
                  end;

     (* Types needed for putting together GP *)
     columntype=array[ktype] of nodedatatype;
     platetype=array[jtype] of ^columntype;
     gridtype=array[itype] of ^platetype;

     (* Fields *)
     fieldcolumntype=array[ktype] of record
                                       c,
                                       ap_old_dt:datatype;
                                     end;
     fieldplatetype=array[jtype] of ^fieldcolumntype;
     fieldtype=record
                 imax:itype;
                 jmax:jtype;
                 kmax:ktype;
                 GP:array[itype] of ^fieldplatetype;
               end;

     (* Flow fields *)
     flowfieldcolumntype=array[ktype] of record
                                           qE,qW,qN,qS,qT,qB:datatype;
                                         end;
     flowfieldplatetype=array[jtype] of ^flowfieldcolumntype;
     flowfieldtype=record
                     imax:itype;
                     jmax:jtype;
                     kmax:ktype;
                     q:array[itype] of ^flowfieldplatetype;
                   end;

     (* Types of functions and procedures to which pointers should point *)
     Afunctype=function(dir:dirtype;i:itype;j:jtype;k:ktype):datatype;
     Bfunctype=function(i:itype;j:jtype;k:ktype):datatype;
     Cfunctype=function(i:itype;j:jtype;k:ktype):mattype;
     Dfunctype=function:boolean;
     Aproctype=procedure;
     aijkproctype=procedure(i:itype;j:jtype;k:ktype);

     {$IFDEF Delphi}
     timetype = TDateTime;
     {$ELSE}
     timetype = datatype;
     {$ENDIF}

     
     vectype = array[htype] of datatype;


var imax:itype; (* Current max of nodes on x-axis *)
    jmax:jtype; (* Current max of nodes on y-axis *)
    kmax:ktype; (* Current max of nodes on z-axis *)

    tim:datatype; (* Global time *)

    iter:longint; (* Number of iterations *)

    GP:gridtype;  (* Main data structure *)

    (* Buffers for GP *)
    cBUF1v,
    cBUF2v:fieldtype;
    cBUF1_has_been_created,
    cBUF2_has_been_created:boolean;

    (* Buffer for the flow field *)
    qBUF:flowfieldtype;
    qBUF_has_been_created:boolean;

    (* Fix points *)
    wFixVal:array[FixType] of record
                                 defined:boolean; (* Flag                 *)
                                 w:datatype; (* Physical coordinate (xyz) *)
                                 h:htype;    (* Index coordinate (ijk)    *)
                               end;

    (* Variables related to the grid *)
    number_of_nodes:longint;
    x,y,z,
    dx,dy,dz,
    dcdxMax,dcdyMax,dcdzMax:axistype;
    dcdwMaxAll:datatype;

    (* Array with fixed-value boundary conditions *)
    cBC:array[nodetypType] of datatype;

    (* Flux probes *)
    flxval,
    flxval_old,
    flxval_change:array[flxtype] of record j,q:datatype end;

    (* Field-value probes *)
    Obsval,
    Obsval_old,
    Obsval_change:array[obstype] of datatype;

    (* Variables related to residuals *)
    residual_b_sum,
    residual_sum,
    residual_sum_change,
    residual_max,
    residual_max_change:datatype;
    i_residual_max:itype;
    j_residual_max:jtype;
    k_residual_max:ktype;

    (* Variables that keeps track of the computations *)
    first_run:boolean;
    unsteady_has_been_started:boolean;
    max_issued_warning_priority:warningtype;
    t1,t2:timetype;
    grid_def_last:aproctype;
    geometry_original:geometrytype;

    (* Warnings *)
    warnings_were_issued:boolean;
    warning_table:array[warningtype] of longint;

    (* Variables related to files *)
    LOG,RES:text; (* Main file variables *)
    LOG_file_is_open,
    RES_file_is_open:boolean;
    log_name,
    res_name,
    temp_field_name:string;
    pltfileset:set of pltfiletype;
    pltfiledirs:array[pltfiletype] of dirtype;

    (* Block matrices needed by solver *)
    AA,BB,CC,DD:vectype;

    (* Control variables in same order as in the Users Guide *)
    runid:runidtype;
    runtitle:runtitletype;
    solution:solutiontype;
    geometry:geometrytype;
    Ly:datatype;
    grid_def:aproctype;
    force_new_grid_in_every_run:boolean;
    boundary_conditions_def:aijkproctype;
    flux_def:aijkproctype;
    probe_def:aproctype;
    materials_def:Cfunctype;
    e_def:Bfunctype;
    beta_def:Bfunctype;
    G_def:Bfunctype;
    lambda_def:Bfunctype;
    D_def:Afunctype;
    initialfield_def:Bfunctype;
    import_initialfield:boolean;
    import_finalfield_guess:boolean;
    export_field:boolean;
    use_fieldbuffer:use_fieldbuffer_type;
    flowfield:filehandletype;
    flowfactor:datatype;
    import_field_name:string;
    export_field_name:string;
    flowfield_name:string;
    plotfiles_def:aijkproctype;
    user_procedure_each_iter_def:aproctype;
    wr_details:boolean;
    wr_main_procedure_id:boolean;
    wr_all_procedure_id:boolean;
    wr_iteration_line_log:boolean;
    wr_iteration_line_screen:boolean;
    wr_residual_during_calc_log:boolean;
    wr_residual_during_calc_screen:boolean;
    wr_flux_during_calc_log:boolean;
    wr_flux_during_calc_screen:boolean;
    wr_probes_during_calc_log:boolean;
    wr_probes_during_calc_screen:boolean;
    wr_final_results_log:boolean;
    wr_final_results_screen:boolean;
    wr_axes:boolean;
    wr_nodes:boolean;
    wr_node_numbers:boolean;
    wr_node_sizes:boolean;
    wr_coefficients:boolean;
    wr_materials_volumes:boolean;
    warning_priority_log:warningtype;
    warning_priority_screen:warningtype;
    solver_def:aproctype;
    scheme:schemetype;
    relax_factor:datatype;
    flux_convset:set of flxtype;
    probe_convset:set of obstype;
    conv_evaluation_period:longint;
    min_iterations:longint;
    max_iterations:longint;
    max_time:datatype;
    max_change:datatype;
    max_residual_sum:datatype;
    dtim :datatype;
    BC_running:boolean;
    BC_running_update_of_cBCs_def:AProcType;
    BC_running_min_iterations:longint;
    BC_running_max_residual_sum_before_new_BC:datatype;
    BC_running_convergence_def:DFuncType;
    wr_BC_running_messages_log:boolean;
    wr_BC_running_messages_screen:boolean;
    press_enter_wanted:boolean;

const FixedBCs=[NOP..NodX]-[NOP,free,NodX];

function FixDir(FV:FixType):dirtype;
function FixName(FV:FixType):string;
procedure check_wFix(IDst:string; wFix:FixType);
procedure check_fix_region_param(id0:string;xFixA,xFixB,yFixA,yFixB,zFixA,zFixB:FixType);

function in_interval_ijk(h,hA,hB:htype;regset:regsettype):boolean;
function in_interval(h:htype;wFixA,wFixB:FixType;wreg:regsettype):boolean;

function in_region_ijk(i,iA,iB:itype; ireg:regsettype;
                       j,jA,jB:jtype; jreg:regsettype;
                       k,kA,kB:ktype; kreg:regsettype):boolean;
function in_region(i:itype;xFixA,xFixB:FixType;xreg:regsettype;
                   j:jtype;yFixA,yFixB:FixType;yreg:regsettype;
                   k:ktype;zFixA,zFixB:FixType;zreg:regsettype):boolean;
function in_plane_ijk(reg:regsettype;
                      i,iA,iB:itype;
                      j,jA,jB:jtype;
                      k,kA,kB:ktype):boolean;
function in_plane(reg:regsettype;
                  i:itype;xFixA,xFixB:FixType;
                  j:jtype;yFixA,yFixB:FixType;
                  k:ktype;zFixA,zFixB:FixType):boolean;
function in_cube_ijk(reg:regsettype;
                      i,iA,iB:itype;
                      j,jA,jB:jtype;
                      k,kA,kB:ktype):boolean;
function in_cube(reg:regsettype;
                 i:itype;xFixA,xFixB:FixType;
                 j:jtype;yFixA,yFixB:FixType;
                 k:ktype;zFixA,zFixB:FixType):boolean;

procedure estimate_fieldvalues_at_lost_NOPs;

procedure set_control_variables_to_defaults;
procedure default_problem;
procedure grid_default;
procedure boundary_conditions_default(i:itype;j:jtype;k:ktype);
procedure fluxes_default(i:itype;j:jtype;k:ktype);
procedure probes_default;
function  e_default(i:itype;j:jtype;k:ktype):datatype;
function  G_default(i:itype;j:jtype;k:ktype):datatype;
function  Lambda_default(i:itype;j:jtype;k:ktype):datatype;

procedure wr_line(var OM:text);
procedure wr_header(var OM:text);
function file_exists(fname:string):boolean;
procedure wr_warning_table(var OM:text);

procedure set_FixVal(wFix:Fixtype;wval:datatype);
procedure set_axis_single(wFixA,wFixB:FixType;
                          hAB:htype;
                          f:funcidtype;
                          pow:datatype);
procedure set_axis_double(wFixA,wFixB:FixType;
                          hAM,hMB:htype;
                          fA,fB:funcidtype;
                          powA,powB,wdiv:datatype);
procedure set_axis_triple(wFixA,wFixB:FixType;
                          hAM,hMM,hMB:htype;
                          fA,fM,fB:funcidtype;
                          powA,powM,powB,wdivA,wdivB:datatype);
procedure check_axes;
function DirName(Dir:DirType):string;

function maxi(i1,i2:integer):integer;
function max(x1,x2:datatype):datatype;
function min(x1,x2:datatype):datatype;
function pw(x,y:datatype):datatype;
function cosech(x:datatype):datatype;



implementation

uses R3Main03,R3Writ03;

procedure check_wFix(IDst:string; wFix:FixType);
const ID='check_wFix';
      Syntax='wFixType=xFix1,..yFix2,..zFix1..wFixLast';
var err:word;
begin
err:=0;
if not( (ord(xFix1)<ord(yFix1)) and (ord(yFix1)<ord(zFix1)) ) then err:=1;
if (ord(wFix)<ord(xFix1)) then err:=2;
if ord(wFix)>=ord(wFixLast) then err:=3;
if err<>0 then
  begin
    writeln('* Error code = ',err:3);
    writeln('* Correct syntax :',syntax);
    writeln('* Check was invoked through ',idst);
    error_std(ID,'wFixType decleration is incorrect');
  end;
end;


function FixDir(FV:FixType):dirtype;
const ID='FixDir';
var dir:dirtype;
begin
dir:=nodir;
check_wFix(id,FV);
if (ord(FV)<ord(yFix1)) then dir:=xdir;
if((ord(yFix1)<=ord(FV)) and (ord(FV)<ord(zFix1))) then dir:=ydir;
if ord(zFix1)<=ord(FV) then dir:=zdir;
FixDir:=dir;
end;

function FixName(FV:FixType):string;
const ID='FixValName';
var st,nost:string;
    no:integer;
begin
no:=0;
check_wFix(id,FV);
case FixDir(FV) of
  xdir: begin st:='xFix'; no:=ord(FV)-ord(xFix1)+1 end;
  ydir: begin st:='yFix'; no:=ord(FV)-ord(yFix1)+1 end;
  zdir: begin st:='zFix'; no:=ord(FV)-ord(zFix1)+1 end;
else
  error_std(ID,'Unknown FixDir(FV)');
end;
str(no,nost);
FixName:=st+nost;
end;


procedure check_fix_region_param(id0:string;
                                 xFixA,xFixB,yFixA,yFixB,zFixA,zFixB:FixType);
const ID='check_fix_region_param';
begin
if not wFixVal[xFixA].defined then
  error_std(ID,'Undefined fixpoint '+Fixname(xFixA));
if not wFixVal[xFixB].defined then
  error_std(ID,'Undefined fixpoint '+Fixname(xFixB));
if not wFixVal[yFixA].defined then
  error_std(ID,'Undefined fixpoint '+Fixname(yFixA));
if not wFixVal[yFixB].defined then
  error_std(ID,'Undefined fixpoint '+Fixname(yFixB));
if not wFixVal[zFixA].defined then
  error_std(ID,'Undefined fixpoint '+Fixname(zFixA));
if not wFixVal[zFixB].defined then
  error_std(ID,'Undefined fixpoint '+Fixname(zFixB));

if not ((FixDir(xFixA)=xdir) and (FixDir(xFixB)=xdir)) then
  error_std(ID,'Not x-axis fixpoints : '+Fixname(xFixA)+' '+Fixname(xFixB));

if not ((FixDir(yFixA)=ydir) and (FixDir(yFixB)=ydir)) then
  error_std(ID,'Not y-axis fixpoints : '+Fixname(yFixA)+' '+Fixname(yFixB));

if not ((FixDir(zFixA)=zdir) and (FixDir(zFixB)=zdir)) then
  error_std(ID,'Not z-axis fixpoints : '+Fixname(zFixA)+' '+Fixname(zFixB));
end; (* check_fix_region_param*)


function in_interval_ijk(h,hA,hB:htype;regset:regsettype):boolean;
const ID='in_interval_ijk';
var res:boolean;
begin
if (hA>hB) then
begin
  writeln('hA = ',hA,'  hB = ',hB);
  writeln(LOG,'hA = ',hA,'  hB = ',hB);
error_std(ID,'hA>hB!');
end;
if not (regset-[inside,outside,eqA,eqB,eqAB]=[]) then
  error_std(ID,'Unknown reg (region)');
res:=false;
if (inside in regset) and (hA<h) and (h<hB) then res:=true;
if (outside in regset) and ((h<hA) or (hB<h)) then res:=true;
if ((eqA in regset) or (eqAB in regset)) and (h=hA) then res:=true;
if ((eqB in regset) or (eqAB in regset)) and (h=hB) then res:=true;
in_interval_ijk:=res;
end;

function in_interval(h:htype;wFixA,wFixB:FixType;wreg:regsettype):boolean;
const ID='in_interval';
var hA,hB:htype;
begin
if not wFixVal[wFixA].defined then
  error_std(ID,'Undefined fixpoint '+Fixname(wFixA));
if not wFixVal[wFixB].defined then
  error_std(ID,'Undefined fixpoint '+Fixname(wFixB));
if not (FixDir(wFixA)=FixDir(wFixB)) then
  error_std(ID,'Fixpoints for same axis: '+Fixname(wFixA)+' '+Fixname(wFixB));
hA:=wFixVal[wFixA].h; hB:=wFixVal[wFixB].h;
in_interval:=in_interval_ijk(h,hA,hB,wreg)
end; (* in_interval*)

function in_region_ijk(i,iA,iB:itype; ireg:regsettype;
                       j,jA,jB:jtype; jreg:regsettype;
                       k,kA,kB:ktype; kreg:regsettype):boolean;
const ID = 'in_region_ijk';
begin
in_region_ijk:=in_interval_ijk(i,iA,iB,ireg) and
               in_interval_ijk(j,jA,jB,jreg) and
               in_interval_ijk(k,kA,kB,kreg);
end; (* in_region_ijk*)

function in_region(i:itype;xFixA,xFixB:FixType;xreg:regsettype;
                       j:jtype;yFixA,yFixB:FixType;yreg:regsettype;
                       k:ktype;zFixA,zFixB:FixType;zreg:regsettype):boolean;
const ID='in_region';
var iA,iB:itype;
    jA,jB:jtype;
    kA,kB:ktype;
begin
check_fix_region_param(ID,xFixA,xFixB,yFixA,yFixB,zFixA,zFixB);
iA:=wFixVal[xFixA].h; iB:=wFixVal[xFixB].h;
jA:=wFixVal[yFixA].h; jB:=wFixVal[yFixB].h;
kA:=wFixVal[zFixA].h; kB:=wFixVal[zFixB].h;
in_region:=in_region_ijk(i,iA,iB,xreg,j,jA,jB,yreg,k,kA,kB,zreg)
end; (* in_region *)

function in_plane_ijk(reg:regsettype;
                      i,iA,iB:itype;
                      j,jA,jB:jtype;
                      k,kA,kB:ktype):boolean;
const ID='in_plane_ijk';
var h,hAB,m,mA,mB,n,nA,nB:htype;
    res_inside,res_outside,res_eq,res_eqm,res_eqn:boolean;
begin
h  :=0; (* Arbitary initialization just to make compiler happy ! *)
hAB:=0;
m  :=0;
mA :=0;
mB :=0;
n  :=0;
nA :=0;
nB :=0;
if not ((iA=iB) or (jA=jB) or (kA=kB)) then
  error_std(ID,'One axis must collapse (iA=iB) or (jA=jB) or (kA=kB)!');
if ((iA=iB) and (jA=jB)) or ((iA=iB) and (kA=kB)) or ((jA=jB) and (kA=kB)) then
  error_std(ID,'Only one axis can collapse '+
               '(iA=iB) and ( (jA<>jB) and (kA<>kB) ) etc');
if ((eqA in reg) or (eqB in reg)) then
  error_std(ID,'Cannot use eqA or eqB in this connection.');
if iA=iB then
  begin
    h:=i; hAB:=iA;
    m:=j; mA:=jA; mB:=jB;
    n:=k; nA:=kA; nB:=kB
  end;
if jA=jB then
  begin
    h:=j; hAB:=jA;
    m:=i; mA:=iA; mB:=iB;
    n:=k; nA:=kA; nB:=kB
  end;
if kA=kB then
  begin
    h:=k; hAB:=kA;
    m:=i; mA:=iA; mB:=iB;
    n:=j; nA:=jA; nB:=jB
  end;
res_inside:=(inside in reg) and
            in_interval_ijk(m,mA,mB,[inside]) and
            in_interval_ijk(n,nA,nB,[inside]);
res_outside:=(outside in reg) and
            (in_interval_ijk(m,mA,mB,[outside]) or
             in_interval_ijk(n,nA,nB,[outside]));
res_eqm    := ((m=mA) or (m=mB)) and in_interval_ijk(n,nA,nB,[inside,eqAB]);
res_eqn    := ((n=nA) or (n=nB)) and in_interval_ijk(m,mA,mB,[inside,eqAB]);
res_eq     :=(eqAB in reg) and (res_eqm or res_eqn);
in_plane_ijk:= (h=hAB) and (res_inside or res_outside or res_eq);
end; (* in_plane_ijk *)

function in_plane(reg:regsettype;
                  i:itype;xFixA,xFixB:FixType;
                  j:jtype;yFixA,yFixB:FixType;
                  k:ktype;zFixA,zFixB:FixType):boolean;
const ID='in_plane';
var iA,iB:itype;
    jA,jB:jtype;
    kA,kB:ktype;
begin
check_fix_region_param(ID,xFixA,xFixB,yFixA,yFixB,zFixA,zFixB);
iA:=wFixVal[xFixA].h; iB:=wFixVal[xFixB].h;
jA:=wFixVal[yFixA].h; jB:=wFixVal[yFixB].h;
kA:=wFixVal[zFixA].h; kB:=wFixVal[zFixB].h;
in_plane:=in_plane_ijk(reg,i,iA,iB,j,jA,jB,k,kA,kB)
end; (* in_plane *)

function in_cube_ijk(reg:regsettype;
                      i,iA,iB:itype;
                      j,jA,jB:jtype;
                      k,kA,kB:ktype):boolean;
const ID='in_cube_ijk';
var res_inside,res_outside,res_eq,res_eqx,res_eqy,res_eqz:boolean;
begin
if not ((iA<>iB) and (jA<>JB) or (kA<>kB)) then
  error_std(ID,'It is illegal to have (iA=iB) or (jA=jB) or (kA=kB)!');
if ((eqA in reg) or (eqB in reg)) then
  error_std(ID,'Cannot use eqA or eqB in this connection.');
res_inside:=(inside in reg) and
            in_interval_ijk(i,iA,iB,[inside]) and
            in_interval_ijk(j,jA,jB,[inside]) and
            in_interval_ijk(k,kA,kB,[inside]);
res_outside:=(outside in reg) and
            (in_interval_ijk(i,iA,iB,[outside]) or
            in_interval_ijk(j,jA,jB,[outside]) or
            in_interval_ijk(k,kA,kB,[outside]));
res_eqx    :=((i=iA) and in_plane_ijk([inside,eqAB],i,iA,iA,j,jA,jB,k,kA,kB)) or
             ((i=iB) and in_plane_ijk([inside,eqAB],i,iB,iB,j,jA,jB,k,kA,kB));
res_eqy    :=((j=jA) and in_plane_ijk([inside,eqAB],i,iA,iB,j,jA,jA,k,kA,kB)) or
             ((j=jB) and in_plane_ijk([inside,eqAB],i,iA,iB,j,jB,jB,k,kA,kB));
res_eqz    :=((k=kA) and in_plane_ijk([inside,eqAB],i,iA,iB,j,jA,jB,k,kA,kA)) or
             ((k=kB) and in_plane_ijk([inside,eqAB],i,iA,iB,j,jA,jB,k,kB,kB));
res_eq     :=(eqAB in reg) and (res_eqx or res_eqy or res_eqz);
in_cube_ijk:= (res_inside or res_outside or res_eq);
end; (* in_plane_ijk *)

function in_cube(reg:regsettype;
                 i:itype;xFixA,xFixB:FixType;
                 j:jtype;yFixA,yFixB:FixType;
                 k:ktype;zFixA,zFixB:FixType):boolean;
const ID='in_cube';
var iA,iB:itype;
    jA,jB:jtype;
    kA,kB:ktype;
begin
check_fix_region_param(ID,xFixA,xFixB,yFixA,yFixB,zFixA,zFixB);
iA:=wFixVal[xFixA].h; iB:=wFixVal[xFixB].h;
jA:=wFixVal[yFixA].h; jB:=wFixVal[yFixB].h;
kA:=wFixVal[zFixA].h; kB:=wFixVal[zFixB].h;
in_cube:=in_cube_ijk(reg,i,iA,iB,j,jA,jB,k,kA,kB)
end; (* in_cube *)


procedure estimate_fieldvalues_at_lost_NOPs;
const ID='estimate_fieldvalues_at_lost_NOPs';
var i:itype;j:jtype; k:ktype;
    dum,dV:datatype;
    cW,dW,
    cE,dE,
    cS,dS,
    cN,dN,
    cB,dB,
    cT,dT,
    csum,dsum:datatype;
    nsum:integer;
    vW,vE,vS,vN,vB,vT:boolean;
begin
if wr_all_procedure_id then writeln(ID,'...');
if wr_all_procedure_id then writeln(LOG,ID,'...');
dW:=0; (* Arbitary initialization just to keep compiler happy! *)
dE:=0; (* Arbitary initialization just to keep compiler happy! *)
dB:=0; (* Arbitary initialization just to keep compiler happy! *)
dT:=0; (* Arbitary initialization just to keep compiler happy! *)
dS:=0; (* Arbitary initialization just to keep compiler happy! *)
dN:=0; (* Arbitary initialization just to keep compiler happy! *)

for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      begin
        if ((GP[i]^[j]^[k].nodetyp=NOP)) then
         GP[i]^[j]^[k].valid_fieldvalue:=false; (* claus !!!!!!!!!!!!!!! *)
      end;

for i:=1 to imax do
  for j:=1 to jmax do
    for k:=1 to kmax do
      begin
        set_cvsize(i,j,k,dum,dum,dum,dum,dum,dum,dV);

      if (GP[i]^[j]^[k].nodetyp=NOP) and
         (dV=0) and
         (GP[i]^[j]^[k].valid_fieldvalue=false) then
        begin (* It may be possible to estimate a fieldvalue for this NOP *)
          vW:=false;
          vE:=false;
          vN:=false;
          vS:=false;
          vB:=false;
          vT:=false;
          if (i>1) and (i<=imax) then
            begin
              if (GP[i-1]^[j]^[k].valid_fieldvalue) then
                begin
                  cW:=GP[i-1]^[j]^[k].c;
                  dW:=xnod(i)-xnod(i-1);
                  vW:=true;
                end;
            end;
         if (i>=1) and (i<imax) then
            begin
              if (GP[i+1]^[j]^[k].valid_fieldvalue) then
                begin
                  cE:=GP[i+1]^[j]^[k].c;
                  dE:=xnod(i+1)-xnod(i);
                  vE:=true;
                end;
            end;

         if (geometry=cartesian3D) then
           begin
             if (j>1) and (j<=jmax) then
               begin
                if (GP[i]^[j-1]^[k].valid_fieldvalue) then
                  begin
                    cS:=GP[i]^[j-1]^[k].c;
                    dS:=ynod(j)-ynod(j-1);
                    vS:=true;
                  end;
              end;
           if (j>=1) and (j<jmax) then
              begin
                if (GP[i]^[j+1]^[k].valid_fieldvalue) then
                  begin
                    cN:=GP[i]^[j+1]^[k].c;
                    dN:=ynod(j+1)-ynod(j);
                    vN:=true;
                  end;
              end;
          end; (* not 2D geometry *)

         if (k>1) and (k<=kmax) then
            begin
              if (GP[i]^[j]^[k-1].valid_fieldvalue) then
                begin
                  cB:=GP[i]^[j]^[k-1].c;
                  dB:=znod(k)-znod(k-1);
                  vB:=true;
                end;
            end;
         if (k>=1) and (k<kmax) then
            begin
              if (GP[i]^[j]^[k+1].valid_fieldvalue) then
                begin
                  cT:=GP[i]^[j]^[k+1].c;
                  dT:=znod(k+1)-znod(k);
                  vT:=true;
                end;
            end;

          if ((vW) and (dW=0)) or
             ((vE) and (dE=0)) or
             ((vN) and (dN=0)) or
             ((vS) and (dS=0)) or
             ((vB) and (dB=0)) or
             ((vT) and (dT=0)) then
               error_node(ID,i,j,k,'Cannot repair this fieldvalue');
         csum:=0;
         dsum:=0;
         nsum:=0;

         if vW then begin csum:=csum+cW/dW; dsum:=dsum+1/dW; inc(nsum) end;
         if vE then begin csum:=csum+cE/dE; dsum:=dsum+1/dE; inc(nsum) end;
         if vN then begin csum:=csum+cN/dN; dsum:=dsum+1/dN; inc(nsum) end;
         if vS then begin csum:=csum+cS/dS; dsum:=dsum+1/dS; inc(nsum) end;
         if vB then begin csum:=csum+cB/dB; dsum:=dsum+1/dB; inc(nsum) end;
         if vT then begin csum:=csum+cT/dT; dsum:=dsum+1/dT; inc(nsum) end;
         if nsum>=2 then
           begin
             if (dsum=0) then error_std(ID,'dsum=0!?');
             GP[i]^[j]^[k].c:=csum/dsum;
             GP[i]^[j]^[k].valid_fieldvalue:=true;
           end;
         if (nsum=1) and wr_details then
           begin
             write(LOG,ID,' ');
             write(LOG,' i = ',i:5);
             write(LOG,' j = ',j:5);
             write(LOG,' k = ',k:5);
             if vW then write(LOG,' vW ');
             if vE then write(LOG,' vE ');
             if vS then write(LOG,' vS ');
             if vN then write(LOG,' vN ');
             if vB then write(LOG,' vB ');
             if vT then write(LOG,' vT ');
             writeln(LOG);
           end;
       end;
    end; (* grand loop *)
end;


(* Default *)

procedure grid_default;
begin
(* x-axis *)
set_FixVal(xFix1,0.0);
set_FixVal(xFix2,1.0);
set_axis_single(xFix1,xFix2,3,FocusA,1.0);

(* y-axis *)
set_FixVal(yFix1,0.0);
set_FixVal(yFix2,1.0);
set_axis_single(yFix1,yFix2,3,FocusA,1.0);

(* z-axis *)
set_FixVal(zFix1,0.0);
set_FixVal(zFix2,1.0);
set_axis_single(zFix1,zFix2,3,FocusA,1.0);
end; (* set_grid_definitions *)

procedure boundary_conditions_default(i:itype;j:jtype;k:ktype);
begin
cBC[fixed1]:=1.0;
cBC[fixed2]:=2.0;
if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix1,zFix1) then set_node(i,j,k,fixed1);

if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix2,zFix2) then set_node(i,j,k,fixed2);
end; (* set_boundary_condition *)


procedure fluxes_default(i:itype;j:jtype;k:ktype);
begin
if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix1,zFix1) then update_flxval(Flx1,top,i,j,k,plus);
if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix2,zFix2) then update_flxval(Flx2,bottom,i,j,k,plus);
end; (* fluxes *)

procedure probes_default;
var c1,dc1:datatype;
    valid1:boolean;
begin
get_fieldvalue(0.5,0.001,0.5,0.001,0.5,0.001,c1,dc1,valid1);
obsval[obs1]:=c1;
end; (* probes *)

function materials_default(i:itype;j:jtype;k:ktype):mattype;
begin
materials_default:=mat1;
end; (* materials *)

function e_default(i:itype;j:jtype;k:ktype):datatype;
begin
e_default:=1;
end; (* e *)

function beta_default(i:itype;j:jtype;k:ktype):datatype;
begin
beta_default:=e_default(i,j,k);
end; (* beta *)

function D_default(dir:dirtype;i:itype;j:jtype;k:ktype):datatype;
begin
D_default:=1e-6;
end; (* D *)

function G_default(i:itype;j:jtype;k:ktype):datatype;
begin
G_default:=0;
end; (* G *)

function Lambda_default(i:itype;j:jtype;k:ktype):datatype;
begin
Lambda_default:=0.0;
end; (* Lambda *)

procedure default_problem;
begin
runid                          := '0000';
runtitle                       := 'Default problem';
solution                       := steady;
geometry                       := cartesian3d;
Ly                             := 1.0;
grid_def                       := grid_default;
force_new_grid_in_every_run    := false;
boundary_conditions_def        := boundary_conditions_default;
flux_def                       := fluxes_default;
probe_def                      := probes_default;
materials_def                  := materials_default;
e_def                          := e_default;
beta_def                       := beta_default;
G_def                          := G_default;
lambda_def                     := lambda_default;
D_def                          := D_default;
initialfield_def               := nil;
import_initialfield            := false;
import_finalfield_guess        := false;
export_field                   := false;
use_fieldbuffer                := no_cBUF;
flowfield                      := none;
flowfactor                     := 1.0;
import_field_name              := '';
export_field_name              := '';
flowfield_name                 := '';
plotfiles_def                  := nil;
user_procedure_each_iter_def   := nil;
wr_details                     := false;
wr_main_procedure_id           := false;
wr_all_procedure_id            := false;
wr_iteration_line_log          := false;
wr_iteration_line_screen       := true;
wr_residual_during_calc_log    := false;
wr_residual_during_calc_screen := false;
wr_flux_during_calc_log        := false;
wr_flux_during_calc_screen     := false;
wr_probes_during_calc_log      := false;
wr_probes_during_calc_screen   := false;
wr_final_results_log           := true;
wr_final_results_screen        := true;
wr_axes                        := true;
wr_nodes                       := false;
wr_node_numbers                := true;
wr_node_sizes                  := false;
wr_coefficients                := false;
wr_materials_volumes           := true;
warning_priority_log           := war_other;
warning_priority_screen        := war_other;
solver_def                     := Find_better_field_thomas;
scheme                         := exact;
relax_factor                   := 1.0;
flux_convset                   := [flx1];
probe_convset                  := [];
conv_evaluation_period         := 100;
min_iterations                 := 5;
max_iterations                 := 1000;
max_time                       := 3*60;
max_change                     := 1e-6;
max_residual_sum               := 1e-4;
dtim                           := 0;
BC_running                     := false;
BC_running_update_of_cBCs_def  := nil;
BC_running_min_iterations      := 100;
BC_running_max_residual_sum_before_new_BC := 1e-9;
BC_running_convergence_def     := nil;
wr_BC_running_messages_log     := false;
wr_BC_running_messages_screen  := false;
press_enter_wanted             := true;
end; (* default problem *)

procedure set_control_variables_to_defaults;
begin
runid                          := '0000';
runtitle                       := 'Default control variables';
solution                       := steady;
geometry                       := cartesian3d;
Ly                             := 1.0;
grid_def                       := nil;
force_new_grid_in_every_run    := false;
boundary_conditions_def        := nil;
flux_def                       := nil;
probe_def                      := nil;
materials_def                  := nil;
e_def                          := nil;
beta_def                       := nil;
G_def                          := nil;
lambda_def                     := nil;
D_def                          := nil;
initialfield_def               := nil;
import_initialfield            := false;
import_finalfield_guess        := false;
export_field                   := false;
use_fieldbuffer                := no_cBUF;
flowfield                      := none;
flowfactor                     := 1.0;
import_field_name              := '';
export_field_name              := '';
flowfield_name                 := '';
plotfiles_def                  := nil;
user_procedure_each_iter_def   := nil;
wr_details                     := false;
wr_main_procedure_id           := false;
wr_all_procedure_id            := false;
wr_iteration_line_log          := false;
wr_iteration_line_screen       := true;
wr_residual_during_calc_log    := false;
wr_residual_during_calc_screen := false;
wr_flux_during_calc_log        := false;
wr_flux_during_calc_screen     := false;
wr_probes_during_calc_log      := false;
wr_probes_during_calc_screen   := false;
wr_final_results_log           := true;
wr_final_results_screen        := true;
wr_axes                        := true;
wr_nodes                       := false;
wr_node_numbers                := true;
wr_node_sizes                  := false;
wr_coefficients                := false;
wr_materials_volumes           := true;
warning_priority_log           := war_other;
warning_priority_screen        := war_other;
solver_def                     := Find_better_field_thomas;
scheme                         := exact;
relax_factor                   := 1.0;
flux_convset                   := [];
probe_convset                  := [];
conv_evaluation_period         := 50;
min_iterations                 := 5;
max_iterations                 := 500;
max_time                       := 3*60;
max_change                     := 1e-6;
max_residual_sum               := 1e-4;
dtim                           := 0;
BC_running                     := false;
BC_running_update_of_cBCs_def  := nil;
BC_running_min_iterations      := 100;
BC_running_max_residual_sum_before_new_BC := 1e-9;
BC_running_convergence_def     := nil;
wr_BC_running_messages_log     := false;
wr_BC_running_messages_screen  := false;
press_enter_wanted             := true;
end; (* set_control_variables_to_defaults *)

procedure wr_line(var OM:text);
var i:integer;
begin
for i:=1 to 75 do write(OM,'-');
writeln(OM);
end;

procedure wr_header(var OM:text);
begin
wr_line(OM);
writeln(OM,'Description   : ',description);
writeln(OM,'Program name  : ',programid+' ('+copyright+')');
writeln(OM,'Version       : ',versionid);
writeln(OM,'Documentation : ',documentation);
(*writeln(OM,'Programmer    : ',programmer1);*)
(*writeln(OM,'                ',programmer2);*)
wr_line(OM);
end;


function file_exists(fname:string):boolean;
var f:text;
begin
{$I-}
assign(F,fname);
reset(F);
close(F);
{I+}
file_exists:=(IOResult=0) and (fname<>'');
end;

procedure wr_warning_table(var OM:text);
var war:warningtype;
begin
for war:=succ(war_none) to pred(war_last) do
  if (warning_table[war]<>0) then
    writeln(OM,'Warning: ',war_string(war):30,' was issued ',warning_table[war]:4,' times');
end;


procedure set_FixVal(wFix:Fixtype;wval:datatype);
const ID='set_FixVal';
begin
if wr_details or wr_main_procedure_id then writeln(ID,' (',FixName(wFix),'= ',wval,')');
if wr_details or wr_main_procedure_id then writeln(LOG,ID,' (',FixName(wFix),'= ',wval,')');
check_wFix(id,wFix);
if wFixval[wFix].defined then
  error_std(ID,'This wFix was already defined: '+FixName(wFix));
wFixval[wFix].w:=wval;
wFixval[wFix].defined:=true;
end;


procedure set_axis_subdiv(var w,dw:axistype;
                              h1,h2:htype;
                              w1,w2:datatype;
                              f:funcidtype;
                              pow:datatype;
                              wdir:dirtype);
const ID='set_axis_subdiv';
var h,h12:htype;
    h1st,h2st,hmaxst,st0,st:string;
    wnorm,hnorm,hmax:datatype;
begin
if wr_all_procedure_id then writeln(ID,' ',DirName(wdir));
if wr_all_procedure_id then writeln(LOG,ID,' ',DirName(wdir));
wnorm:=0; (* Arbitary initialization just to keep compiler happy! *)
hmax:=0;  (* Arbitary initialization just to keep compiler happy! *)
w[h1]:=w1;
w[h2]:=w2;
h12:=h2-h1;
case wdir of
  xdir: begin hmax:=imaxTot; st0:='imaxTot'; end;
  ydir: begin hmax:=jmaxTot; st0:='jmaxTot'; end;
  zdir: begin hmax:=kmaxTot; st0:='kmaxTot'; end;
else
  error_std(ID,'Illegal direction : '+dirname(wdir));
end;
str(h1,h1st);
str(h2,h2st);
str(hmax:12:1,hmaxst);
st:=' h1='+h1st+' h2='+h2st+' '+st0+'(hmax)='+hmaxst;
if h1>hmax then error_std(ID,'h1>'+st0+' (hmax) '+dirname(wdir)+st+cr+'Increase '+st0);
if h2>hmax then error_std(ID,'h2>'+st0+' (hmax) '+dirname(wdir)+st+cr+'Increase '+st0);
if h1>h2 then error_std(ID,'h1>h2 '+dirname(wdir)+st);
if (h12<=0) then error_std(ID,'h12<=0 '+dirname(wdir)+st);
if (w1>=w2) then error_std(ID,'w1>=w2 '+dirname(wdir)+st);
for h:=h1 to h2-1 do
  begin
    hnorm:=(h-h1)/(h2-h1);
    case f of
      FocusA: wnorm:=pw(hnorm,pow);
      FocusB: wnorm:=1-pw(1-hnorm,pow)
    else
      error_std(ID,'Unknown focus function. '+dirname(wdir)+st);
    end;
    w[h]:=w1+wnorm*(w2-w1);
  end;
for h:=h1 to h2-1 do
  begin
    dw[h]:=w[h+1]-w[h];
    if dw[h]<0 then error_std(ID,'dw[h]<=0');
  end;
case wdir of
  xdir: begin
          imax:=maxi(imax,h2);
          if imax>imaxtot then error_std(ID,'Too many x-nodes.'+st);
          x:=w;
          dx:=dw
        end;
  ydir: begin
          jmax:=maxi(jmax,h2);
          if jmax>jmaxtot then error_std(ID,'Too many y-nodes.'+st);
          y:=w;
          dy:=dw
          end;
  zdir: begin
          kmax:=maxi(kmax,h2);
          if kmax>kmaxtot then error_std(ID,'Too many z-nodes.'+st);
          z:=w;
          dz:=dw
        end;
else
  error_std(ID,'Illegal direction '+dirname(wdir)+st);
end; (* case *)

end; (* set_axis_subdiv *)

procedure check_Fix_pair(IDcall:string;wFixA,wFixB:FixType);
const ID='check_fix_pair';
var err:integer;
begin
err:=0;
if (FixDir(wFixA)<>FixDir(wFixB)) then err:=1;
if not wFixVal[wFixA].defined     then err:=err+10;
if not wFixVal[wFixb].defined     then err:=err+20;
if (ord(wFixA)>=ord(wFixB))       then err:=err+100;
if (ord(wFixB)-ord(wFixA))>1      then err:=err+1000;
if err<>0 then
  begin
    writeln('* Message from ',ID,' called by ',IDcall);
    writeln('* Error code = ',err);
    writeln('* Code    1 : wFixA and wFixB are not in the same direction');
    writeln('* Code   10 : wFixA is undefined');
    writeln('* Code   20 : wFixB is undefined');
    writeln('* Code  100 : ord(wFixA)>=ord(wFixB) !');
    writeln('* Code 1000 : (ord(wFixA)-ord(wFixB))>1 !');
    writeln('* wFixA = ',FixName(wFixA));
    writeln('* wFixB = ',FixName(wFixB));
    error_std(IDcall,'Some problem with wFixA and/or wFixB (see code).');
  end;
if wFixVal[wFixA].h = 0 then
  begin
    writeln('* Fix-point = ',FixName(wFixA));
    error_std(ID,'h is undefined for this Fix-point.'+cr+
    'Have you treated all the previous fix points on this axis?');
  end;

end;

procedure set_axis_single(wFixA,wFixB:FixType;
                          hAB:htype;
                          f:funcidtype;
                          pow:datatype);
const id='set_axis_single';
var w,dw:axistype;
    hA,hB,hmax:htype;
    wdir:dirtype;
begin
hmax:=0; (* Arbitary initialization just to keep compiler happy !*)
wdir:=FixDir(wFixA);
if wr_details then writeln(ID,': ',FixName(wFixA),' - ',FixName(wFixB));
if wr_details then writeln(LOG,ID,': ',FixName(wFixA),' - ',FixName(wFixB));
if hAB<=0 then
  begin
    writeln('* hAB = ',hAB:5);
    writeln('* wFixA = ',FixName(wFixA));
    writeln('* wFixB = ',FixName(wFixB));
    error_std(ID,'hAB<=0');
  end;
case FixDir(wFixA) of
  xdir: begin hmax:=imaxTot; w:=x; dw:=dx end;
  ydir: begin hmax:=jmaxTot; w:=y; dw:=dy end;
  zdir: begin hmax:=kmaxTot; w:=z; dw:=dz end;
else
  error_std(ID,'Unknown direction of '+FixName(wFixA));
end; (* case *)
if (wFixA=xFix1) or
   (wFixA=yFix1) or
   (wFixA=zFix1) then wFixVal[wFixA].h:=1; (* Initial condition *)
check_fix_pair(ID,wFixA,wFixB);
(* Set w- and dw-values for hA and hb (plus w[hB+1]) *)
hA:=wFixVal[wFixA].h;
hB:=hA+hAB+1;
if (hB>hmax) then
begin
  wr_memory_status(LOG,'');
  wr_memory_status(output,'');
  writeln('Number of nodes requested = ',hAB);
  writeln(LOG,'Number of nodes requested = ',hAB);
  error_std(ID,'Too many nodes between '+FixName(wFixA)+' and '+FixName(wFixB)); (* hB too large *)
end;
wFixVal[wFixB].h:=hB;
w[hA]:=wFixVal[wFixA].w;
w[hB]:=wFixVal[wFixB].w;
w[hB+1]:=wFixVal[wFixB].w;
dw[hA]:=0;
dw[hB]:=0;
(* set values for hA+1 to hB-1 *)
set_axis_subdiv(w,dw,hA+1,hB,wFixVal[wFixA].w,wFixVal[wFixB].w,f,pow,wdir);
end; (* set_axis_single *)

procedure set_axis_double(wFixA,wFixB:FixType;
                          hAM,hMB:htype;
                          fA,fB:funcidtype;
                          powA,powB,wdiv:datatype);
const id='set_axis_double';
var w,dw:axistype;
    hA,hM,hB,hmax:htype;
    wdir:dirtype;
begin
hmax:=0; (* Arbitary initialization just to keep compiler happy !*)
wdir:=FixDir(wFixA);
if wr_details then writeln(ID,': ',FixName(wFixA),' - ',FixName(wFixB));
if wr_details then writeln(LOG,ID,': ',FixName(wFixA),' - ',FixName(wFixB));
if (hAM<=0) or (hMB<=0) then
  begin
    writeln('* hAM = ',hAM:5,' hMB = ',hMB:5);
    writeln('* wFixA = ',FixName(wFixA));
    writeln('* wFixB = ',FixName(wFixB));
    error_std(ID,'hAM or hMB too small (<=0)');
  end;
if (wdiv<=0) or (wdiv>=1) then
  begin
    writeln('* wdiv = ',wdiv:12:6);
    error_std(ID,'wdiv should be in the range 0+ to 1-.');
  end;
case FixDir(wFixA) of
  xdir: begin hmax:=imaxTot; w:=x; dw:=dx end;
  ydir: begin hmax:=jmaxTot; w:=y; dw:=dy end;
  zdir: begin hmax:=kmaxTot; w:=z; dw:=dz end;
else
  error_std(ID,'Unknown direction of '+FixName(wFixA));
end; (* case *)
if (wFixA=xFix1) or
   (wFixA=yFix1) or
   (wFixA=zFix1) then wFixVal[wFixA].h:=1; (* Initial condition *)
check_fix_pair(ID,wFixA,wFixB);
(* Set w- and dw-values for hA and hb (plus w[hB+1]) *)
hA:=wFixVal[wFixA].h;
hM:=hA+hAM+1;
hB:=hM+hMB+1;
if (hM>hmax) or (hB>hmax) then
begin
  wr_memory_status(LOG,'');
  wr_memory_status(output,'');
  writeln('Number of nodes requested = ',hAM,' + ',hMB,' = ',hAM+hMB);
  writeln(LOG,'Number of nodes requested = ',hAM,' + ',hMB,' = ',hAM+hMB);
  error_std(ID,'Too many nodes between '+FixName(wFixA)+' and '+FixName(wFixB)); (* hB too large *)
end;
wFixVal[wFixB].h:=hB;
w[hA]:=wFixVal[wFixA].w;
w[hB]:=wFixVal[wFixB].w;
w[hM]:=w[hA]+wdiv*(w[hB]-w[hA]);
w[hB+1]:=wFixVal[wFixB].w;
dw[hA]:=0;
dw[hB]:=0;
(* set values for hA+1 to hB-1 *)
set_axis_subdiv(w,dw,hA+1,hM,w[hA],w[hM],fA,powA,wdir);
set_axis_subdiv(w,dw,hM,hB,w[hM],w[hB],fB,powB,wdir);
end; (* set_axis_double *)

procedure set_axis_triple(wFixA,wFixB:FixType;
                          hAM,hMM,hMB:htype;
                          fA,fM,fB:funcidtype;
                          powA,powM,powB,wdivA,wdivB:datatype);
const id='set_axis_triple';
var w,dw:axistype;
    hA,hMaa,hMbb,hB,hmax:htype;
    wdir:dirtype;
begin
hmax:=0; (* Arbitary initialization just to keep compiler happy !*)
wdir:=FixDir(wFixA);
if wr_details then writeln(ID,': ',FixName(wFixA),' - ',FixName(wFixB));
if wr_details then writeln(LOG,ID,': ',FixName(wFixA),' - ',FixName(wFixB));
if (hAM<=0) or (hMM<=0) or (hMB<=0) then
  begin
    writeln('* hAM = ',hAM:5,' hMM = ',hMM:5,' hMB = ',hMB:5);
    writeln('* wFixA = ',FixName(wFixA));
    writeln('* wFixB = ',FixName(wFixB));
    error_std(ID,'hAM, hMM or hMB too small (<=0)');
  end;
if (wdivA<=0) or (wdivA>=1) or
   (wdivB<=0) or (wdivB>=1) or
   (wdivA>=wdivB) then
  begin
    writeln('* wdivA = ',wdivA:12:6);
    writeln('* wdivB = ',wdivB:12:6);
    error_std(ID,'wdivA and wdivB should each be in the range 0+ to 1-.'+cr+
                 'Furthermore it is required that wdivA<wdivB.');
  end;

case FixDir(wFixA) of
  xdir: begin hmax:=imaxTot; w:=x; dw:=dx end;
  ydir: begin hmax:=jmaxTot; w:=y; dw:=dy end;
  zdir: begin hmax:=kmaxTot; w:=z; dw:=dz end;
else
  error_std(ID,'Unknown direction of '+FixName(wFixA));
end; (* case *)
if (wFixA=xFix1) or
   (wFixA=yFix1) or
   (wFixA=zFix1) then wFixVal[wFixA].h:=1; (* Initial condition *)
check_fix_pair(ID,wFixA,wFixB);
(* Set w- and dw-values for hA and hb (plus w[hB+1]) *)
hA  :=wFixVal[wFixA].h;
hMaa:=hA+hAM+1;
hMbb:=hMaa+hMM+1;
hB:=hMbb+hMB+1;
if hB>hmax then error_std(ID,'hB too large !');
wFixVal[wFixB].h:=hB;
w[hA]:=wFixVal[wFixA].w;
w[hB]:=wFixVal[wFixB].w;
w[hMaa]:=w[hA]+wdivA*(w[hB]-w[hA]);
w[hMbb]:=w[hA]+wdivB*(w[hB]-w[hA]);
if (w[hMaa]>w[hMbb]) then
  error_std(ID,'Probably wrong axis breakdown');
w[hB+1]:=wFixVal[wFixB].w;
dw[hA]:=0;
dw[hB]:=0;
(* set values for hA+1 to hB-1 *)
set_axis_subdiv(w,dw,hA+1,hMaa,w[hA],w[hMaa],fA,powA,wdir);
set_axis_subdiv(w,dw,hMaa,hMbb,w[hMaa],w[hMbb],fM,powM,wdir);
set_axis_subdiv(w,dw,hMbb,hB,w[hMbb],w[hB],fB,powB,wdir);
end; (* set_axis_triple*)

procedure check_axes;
const ID='check_axes';
var err,err1,err2,err3,err10,err20,err30:integer; h:htype;
begin
if wr_main_procedure_id then writeln(ID,'...');
if wr_main_procedure_id then writeln(LOG,ID,'...');
err:=0; err1:=0; err2:=0; err3:=0; err10:=0; err20:=0; err30:=0;
if not (wFixVal[xFix1].defined) then err:=err+1;
if not (wFixVal[xFix2].defined) then err:=err+2;
if not (wFixVal[yFix1].defined) then err:=err+10;
if not (wFixVal[yFix2].defined) then err:=err+20;
if not (wFixVal[zFix1].defined) then err:=err+100;
if not (wFixVal[zFix2].defined) then err:=err+200;
for h:=2 to imax do
 if (x[h-1]>x[h]) then err1:=1000;
for h:=2 to jmax do
 if (y[h-1]>y[h]) then err2:=2000;
for h:=2 to kmax do
 if (z[h-1]>z[h]) then err3:=3000;
for h:=1 to imax do
 if (dx[h]<0) then err10:=10000;
for h:=1 to jmax do
 if (dy[h]<0) then err20:=30000;
for h:=1 to kmax do
 if (dz[h]<0) then err30:=20000;
err:=err+err1+err2+err3+err10+err20+err30;
if (err<>0) then
  begin
    writeln('* Error code = ',err:7);
    writeln('*    1: xFix1 is undefined.');
    writeln('*    2: xFix2 is undefined.');
    writeln('*   10: yFix1 is undefined.');
    writeln('*   20: yFix2 is undefined.');
    writeln('*  100: zFix1 is undefined.');
    writeln('*  200: zFix2 is undefined.');
    writeln('* 1000,2000,3000: x,y,z-arrayx[i] is unsorted !');
    writeln('* 10000,20000,30000: d[i]x,dy[j],dz[k] <0 !');
    error_std(ID,'Some problems (see code).');
  end;
end;

function DirName(Dir:DirType):string;
const ID='DirName';
var st:string;
begin
case Dir of
  east    : st:='east  ';
  west    : st:='west  ';
  north   : st:='north ';
  south   : st:='south ';
  top     : st:='top   ';
  bottom  : st:='bottom';
  xdir    : st:='xdir  ';
  ydir    : st:='ydir  ';
  zdir    : st:='zdir  ';
  nodir   : st:='NODIR '
else
  error_std(ID,'Unknown dir.');
end;
DirName:=st;
end; (* DirName *)


function maxi(i1,i2:integer):integer;
begin
if i1>i2 then maxi:=i1 else maxi:=i2
end;

function max(x1,x2:datatype):datatype;
begin
if x1>x2 then max:=x1 else max:=x2
end;

function min(x1,x2:datatype):datatype;
begin
if x1<x2 then min:=x1 else min:=x2
end;

function pw(x,y:datatype):datatype;
const ID='pw';
      zero=1e-139;
var pw2:datatype;
begin
pw2:=0;
if x<0 then Error_std(ID,'pw: x<0');
if (abs(x)<zero) then pw2:=0;
if (abs(x)>zero) and (abs(y)<zero) then pw2:=1;
if x>zero then pw2:=exp(y*ln(x));
pw:=pw2;
end;


function cosech(x:datatype):datatype;
const ID='cosech';
begin
if x=0 then error_std(ID,'x=0');
if (x>10000) or (x<-10000) then
  cosech:=0.0
else
  cosech:=2/(exp(x)-exp(-x));
end;



begin
end.