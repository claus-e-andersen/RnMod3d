program f0102prg;
(* --------------------- RnMod3d jobfile ---------------------- *)
(* Project: User guide example:                                 *)
(*          Steady radon diffusion + advection, 1D.             *)
(* Created: May 24, 1999                                        *)
(* Revised: July 17, 2000                                       *)

{$I R3dirs03}
uses R3Defi03,R3Main03,R3Writ03;

const lambda_use = 2.09838e-6;
      mu         = 17.5e-6;

var ksoil,cS,dP,velocity,Lz,Dsoil,esoil,Gsoil:datatype;

procedure grid;
begin
set_FixVal(xFix1,0.0);
set_FixVal(xFix2,1.0);
set_axis_single(xFix1,xFix2,1,FocusA,1.0);

set_FixVal(yFix1,0.0);
set_FixVal(yFix2,1.0);
set_axis_single(yFix1,yFix2,1,FocusA,1.0);

set_FixVal(zFix1, 0.0);
set_FixVal(zFix2,  Lz);
set_axis_double(zFix1,zFix2,30,30,FocusA,FocusB,2,2,0.5);
end;

procedure boundary_conditions_Soilgas(i:itype;j:jtype;k:ktype);
begin
cBC[fixed1]:=dP;
if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix1,zFix1) then set_node(i,j,k,fixed1);
cBC[fixed2]:=0.0;
if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix2,zFix2) then set_node(i,j,k,fixed2);
end;

procedure boundary_conditions_Rn(i:itype;j:jtype;k:ktype);
begin
cBC[fixed1]:=cS;
if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix1,zFix1) then set_node(i,j,k,fixed1);
cBC[fixed2]:=0.0;
if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix2,zFix2) then set_node(i,j,k,fixed2);
end;

procedure fluxes(i:itype;j:jtype;k:ktype);
begin
if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix1,zFix1) then update_flxval(Flx1,top,i,j,k,plus);
if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix2,zFix2) then update_flxval(Flx2,bottom,i,j,k,plus);
end;

procedure probes;
var cc:datatype; valid:boolean;
begin
cc:=fieldvalue(0.5,0.5,Lz/2,valid);
if not valid then cc:=0.0;
obsval[obs1]:=cc;
end;

function materials(i:itype;j:jtype;k:ktype):mattype;
begin
materials:=mat1;
end;

function e_soilgas(i:itype;j:jtype;k:ktype):datatype;
begin
e_soilgas:=0;
end;

function beta_soilgas(i:itype;j:jtype;k:ktype):datatype;
begin
beta_soilgas:=0;
end;

function D_soilgas(dir:dirtype;i:itype;j:jtype;k:ktype):datatype;
begin
D_soilgas:=ksoil/mu;
end;

function G_soilgas(i:itype;j:jtype;k:ktype):datatype;
begin
G_soilgas:=0;
end;

function lambda_soilgas(i:itype;j:jtype;k:ktype):datatype;
begin
lambda_soilgas:=0;
end;

function e_Rn(i:itype;j:jtype;k:ktype):datatype;
begin
e_Rn:=esoil;
end;

function beta_Rn(i:itype;j:jtype;k:ktype):datatype;
begin
beta_Rn:=e_Rn(i,j,k);
end;

function D_Rn(dir:dirtype;i:itype;j:jtype;k:ktype):datatype;
begin
D_Rn:=Dsoil;
end;

function G_Rn(i:itype;j:jtype;k:ktype):datatype;
begin
G_Rn:=Gsoil;
end;

function lambda_Rn(i:itype;j:jtype;k:ktype):datatype;
begin
lambda_Rn:=lambda_use;
end;

function sinh(x:datatype):datatype;
var z:datatype;
begin
z:=exp(x);
sinh:=(z-1/z)/2
end;

function cosh(x:datatype):datatype;
var z:datatype;
begin
z:=exp(x);
cosh:=(z+1/z)/2
end;

function c_exact(z:datatype):datatype;
var v,D,cinf,s,alpha,Ld:datatype;
(* See NBS technical note 1139, p. 26 *)
begin
D:=Dsoil;
v:=velocity;
alpha:=v/2/D;
Ld:=sqrt(D/esoil/lambda_use);
s:=1/sqrt(sqr(alpha) + sqr(1/Ld));
cinf:=Gsoil/lambda_use;
c_exact:=cinf*(1-1/sinh(Lz/s)*(exp(v*z/2/D)*sinh((Lz-z)/s) + exp(-v*(Lz-z)/2/D)*sinh(z/s)))+
         cS*exp(v*z/2/D)*sinh((Lz-z)/s) / sinh(Lz/s);
end;

function j_exact:datatype;
var v,D,cinf,s,alpha,Ld:datatype;
(* See NBS technical note 1139, p. 26 *)
begin
D:=Dsoil;
v:=velocity;
alpha:=v/2/D;
Ld:=sqrt(D/esoil/lambda_use);
s:=1/sqrt(sqr(alpha) + sqr(1/Ld));
cinf:=Gsoil/lambda_use;
j_exact:=cinf*(V/2 + D/s/sinh(Lz/s)*(cosh(Lz/s)-exp(v*Lz/2/D)))+
         cS*D/s*exp(v*Lz/2/D)/sinh(Lz/s);
end;


procedure wr_flux;
begin
writeln(LOG,' dP = ',dP:6:2);
writeln(LOG,'RnMod3d Rn flux at z=0: ',FlxVal[flx2].j:16,' Bq/m2/s');
writeln(LOG,'Exact Rn flux at z=0:   ',j_exact:16,' Bq/m2/s');
writeln(LOG,'Deviation:              ',100*(FlxVal[flx2].j-j_exact)/j_exact:16:4,' %');
end;

procedure wr_profile;
var Nsteps,zstart,zstop,dzz,zz,cc:datatype;
    valid:boolean;
begin
(* This procedure finds the field at (x,y,z) where x=0.5m           *)
(* and y=0.5m, and z is looped through the values from top to       *)
(* bottom.                                                          *)
Nsteps:=800;
if not (wFixVal[zFix1].defined and wFixVal[zFix2].defined) then
  error_std('wr_profile','Undefined fixpoints!');
zstart:=wFixVal[zFix1].w;
zstop :=wFixVal[zFix2].w;
dzz:=(zstop-zstart)/Nsteps;
zz :=zstart;
writeln(RES,'z':12,',','c':12,',','cexact':12);
while (zz<zstop) do
  begin
    cc:=fieldvalue(0.5,0.5,zz,valid);
    if valid then
      writeln(RES,zz:12:6,',',cc:12:6,',',c_exact(zz):12:6);
    zz:=zz+dzz;
  end;
end;

begin (* main *)
runid                          := '0102';
runtitle                       := 'User guide example: Steady Rn diff. and adv.';
solution                       := steady;
geometry                       := cartesian3d;
Ly                             := 1.0;
grid_def                       := grid;
flux_def                       := fluxes;
probe_def                      := probes;
materials_def                  := materials;
flux_convset                   := [flx1,flx2];
probe_convset                  := [obs1];
conv_evaluation_period         := 200;
min_iterations                 := 100;
max_iterations                 := 5000;
wr_axes                        := false;
wr_node_numbers                := false;
wr_materials_volumes           := false;

(* User-defined constants *)
Lz    := 5;      (* Column depth        *)
ksoil := 1e-11;  (* Soil permeability   *)
cS    := 5000;   (* Radon conc. at z=0  *)
dP    := -100;   (* Pressure difference *)
Dsoil := 1e-6;   (* Diffusivity         *)
esoil := 0.3;    (* Porosity            *)
Gsoil := 10000*lambda_use;  (* Generation rate *)

velocity:=ksoil/mu*dP/Lz;

(* First, the soil gas problem *)
boundary_conditions_def      := boundary_conditions_soilgas;
D_def                        := D_soilgas;
e_def                        := e_soilgas;
beta_def                     := beta_soilgas;
G_def                        := G_soilgas;
lambda_def                   := lambda_soilgas;
flowfield                    := export_to_qBUF;
relax_factor                 := 1.9;
max_change                   := 1e-12;
max_residual_sum             := 3e-16;
run_model;

(* Second, the radon problem *)
flowfield                    := import_from_qBUF;
boundary_conditions_def      := boundary_conditions_Rn;
D_def                        := D_Rn;
e_def                        := e_Rn;
beta_def                     := beta_Rn;
G_def                        := G_Rn;
lambda_def                   := lambda_Rn;
relax_factor                 := 1.0;
max_change                   := 1e-12;
max_residual_sum             := 3e-16;
run_model;

wr_flux;
wr_profile;
close_model;
end.