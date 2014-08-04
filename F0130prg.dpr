program F0130prg;
(* --------------------- RnMod3d jobfile ---------------------- *)
(* Project: User guide example:                                 *)
(*          House simulation (slab-on-grade)                    *)
(* Created: September 26, 1998                                  *)
(* Revised: July 18, 2000                                       *)

{$I R3dirs03}
uses R3Defi03,R3Main03,R3Writ03;

const
LambdaRn222    = 2.09838e-6; (* 1/s *)
mu             = 18.0e-6;    (* Pa s *)
rho_g          = 2.7e3;      (* kg/m3 *)
LOstwald       = 0.30;       (* water/gas partitioning *)
deltaP         = -1.0;       (* Pa *)

(* Horizontal (x) dimensions, m *)
Lx_soil        =  20.00;
Lx_slab        =  5.6419;
Lx_footer      =  0.300;
Lx_gap         =  0.003;

(* Vertical (z) dimensions, m *)
Lz_soil        = 10.00;
Lz_slab        =  0.10;
Lz_gravel      =  0.15;
Lz_footer      =  0.80;

(* Radium-226 concentration, Bq/kg *)
ARa_soil       = 40;
ARa_slab       = 50;
ARa_gravel     = 40;
ARa_footing    = 0;
ARa_gap        = 0;

(* Fraction of emanation, - *)
f_soil         = 0.2;
f_slab         = 0.1;
f_gravel       = 0.2;
f_footing      = 0;
f_gap          = 0;

(* Porosity, - *)
etot_soil      = 0.25;
etot_slab      = 0.20;
etot_gravel    = 0.40;
etot_footing   = 0.20;
etot_gap       = 1.00;

(* Volumetric water content, - *)
msat_soil      = 0.20;
msat_slab      = 0.0;
msat_gravel    = 0.0;
msat_footing   = 0.0;
msat_gap       = 0.0;

(* Bulk diffusivity, m2/s *)
D_soil         = 4.3e-7;
D_slab         = 2.0e-8;
D_gravel       = 1.8e-6;
D_footing      = 1.0e-10;
D_gap          = 1.2e-5;

(* Gas permeability, m2 *)
k_soil         = 1e-11;
k_slab         = 1e-15;
k_gravel       = 5e-9;
k_footing      = 1e-15;
k_gap          = 7.5e-7;

procedure grid;
begin
(* x-axis *)
set_FixVal(xFix1,0.000);
set_FixVal(xFix2,Lx_slab-Lx_gap);
set_FixVal(xFix3,Lx_slab);
set_FixVal(xFix4,Lx_slab+Lx_footer);
set_FixVal(xFix5,Lx_soil);
set_axis_double(xFix1,xFix2,15,15,FocusB,FocusB,2.1,3.0,0.97);
set_axis_single(xFix2,xFix3,5,FocusA,1.5);
set_axis_double(xFix3,xFix4,4,4,FocusA,FocusB,2.0,2.0,0.5);
set_axis_single(xFix4,xFix5,20,FocusA,2.5);

(* z-axis *)
set_FixVal(zFix1,-Lz_soil);
set_FixVal(zFix2,-Lz_footer);
set_FixVal(zFix3,-Lz_slab-Lz_gravel);
set_FixVal(zFix4,-Lz_slab);
set_FixVal(zFix5, 0.00);
set_axis_double(zFix1,zFix2,6,14,FocusA,FocusB,2.0,2.0,0.5);
set_axis_double(zFix2,zFix3,6,8,FocusA,FocusB,1.8,1.8,0.5);
set_axis_double(zFix3,zFix4,15,5,FocusB,FocusB,2.0,2,0.95);
set_axis_single(zFix4,zFix5,4,FocusB,3.0);
end;

procedure boundary_conditions_soilgas(i:itype;j:jtype;k:ktype);
begin
cBC[fixed1]:=deltaP;
cBC[fixed2]:=0;
if in_plane([inside,eqAB], (* Observe: Full slab, not just the gap *)
            i,xFix1,xFix3,
            j,yFix1,yFix2,
            k,zFix5,zFix5) then
  change_node(i,j,k,fixed1,ConX,ConX,ConX,ConX,ConX,ConX);
if in_plane([inside,eqAB], (* Atmospheric surface *)
            i,xFix4,xFix5,
            j,yFix1,yFix2,
            k,zFix5,zFix5) then
  change_node(i,j,k,fixed2,ConX,ConX,ConX,ConX,ConX,ConX);
end;

procedure boundary_conditions_radon(i:itype;j:jtype;k:ktype);
begin
cBC[fixed1]:=0;
cBC[fixed2]:=0;
if in_plane([inside,eqAB], (* Observe: Full slab, not just the gap *)
            i,xFix1,xFix3,
            j,yFix1,yFix2,
            k,zFix5,zFix5) then
  change_node(i,j,k,fixed1,ConX,ConX,ConX,ConX,ConX,ConX);
if in_plane([inside,eqAB], (* Atmospheric surface *)
            i,xFix4,xFix5,
            j,yFix1,yFix2,
            k,zFix5,zFix5) then
  change_node(i,j,k,fixed2,ConX,ConX,ConX,ConX,ConX,ConX);
end;

procedure fluxes(i:itype;j:jtype;k:ktype);
begin
if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix5,zFix5) then
  begin
    update_flxval(Flx1,bottom,i,j,k,plus); (* slab *)
    update_flxval(Flx4,bottom,i,j,k,plus); (* add to total house entry *)
  end;
if in_plane([inside,eqAB],
            i,xFix2,xFix3,
            j,yFix1,yFix2,
            k,zFix5,zFix5) then
  begin
    update_flxval(Flx2,bottom,i,j,k,plus); (* gap *)
    update_flxval(Flx4,bottom,i,j,k,plus); (* add to total house entry *)
  end;
if in_plane([inside,eqAB],
            i,xFix4,xFix5,
            j,yFix1,yFix2,
            k,zFix5,zFix5) then
  update_flxval(Flx3,bottom,i,j,k,plus); (* atm. surface *)
end; (* fluxes *)

procedure probes;
var c1,dc1:datatype;
    valid1:boolean;
begin
get_fieldvalue2d((wFixVal[xfix2].w+wFixVal[xfix3].w)/2,0.0005,wFixVal[zfix4].w,0.000001,c1,dc1,valid1);
obsval[obs1]:=c1;
get_fieldvalue2d((wFixVal[xfix3].w+wFixVal[xfix4].w)/2,0.001,wFixVal[zfix2].w-0.05,0.02,c1,dc1,valid1);
obsval[obs2]:=c1;
get_fieldvalue2d(wFixVal[xfix1].w+0.3,0.001,wFixVal[zfix1].w+0.3,0.02,c1,dc1,valid1);
obsval[obs3]:=c1;
get_fieldvalue2d(wFixVal[xfix5].w-0.3,0.001,wFixVal[zfix1].w+0.3,0.02,c1,dc1,valid1);
obsval[obs4]:=c1;
get_fieldvalue2d(wFixVal[xfix5].w-0.3,0.001,wFixVal[zfix5].w-0.3,0.02,c1,dc1,valid1);
obsval[obs5]:=c1;
end; (* probes *)

function materials(i:itype;j:jtype;k:ktype):mattype;
var mat:mattype;
begin
mat:=mat1; (* soil *)
if in_region(i,xFix1,xFix2,[inside,eqab],
             j,yFix1,yFix2,[inside,eqab],
             k,zFix4,zFix5,[inside,eqab]) then mat:=mat2; (* slab *)
if in_region(i,xFix1,xFix3,[inside,eqab],
             j,yFix1,yFix2,[inside,eqab],
             k,zFix3,zFix4,[inside,eqab]) then mat:=mat3; (* gravel *)
if in_region(i,xFix3,xFix4,[inside,eqab],
             j,yFix1,yFix2,[inside,eqab],
             k,zFix2,zFix5,[inside,eqab]) then mat:=mat4; (* footing *)
if in_region(i,xFix2,xFix3,[inside,eqab],
             j,yFix1,yFix2,[inside,eqab],
             k,zFix4,zFix5,[inside,eqab]) then mat:=mat5; (* gap *)
materials:=mat;
end; (* materials *)

function m(i:itype;j:jtype;k:ktype):datatype;
var mm:datatype;
begin
mm:=0;
case materials(i,j,k) of
  mat1: mm:=msat_soil;
  mat2: mm:=msat_slab;
  mat3: mm:=msat_gravel;
  mat4: mm:=msat_footing;
  mat5: mm:=msat_gap;
else
  error_std('m','Unknown material');
end;
m:=mm;
end;

function e_soilgas(i:itype;j:jtype;k:ktype):datatype;
begin
e_soilgas:=0;
end;

function beta_soilgas(i:itype;j:jtype;k:ktype):datatype;
begin
beta_soilgas:=0;
end;

function G_soilgas(i:itype;j:jtype;k:ktype):datatype;
begin
G_soilgas:=0;
end;

function Lambda_soilgas(i:itype;j:jtype;k:ktype):datatype;
begin
Lambda_soilgas:=0;
end;

function D_soilgas(dir:dirtype;i:itype;j:jtype;k:ktype):datatype;
var kk:datatype;
begin
kk:=0;
case materials(i,j,k) of
  mat1: kk:=k_soil;
  mat2: kk:=k_slab;
  mat3: kk:=k_gravel;
  mat4: kk:=k_footing;
  mat5: kk:=k_gap;
else
    error_std('D_soilgas','Unknown material');
end;
D_soilgas:=kk/mu
end;

function e_radon(i:itype;j:jtype;k:ktype):datatype;
var ee:datatype;
begin
ee:=0;
case materials(i,j,k) of
  mat1: ee:=etot_soil;
  mat2: ee:=etot_slab;
  mat3: ee:=etot_gravel;
  mat4: ee:=etot_footing;
  mat5: ee:=etot_gap;
else
    error_std('e','Unknown material');
end;
e_radon:=ee;
end;

function beta_radon(i:itype;j:jtype;k:ktype):datatype;
var ea,ew:datatype;
begin
ew:=m(i,j,k)*e_radon(i,j,k);
ea:=e_radon(i,j,k)-ew;
beta_radon:=ea+LOstwald*ew;
end;

function G_radon(i:itype;j:jtype;k:ktype):datatype;
var GG,ee,lam:datatype;
begin
GG:=0;
ee:=e_radon(i,j,k);
lam:=lambdaRn222;
case materials(i,j,k) of
  mat1: GG:=rho_g*(1-ee)/ee*lam*f_soil    * ARa_soil;
  mat2: GG:=rho_g*(1-ee)/ee*lam*f_slab    * ARa_slab;
  mat3: GG:=rho_g*(1-ee)/ee*lam*f_gravel  * ARa_gravel;
  mat4: GG:=rho_g*(1-ee)/ee*lam*f_footing * ARa_footing;
  mat5: GG:=rho_g*(1-ee)/ee*lam*f_gap     * ARa_gap;
else
  error_std('G','Unknown material');
end;
G_radon:=GG;
end;

function Lambda_radon(i:itype;j:jtype;k:ktype):datatype;
begin
Lambda_radon:=LambdaRn222;
end;

function D_radon(dir:dirtype;i:itype;j:jtype;k:ktype):datatype;
var DD:datatype;
begin
DD:=0;
case materials(i,j,k) of
  mat1: DD:=D_soil;
  mat2: DD:=D_slab;
  mat3: DD:=D_gravel;
  mat4: DD:=D_footing;
  mat5: DD:=D_gap;
else
    error_std('D_radon','Unknown material');
end;
D_radon:=DD;
end;

begin (* main *)
runid                        := '0130';
solution                     := steady;
geometry                     := cylindrical2d;
grid_def                     := grid;
flux_def                     := fluxes;
probe_def                    := probes;
materials_def                := materials;
wr_iteration_line_screen     := true;
wr_flux_during_calc_screen   := true;
wr_axes                      := false;

(* First do the soil-gas simulation *)
runtitle                     := 'Slab-on-grade house (pressure)';
boundary_conditions_def      := boundary_conditions_soilgas;
e_def                        := e_soilgas;
beta_def                     := beta_soilgas;
G_def                        := G_soilgas;
lambda_def                   := lambda_soilgas;
D_def                        := D_soilgas;
import_finalfield_guess      := true;
export_field                 := true;
flowfield                    := export;
import_field_name            := 'PRES00.dat';
export_field_name            := import_field_name;
relax_factor                 := 1.98;
flux_convset                 := [flx1..flx3];
probe_convset                := [obs1..obs4];
conv_evaluation_period       := 300;
min_iterations               := 150;
max_iterations               := 10000;
max_time                     := 60*60;
max_change                   := 1e-10;
max_residual_sum             := 1e-8;
run_model; (* Soil gas run *)

(* Then do the radon simulation *)
runtitle                     := 'Slab-on-grade house (radon)';
boundary_conditions_def      := boundary_conditions_radon;
e_def                        := e_radon;
beta_def                     := beta_radon;
G_def                        := G_radon;
lambda_def                   := lambda_radon;
D_def                        := D_radon;
import_finalfield_guess      := true;
export_field                 := true;
flowfield                    := import;
import_field_name            := 'Rn0000.dat';
export_field_name            := import_field_name;
relax_factor                 := 1.0;
flux_convset                 := [flx1..flx3];
probe_convset                := [obs1..obs4];
conv_evaluation_period       := 300;
min_iterations               := 150;
max_iterations               := 20000;
max_time                     := 60*60;
max_change                   := 1e-10;
max_residual_sum             := 1e-8;
run_model; (* Radon run *)

writeln(LOG,'The total soil-gas entry into the house is (Flx4)      = ',FlxVal[Flx4].Q:16,' m3/s');
writeln(LOG,'The total radon entry into the house is (Flx4)         = ',FlxVal[Flx4].J:16,' Bq/s');
close_model;
end.