program f0103prg;
(* --------------------- RnMod3d jobfile ---------------------- *)
(* Project: User guide example: Transient gas flow (1D)         *)
(* Created: October 12, 1999                                    *)
(* Revised: July 17, 2000                                       *)

{$I R3dirs03}
uses R3Defi03,R3Main03,R3Writ03;

const mu     = 17.5e-6;
      easoil = 0.2;
      ksoil  = 1e-14;
      Tper   = 10*3600;
      phi    = pi/2;
      p1     = 3;
      P0     = 100000;
      omega  = 2*pi/Tper;
      Dp     = ksoil*P0/easoil/mu;
      Lz     = 5.0;

      z_obs1 = 0.2; (* probe locations *)
      z_obs2 = 1.0;
      z_obs3 = 2.5;
      z_obs4 = 4.0;
      z_obs5 = 4.8;

procedure grid;
begin
set_FixVal(xFix1,0.0);
set_FixVal(xFix2,1.0);
set_axis_single(xFix1,xFix2,1,FocusA,1.0);

set_FixVal(yFix1,0.0);
set_FixVal(yFix2,1.0);
set_axis_single(yFix1,yFix2,1,FocusA,1.0);

set_FixVal(zFix1,0.0);
set_FixVal(zFix2,Lz);
set_axis_double(zFix1,zFix2,10,10,FocusA,FocusB,2,2,0.5);
end;

procedure boundary_conditions_Soilgas(i:itype;j:jtype;k:ktype);
begin
cBC[fixed1]:=0;
if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix1,zFix1) then set_node(i,j,k,fixed1);

cBC[fixed2]:=0;
if tim>0 then
  cBC[fixed2]:=p1 * sin(omega*tim + phi);
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
var valid:boolean;
begin
obsval[obs1]:=fieldvalue(0.5,0.5, z_obs1,valid);
obsval[obs2]:=fieldvalue(0.5,0.5, z_obs2,valid);
obsval[obs3]:=fieldvalue(0.5,0.5, z_obs3,valid);
obsval[obs4]:=fieldvalue(0.5,0.5, z_obs4,valid);
obsval[obs5]:=fieldvalue(0.5,0.5, z_obs5,valid);
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
beta_soilgas:=easoil/P0;
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

begin (* main *)
runid                        := '0103';
runtitle                     := 'User guide example: Transient gas flow in slab';
solution                     := steady;
geometry                     := cartesian3d;
grid_def                     := grid;
force_new_grid_in_every_run  := false;
boundary_conditions_def      := boundary_conditions_soilgas;
flux_def                     := fluxes;
probe_def                    := probes;

materials_def                := materials;
D_def                        := D_soilgas;
e_def                        := e_soilgas;
beta_def                     := beta_soilgas;
G_def                        := G_soilgas;
lambda_def                   := lambda_soilgas;

flux_convset                 := [flx2];
probe_convset                := [obs1,obs2,obs3,obs4,obs5];
conv_evaluation_period       := 400;
min_iterations               := 70;
max_iterations               := 10000;
max_time                     := 5*60;
max_change                   := 1e-12;
max_residual_sum             := 3e-9;

wr_iteration_line_screen     := false;
wr_final_results_screen      := false;
wr_axes                      := false;
wr_node_numbers              := false;
wr_materials_volumes         := false;


(* First do steady-state for t=0 *)
tim:=0;
run_model;

(* Then do the unsteady part *)
solution:=unsteady;
dtim:=Tper/500;

(* Write header w. labels *)
writeln('tim/Tper':16,' ','dtim/Tper':16,' ','cBC[Fixed2]':16,' ','obsval[obs5]':16);
writeln(RES,'tim':16,', ',
            'hr':16,', ',
            'Patm':16,', ',
            'Q1':16,', ',
            'Q*2':16,', ',
            'P1':16,', ',
            'P2':16,', ',
            'P3':16,', ',
            'P4':16,', ',
            'P5':16);

repeat
  writeln(tim/Tper:16:4,' ',dtim/Tper:16:4,' ',cBC[Fixed2]:16:4,' ',obsval[obs5]:16:4);
  writeln(RES,tim:16,', ',
              tim/3600:16,', ',
              cBC[Fixed1]:16,', ',
              FlxVal[Flx1].j:16,', ',
              FlxVal[Flx2].j:16,', ',
              obsval[obs1]:16,', ',
              obsval[obs2]:16,', ',
              obsval[obs3]:16,', ',
              obsval[obs4]:16,', ',
              obsval[obs5]:16);
  tim:=tim+dtim;
  run_model;
until tim>4*Tper;

wr_gridfiles;
close_model;
end.
