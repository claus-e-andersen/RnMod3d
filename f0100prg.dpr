program F0100prg;
(* --------------------- RnMod3d jobfile ---------------------- *)
(* Project: User guide example: Steady soil-gas flow, 1D        *)
(* Created: May 24, 1999                                        *)
(* Revised: July 17, 2000                                       *)

{$I R3dirs03}
uses R3Defi03,R3Main03,R3Writ03;

procedure grid;
begin
set_FixVal(xFix1,0.0);
set_FixVal(xFix2,1.0);
set_axis_single(xFix1,xFix2,1,FocusA,1.0);

set_FixVal(yFix1,0.0);
set_FixVal(yFix2,1.0);
set_axis_single(yFix1,yFix2,1,FocusA,1.0);

set_FixVal(zFix1,-3.0);
set_FixVal(zFix2, 0.0);
set_axis_double(zFix1,zFix2,30,30,FocusA,FocusB,1.1,1.1,0.5);
end;

procedure boundary_conditions(i:itype;j:jtype;k:ktype);
begin
cBC[fixed1]:=0;
cBC[fixed2]:=-3.0;
if in_plane([inside,eqAB],
            i,xFix1,xFix2,
            j,yFix1,yFix2,
            k,zFix1,zFix1) then set_node(i,j,k,fixed1);
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
cc:=fieldvalue(0.5,0.5,-1.5,valid);
if not valid then cc:=0.0;
obsval[obs1]:=cc;
end;

function materials(i:itype;j:jtype;k:ktype):mattype;
begin
materials:=mat1;
end;

function e(i:itype;j:jtype;k:ktype):datatype;
begin
e:=0;
end;

function beta(i:itype;j:jtype;k:ktype):datatype;
begin
beta:=0;
end;

function D(dir:dirtype;i:itype;j:jtype;k:ktype):datatype;
var mu:datatype;
begin
mu:=17.5e-6;
D:=2e-10/mu;
end;

function G(i:itype;j:jtype;k:ktype):datatype;
begin
G:=0;
end;

function lambda(i:itype;j:jtype;k:ktype):datatype;
begin
lambda:=0;
end;

begin (* main *)
runid                          := '0100';
runtitle                       := 'User guide example: Steady soil-gas flow, 1D';
solution                       := steady;
geometry                       := cartesian3d;
Ly                             := 1.0;
grid_def                       := grid;
force_new_grid_in_every_run    := false;
boundary_conditions_def        := boundary_conditions;
flux_def                       := fluxes;
probe_def                      := probes;
materials_def                  := materials;
e_def                          := e;
beta_def                       := beta;
G_def                          := G;
lambda_def                     := lambda;
D_def                          := D;
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
flux_convset                   := [flx1,flx2];
probe_convset                  := [obs1];
conv_evaluation_period         := 100;
min_iterations                 := 50;
max_iterations                 := 5000;
max_time                       := 5*60;
max_change                     := 1e-9;
max_residual_sum               := 3e-20;
dtim                           := 0;
BC_running                     := false;
BC_running_update_of_cBCs_def  := nil;
BC_running_min_iterations      := 100;
BC_running_max_residual_sum_before_new_BC := 1e-9;
BC_running_convergence_def     := nil;
wr_BC_running_messages_log     := false;
wr_BC_running_messages_screen  := false;
press_enter_wanted             := true;

run_model;
close_model;
end.