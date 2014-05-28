open Bigarray
open Ida
open Ida_utils

(* constants taken from c-#defines *)
let nout = 11
let mm = 10
let neq = mm * mm
let mmf = float_of_int mm
let dx = 1. /. (mmf -. 1.)
let coeff = 1. /. ( dx *. dx )
let bval = 0.1

let heatres {t=tres; yp=up; yy=uu; res=resval} = 
  (* Initialize resval to uu, to take care of boundary equations. *)
  Array1.blit uu resval ;

  (* Loop over interior points; set res = up - (central difference). *)
  for j = 1 to mm - 2 do
    let offset = mm*j in
    for i = 1 to mm - 2 do      
      let loc = offset + i in
      let f : float = 
	(up.{loc} -. coeff *. ( 
		    uu.{loc-1} +. uu.{loc+1} +. 
		      uu.{loc-mm} +. uu.{loc+mm} -. 4. *. uu.{loc} )) in
      resval.{loc} <- f (*;
  
      Printf.printf "upv[%d] == %f  " loc (up.{loc}) ;
      Printf.printf "resv[%d] = %f\n" loc (resval.{loc})*)
    done;
  done;
  0

   
(* 
 * Print first lines of output (problem description)
 *)
let print_header rtol atol =
  Printf.printf "\nidaHeat2D_bnd: Heat equation, serial example problem for IDA\n";
  Printf.printf "          Discretized heat equation on 2D unit square.\n";
  Printf.printf "          Zero boundary conditions,";
  Printf.printf " polynomial initial conditions.\n";
  Printf.printf "          Mesh dimensions: %d x %d" mm mm;
  Printf.printf "        Total system size: %d\n\n"  neq;
  Printf.printf "Tolerance parameters:  rtol = %g   atol = %g\n"  rtol atol;
  Printf.printf "Constraints set to force all solution components >= 0. \n";
  Printf.printf "Linear solver: IDABAND, banded direct solver \n";
  Printf.printf "       difference quotient Jacobian, half-bandwidths = %d \n" mm;
  Printf.printf "IDACalcIC called with input boundary values = %g \n" bval;
  (* Print output table heading and initial line of table. *)
  Printf.printf "\n   Output Summary (umax = max-norm of solution) \n\n";
  Printf.printf "  time       umax     k  nst  nni  nje   nre   nreLS    h      \n" ;
  Printf.printf " .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . \n" 

(*
 * Print Output
 *)
let print_output ida t uu = 
  let umax = fvector_max_norm uu in
  let kused = ref 0 in
  check_flag "IDAGetLastOrder" (ida_get_last_order ida kused) ;

  let nst = ref 0 in
  check_flag "IDAGetNumSteps" (ida_get_num_steps ida nst) ;

  let nni = ref 0 in
  check_flag "IDAGetNumNonlinSolvIters" (ida_get_num_nonlin_solv_iters ida nni) ;

  let nre = ref 0 in
  check_flag "IDAGetNumResEvals" (ida_get_num_res_evals ida nre) ;

  let hused = ref 0. in
  check_flag "IDAGetLastStep" (ida_get_last_step ida hused) ;

  let nje = ref 0 in
  check_flag "IDAGetNumJacEvals" (ida_dls_get_num_jac_evals ida nje) ;

  let nreLS = ref 0 in
  check_flag "IDADlsGetNumResEvals" (ida_dls_get_num_res_evals ida nreLS) ;
  Printf.printf " %5.2f %13.5e  %d  %3d  %3d  %3d  %4d  %4d  %9.2e \n"  t umax !kused !nst !nni !nje !nre !nreLS !hused
		

let set_initial_profile uu up id res = 
  (* Initialize id to 1's. *)
  Array1.fill id 1. ;

  (* Initialize uu on all grid points. *)
  for j = 0 to mm - 1 do
    let yfact = dx *. (float_of_int j) in
    let offset = mm * j in
    for i = 0 to mm - 1 do
      let xfact = dx *. (float_of_int i) in
      let loc = offset + i in
      uu.{loc} <- (16. *. xfact *. (1. -. xfact) *. yfact *. (1. -. yfact))
    done      
  done ;

  (* heatres sets res to negative of ODE RHS values at interior points. *)
  ignore ( heatres {t=0.; yp=up; yy=uu; res=res} ) ;
  
  (* Copy -res into up to get correct interior initial up values. *)
  fvector_scale (-1.) res up;

  (* Finally, set values of u, up, and id at boundary points. *)
  for j = 0 to mm - 1 do
    let offset = mm * j in
    for i = 0 to mm - 1 do
      let loc = offset + i in
      if j = 0 || j = (mm - 1) || i = 0 || i = (mm - 1) then (
	uu.{loc} <- bval; 
        up.{loc} <- 0.;
        id.{loc} <- 0. ) 
    done
  done

let () = 
  let uu = fvector neq in
  let up = fvector neq in
  let id = fvector neq in
  let res = fvector neq in
  let constraints = fvector neq in

  (* Initialize uu, up, id. *)
  set_initial_profile uu up id res;

  (* Set constraints to all 1's for nonnegative solution values. *)
  Array1.fill constraints 1. ;

  (* Set remaining input parameters. *)
  let t0 = 0. in
  let t1 = 0.01 in
  let rtol = 0. in
  let atol = 1.0e-3 in

  (* Call IDACreate and IDAMalloc to initialize solution *)
  let ida = ida_create () in

  check_flag "IDASetId" (ida_set_id ida id) ;

  check_flag "IDASetConstraints" (ida_set_constraints ida constraints) ;

  let num_state = { t = 0.0; yp = up; yy = uu; res = res } in

  let event_state = { e_t = 0.0; e_y = fvector neq ; e_yp = fvector neq ; e_gi = fvector 0 } in
  
  let event_roots state = 0 in

  let ida_ctxt = { residual = heatres ; state = num_state ; event_roots ; event_state } in

  check_flag "IDAInit" (ida_init ida ida_ctxt) ;

  check_flag "IDASStolerances" (ida_ss_tolerances ida rtol atol);

  (* Call IDABand to specify the linear solver. *)
  let mu = mm in
  let ml = mm in

  check_flag "IDABand" (ida_band ida neq mu ml) ;
 
  (* Call IDACalcIC to correct the initial values. *)
  check_flag "IDACalcIC" (ida_calc_ic ida ida_ya_ydp_init t1) ;

  (* Print output heading. *)
  print_header rtol atol;
  
  print_output ida t0 uu;

  (* Loop over output times, call IDASolve, and print results. *)
  let rec sim_loop step = function
    | iout when iout = nout + 1 -> ()
    | iout -> check_flag "IDASolve" (ida_solve ida step ida_normal) ; 
	      print_output ida step.tret uu ; 
	      sim_loop { step with tout = step.tout *. 2. } (iout + 1)
  in
  sim_loop { tout = t1 ; tret = 0.} 1 ;
    
  (* Print remaining counters *)
  let netf = ref 0 in
  check_flag "IDAGetNumErrTestFails" (ida_get_num_err_test_fails ida netf);

  let ncfn = ref 0 in
  check_flag "IDAGetNumNonlinSolvConvFails" (ida_get_num_nonlin_solv_conv_fails ida ncfn);

  Printf.printf "\n netf = %d, ncfn = %d \n" (!netf) (!ncfn)
