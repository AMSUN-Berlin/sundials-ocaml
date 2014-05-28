(*
 * Copyright (c) 2014, TU Berlin
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *   * Neither the name of the TU Berlin nor the
 *     names of its contributors may be used to endorse or promote products
 *     derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL TU Berlin BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *)

open Bigarray
open Batteries

type fvector = (float, float64_elt, c_layout) Array1.t

type ida_solver

external ida_create : unit -> ida_solver = "sundials_ml_ida_create"

type numeric_state = {
  t : float;
  yp : fvector ;
  yy : fvector;
  res : fvector;
}

type event_state = {
  e_t : float ;
  e_y : fvector ;
  e_yp : fvector ;
  e_gi : fvector ;
}

type residual = numeric_state -> int

type events = event_state -> int

type ida_ctxt = {
  residual : residual;
  state : numeric_state;
  event_roots : events ;
  event_state : event_state
}

type solution_step = {
  tout : float;
  tret : float;
}

external ida_init : ida_solver -> ida_ctxt -> int = "sundials_ml_ida_init"

external ida_reinit : ida_solver -> float -> ida_ctxt -> int = "sundials_ml_ida_reinit"

external ida_get_ctxt : ida_solver -> ida_ctxt = "sundials_ml_ida_ctxt"

(* flags *)
external _ida_normal : unit -> int = "sundials_ml_ida_normal"
external _ida_one_step : unit -> int = "sundials_ml_ida_one_step"
external _ida_ya_ydp_init : unit -> int = "sundials_ml_ida_ya_ydp_init"
external _ida_y_init : unit -> int = "sundials_ml_ida_y_init"
external _ida_root_return : unit -> int = "sundials_ml_ida_root_return"

let ida_normal = _ida_normal ()
let ida_one_step = _ida_one_step ()
let ida_ya_ydp_init = _ida_ya_ydp_init ()
let ida_y_init = _ida_y_init ()
let ida_root_return = _ida_root_return ()

external ida_solve : ida_solver -> solution_step -> int -> int = "sundials_ml_ida_solve"

external ida_set_id : ida_solver -> fvector -> int = "sundials_ml_ida_set_id"

external ida_set_constraints : ida_solver -> fvector -> int = "sundials_ml_ida_set_constraints"

external ida_ss_tolerances : ida_solver -> float -> float -> int = "sundials_ml_ida_ss_tolerances"

external ida_band : ida_solver -> int -> int -> int -> int = "sundials_ml_ida_band"

external ida_dense : ida_solver -> int -> int = "sundials_ml_ida_dense"

external ida_calc_ic : ida_solver -> int -> float -> int = "sundials_ml_ida_calc_ic"

external ida_get_last_order : ida_solver -> int ref -> int = "sundials_ml_ida_get_last_order"

external ida_get_num_steps : ida_solver -> int ref -> int = "sundials_ml_ida_get_num_steps"

external ida_get_num_steps : ida_solver -> int ref -> int = "sundials_ml_ida_get_num_steps"

external ida_get_num_nonlin_solv_iters : ida_solver -> int ref -> int = "sundials_ml_ida_get_num_nonlin_solv_iters"

external ida_get_num_res_evals : ida_solver -> int ref -> int = "sundials_ml_ida_get_num_res_evals"

external ida_get_last_step : ida_solver -> float ref -> int = "sundials_ml_ida_get_last_step"

external ida_dls_get_num_jac_evals : ida_solver -> int ref -> int = "sundials_ml_ida_dls_get_num_jac_evals"

external ida_dls_get_num_res_evals : ida_solver -> int ref -> int = "sundials_ml_ida_dls_get_num_res_evals"

external ida_get_num_err_test_fails : ida_solver -> int ref -> int = "sundials_ml_ida_get_num_err_test_fails"

external ida_get_num_nonlin_solv_conv_fails : ida_solver -> int ref -> int = "sundials_ml_ida_get_num_nonlin_solv_conv_fails"

external ida_get_root_info : ida_solver -> int array -> int = "sundials_ml_ida_get_root_info"

let check_flag call_msg = function 
    err when err < 0 -> raise (Failure (Printf.sprintf "SUNDIALS_ERROR: %s() failed with flag = %d" call_msg err))
  | _ -> ()
    
