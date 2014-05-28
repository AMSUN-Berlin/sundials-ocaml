/*
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
 */

#include <stdio.h>
#include <assert.h>

#include <caml/alloc.h>
#include <caml/custom.h>
#include <caml/memory.h>
#include <caml/fail.h>
#include <caml/callback.h>
#include <caml/bigarray.h>

#include <ida/ida.h>
#include <ida/ida_dense.h>
//#include <ida/ida_lapack.h> //TODO: This file is not installed?
#include <ida/ida_band.h>
#include <ida/ida_spgmr.h>
#include <ida/ida_spbcgs.h>
#include <ida/ida_sptfqmr.h>

#include <sundials/sundials_nvector.h>
#include <nvector/nvector_serial.h>

struct sundials_ml_ida_ctxt {
  value* root; //root of globally stored OCaml memory
  void* mem; //IDA memory
};

/* Ocaml def. of the root memory:
type numeric_state = {
  t0 : float;
  yp : fvector ;
  yy : fvector;
  res : fvector;
}

type event_state = {
  t : float ;
  y : fvector ;
  yp : fvector ;
  gi : fvector ;
}

type residual = numeric_state -> int

type events = event_state -> int

type ida_ctxt = {
  residual : residual;
  state : numeric_state;
  events : events ;
  event_state : event_state
}
*/

#define IDA_MEM(val) ( ( (struct sundials_ml_ida_ctxt*) Data_custom_val(val) ) -> mem )

#define ROOT(val) ( ( (struct sundials_ml_ida_ctxt*) Data_custom_val(val) ) -> root )

#define IDA_CTXT(val) ( *ROOT(val) ) /* just dereference the root-pointer to get the ocaml value */

#define RESIDUAL(val) ( Field(IDA_CTXT(val), 0 ) )
#define NUMSTATE(val) ( Field(IDA_CTXT(val), 1 ) )
#define EVENTS(val) ( Field(IDA_CTXT(val), 2 ) )
#define EVENTSTATE(val) ( Field(IDA_CTXT(val), 3 ) )

#define NUMSTATE_T0(val) ( Field(NUMSTATE(val), 0 ) )
#define NUMSTATE_YP(val) ( Field(NUMSTATE(val), 1 ) )
#define NUMSTATE_YY(val) ( Field(NUMSTATE(val), 2 ) ) 
#define NUMSTATE_RES(val) ( Field(NUMSTATE(val), 3 ) ) 


/* VERY big magic to avoid allocations
   The theory here is, that:
    1. all N_Vector operations are uniform 
       (no operation impl queries the ops field)
    2. no (serial) operation actually expects heap data
       (I checked for version 2.5.0, but anything else
        would be silly)
    3. Hence we can hand over stack-pointers
    4. We can avoid _any_ allocation that way and still
       mask any OCaml array as a N_Vector_Serial     
 */
#define BA_STACK_NVECTOR(X, NV)                             \
  struct caml_ba_array* X##__BA = Caml_ba_array_val(X);     \
  double* X##__data = X##__BA->data ;                       \
  struct _N_VectorContent_Serial X##_serial_wrap = {        \
    X##__BA->dim[0], /* length */                           \
    FALSE, /* own_data */                                   \
    (realtype*)X##__data, /* data */                        \
  };                                                        \
  struct _generic_N_Vector NV = {         	            \
    &X##_serial_wrap,				            \
    &static_serial_ops,                                     \
  };                                                        \

void static_vector_destroy(N_Vector v) { }

static struct _generic_N_Vector_Ops static_serial_ops = {
  N_VClone_Serial,
  N_VCloneEmpty_Serial,
  static_vector_destroy,
  N_VSpace_Serial,
  N_VGetArrayPointer_Serial,
  N_VSetArrayPointer_Serial,
  N_VLinearSum_Serial,
  N_VConst_Serial,
  N_VProd_Serial,
  N_VDiv_Serial,
  N_VScale_Serial,
  N_VAbs_Serial,
  N_VInv_Serial,
  N_VAddConst_Serial,
  N_VDotProd_Serial,
  N_VMaxNorm_Serial,
  N_VWrmsNorm_Serial,
  N_VWrmsNormMask_Serial,
  N_VMin_Serial,
  N_VWL2Norm_Serial,
  N_VL1Norm_Serial,
  N_VCompare_Serial,
  N_VInvTest_Serial,
  N_VConstrMask_Serial,
  N_VMinQuotient_Serial,
};

CAMLprim value sundials_ml_fvector_scale(value s, value x, value z) {
  CAMLparam3(s,x,z);
  const double ds = Double_val(s);
  struct caml_ba_array* ba_x = Caml_ba_array_val(x);
  struct caml_ba_array* ba_z = Caml_ba_array_val(z);
  double* dx = (double*) ba_x -> data;
  double* dz = (double*) ba_z -> data;

  for(int i = 0; i < ba_x->dim[0]; i++)
    dz[i] = dx[i] * ds;

  CAMLreturn(Val_unit);
}

CAMLprim value sundials_ml_fvector_get(value a, value i) {
  CAMLparam2(a, i);
  double* d = Caml_ba_array_val(a)->data;
  CAMLreturn(caml_copy_double(d[Int_val(i)]));
}

CAMLprim value sundials_ml_fvector_maxnorm(value fvector) {
  CAMLparam1(fvector);
  BA_STACK_NVECTOR(fvector, nv_f);
  const double n = N_VMaxNorm_Serial(&nv_f);
  CAMLreturn(caml_copy_double(n));
}

/**
 * Garbage collector callback
 * this is actually unlikely to ever get called, since the integrator will probably
 * be in use for the whole runtime and ocaml does not seem to run any finalizers on
 * program termination 
 */
void sundials_ml_ida_solver_finalize(value v) {
  printf("Finalizing ida solver...\n");

  caml_remove_global_root(ROOT(v));  

  free(ROOT(v));
  IDAFree(IDA_MEM(v));
}

CAMLprim value sundials_ml_ida_create(value unit) {

  static struct custom_operations ida_ctxt_ops = {
    "sundials_ml_ida_solver",
    sundials_ml_ida_solver_finalize,
    custom_compare_default,
    custom_hash_default,
    custom_serialize_default,
    custom_deserialize_default,
  };

  CAMLparam0 ();
  CAMLlocal1 (block);
 
  block = caml_alloc_custom(&ida_ctxt_ops, sizeof(struct sundials_ml_ida_ctxt), 1, 10);

  ROOT(block) = malloc(sizeof(value));
  IDA_MEM(block) = IDACreate();

  IDASetUserData(IDA_MEM(block), ROOT(block));

  CAMLreturn (block);
}

int sundials_ml_event_wrapper(realtype tt, N_Vector yy, N_Vector yp, realtype *gout, void* user_data) {
  value ev = Field(*(value*)user_data, 2);
  value ev_state = Field(*(value*)user_data, 3);

  double* t = (double*)Field(ev_state, 0);
  *t = tt;

  double* old_y  = Caml_ba_array_val(Field(ev_state, 1))->data;
  double* old_yp = Caml_ba_array_val(Field(ev_state, 2))->data;
  double* old_gi = Caml_ba_array_val(Field(ev_state, 3))->data;
  
  double* new_y  =  NV_DATA_S(yy);
  double* new_yp =  NV_DATA_S(yp);

  Caml_ba_array_val(Field(ev_state, 1))->data = new_y;
  Caml_ba_array_val(Field(ev_state, 2))->data = new_yp;
  Caml_ba_array_val(Field(ev_state, 3))->data = gout;

  value ret = caml_callback(ev, ev_state);

  /* because we might have triggered a GC cycle, num_state can be invalid */
  ev_state = Field(*(value*)user_data, 3);

  Caml_ba_array_val(Field(ev_state, 1))->data = old_y;
  Caml_ba_array_val(Field(ev_state, 2))->data = old_yp;
  Caml_ba_array_val(Field(ev_state, 3))->data = old_gi;

  return Int_val (ret);  
} 

int sundials_ml_residual_wrapper(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data) {

  value res = Field(*(value*)user_data, 0);
  value num_state = Field(*(value*)user_data, 1);

  double* t = (double*)Field(num_state, 0);
  *t = tt;

  double* old_yp = Caml_ba_array_val(Field(num_state, 1))->data;
  double* old_yy = Caml_ba_array_val(Field(num_state, 2))->data;
  double* old_rr = Caml_ba_array_val(Field(num_state, 3))->data;

  double* new_yy =  NV_DATA_S(yy);
  double* new_yp =  NV_DATA_S(yp);
  double* new_rr =  NV_DATA_S(rr);

  Caml_ba_array_val(Field(num_state, 1))->data = new_yp;
  Caml_ba_array_val(Field(num_state, 2))->data = new_yy;
  Caml_ba_array_val(Field(num_state, 3))->data = new_rr;

  value ret = caml_callback(res, num_state);
  /* because we might have triggered a GC cycle, num_state can be invalid */
  num_state = Field(*(value*)user_data, 1);

  Caml_ba_array_val(Field(num_state, 1))->data = old_yp;
  Caml_ba_array_val(Field(num_state, 2))->data = old_yy;
  Caml_ba_array_val(Field(num_state, 3))->data = old_rr;

  return Int_val (ret);
}

CAMLprim value sundials_ml_ida_init(value ida_solver, value ida_ctxt) {
  CAMLparam2(ida_solver, ida_ctxt);

  assert (Tag_val(ida_ctxt) == 0);
  assert (Tag_val(Field(ida_ctxt, 0)) == Closure_tag);
  assert (Tag_val(Field(ida_ctxt, 1)) == 0 );
  assert (Tag_val(Field(Field(ida_ctxt, 1), 0)) == Double_tag );

  IDA_CTXT(ida_solver) = ida_ctxt;
  caml_register_global_root(&IDA_CTXT(ida_solver));  

  const realtype rt_t0 = Double_val(NUMSTATE_T0(ida_solver));
  value y0 = NUMSTATE_YY(ida_solver);
  value yp0 = NUMSTATE_YP(ida_solver);

  BA_STACK_NVECTOR(y0, nv_y0);
  BA_STACK_NVECTOR(yp0, nv_yp0);

  value gi = Field(EVENTSTATE(ida_solver), 3);
  const intnat ev_len = Caml_ba_array_val(gi)->dim[0];
  
  const int ret = IDAInit(IDA_MEM(ida_solver), &sundials_ml_residual_wrapper, rt_t0, &nv_y0, &nv_yp0);

  if (ev_len > 0) {
    IDARootInit(IDA_MEM(ida_solver), ev_len, sundials_ml_event_wrapper);
  }

  CAMLreturn(Val_int(ret));   
}

CAMLprim value sundials_ml_ida_reinit(value ida_solver, value t, value ida_ctxt) {
  CAMLparam2(ida_solver, ida_ctxt);

  const realtype rt_t0 = Double_val(t);
  value y0 = NUMSTATE_YY(ida_solver);
  value yp0 = NUMSTATE_YP(ida_solver);

  BA_STACK_NVECTOR(y0, nv_y0);
  BA_STACK_NVECTOR(yp0, nv_yp0);

  const int ret = IDAReInit(IDA_MEM(ida_solver), rt_t0, &nv_y0, &nv_yp0);
  CAMLreturn(Val_int(ret));  
}

CAMLprim value sundials_ml_ida_ctxt(value ida) {
  CAMLparam1(ida);
  CAMLreturn(IDA_CTXT(ida));
}

/*
 * Tolerances
 */
CAMLprim value sundials_ml_ida_ss_tolerances(value ida_solver, value reltol, value abstol) {
  CAMLparam3(ida_solver, reltol, abstol);
  const int ret = IDASStolerances(IDA_MEM(ida_solver), Double_val(reltol), Double_val(abstol));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_sv_tolerances(value ida_solver, value reltol, value abstol) {
  CAMLparam3(ida_solver, reltol, abstol);
  BA_STACK_NVECTOR(abstol, nv_abstol);
  const int ret = IDASVtolerances(IDA_MEM(ida_solver), Double_val(reltol), &nv_abstol);
  CAMLreturn(Val_int(ret));  
}

//TODO: CAMLprim value sundials_ml_ida_wf_tolerances(value ida_solver, value efun);

/*
 * Linear Solvers
 */
CAMLprim value sundials_ml_ida_dense(value ida_solver, value N) {
  CAMLparam2(ida_solver, N);
  const int ret = IDADense(IDA_MEM(ida_solver), Int_val(N));
  CAMLreturn(Val_int(ret));
}

/*
TODO: fix lapack support
CAMLprim value sundials_ml_ida_lapack_dense(value ida_solver, value N) {
  CAMLparam2(ida_solver, N);
  const int ret = IDALapackDense(IDA_MEM(ida_solver), Int_val(N));
  CAMLreturn(Val_int(ret));
}
*/

CAMLprim value sundials_ml_ida_band(value ida_solver, value N, value mupper, value mlower) {
  CAMLparam4(ida_solver, N, mupper, mlower);
  const int ret = IDABand(IDA_MEM(ida_solver), Int_val(N), Int_val(mupper), Int_val(mlower));
  CAMLreturn(Val_int(ret));
}

/*
TODO: fix lapack support
CAMLprim value sundials_ml_ida_lapack_band(value ida_solver, value N, value mupper, value mlower) {
  CAMLparam4(ida_solver, N, mupper, mlower);
  const int ret = IDALapackBand(IDA_MEM(ida_solver), Int_val(N), Int_val(mupper), Int_val(mlower));
  CAMLreturn(Val_int(ret));
}
*/

CAMLprim value sundials_ml_ida_spmgr(value ida_solver, value maxl) {
  CAMLparam2(ida_solver, maxl);
  const int ret = IDASpgmr(IDA_MEM(ida_solver), Int_val(maxl));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_spbcg(value ida_solver, value maxl) {
  CAMLparam2(ida_solver, maxl);
  const int ret = IDASpbcg(IDA_MEM(ida_solver), Int_val(maxl));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_sptfqmr(value ida_solver, value maxl) {
  CAMLparam2(ida_solver, maxl);
  const int ret = IDASptfqmr(IDA_MEM(ida_solver), Int_val(maxl));
  CAMLreturn(Val_int(ret));
}

/*
 * Initial condition
 */
CAMLprim value sundials_ml_ida_calc_ic(value ida_solver, value icopt, value tout1) {
  CAMLparam3(ida_solver, icopt, tout1);
  const int ret = IDACalcIC(IDA_MEM(ida_solver), Int_val(icopt), Double_val(tout1));
  CAMLreturn(Val_int(ret));
}

/*
 * Solver
 */
CAMLprim value sundials_ml_ida_solve(value ida_solver, value solution_step, value itask) {
  CAMLparam3(ida_solver, solution_step, itask);
  
  /* solution_step contains only floats, so we use the OCaml optimized data layout */
  realtype tout = ((realtype*) solution_step)[0];
  realtype* tret = ((realtype*) solution_step) + 1;

  value y = NUMSTATE_YY(ida_solver);
  value yp = NUMSTATE_YP(ida_solver);

  BA_STACK_NVECTOR(y, nv_yret);
  BA_STACK_NVECTOR(yp, nv_ypret);

  const int ret = IDASolve(IDA_MEM(ida_solver), tout, tret, &nv_yret, &nv_ypret, Int_val(itask));
  CAMLreturn(Val_int(ret));
}

/*
 * Flags
 */
CAMLprim value sundials_ml_ida_root_return(value unit) {
  CAMLparam1(unit);
  int flag = IDA_ROOT_RETURN;
  CAMLreturn(Val_int(flag));
}

CAMLprim value sundials_ml_ida_normal(value unit) {
  CAMLparam1(unit);
  int flag = IDA_NORMAL;
  CAMLreturn(Val_int(flag));
}

CAMLprim value sundials_ml_ida_one_step(value unit) {
  CAMLparam1(unit);
  CAMLreturn(Val_int(IDA_ONE_STEP));
}

CAMLprim value sundials_ml_ida_ya_ydp_init(value unit) {
  CAMLparam1(unit);
  CAMLreturn(Val_int(IDA_YA_YDP_INIT));
}

CAMLprim value sundials_ml_ida_y_init(value unit) {
  CAMLparam1(unit);
  CAMLreturn(Val_int(IDA_Y_INIT));
}

/*
 * Main solver optional input functions
 * TODO: complete
 */
CAMLprim value sundials_ml_ida_set_max_ord(value ida_solver, value maxord) {
  CAMLparam2(ida_solver, maxord);
  const int ret = IDASetMaxOrd(IDA_MEM(ida_solver), Int_val(maxord));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_set_max_num_steps(value ida_solver, value maxsteps) {
  CAMLparam2(ida_solver, maxsteps);
  const int ret = IDASetMaxNumSteps(IDA_MEM(ida_solver), Int_val(maxsteps));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_set_init_step(value ida_solver, value hin) {
  CAMLparam2(ida_solver, hin);
  const int ret = IDASetInitStep(IDA_MEM(ida_solver), Double_val(hin));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_set_max_step(value ida_solver, value hmax) {
  CAMLparam2(ida_solver, hmax);
  const int ret = IDASetMaxStep(IDA_MEM(ida_solver), Double_val(hmax));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_set_stop_time(value ida_solver, value tstop) {
  CAMLparam2(ida_solver, tstop);
  const int ret = IDASetStopTime(IDA_MEM(ida_solver), Double_val(tstop));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_set_max_err_test_fails(value ida_solver, value maxnef) {
  CAMLparam2(ida_solver, maxnef);
  const int ret = IDASetMaxErrTestFails(IDA_MEM(ida_solver), Int_val(maxnef));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_set_id(value ida_solver, value id)  {
  CAMLparam2(ida_solver, id);
  BA_STACK_NVECTOR(id, nv_id);
  const int ret = IDASetId(IDA_MEM(ida_solver), &nv_id);
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_set_constraints(value ida_solver, value constraints)  {
  CAMLparam2(ida_solver, constraints);
  BA_STACK_NVECTOR(constraints, nv_constraints);
  const int ret = IDASetConstraints(IDA_MEM(ida_solver), &nv_constraints);
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_get_root_info(value ida_solver, value rootsfound) {
  CAMLparam2(ida_solver, rootsfound);
  int len = Wosize_val(rootsfound);
  int roots[len];
  /* totally stupid API from the IDA side, they copy the array and then we copy ... */
  const int ret = IDAGetRootInfo(IDA_MEM(ida_solver), roots);

  for (int i = 0; i < len; i++)
    Field(rootsfound, i) = Val_int(roots[i]);

  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_get_last_order(value ida_solver, value klast) {
  CAMLparam2(ida_solver, klast);
  int _klast;
  const int ret = IDAGetLastOrder(IDA_MEM(ida_solver), &_klast);
  Store_field(klast, 0, Val_int(_klast));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_get_num_steps(value ida_solver, value nsteps) {
  CAMLparam2(ida_solver, nsteps);
  long int _nsteps;
  const int ret = IDAGetNumSteps(IDA_MEM(ida_solver), &_nsteps);
  Store_field(nsteps, 0, Val_int(_nsteps));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_get_num_nonlin_solv_iters(value ida_solver, value nniters) {
  CAMLparam2(ida_solver, nniters);
  long int _nniters;
  const int ret = IDAGetNumNonlinSolvIters(IDA_MEM(ida_solver), &_nniters);
  Store_field(nniters, 0, Val_int(_nniters));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_get_num_nonlin_solv_conv_fails(value ida_solver, value nncfails) {
  CAMLparam2(ida_solver, nncfails);
  long int _nncfails;
  const int ret = IDAGetNumNonlinSolvConvFails(IDA_MEM(ida_solver), &_nncfails);
  Store_field(nncfails, 0, Val_int(_nncfails));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_get_last_step(value ida_solver, value hlast) {
  CAMLparam2(ida_solver, hlast);
  double* _hlast = (double*)Field(hlast, 0);
  const int ret = IDAGetLastStep(IDA_MEM(ida_solver), _hlast);
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_get_num_res_evals(value ida_solver, value nrevals) {
  CAMLparam2(ida_solver, nrevals);
  long int _nrevals;
  const int ret = IDAGetNumResEvals(IDA_MEM(ida_solver), &_nrevals);
  Store_field(nrevals, 0, Val_int(_nrevals));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_dls_get_num_res_evals(value ida_solver, value nrevalsLS) {
  CAMLparam2(ida_solver, nrevalsLS);
  long int _nrevalsLS;
  const int ret = IDADlsGetNumResEvals(IDA_MEM(ida_solver), &_nrevalsLS);
  Store_field(nrevalsLS, 0, Val_int(_nrevalsLS));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_dls_get_num_jac_evals(value ida_solver, value njevals) {
  CAMLparam2(ida_solver, njevals);
  long int _njevals;
  const int ret = IDADlsGetNumJacEvals(IDA_MEM(ida_solver), &_njevals);
  Store_field(njevals, 0, Val_int(_njevals));
  CAMLreturn(Val_int(ret));
}

CAMLprim value sundials_ml_ida_get_num_err_test_fails(value ida_solver, value netfails) {
  CAMLparam2(ida_solver, netfails);
  long int _netfails;
  const int ret = IDAGetNumErrTestFails(IDA_MEM(ida_solver), &_netfails);
  Store_field(netfails, 0, Val_int(_netfails));
  CAMLreturn(Val_int(ret));
}




