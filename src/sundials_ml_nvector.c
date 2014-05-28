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

#include "sundials_ml_nvector.h"

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
  
    This macro takes an identifier (i.e. a function 
    parameter) of type N_Vector and "casts" it to
    a N_Vector_Serial, when the ops field is the 
    managed ops structure.   
 */
#define SERIAL_CAST(X, NV)                            \
  NV = X;                                             \
  value NV##__v = (value) X;                          \
  value NV##__array = Field(NV##__v, 0);              \
  mlsize_t NV##__l = caml_array_length(NV##__array);  \
  struct _N_VectorContent_Serial X##_serial_wrap = {  \
    NV##__l, /* length */                             \
    FALSE, /* own_data */                             \
    (realtype*)NV##__array, /* data */                \
  };                                                  \
                                                      \
  struct _generic_N_Vector X##_wrap = {  	      \
    &X##_serial_wrap,				      \
    NULL                                              \
  };                                                  \
  if (NV -> ops == &managed_ops) {		      \
    NV = &X##_wrap;                                   \
  }                                                   \

#define SERIAL_WRAP(X) SERIAL_CAST(X, X)
#define DECL_SERIAL_CAST(X, NV) N_Vector NV; SERIAL_CAST(X, NV)

static struct _generic_N_Vector_Ops managed_ops = {
  managed_clone,           // N_VClone will return an _unmanaged_ vector
  managed_clone_empty,     // N_VCloneEmpty will return an _unmanaged_ vector
  managed_destroy,         // N_VDestroy is essentially a no-op (since this vector is managed)
  managed_space,           // N_VSpace needs to be a rough estimate only (used for diag)
  managed_get_ptr,
  NULL,                    // N_VSetArrayPointer is not supported for obvious reasons
  managed_linearsum,       // N_VLinearSum
  managed_const,           // N_VConst
  managed_prod,            // N_VProd
  managed_div,             // N_VDiv
  managed_scale,           // N_VScale
  managed_abs             , //N_VAbs
  managed_inv             , //N_VInv
  managed_addconst        , //N_VAddConst
  managed_dotprod         , //N_VDotProd
  managed_max_norm        , //N_VMaxNorm
  managed_wrms_norm       , //N_VWrmsNorm
  managed_wrms_norm_mask  , //N_VWrmsNormMask
  managed_min             , //N_VMin
  managed_wl2_norm        , //N_VWL2Norm
  managed_l1_norm         , //N_VL1Norm
  managed_compare         , //N_VCompare
  managed_inv_test        , //N_VInvTest
  managed_constr_mask     , //N_VConstrMask
  managed_min_quotient    , //N_VMinQuotient
};

void managed_destroy(N_Vector v) { }

N_Vector managed_clone(N_Vector w) {
  value v = (value) w;
  value array = Field(v, 0);
  mlsize_t l = caml_array_length(array);
  N_Vector ret = N_VNew_Serial(l);
  memcpy((void*)array, NV_DATA_S(ret), l);
  //printf("Cloned managed vector with %d elements to %x, ops:%x \n", l, (int)ret, (int) (ret->ops));
  return ret;
}

N_Vector managed_clone_empty(N_Vector w) {
  value v = (value) w;
  value array = Field(v, 0);
  mlsize_t l = caml_array_length(array);
  N_Vector ret = N_VNewEmpty_Serial(l);
  //printf("Created unmanaged vector with %d elements in %x\n", l, (int)ret);
  return ret;
}

void managed_space(N_Vector v, long int *lrw, long int *liw) {
  SERIAL_WRAP(v);
  return N_VSpace_Serial(v, lrw, liw);
}

realtype* managed_get_ptr(N_Vector w) {
  value v = (value) w;
  value array = Field(v, 0);
  return (realtype*)array;
}

void managed_linearsum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z) {
  SERIAL_WRAP(x); SERIAL_WRAP(y); SERIAL_WRAP(z);
  return N_VLinearSum_Serial(a, x, b, y, z);
}

void managed_const(realtype c, N_Vector z) {
  SERIAL_WRAP(z);
  return N_VConst_Serial(c, z);
}

void managed_prod(N_Vector x, N_Vector y, N_Vector z) {
  SERIAL_WRAP(x); SERIAL_WRAP(y); SERIAL_WRAP(z);
  return N_VProd_Serial(x,y,z);
}

void managed_div(N_Vector x, N_Vector y, N_Vector z) {
  SERIAL_WRAP(x); SERIAL_WRAP(y); SERIAL_WRAP(z);
  return N_VDiv_Serial(x,y,z);
}

void managed_scale(realtype c, N_Vector x, N_Vector z) {
  SERIAL_WRAP(x); SERIAL_WRAP(z);
  return N_VScale_Serial(c,x,z);
}

void managed_abs(N_Vector x, N_Vector z) {
  SERIAL_WRAP(x); SERIAL_WRAP(z);
  return N_VAbs_Serial(x,z);
}

void managed_inv(N_Vector x, N_Vector z) {
  SERIAL_WRAP(x); SERIAL_WRAP(z);
  return N_VInv_Serial(x,z);  
}

void managed_addconst(N_Vector x, realtype b, N_Vector z) {
  SERIAL_WRAP(x); SERIAL_WRAP(z);
  return N_VAddConst_Serial(x,b,z);  
}

realtype managed_dotprod(N_Vector x, N_Vector y) {
  SERIAL_WRAP(x); SERIAL_WRAP(y);
  return N_VDotProd_Serial(x,y);  
}

realtype managed_max_norm(N_Vector x) {
  SERIAL_WRAP(x);
  return N_VMaxNorm_Serial(x);  
}

realtype managed_wrms_norm(N_Vector x, N_Vector w) {
  SERIAL_WRAP(x); SERIAL_WRAP(w);
  return N_VWrmsNorm(x, w);  
}

realtype managed_wrms_norm_mask(N_Vector x, N_Vector w, N_Vector id)  {
  SERIAL_WRAP(x); SERIAL_WRAP(w); SERIAL_WRAP(id);
  return N_VWrmsNormMask(x, w, id);  
}

realtype managed_min(N_Vector x)  {
  SERIAL_WRAP(x);
  return N_VMin(x);  
}

realtype managed_wl2_norm(N_Vector x, N_Vector w) {
  SERIAL_WRAP(x); SERIAL_WRAP(w);
  return N_VWL2Norm_Serial(x, w);  
}

realtype managed_l1_norm(N_Vector x) {
  SERIAL_WRAP(x);
  return N_VL1Norm_Serial(x);
}

void managed_compare(realtype c, N_Vector x, N_Vector z) {
  SERIAL_WRAP(x); SERIAL_WRAP(z);
  return N_VCompare(c, x, z);  
}

booleantype managed_inv_test(N_Vector x, N_Vector z) {
  SERIAL_WRAP(x); SERIAL_WRAP(z);
  return N_VInvTest(x, z);  
}

booleantype managed_constr_mask(N_Vector c, N_Vector x, N_Vector m)  {
  SERIAL_WRAP(c); SERIAL_WRAP(x); SERIAL_WRAP(m);
  return N_VConstrMask_Serial(c, x, m);  
}

realtype managed_min_quotient(N_Vector num, N_Vector denom) {
  SERIAL_WRAP(num); SERIAL_WRAP(denom);
  return N_VInvTest(num, denom); 
}

CAMLprim value sundials_ml_array_to_nvector(value array) {
  CAMLparam1 (array);

  CAMLlocal1 (block);
  block = caml_alloc(2, 0);
  //printf("Allocating vector %x from array.\n", block);
  
  Store_field(block, 0, array);
  Store_field(block, 1, &managed_ops); 

  CAMLreturn (block);     
}

CAMLprim value sundials_ml_nvector_to_array(value v) {
  CAMLparam1 (v);
  CAMLreturn (Field(v, 0));     
}

CAMLprim value sundials_ml_nvector_new(value len) { 
  CAMLparam1 (len);
  size_t length = Int_val(len);
  //printf("Allocating vector of %d elements\n", length);

  CAMLlocal1 (array);
  array = caml_alloc(length, Double_array_tag);

  CAMLreturn ( sundials_ml_array_to_nvector(array) );
}

CAMLprim value sundials_ml_nvector_linear_sum (value a, value x, value b, value y) {
  CAMLparam4(a, x, b, y);

  DECL_SERIAL_CAST(x, nv_x); 
  DECL_SERIAL_CAST(y, nv_y);

  size_t len = NV_LENGTH_S(nv_x);
  CAMLlocal1 (z) ;

  z = sundials_ml_nvector_new (Val_int(len)) ;

  DECL_SERIAL_CAST(z, nv_z);

  N_VLinearSum_Serial(Double_val(a), nv_x, Double_val(b), nv_y, nv_z);

  CAMLreturn(z);
}

CAMLprim value sundials_ml_nvector_get(value idx, value nvector) { //not size_t because of nvector API
  CAMLparam2 (idx, nvector);
  if ( ((N_Vector)nvector)->ops == &managed_ops) {    
    CAMLreturn(caml_copy_double(  Double_field(Field(nvector, 0), Int_val(idx))  ));
  } else {
    CAMLreturn(caml_copy_double(  NV_Ith_S( ((N_Vector)nvector) , Int_val(idx))  ));
  }
}
 
CAMLprim value sundials_ml_nvector_set(value idx, value v, value nvector) { //not size_t because of nvector API
  CAMLparam3 (idx, v, nvector);
  if ( ((N_Vector)nvector)->ops == &managed_ops) {
    Store_double_field(Field(nvector, 0), Int_val(idx), Double_val(v) );
  } else {
    NV_Ith_S( ((N_Vector)nvector) , Int_val(idx)) = Double_val(v);
  }

  CAMLreturn (Val_unit);
}

CAMLprim value sundials_ml_nvector_scale(value c, value x, value z) {
  CAMLparam3(c, x, z);
  DECL_SERIAL_CAST(x, nv_x);
  DECL_SERIAL_CAST(z, nv_z);
  
  N_VScale_Serial(Double_val(c), nv_x, nv_z);

  CAMLreturn(Val_unit);
}

CAMLprim value sundials_ml_nvector_const(value c, value z) {
  CAMLparam2(c, z);
  DECL_SERIAL_CAST(z, nv_z);
  
  N_VConst_Serial(Double_val(c), nv_z);

  CAMLreturn(Val_unit);
}

CAMLprim value sundials_ml_nvector_max_norm(value z) {
  CAMLparam1(z);
  DECL_SERIAL_CAST(z, nv_z);
  const double max = N_VMaxNorm_Serial(nv_z);
  CAMLreturn(caml_copy_double(max));
}
