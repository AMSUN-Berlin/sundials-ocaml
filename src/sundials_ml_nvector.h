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

#ifndef SUNDIALS_ML_NVECTOR_H
#define SUNDIALS_ML_NVECTOR_H

#include <stdio.h>
#include <string.h>

#include <caml/config.h>

#ifdef ARCH_ALIGN_DOUBLE 
/* In case of some architectures, OCaml uses two words to store each double, since we want
   to pass a float array directly to sundials, we cannot run on these architectures.
*/
#error "ARCH_ALIGN_DOUBLE is set in the OCaml config.h, this is not supported by this library. Sorry."
#endif

//TODO: check if sundials realtype == double

#include <caml/alloc.h>
#include <caml/custom.h>
#include <caml/memory.h>
#include <caml/fail.h>
#include <caml/callback.h>

#include <sundials/sundials_nvector.h>
#include <nvector/nvector_serial.h>

N_Vector managed_clone(N_Vector w);
N_Vector managed_clone_empty(N_Vector w);
void managed_destroy(N_Vector v);
void managed_space(N_Vector v, long int *lrw, long int *liw);
realtype* managed_get_ptr(N_Vector v);
void managed_linearsum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
void managed_const(realtype c, N_Vector z);
void managed_prod(N_Vector x, N_Vector y, N_Vector z);
void managed_div(N_Vector x, N_Vector y, N_Vector z);
void managed_scale(realtype c, N_Vector x, N_Vector z);
void managed_abs(N_Vector x, N_Vector z);
void managed_inv(N_Vector x, N_Vector z);
void managed_addconst(N_Vector x, realtype b, N_Vector z);
realtype managed_dotprod(N_Vector x, N_Vector y);
realtype managed_max_norm(N_Vector x);
realtype managed_wrms_norm(N_Vector x, N_Vector w);
realtype managed_wrms_norm_mask(N_Vector x, N_Vector w, N_Vector id);
realtype managed_min(N_Vector x);
realtype managed_wl2_norm(N_Vector x, N_Vector w);
realtype managed_l1_norm(N_Vector x);
void managed_compare(realtype c, N_Vector x, N_Vector z);
booleantype managed_inv_test(N_Vector x, N_Vector z);
booleantype managed_constr_mask(N_Vector c, N_Vector x, N_Vector m);
realtype managed_min_quotient(N_Vector num, N_Vector denom);

static struct _generic_N_Vector_Ops managed_ops;

#endif
