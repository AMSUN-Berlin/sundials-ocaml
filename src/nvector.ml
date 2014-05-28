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

type nvector = {
	 nv_data : float array ;
	 nv_ops : int
       }

external nvector_create : int -> nvector = "sundials_ml_nvector_new"

external nvector_from_array : float array -> nvector = "sundials_ml_array_to_nvector"

(**
 * create an array from an nvector
 * This will either create a copy (if the nvector is not managed by ocaml) or return
 * the embedded array - thus care needs to be taken when changing the content.
 *)
external array_from_nvector : nvector -> float array = "sundials_ml_nvector_to_array"

external nvector_linear_sum : float -> nvector -> float -> nvector -> nvector = "sundials_ml_nvector_linear_sum"

external nvector_get : int -> nvector -> float = "sundials_ml_nvector_get"

external nvector_set : int -> float -> nvector -> unit = "sundials_ml_nvector_set"

external nvector_scale : float -> nvector -> nvector -> unit = "sundials_ml_nvector_scale"

external nvector_max_norm : nvector -> float = "sundials_ml_nvector_max_norm"

external nvector_const : float -> nvector -> unit = "sundials_ml_nvector_const"
