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


(** an unknown is a entry of two dimensions: derivation and number*)
type unknown = {
  u_idx : int; 
  u_der : int;

  u_names : string list; (** for pretty printing and results *)
}

let der u = {u with u_der = u.u_der + 1}

(** flat unknowns are indices in one of the two state vectors *)
type flatU = Yy of int (** state or algebraic variable *)
	   | Yp of int (** highest derivative *)

(** a flat layout maps unknowns to a flat representation *)
type flat_layout = {
  u2flat : unknown -> flatU ;
  flat2u : flatU -> unknown ;
}

module Monad = struct   

  type 'a m = (int * unknown list) -> ((int * unknown list) * 'a)

  let return a s = (s, a)

  let bind : 'a m -> ('a -> 'b m) -> 'b m = 
    function am -> function f -> function (n, us) -> 
					  let ((n', us'), a) = am (n, us) in
					  (f a) (n', us')
  
  let make name (n, us) = let u' = {u_idx = n; u_der = 0; u_names = [name]} in
			  ((n+1, u'::us), u')

  let unknowns (n, us) = ((n, us), us)
		       			 
end
