
open Ida
open Nvector

open Core.Std
open Core_bench.Std

let rec range a b =
  if a > b then []
  else a :: range (a + 1) b

let make_array n = Array.create n 1.

let make_vector n = nvector_create n

let foo n = (make_array n).(n-1)

let vector_create n = Bench.Test.create (fun() -> ignore ( make_vector n )) ~name:(Printf.sprintf "nvector allocation (size %d)" n)

let array_create n = Bench.Test.create (fun() -> ignore ( foo n )) ~name:(Printf.sprintf "array allocation (size %d)" n)

(*let vector_add = Bench.Test.create (fun() -> ignore ( nvector_add_const n )) ~name:(Printf.sprintf "float array to nvector and const addition (size %d)" n)*)

let tests : Bench.Test.t list  = 
  (List.map (range 1 20) vector_create) @ (List.map (range 1 20) array_create)


let () = 
  tests
  |> Bench.make_command 
  |> Command.run

