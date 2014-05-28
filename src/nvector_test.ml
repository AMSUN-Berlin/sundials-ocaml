open Ida
open Nvector

open Core.Std
open Batteries
open Kaputt.Abbreviations

let print_nvector v = IO.to_string (Array.print Float.print) (array_from_nvector v)

let () =
  Printf.printf("Running...\n") ;

  Test.add_simple_test
    ~title:"Identity"
    (fun () -> Assert.equal ~prn:print_nvector (nvector_from_array [| 0. ; 1. ; 2. |]) (nvector_from_array [| 0. ; 1. ; 2. |]) );

  Test.add_simple_test
    ~title:"Linear Sum"
    (fun () -> Assert.equal ~prn:print_nvector (nvector_from_array [| 0. ; 1. ; 2. |]) (
			      nvector_linear_sum 1. ( nvector_from_array [| 0. ; 1. ; 0. |] )
						 2. ( nvector_from_array [| 0. ; 0. ; 1. |] )
			    )
    );


  Test.launch_tests ()
		    
