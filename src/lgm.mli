open Owl

(** negloglik data w sigma *)
val negloglik : Mat.mat -> Mat.mat -> Mat.mat -> float

(** infer ?verbose=true ?(tol=1E-5) ?(every=10) ?(max_steps=1000) ?(model=`fa) data d *)
val infer
  :  ?verbose:bool
  -> ?tol:float
  -> ?every:int
  -> ?max_steps:int
  -> ?model:[ `fa | `ppca ]
  -> Mat.mat
  -> int
  -> Mat.mat * Mat.mat * Mat.mat * Mat.mat * float
