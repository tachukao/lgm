open Owl

let negloglik data w sigma =
  let n, n_samples = Mat.shape data in
  let xs = Mat.(data - mean ~axis:1 data) in
  let d = Mat.col_num w in
  let inv_gamma = 
    (* use matrix inversion lemma to speed up computations *)
    let inv_sigma = Mat.(1. $/ sigma) in
    Mat.(inv_sigma * ((eye n) - w *@ (inv ((eye d) + (transpose w) *@ (inv_sigma * w))) *@ (transpose w) * (transpose inv_sigma))) in
  let nll = 
    Mat.(sum' (((transpose xs) *@ inv_gamma) * (transpose xs))) /. (float n_samples) 
    -. (Linalg.D.logdet inv_gamma) in 
  nll 

let infer ?(verbose=true) ?(tol=1E-5) ?(every=10) ?(max_steps=1000) ?(model=`fa) data d = 
  let n, n_samples = Mat.shape data in
  let mu = Mat.(mean ~axis:1 data) in
  (* remove mean *)
  let xs = Mat.(data - mu) in
  (* initialise parameters *)
  let w = Mat.gaussian n d in
  let sigma = Mat.ones n 1 in
  let rec iterate nll_ step w sigma =
    let inv_gamma = 
      (* use matrix inversion lemma to speed up computations *)
      let inv_sigma = Mat.(1. $/ sigma) in
      Mat.(inv_sigma * ((eye n) - w *@ (inv ((eye d) + (transpose w) *@ (inv_sigma * w))) *@ (transpose w) * (transpose inv_sigma))) in
    let nll = 
      Mat.(sum' (((transpose xs) *@ inv_gamma) * (transpose xs))) /. (float n_samples) 
      -. (Linalg.D.logdet inv_gamma) in 
    let pct_change = (nll_ -. nll) /. nll_ in
    if step mod every = 0 && verbose then 
      Printf.printf "step %i | nll: %f | pct_change: %f \n%!" step nll pct_change;
    let zs = Mat.((transpose w) *@ inv_gamma *@ xs) in
    if step < max_steps && pct_change > tol then begin
      (* calculate the posterio over latents *)
      let var_zs = 
        zs 
        |> Mat.map_cols (fun z -> Arr.expand Mat.((z *@ transpose z) + (eye d) - (transpose w) *@ inv_gamma *@ w) 3) 
        |> Arr.concatenate ~axis:0 |> Arr.sum ~axis:0 |> Arr.squeeze in
      (* update parameters *)
      let w = 
        let inv_p = Mat.(inv var_zs) in
        let xz = 
          zs
          |> Mat.mapi_cols (fun i z -> Arr.expand Mat.(Mat.(col xs i) *@ transpose z) 3) 
          |> Arr.concatenate ~axis:0 |> Arr.sum ~axis:0 |> Arr.squeeze in
        Mat.(xz *@ inv_p) in
      let sigma = 
        let s = 
          zs
          |> Mat.mapi_cols (fun i z -> 
              let x = Mat.col xs i in
              Arr.expand Mat.((x *@ transpose x) - (w *@ z *@ transpose x)) 3) 
          |> Arr.concatenate ~axis:0 |> Arr.sum ~axis:0 |> Arr.squeeze in 
        let s = Mat.(s /$ float n_samples) in 
        match model with
        | `ppca -> 
          let s = Mat.(trace s) /. float n  in
          Mat.((ones n 1) *$ s)
        | `fa -> Mat.(transpose (diag s)) in
      iterate nll (succ step) w sigma 
    end else mu, w, sigma, zs, nll  in
  iterate (1E10) 0 w sigma 

