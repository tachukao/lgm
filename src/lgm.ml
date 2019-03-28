open Owl

let negloglik data w sigma =
  let n, n_samples = Mat.shape data in
  let xs = Mat.(data - mean ~axis:1 data) in
  let d = Mat.col_num w in
  let inv_gamma =
    (* use matrix inversion lemma to speed up computations *)
    let inv_sigma = Mat.(1. $/ sigma) in
    Mat.(
      inv_sigma
      * ( eye n
        - w
          *@ inv (eye d + (transpose w *@ (inv_sigma * w)))
          *@ transpose w * transpose inv_sigma ))
  in
  let nll =
    (Mat.(sum' (transpose xs *@ inv_gamma * transpose xs)) /. float n_samples)
    -. Linalg.D.logdet inv_gamma
  in
  nll

let infer ?(verbose = true) ?(tol = 1E-5) ?(every = 10) ?(max_steps = 1000)
    ?(model = `fa) data d =
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
      Mat.(
        inv_sigma
        * ( eye n
          - w
            *@ inv (eye d + (transpose w *@ (inv_sigma * w)))
            *@ transpose w * transpose inv_sigma ))
    in
    let nll =
      (Mat.(sum' (transpose xs *@ inv_gamma * transpose xs)) /. float n_samples)
      -. Linalg.D.logdet inv_gamma
    in
    let pct_change = (nll_ -. nll) /. nll_ in
    if step mod every = 0 && verbose then
      Printf.printf "step %i | nll: %f | pct_change: %f \n%!" step nll
        pct_change ;
    (* calculate the posterio over latents *)
    let zs = Mat.(transpose w *@ inv_gamma *@ xs) in
    if step < max_steps && pct_change > tol then
      let zz, xz, xx =
        let rec loop_add i zz xz xx =
          if i < n_samples then
            let z = Mat.col zs i in
            let x = Mat.col xs i in
            let zz = Mat.(zz + (z *@ transpose z)) in
            let xz = Mat.(xz + (x *@ transpose z)) in
            let xx = Mat.(xx + (x *@ transpose x)) in
            loop_add (succ i) zz xz xx
          else (zz, xz, xx)
        in
        loop_add 0 Mat.(zeros d d) Mat.(zeros n d) Mat.(zeros n n)
      in
      let var_zs =
        Mat.(
          zz + ((eye d - (transpose w *@ inv_gamma *@ w)) *$ float n_samples))
      in
      (* update parameters *)
      let w =
        let inv_p = Mat.(inv var_zs) in
        Mat.(xz *@ inv_p)
      in
      let sigma =
        let s = Mat.(xx - (w *@ (transpose xz))) in
        let s = Mat.(s /$ float n_samples) in
        match model with
        | `ppca ->
            let s = Mat.(trace s) /. float n in
            Mat.(ones n 1 *$ s)
        | `fa -> Mat.(transpose (diag s))
      in
      iterate nll (succ step) w sigma
    else (mu, w, sigma, zs, nll)
  in
  iterate 1E10 0 w sigma
