open Owl
open Printf
open Lgm

let model = `ppca

module P = struct
  let n = 20

  let d = 10

  let n_samples = 1500

  let w = Mat.(3. $* gaussian n d)

  let sigma =
    match model with
    | `ppca -> Mat.(sqr (gaussian 1 1) * ones n 1)
    | `fa   -> Mat.(sqr (gaussian n 1))


  let xs =
    let z = Mat.gaussian d n_samples in
    Mat.(w *@ z)


  let data =
    let noise = Mat.gaussian n n_samples in
    Mat.(xs + (sqrt sigma * noise))


  let test_zs = Mat.gaussian d n_samples

  let test_nll =
    let test_xs = Mat.(w *@ test_zs) in
    let noise = Mat.gaussian n n_samples in
    let test_data = Mat.(test_xs + (sqrt sigma * noise)) in
    fun w sigma -> negloglik test_data w sigma


  let test_rms_error wp =
    let prediction = Mat.(wp *@ test_zs) in
    let test_xs = Mat.(w *@ test_zs) in
    Mat.(l2norm' (prediction - test_xs)) /. Mat.(l2norm' test_xs)


  let train_rms_error wp zs =
    let prediction = Mat.(wp *@ zs) in
    Mat.(l2norm' (xs - prediction)) /. Mat.(l2norm' xs)
end

let _ =
  let data = P.data in
  let _, w, sigma, zs, train_nll = infer ~model data 10 in
  let train_rms_err = P.train_rms_error w zs in
  let test_nll = P.test_nll w sigma in
  printf "TRAIN negative loglik: %f | TEST NLL: %f\n" train_nll test_nll;
  printf "TRAIN RMS ERROR: %f\n" train_rms_err
