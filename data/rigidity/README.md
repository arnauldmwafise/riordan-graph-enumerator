# data/rigidity — class sizes of random Riordan graphs with f'(0) = 1

Each file records one run: random labeled Riordan graphs `H` of order n with `f'(0) = 1`, drawn
uniformly, together with the exact class size `|C(H)|`. `|C(H)|` is the number of labeled Riordan
graphs isomorphic to `H`, found by enumerating all Riordan labelings of `H`. Searches that exceeded
2·10⁷ nodes were stopped and are reported as "aborted".

Line format: `n samples | R0: … | R1: alone A, rev-pair P, other O, aborted X | max class M | mean nodes K`.
"rev-pair" means `C(H) = {H, rev H}`, where `rev` reverses the vertex order. "alone" means
`|C(H)| = 1`, and "other" means `|C(H)| > 2`.

These files were produced by an earlier stand-alone program. The same experiment is now run with
the sampler, which reports `class_size` and `reversal_pair` for every sample:

    build/rgsample 20 --samples 400 --fprime 1 --class-size --node-limit 20000000 --format csv --seed 11
