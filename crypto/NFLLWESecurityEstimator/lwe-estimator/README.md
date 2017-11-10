# Estimator for the Bit Security of LWE Instances

This [Sage](http://sagemath.org) module provides functions for estimating the bit security of [Learning with Errors](https://en.wikipedia.org/wiki/Learning_with_errors) instances.

The main intend of this estimator is to give designers an easy way to choose parameters resisting known attacks and to enable cryptanalysts to compare their results and ideas with other techniques known in the literature.

## Usage Examples ##

    sage: load("https://bitbucket.org/malb/lwe-estimator/raw/HEAD/estimator.py")
    sage: n, alpha, q = unpack_lwe(Regev(128))
    sage: set_verbose(1)
    sage: costs = estimate_lwe(n, alpha, q, skip="arora-gb")
    sage: print cost_str(costs["dec"])
    bop:   ≈2^56.1,  oracle:   ≈2^14.5,  δ_0: 1.0097372,  bkz2:   ≈2^55.0,  k:        90,  fplll:  ≈2^126.0,  sieve:   ≈2^63.8,  enum:   ≈2^40.2,  enumop:   ≈2^55.3,  ε: 0.0625000

## Try it online ##

You can [try the estimator online](http://aleph.sagemath.org/?z=eJxNjcEKwjAQBe-F_kPoqYXYjZWkKHgQFPyLkOhii6mJyWrx782hiO84MPOcN9e6GohC2gHYkezrckdqfbzBZJwFN-MKE42TIR8hmhnOp8MRfqgNn6opiwdnxoXBcPZke9ZJxZlohRDbXknVSbGMMyXlpi-LhKTfGK1PWK-zr7O1NFHnz_ov2HwBPwsyhw==&lang=sage) using the [Sage Math Cell](http://aleph.sagemath.org/) server. 

## Coverage ##

At present the following algorithms are covered by this estimator.

- exhaustive search + a meet-in-the-middle variant
- Coded-BKW [C:GuoJohSta15]
- lattice-reduction based distinguishing
- lattice-reduction based decoding [RSA:LinPei11]
- solving BDD via Kannan embedding as described in [ICISC:AlbFitGop13]
- using Gröbner bases to solve LWE [EPRINT:ACFP14]

For small secrets, the following variants are covered:

- exhaustive search + a meet-in-the-middle variant for small secrets
- modulus switching in combination with any of the algorithms mentioned above
- small-secret embedding [ACISP:BaiGal14]
- small-secret BKW [PKC:AFFP14]
- small, sparse secret dual lattice attacks [EPRINT:Albrecht17]

Above, we use [cryptobib](http://cryptobib.di.ens.fr)-style bibtex keys as references.

## Evolution ##

This code is evolving, new results are added and bugs are fixed. Hence, estimations from earlier versions might not match current estimations. This is annoying but unavoidable at present. We recommend to also state the commit that was used when referencing this project.

We also encourage authors to let us know if their paper uses this code. In particular, we thrive to tag commits with those cryptobib ePrint references that use it. For example, [this commit](https://bitbucket.org/malb/lwe-estimator/src/6295aa59048daa5d9598378386cb61887a1fe949/?at=EPRINT_Albrecht17) corresponds to this [ePrint entry](https://ia.cr/2017/047).

## Contributions ##

Our intend is for this estimator to be maintained by the research community. For example, we encourage algorithm designers to add their own algorithms to this estimator and we are happy to help with that process.

More generally, all contributions such as bugfixes, documentation and tests are welcome. Please go ahead and submit your pull requests. Also, don’t forget to add yourself to the list of contributors below in your pull requests.

At present, this estimator is maintained by Martin Albrecht. Contributors are:

- Martin Albrecht
- Florian Göpfert
- Cedric Lefebvre
- Rachel Player
- Markus Schmidt
- Sam Scott

Please follow [PEP8](https://www.python.org/dev/peps/pep-0008/) in your submissions. You can use [flake8](http://flake8.pycqa.org/en/latest/) to check for compliance. We use the following flake8 configuration (to allow longer line numbers and more complex functions):

```
[flake8]
max-line-length = 120
max-complexity = 16
ignore = E22,E241
```

## Bugs ##

If you run into a bug, please open an [issue on bitbucket](https://bitbucket.org/malb/lwe-estimator/issues?status=new&status=open). Also, please check there first if the issue has already been reported.


## Citing ##

If you use this estimator in your work, please cite

> Martin R. Albrecht, Rachel Player and Sam Scott.  
> *On the concrete hardness of Learning with Errors*.  
> Journal of Mathematical Cryptology.  
> Volume 9, Issue 3, Pages 169–203,  
> ISSN (Online) 1862-2984, ISSN (Print) 1862-2976  
> DOI: 10.1515/jmc-2015-0016, October 2015

A pre-print is available as

> Cryptology ePrint Archive, Report 2015/046, 2015.
> https://eprint.iacr.org/2015/046

A high-level overview of that work is given, for instance, in this [talk](https://martinralbrecht.files.wordpress.com/2015/05/20150507-lwe-survey-london.pdf).

