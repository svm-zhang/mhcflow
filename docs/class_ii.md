## Extend to Class II typing

The original `polysolver` algorithm has been well-known for genotyping Class
I alleles. However, in theory it should also be able to apply to the Class II
case, with certain modification as well as a set of Class II references.

To type the Class II alleles, you only need to swap in the new reference
data, and the CLI is the same as we have shown for the Class I case.

I have done some preliminary benchmark on Class II typing using some samples from
1000 genome project. The result is suprisingly not too shady and can be found
[here](https://github.com/svm-zhang/hla_benchmark).

You can also find all Class II-related reference data within the `reference`
folder in this repo.
