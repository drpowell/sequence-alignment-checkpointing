# Sequence alignment

This repository contains a number of related programs for sequence alignment.  These implementations explore the idea of reducing space complexity by using a check-pointing technique presented in the [papers](#papers) below.

  * Simple edit costs for aligning 2 sequences : [align2str_checkp](align2str_checkp/README)
  * Affine gap costs for 2 sequences : [align2str_linear_checkp](align2str_linear_checkp/README)
  * Affine gap costs for 3 sequences : [align3str_checkp](align3str_checkp/README)

There is also an illustrative implementation of the 2 sequences, affine costs with checkpointing in python : [dpa_lcheckp_2str.py](dpa_lcheckp_2str.py)


### <a name="papers"></a>Papers

    D. R. Powell, L. Allison, T. I. Dix.
    Fast, optimal alignment of three sequences using linear gap costs.
    J. Theor. Biol., Vol.207(3), pp.325-336, 2000. [doi:10.1006/jtbi.2000.2177]

    L. Allison, D. Powell, T. I. Dix.
    Modelling is more versatile than shuffling.
    School of Computer Science and Software Engineering, Monash University, Australia 3800, TR#2000/83, 2000.

    D. R. Powell, L. Allison, T. I. Dix.
    A versatile divide and conquer technique for optimal string alignment.
    Inf. Proc. Lett., Vol.70(3),, pp.127-139, 1999, [doi:10.1016/S0020-0190(99)00053-8].

    L. Allison, D. Powell, T. I. Dix.
    Compression and approximate matching.
    Computer J., Vol.42(1) pp.1-10, 1999.

    D. R. Powell, L. Allison, T. I. Dix, D. L. Dowe.
    Alignment of low information sequences.
    Australasian Computer Science Theory Symposium, CATS'98, pp.215-230, NUS, isbn:9813083921, 1998.

    D. R. Powell, D. L. Dowe, L. Allison, T. I. Dix.
    Discovering simple DNA sequences by compression.
    Pacific Symposium on Biocomputing '98, pp.597-608, 1998.
