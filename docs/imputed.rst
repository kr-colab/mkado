Imputed MK Test
===============

MKado implements the imputed MK test from `Murga-Moreno et al. (2022)`_, which corrects for slightly deleterious mutations by imputing their count from the synonymous frequency spectrum rather than discarding low-frequency polymorphisms.

Background
----------

The standard MK test is biased by **slightly deleterious mutations** that inflate nonsynonymous polymorphism (Pn) without contributing to divergence (Dn). Several approaches exist to address this:

- **Asymptotic MK test**: Fits a curve across the full frequency spectrum and extrapolates to x=1
- **Standard MK with frequency filtering**: Discards all low-frequency polymorphisms below a cutoff (e.g., ``--min-freq``, ``--no-singletons``)

The imputed MK test takes a different approach: instead of discarding low-frequency data, it **estimates how many** low-frequency nonsynonymous polymorphisms are slightly deleterious, using the synonymous frequency spectrum as a neutral reference.

The Key Insight
---------------

Under neutrality, nonsynonymous and synonymous polymorphisms should have the same ratio of low-frequency to high-frequency variants. An excess of low-frequency nonsynonymous polymorphisms (relative to synonymous) indicates segregating slightly deleterious mutations.

By imputing this excess, the test:

- Retains more data than methods that discard all low-frequency polymorphisms
- Increases statistical power at the gene level
- Decomposes the distribution of fitness effects (DFE) into interpretable fractions

The Imputation Formula
----------------------

Polymorphisms are split at a derived allele frequency (DAF) cutoff (default 15%):

1. Count low-frequency (DAF <= cutoff) and high-frequency (DAF > cutoff) variants separately for nonsynonymous (P\ :sub:`n`) and synonymous (P\ :sub:`s`) classes

2. Compute the neutral ratio from synonymous polymorphisms:

   .. math::

      r = \frac{P_{s,low}}{P_{s,high}}

3. Impute the number of weakly deleterious nonsynonymous polymorphisms:

   .. math::

      P_{wd} = P_{n,low} - P_{n,high} \times r

   This is clamped to >= 0.

4. Compute neutral nonsynonymous polymorphisms:

   .. math::

      P_{n,neutral} = P_n - P_{wd}

5. Calculate corrected alpha:

   .. math::

      \alpha = 1 - \frac{P_{n,neutral}}{P_s} \times \frac{D_s}{D_n}

6. Significance is assessed with Fisher's exact test on the corrected 2x2 table.

DFE Fractions
-------------

When the number of synonymous (m\ :sub:`0`) and nonsynonymous (m\ :sub:`i`) sites are provided, the test decomposes the DFE into four fractions:

- **alpha (a)**: Fraction of adaptive substitutions
- **f**: Fraction of effectively neutral nonsynonymous mutations
- **b**: Fraction of weakly deleterious mutations
- **d**: Fraction of strongly deleterious mutations (d = 1 - f - b)

These fractions sum to 1 and describe the full distribution of fitness effects for new nonsynonymous mutations.

Usage
-----

Single Gene Analysis
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Basic imputed MK test (default 15% DAF cutoff)
   mkado test alignment.fa -i ingroup -o outgroup --imputed

   # With custom DAF cutoff (10%)
   mkado test alignment.fa -i ingroup -o outgroup --imputed --min-freq 0.10

   # Separate ingroup/outgroup files
   mkado test ingroup.fa outgroup.fa --imputed

.. note::

   The ``--min-freq`` option is reused as the DAF cutoff when ``--imputed`` is set. If ``--min-freq`` is not specified, the default cutoff of 0.15 (15%) is used.

Batch Analysis (Aggregated)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For multi-gene analyses, pooling data across genes increases power:

.. code-block:: bash

   # Pool polymorphisms and divergence across all genes
   mkado batch alignments/ -i ingroup -o outgroup --imputed

Batch Analysis (Per-Gene)
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Run imputed test separately for each gene
   mkado batch alignments/ -i ingroup -o outgroup --imputed --per-gene

Bootstrap Confidence Intervals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``--bootstrap N`` is set (with N > 0), the imputed test runs an
additional case-resampling bootstrap to produce a 95% CI on alpha. Each
replicate resamples the polymorphism list with replacement and re-runs
the imputed MK algorithm; the 2.5/97.5 percentiles of the resulting
alpha distribution are reported as ``alpha_CI_low`` / ``alpha_CI_high``.
The omega decomposition CIs are derived from the alpha CI by analytical
scaling (``omega_a_ci = alpha_ci * omega``), mirroring the asymptotic
test.

.. code-block:: bash

   # Imputed test with 500-replicate bootstrap
   mkado test alignment.fa -i ingroup -o outgroup --imputed --bootstrap 500

   # Aggregated imputed batch with bootstrap
   mkado batch alignments/ -i ingroup -o outgroup --imputed --bootstrap 200

The legacy default (``--bootstrap 100``) computes the bootstrap; pass
``--bootstrap 0`` to disable CI computation entirely.

Output
------

The imputed test reports:

- **alpha**: Corrected proportion of adaptive substitutions
- **p_value**: Fisher's exact test on the corrected contingency table
- **Pwd**: Imputed count of weakly deleterious nonsynonymous polymorphisms
- **Pn_neutral**: Nonsynonymous polymorphisms after removing imputed slightly deleterious mutations
- **Dn, Ds, Pn, Ps**: Raw counts
- **cutoff**: The DAF cutoff used
- **Ln, Ls**: Nei-Gojobori non-synonymous and synonymous site totals
- **omega, omega_a, omega_na**: dN/dS and the adaptive/non-adaptive
  decomposition (``omega_a = alpha * omega``) following
  `Gossmann, Keightley & Eyre-Walker 2012`_, applied to MK counts by
  `Coronado-Zamora et al. 2019`_. See :doc:`omega` for the decomposition
  formula.
- **alpha_CI_low, alpha_CI_high**: 95% bootstrap CI on alpha (when
  ``--bootstrap > 0``)
- **omega_a_CI_low/high, omega_na_CI_low/high**: 95% CIs on the omega
  decomposition, derived from the alpha CI scaled by omega
- **ci_method**: ``"bootstrap"`` when CI was computed, ``None`` otherwise

Example output (pretty format):

.. code-block:: text

   Imputed MK Test Results:
     Divergence:    Dn=6, Ds=8
     Polymorphism:  Pn=11, Ps=17
     DAF cutoff:    0.15
     Imputed Pwd:   4.82
     Pn (neutral):  6.18
     Alpha:         0.5154
     p-value:       0.0891
     Alpha 95% CI [bootstrap]: (0.3210, 0.7050)
     Sites:         Ln=573.00, Ls=195.00
     omega:         0.2552 (omega_a=0.1315, omega_na=0.1237)
       omega_a 95% CI:  (0.0819, 0.1799)
       omega_na 95% CI: (0.0753, 0.1733)

Comparison with Other Methods
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Method
     - Approach
     - Best used when
   * - Standard MK
     - No correction
     - Quick assessment; comparing specific genes
   * - Asymptotic MK
     - Curve fitting across frequency spectrum
     - Genome-wide analyses with many polymorphisms
   * - Standard MK with frequency filtering
     - Discards low-frequency polymorphisms (e.g., ``--min-freq``, ``--no-singletons``)
     - Simple correction, but loses data
   * - **Imputed MK**
     - Imputes slightly deleterious count from synonymous spectrum
     - Gene-level analyses; maximizing statistical power

The imputed test is particularly useful when:

- You want gene-level significance (p-values) rather than only genome-wide estimates
- You want to retain as much data as possible
- You want to decompose the DFE into interpretable fractions

When to Use the Imputed MK Test
-------------------------------

**Use imputed MK when:**

- Analyzing individual genes or small gene sets where power matters
- You want per-gene corrected alpha with p-values
- You want DFE decomposition (with site count information)
- The asymptotic test lacks sufficient data for curve fitting

**Consider alternatives when:**

- You have genome-wide data with thousands of polymorphisms (asymptotic MK may be more robust)
- You want frequency-bin visualization of alpha(x) (use asymptotic MK with ``--plot-asymptotic``)
- You need unbiased multi-gene weighting without frequency modeling (use alpha_TG)

Interpreting Results
--------------------

- **α > 0**: Evidence for positive selection (proportion of adaptive substitutions)
- **α ≈ 0**: Consistent with neutral evolution
- **α < 0**: Excess of nonsynonymous polymorphism relative to divergence, suggesting segregating weakly deleterious mutations remain even after correction
- **Pwd > 0**: Low-frequency nonsynonymous excess detected and corrected
- **Pwd = 0**: No evidence for segregating slightly deleterious mutations (or all polymorphisms are high-frequency)

Choosing the DAF Cutoff
-----------------------

The default cutoff of 15% follows the recommendation of `Murga-Moreno et al. (2022)`_. The cutoff defines the boundary between "low-frequency" and "high-frequency" polymorphisms:

- **Lower cutoff** (e.g., 5%): Only the rarest variants are considered low-frequency. More conservative imputation.
- **Higher cutoff** (e.g., 25%): More variants classified as low-frequency. More aggressive imputation.

The optimal cutoff depends on the effective population size and strength of selection in your system.

Reference
---------

.. _Murga-Moreno et al. (2022): https://doi.org/10.1093/g3journal/jkac206
.. _Gossmann, Keightley & Eyre-Walker 2012: https://doi.org/10.1093/gbe/evs027
.. _Coronado-Zamora et al. 2019: https://doi.org/10.1093/gbe/evz046

Murga-Moreno J, Coronado-Zamora M, Casillas S, Barbadilla A (2022) impMKT: the imputed McDonald and Kreitman test, a straightforward correction that significantly increases the evidence of positive selection of the McDonald and Kreitman test at the gene level. *G3: Genes, Genomes, Genetics* 12(10):jkac206. https://doi.org/10.1093/g3journal/jkac206

Gossmann TI, Keightley PD, Eyre-Walker A (2012) The effect of variation in the effective population size on the rate of adaptive molecular evolution in eukaryotes. *Genome Biology and Evolution* 4(5):658-667. https://doi.org/10.1093/gbe/evs027

Coronado-Zamora M, Salvador-Martínez I, Castellano D, Barbadilla A, Salazar-Ciudad I (2019) Adaptation and conservation throughout the *Drosophila melanogaster* life-cycle. *Genome Biology and Evolution* 11(5):1463-1482. https://doi.org/10.1093/gbe/evz046
