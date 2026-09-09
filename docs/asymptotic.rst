Asymptotic MK Test
==================

MKado implements the asymptotic MK test from `Messer & Petrov (2013)`_, which corrects for the bias caused by slightly deleterious mutations that inflate polymorphism counts.

Background
----------

The standard MK test assumes that all non-synonymous polymorphisms are either neutral or adaptive. In reality, many non-synonymous polymorphisms are **weakly deleterious** — harmful enough to eventually be removed by selection but not so harmful that they're immediately eliminated.

These weakly deleterious mutations:

- Contribute to Pn (non-synonymous polymorphism)
- Rarely reach fixation, so they don't contribute to Dn
- Inflate the Pn/Dn ratio
- Cause alpha to be **underestimated** (often yielding negative values)

The asymptotic MK test addresses this bias by examining how alpha varies with **derived allele frequency** and extrapolating to high frequencies where deleterious variants have been purged.

The Key Insight
---------------

Weakly deleterious mutations have a characteristic frequency signature:

- **Common at low frequencies**: Selection hasn't had time to remove them
- **Rare at high frequencies**: Selection has purged most of them

By calculating alpha separately at different frequency classes, we can observe this pattern and extrapolate to what alpha would be if all deleterious polymorphisms were removed.

The α(x) Formula
----------------

At each derived allele frequency bin *x*, alpha is calculated as:

.. math::

   \alpha(x) = 1 - \frac{D_s}{D_n} \times \frac{P_n(x)}{P_s(x)}

Where:

- D\ :sub:`n`, D\ :sub:`s` = total nonsynonymous and synonymous divergence
- P\ :sub:`n`\ (x), P\ :sub:`s`\ (x) = nonsynonymous and synonymous polymorphisms in frequency bin *x*

As *x* approaches 1.0, deleterious polymorphisms are increasingly purged, and α(x) converges to the true proportion of adaptive substitutions.

SFS Construction: At-x vs Above-x
---------------------------------

MKado supports two definitions of P\ :sub:`n`\ (x) and P\ :sub:`s`\ (x) — the per-bin form from `Messer & Petrov 2013`_ and the cumulative form introduced by `Uricchio et al. 2019`_. They share the asymptote at *x* = 1 but differ in finite-sample stability.

**At-x mode** (``--sfs-mode at``, default):

.. math::

   P_n(x), P_s(x) = \text{count of polymorphisms in bin centered at } x

This is the per-bin SFS. At large sample sizes the high-frequency bins become sparse, which inflates the per-bin α(x) noise and destabilizes the curve fit.

**Above-x mode** (``--sfs-mode above``):

.. math::

   P_n(x), P_s(x) = \text{count of polymorphisms with derived frequency } \ge x

This is the inclusive right-tail cumulative SFS. The two modes share the asymptote at *x* = 1 (where both definitions go to zero), but the cumulative form averages out per-bin noise and is more stable as sample size grows.

The default is ``at`` for backward compatibility. The ``above`` form follows the convention used by `MKtest.jl <https://github.com/jmurga/MKtest.jl>`_ and the analyses in Uricchio et al. 2019; pass ``--sfs-mode above`` to switch.

.. code-block:: bash

   mkado batch alignments/ -i ingroup -o outgroup -a --sfs-mode above

Determining Derived Allele Frequency
------------------------------------

To calculate derived frequency, we need to distinguish **ancestral** (original) from **derived** (new mutation) alleles. MKado uses the outgroup:

1. At each polymorphic codon in the ingroup, identify which alleles are present
2. Compare with the outgroup codon at the same position
3. The allele shared between ingroup and outgroup is **ancestral**
4. **Derived frequency** = 1.0 − frequency(ancestral allele)

**Example:**

If a codon position has allele A at 80% and allele B at 20% in the ingroup, and the outgroup has allele A:

- A is ancestral (shared with outgroup)
- B is derived (new mutation)
- Derived allele frequency = 0.20

.. note::

   Sites where no ingroup allele matches the outgroup cannot be polarized and are excluded from frequency spectrum analysis.

Model Fitting and Selection
---------------------------

MKado fits two candidate models to the α(x) curve:

**Exponential model** (3 parameters):

.. math::

   \alpha(x) = a + b \cdot e^{-c \cdot x}

**Linear model** (2 parameters):

.. math::

   \alpha(x) = a + b \cdot x

The exponential model captures the expected shape when weakly deleterious mutations cause alpha to increase with frequency. The linear model provides a robust fallback when data are sparse or noisy.

**Model selection** follows the `asymptoticMK R package <https://github.com/MesserLab/asymptoticMK>`_ convention (`Haller & Messer 2017`_):

1. If exponential fit produces a confidence interval width > 100, use the **linear model** (exponential is unstable)
2. Otherwise, compare both models using **AIC** (Akaike Information Criterion)
3. Select the model with lower AIC
4. If both fits fail, report the highest-frequency bin's alpha as a fallback

Usage
-----

Single Gene Analysis
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Basic asymptotic MK test
   mkado test alignment.fa -i ingroup -o outgroup -a

   # With 20 frequency bins (default is 10)
   mkado test alignment.fa -i ingroup -o outgroup -a -b 20

   # Generate alpha(x) plot
   mkado test alignment.fa -i ingroup -o outgroup -a --plot-asymptotic alpha.png

Batch Analysis (Aggregated)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For multi-gene analyses, aggregating data across genes provides more statistical power:

.. code-block:: bash

   # Pool polymorphisms across all genes
   mkado batch alignments/ -i ingroup -o outgroup -a

   # With custom frequency bins and plot
   mkado batch alignments/ -i ingroup -o outgroup -a -b 20 --plot-asymptotic alpha.png

   # Control which frequency bins are used for curve fitting
   mkado batch alignments/ -i ingroup -o outgroup -a --freq-cutoffs 0.15,0.85

Batch Analysis (Per-Gene)
^^^^^^^^^^^^^^^^^^^^^^^^^

Each gene fits its own curve, but ``--freq-cutoffs`` still restricts
which bins go into that per-gene fit, the same as in aggregated mode.

.. code-block:: bash

   # Separate asymptotic test for each gene
   mkado batch alignments/ -i ingroup -o outgroup -a --per-gene

   # Per-gene fits, restricted to bins between 20% and 80% frequency
   mkado batch alignments/ -i ingroup -o outgroup -a --per-gene --freq-cutoffs 0.2,0.8

Confidence Interval Method
^^^^^^^^^^^^^^^^^^^^^^^^^^

For aggregated analyses, MKado offers two CI methods via ``--ci-method``:

- **monte-carlo** (default): samples curve-fit parameters from a
  multivariate normal of the fit covariance matrix. Fast (the model is
  evaluated, not refit, per draw) but parametric — assumes the fit
  uncertainty is well-described by its covariance.
- **bootstrap**: case-resampling of the pooled polymorphism list with
  replacement; the curve is refit each replicate. More principled for
  small per-bin counts (the typical MK setting); slower because each
  replicate refits.

.. code-block:: bash

   # Default: parametric Monte Carlo CI
   mkado batch alignments/ -i ingroup -o outgroup -a

   # Bootstrap CI (refit per replicate)
   mkado batch alignments/ -i ingroup -o outgroup -a --ci-method bootstrap

   # Bootstrap with more replicates
   mkado batch alignments/ -i ingroup -o outgroup -a --ci-method bootstrap --bootstrap 500

The bootstrap is held to whichever model (exponential or linear) the
point estimate selected via AIC, so the CI is comparable across methods.

The flag is silently accepted but has no effect on per-gene asymptotic
runs (those have always used a polymorphism-list bootstrap natively).

Output
------

The asymptotic test reports:

- **alpha_asymptotic**: Extrapolated alpha at derived frequency = 1.0
- **CI_low, CI_high**: 95% confidence interval (method recorded in ``ci_method``)
- **ci_method**: ``"monte-carlo"`` (default) or ``"bootstrap"`` (case-resampling)
- **model_type**: Selected model (exponential or linear)
- **a, b, c**: Fitted model parameters (c only for exponential)
- **Ln, Ls**: Nei-Gojobori non-synonymous and synonymous site totals
- **omega**: dN/dS ratio ``(Dn/Ds) * (Ls/Ln)``
- **omega_a, omega_na**: Adaptive and non-adaptive substitution rates
  (`Gossmann, Keightley & Eyre-Walker 2012`_; applied to MK counts by
  `Coronado-Zamora et al. 2019`_)
- **omega_a_CI_low/high, omega_na_CI_low/high**: 95% CIs derived analytically by scaling the alpha CI by omega (omega_na percentiles flip since ``(1 - alpha)`` is monotonically decreasing). See :doc:`omega` for the rationale.

Example output (pretty format):

.. code-block:: text

   Asymptotic MK Test Results:
     Asymptotic α: 0.5997 (95% CI [monte-carlo]: 0.5162 - 0.6853)
     Divergence: Dn=18828, Ds=49857
     SFS mode: at
     Polymorphism: Pn=8150, Ps=26641
     Sites: Ln=473618.29, Ls=148023.71
     omega: 0.1180 (omega_a=0.0708, omega_na=0.0472)
       omega_a 95% CI:  (0.0609, 0.0809)
       omega_na 95% CI: (0.0371, 0.0571)
     Genes aggregated: 400
     Fit (linear): α(x) = 0.4336 + 0.1661 * x

The Alpha(x) Plot
-----------------

The ``--plot-asymptotic`` option generates a visualization showing:

- **Scatter points**: Observed α values at each frequency bin
- **Fitted curve**: Exponential or linear fit
- **Horizontal band**: Asymptotic α estimate with 95% CI
- **X-axis**: Derived allele frequency
- **Y-axis**: α(x)

This plot helps assess:

- Whether α(x) increases with frequency (expected pattern)
- The quality of the curve fit
- How much deleterious polymorphism affects lower frequencies

Frequency Cutoffs
-----------------

The ``--freq-cutoffs`` option controls which frequency bins are used for curve fitting.
It applies the same way to ``mkado test -a``, aggregated ``mkado batch -a``, and
per-gene ``mkado batch -a --per-gene`` (and the ``mkado vcf`` equivalents):

.. code-block:: bash

   # Fit using bins between 10% and 90% frequency (default)
   mkado batch alignments/ -i sp1 -o sp2 -a --freq-cutoffs 0.1,0.9

   # Stricter cutoffs for noisy data
   mkado batch alignments/ -i sp1 -o sp2 -a --freq-cutoffs 0.2,0.8

   # Same option on a single alignment
   mkado test alignment.fa -i sp1 -o sp2 -a --freq-cutoffs 0.2,0.8

This does **not** exclude polymorphisms from the total counts — it only affects which bins inform the curve fit. Extreme frequency bins often have sparse data and can destabilize the fit.

Comparison with Other Methods
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Method
     - Corrects for
     - Best used when
   * - Standard α
     - Nothing
     - Quick assessment; comparing specific genes
   * - Imputed MK
     - Weakly deleterious mutations (by imputation)
     - Gene-level analyses; maximizing power with limited data
   * - α_TG
     - Sample size heterogeneity across genes
     - Multi-gene analysis with minimal deleterious load
   * - **Asymptotic α**
     - Weakly deleterious mutations
     - Most genome-wide analyses (recommended default)

**Typical pattern** across methods:

- Standard α is often **negative** due to excess Pn from deleterious mutations
- α_TG is **less biased** by gene sample size but still affected by deleterious load
- Asymptotic α is typically **higher and positive**, revealing adaptive substitutions

When to Use Asymptotic α
------------------------

**Use asymptotic α when:**

- Analyzing genome-wide or multi-gene data
- Weakly deleterious mutations are a concern (most analyses)
- You need the most accurate estimate of adaptive substitution rate
- You have sufficient polymorphism data (hundreds of segregating sites)

**Consider alternatives when:**

- Analyzing a single gene with few polymorphisms (standard test may be more appropriate)
- Polymorphism data is too sparse for reliable frequency binning
- You specifically want an unbiased multi-gene estimate without frequency modeling (use α_TG)

Interpreting Results
--------------------

- **α ≈ 0**: Little evidence for adaptive evolution; most substitutions are neutral
- **α > 0**: Proportion of substitutions driven by positive selection
- **α < 0**: Unusual; may indicate model misspecification or very strong deleterious load

For typical Drosophila genome-wide analyses, asymptotic α is often 0.4–0.6, indicating that 40–60% of amino acid substitutions were adaptive.

References
----------

.. _Haller & Messer 2017: https://doi.org/10.1534/g3.117.039693
.. _Messer & Petrov 2013: https://doi.org/10.1073/pnas.1220835110
.. _Messer & Petrov (2013): https://doi.org/10.1073/pnas.1220835110
.. _Uricchio et al. 2019: https://doi.org/10.1038/s41559-019-0890-6
.. _Gossmann, Keightley & Eyre-Walker 2012: https://doi.org/10.1093/gbe/evs027
.. _Coronado-Zamora et al. 2019: https://doi.org/10.1093/gbe/evz046

Haller BC, Messer PW (2017) asymptoticMK: A web-based tool for the asymptotic McDonald–Kreitman test. *G3: Genes, Genomes, Genetics* 7(5):1569-1575. https://doi.org/10.1534/g3.117.039693

Messer PW, Petrov DA (2013) Frequent adaptation and the McDonald–Kreitman test. *PNAS* 110(21):8615-8620. https://doi.org/10.1073/pnas.1220835110

Uricchio LH, Petrov DA, Enard D (2019) Exploiting selection at linked sites to infer the rate and strength of adaptation. *Nature Ecology & Evolution* 3:977-984. https://doi.org/10.1038/s41559-019-0890-6

Gossmann TI, Keightley PD, Eyre-Walker A (2012) The effect of variation in the effective population size on the rate of adaptive molecular evolution in eukaryotes. *Genome Biology and Evolution* 4(5):658-667. https://doi.org/10.1093/gbe/evs027

Coronado-Zamora M, Salvador-Martínez I, Castellano D, Barbadilla A, Salazar-Ciudad I (2019) Adaptation and conservation throughout the *Drosophila melanogaster* life-cycle. *Genome Biology and Evolution* 11(5):1463-1482. https://doi.org/10.1093/gbe/evz046
