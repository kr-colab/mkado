Omega Decomposition (ω, ω_a, ω_na)
===================================

MKado reports the rate of nonsynonymous substitution per site (ω), and where
appropriate, its decomposition into adaptive (ω_a) and non-adaptive (ω_na)
components. The decomposition was introduced by
`Gossmann, Keightley & Eyre-Walker (2012)`_ and applied in the
McDonald-Kreitman context — the same form MKado uses — by
`Coronado-Zamora et al. (2019)`_ in their analysis of selection across the
*Drosophila* life cycle.

Background
----------

The classical McDonald-Kreitman summary statistics (α, NI, DoS) are
*proportions* — they describe what fraction of substitutions were adaptive,
but not the *rate* at which adaptation occurs. As `Gossmann et al. (2012)`_
argue, comparing α between species or genomic regions is misleading when the
underlying neutral substitution rate varies. The authors recommend reporting
ω_a — the rate of adaptive substitution per site, scaled by the synonymous
substitution rate — as a more appropriate cross-comparison statistic.
`Coronado-Zamora et al. (2019)`_ later applied the same decomposition
directly to MK-style fixed/polymorphic counts (the form MKado computes),
showing how ω_a and ω_na partition the substitution rate into adaptive and
slightly-deleterious components across functional categories.

The decomposition
-----------------

Following `Coronado-Zamora et al. (2019)`_:

.. math::

   \omega    &= \frac{D_n}{D_s} \cdot \frac{L_s}{L_n} \\
   \omega_a  &= \alpha \cdot \omega \\
   \omega_{na} &= (1 - \alpha) \cdot \omega = \omega - \omega_a

Where:

- :math:`D_n, D_s` — nonsynonymous and synonymous fixed differences
- :math:`L_n, L_s` — nonsynonymous and synonymous site totals (Nei-Gojobori
  averaged across analyzed codons in both groups)
- :math:`\alpha` — proportion of adaptive substitutions

ω = (Dn × Ls) / (Ds × Ln) is the standard dN/dS ratio with Nei-Gojobori
site weighting; ω_a is the adaptive-substitution rate per synonymous site;
ω_na is its complement.

Site counting
-------------

MKado computes :math:`L_n` and :math:`L_s` with the **Nei-Gojobori (1986)**
method: for each clean codon, count the fraction of single-nucleotide
substitutions at each site that would be synonymous, sum across the three
codon positions, and average across all clean codons in both ingroup and
outgroup at each analyzed position.

.. note::

   `Gossmann et al. (2012)`_ used the F3×4 codon-frequency model implemented in
   PAML, not Nei-Gojobori. The two methods give numerically different site
   counts, so MKado's ω, ω_a, and ω_na are **not strictly comparable** to
   numbers reported under the F3×4 convention. The decomposition itself
   (α × ω) is identical; the difference is in how :math:`L_n` and :math:`L_s`
   are estimated.

.. _counting-differences:

Counting differences
--------------------

:math:`D_n` and :math:`D_s` are classified one codon pair at a time, as are
:math:`P_n` and :math:`P_s` for FASTA input. A pair differing at a single
position is unambiguous. A pair differing at two or three positions can be
reached by several orderings of the same length, each passing through different
intermediate codons, so each carries its own split of synonymous and
replacement steps.

MKado discards orderings that pass through a stop codon and keeps the one with
the **fewest replacements**. For AAA to CGC, the six orderings carry 1, 2, 3, 3,
3 and 3 replacements, and MKado reports 1.

VCF input does not use this rule for polymorphism. Each variant is classified on
its own against the reference codon, so a codon carrying two SNPs contributes
two independent single-position changes.

.. note::

   This is a parsimony rule, not the `Nei & Gojobori (1986)`_ treatment of
   differences, which averages the counts over the surviving orderings. The two
   agree whenever a pair differs at one position, and often otherwise: under the
   standard code they agree for 79% of pairs differing at two positions and 41%
   of pairs differing at three. Where they differ, parsimony reports fewer
   replacements, never more, so :math:`D_n` runs low and :math:`D_s` high
   relative to the averaged convention, which lowers ω and biases
   :math:`\alpha` downward. The mean shortfall is 0.11 replacements per
   two-position pair and 0.30 per three-position pair.

   MKado's :math:`L_n` and :math:`L_s` **are** Nei-Gojobori, as described above.
   The two conventions meet only in the site counts.

.. warning::

   When every ordering is blocked by a stop codon the pair is dropped from the
   counts rather than classified. Under the vertebrate mitochondrial code this
   silently discards four ordinary sense-codon pairs, AAA and AAG against TGA
   and TGG, all of them Lys against Trp replacements.

When ω_a / ω_na are reported (and when they are not)
-----------------------------------------------------

ω is well-defined and well-behaved per gene: it is the standard dN/dS ratio,
routinely reported per-gene by tools like PAML's codeml. MKado therefore
reports ω on every result type.

ω_a and ω_na, by contrast, **inherit the variance of α**. Per-gene Smith &
Eyre-Walker α is known to be noisy on the small Dn/Ds/Pn/Ps counts typical
of a single gene; multiplying by ω propagates that noise into ω_a. Following
the convention of `Gossmann et al. (2012)`_, `Galtier (2016)`_, and
`Coronado-Zamora et al. (2019)`_, MKado therefore restricts ω_a / ω_na to
result types whose α estimator is appropriate for the decomposition:

.. list-table::
   :header-rows: 1
   :widths: 35 15 15 35

   * - Result type
     - ω
     - ω_a / ω_na
     - α estimator feeding ω_a
   * - ``MKResult`` (per-gene MK)
     - ✓
     - ✗
     - per-gene Smith-Eyre-Walker α
   * - ``PolarizedMKResult`` (per-gene)
     - ✓
     - ✗
     - per-gene polarized α
   * - ``AsymptoticMKResult`` (per-gene)
     - ✓
     - ✓
     - asymptotic α (`Messer & Petrov 2013`_)
   * - ``AsymptoticMKResult`` (aggregated)
     - ✓
     - ✓
     - asymptotic α — **canonical Gossmann-style ω_a**
   * - ``ImputedMKResult``
     - ✓
     - ✓
     - imputed α (`Murga-Moreno et al. 2022`_)
   * - ``AlphaTGResult``
     - ✓
     - ✓
     - Tarone-Greenland weighted α

Rationale: per-gene Smith-Eyre-Walker α is noisy because Dn, Ds, Pn, Ps
counts on a single gene are typically small. Multiplying that noisy α by ω
produces a ω_a whose variance is dominated by α-variance, not by any
biological signal. The asymptotic, Tarone-Greenland, and imputed estimators
were designed for use with sparse counts and are explicitly intended for
this kind of decomposition.

If you need a per-gene ω_a estimate with proper uncertainty quantification,
consider the SnIPRE Bayesian hierarchical estimator (`Eilertson et al. 2012`_)
— not currently implemented in MKado.

Confidence intervals
--------------------

For result types that already perform a sampling procedure to estimate α,
MKado also reports 95% CIs on ω_a and ω_na (and on ω itself where it has a
sampling distribution):

.. list-table::
   :header-rows: 1
   :widths: 30 20 25 25

   * - Result type
     - ω CI
     - ω_a CI
     - ω_na CI
   * - ``AsymptoticMKResult``
     - n/a (constant)
     - ``α_CI × ω``
     - ``(1 − α_CI) × ω`` (flipped)
   * - ``AlphaTGResult``
     - bootstrap
     - bootstrap
     - bootstrap

For ``AsymptoticMKResult`` the asymptotic-fit Monte Carlo (or exponential-fit
bootstrap) only resamples polymorphisms — Dn, Ds, Ln, Ls are constants of the
alignment. ω therefore has no sampling distribution and is reported as a
point estimate. The CIs on ω_a and ω_na are derived analytically by scaling
the α CI (``ω_a_ci_low = α_ci_low × ω``, etc.); the ω_na percentiles flip
because ``(1 − α)`` is monotonically decreasing in α.

For ``AlphaTGResult`` the bootstrap resamples genes, so Dn, Ds, Ln, Ls all
vary per replicate and ω has its own bootstrap distribution. ω, ω_a, and
ω_na are recomputed per replicate and 2.5%/97.5% percentiles taken.

ImputedMKResult does not currently report CIs on ω_a / ω_na — the imputed
test produces a single point estimate. If you need CIs there, consider
running the imputed test on bootstrap replicates of the input alignment.

Edge cases
----------

ω, ω_a, and ω_na return ``None`` (rendered as ``NA`` in TSV, ``null`` in JSON)
when:

- :math:`D_s = 0` — denominator undefined
- :math:`L_s = 0` or :math:`L_n = 0` — site counts undefined
- :math:`L_n` / :math:`L_s` were not computed (e.g. legacy
  ``mk_test_from_counts()`` calls without the ``ln`` / ``ls`` arguments)

Additionally, ω_a and ω_na are ``None`` when α is ``None`` — typically when
:math:`D_n = 0` or :math:`P_s = 0`.

References
----------

.. _Gossmann, Keightley & Eyre-Walker (2012): https://doi.org/10.1093/gbe/evs027
.. _Gossmann et al. (2012): https://doi.org/10.1093/gbe/evs027
.. _Nei & Gojobori (1986): https://doi.org/10.1093/oxfordjournals.molbev.a040410
.. _Galtier (2016): https://doi.org/10.1371/journal.pgen.1005774
.. _Coronado-Zamora et al. (2019): https://doi.org/10.1093/gbe/evz046
.. _Messer & Petrov 2013: https://doi.org/10.1073/pnas.1220835110
.. _Murga-Moreno et al. 2022: https://doi.org/10.1093/g3journal/jkac206
.. _Eilertson et al. 2012: https://doi.org/10.1371/journal.pgen.1002806

Gossmann TI, Keightley PD, Eyre-Walker A (2012). The effect of variation in
the effective population size on the rate of adaptive molecular evolution
in eukaryotes. *Genome Biology and Evolution* 4(5):658-667.
https://doi.org/10.1093/gbe/evs027

Coronado-Zamora M, Salvador-Martínez I, Castellano D, Barbadilla A, Salazar-Ciudad I (2019).
Adaptation and conservation throughout the *Drosophila melanogaster* life-cycle.
*Genome Biology and Evolution* 11(5):1463-1482.
https://doi.org/10.1093/gbe/evz046

Galtier N (2016). Adaptive protein evolution in animals and the effective
population size hypothesis. *PLoS Genetics* 12(1):e1005774.
https://doi.org/10.1371/journal.pgen.1005774

Nei M, Gojobori T (1986). Simple methods for estimating the numbers of
synonymous and nonsynonymous nucleotide substitutions. *Molecular Biology
and Evolution* 3(5):418-426.
