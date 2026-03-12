Getting Started Tutorial
========================

This tutorial walks through using MKado for McDonald-Kreitman analysis, from running your first test to interpreting results.

Prerequisites
-------------

1. Install MKado (see :doc:`installation`)
2. Have aligned coding sequences in FASTA format

The McDonald-Kreitman Test
--------------------------

The MK test (`McDonald & Kreitman 1991`_) compares the ratio of non-synonymous to synonymous changes within species (polymorphism) versus between species (divergence). Under neutral evolution, these ratios should be equal.

The test produces a 2x2 contingency table:

============  ===============  =============
              Non-synonymous   Synonymous
============  ===============  =============
Divergence    Dn               Ds
Polymorphism  Pn               Ps
============  ===============  =============

Your First MK Test
------------------

Let's run a basic MK test using the example data.

**Step 1: Examine your data**

.. code-block:: bash

   # Get information about an alignment file
   mkado info examples/anopheles_batch/AGAP000078.fa

This shows the number of sequences, sequence lengths, and sequence names.

**Step 2: Run the MK test**

.. code-block:: bash

   # Standard MK test: gamb (ingroup) vs afun (outgroup)
   mkado test examples/anopheles_batch/AGAP000078.fa -i gamb -o afun

**Step 3: Interpret the output**

The output shows:

- **Dn, Ds**: Fixed non-synonymous and synonymous differences between species
- **Pn, Ps**: Polymorphic non-synonymous and synonymous sites within species
- **Fisher's p-value**: Statistical significance of the departure from neutrality
- **Neutrality Index (NI)**: (Pn/Ps) / (Dn/Ds) - ratio of ratios
- **Alpha**: Proportion of substitutions driven by positive selection

Interpreting Results
--------------------

Neutrality Index (NI)
^^^^^^^^^^^^^^^^^^^^^

The Neutrality Index (`Rand & Kann 1996`_) measures the direction and degree of departure from neutral evolution:

- **NI = 1**: Consistent with neutral evolution
- **NI > 1**: Excess polymorphism (segregating weakly deleterious variants)
- **NI < 1**: Excess divergence (positive selection)

Alpha
^^^^^

Alpha, the proportion of adaptive substitutions (`Smith & Eyre-Walker 2002`_):

- **Alpha = 0**: No adaptive substitutions
- **Alpha > 0**: Proportion of fixed differences due to positive selection
- **Alpha < 0**: Excess polymorphism relative to divergence

Asymptotic MK Test
------------------

The standard MK test can be biased by slightly deleterious mutations that segregate as polymorphisms but rarely reach fixation. These inflate Pn relative to Dn, causing alpha to be underestimated. The asymptotic MK test (`Messer & Petrov 2013`_) addresses this by examining how alpha varies with **derived allele frequency** in the ingroup.

The key insight is that weakly deleterious mutations are more common at low frequencies (where they haven't yet been purged by selection) and rare at high frequencies. By extrapolating alpha to a derived frequency of 1.0, we estimate what alpha would be if all deleterious polymorphisms were removed.

Determining Derived Allele Frequency
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To calculate derived allele frequency, we need to identify which allele is **ancestral** (the original state) versus **derived** (the new mutation). MKado uses the outgroup to polarize mutations:

1. At each polymorphic codon position in the ingroup, examine which codons are present
2. Compare with the outgroup codon at the same position
3. The allele shared between ingroup and outgroup is assumed to be **ancestral**
4. **Derived frequency** = 1.0 − (frequency of ancestral allele in ingroup)

For example, if a codon position has allele A at 80% and allele B at 20% in the ingroup, and the outgroup has allele A, then A is ancestral and the derived allele frequency is 0.20.

.. note::

   Sites where no ingroup allele matches the outgroup cannot be polarized and are excluded from the frequency spectrum analysis.

The Asymptotic Method
^^^^^^^^^^^^^^^^^^^^^

Once derived frequencies are computed, the test proceeds as follows:

1. **Bin polymorphisms by derived frequency** — Group Pn and Ps into frequency bins (e.g., 0.0–0.05, 0.05–0.10, etc.)
2. **Calculate alpha at each bin** — α(x) = 1 − (Ds/Dn) × (Pn(x)/Ps(x))
3. **Fit a curve** — Fit a model to the per-bin alpha estimates (see Model Selection below)
4. **Extrapolate to x = 1** — The asymptotic alpha is the fitted value at derived frequency = 1.0

Model Selection
^^^^^^^^^^^^^^^

MKado fits two candidate models to the α(x) curve and automatically selects the best one:

- **Exponential model**: α(x) = a + b·exp(−c·x) — 3 parameters
- **Linear model**: α(x) = a + b·x — 2 parameters

The exponential model captures the expected shape when weakly deleterious mutations cause alpha to increase with frequency. However, when data are noisy or sparse, the exponential fit can be unstable.

Model selection follows the `asymptoticMK <https://github.com/MesserLab/asymptoticMK>`_ R package convention (`Haller & Messer 2017`_):

1. If the exponential fit produces a confidence interval width > 100, it is considered unstable and the **linear model** is used
2. Otherwise, both models are compared using **AIC** (Akaike Information Criterion), which balances fit quality against model complexity
3. The model with the lower AIC is selected
4. If both fits fail, the highest-frequency bin's alpha value is reported as a fallback

The output reports which model was selected (exponential or linear) along with the fitted parameters.

.. code-block:: bash

   # Run asymptotic MK test
   mkado test examples/anopheles_batch/AGAP000078.fa -i gamb -o afun -a

   # Customize number of frequency bins
   mkado test examples/anopheles_batch/AGAP000078.fa -i gamb -o afun -a -b 20

The output includes:

- **Alpha asymptotic**: Extrapolated alpha value at derived frequency = 1
- **95% CI**: Confidence interval for alpha asymptotic
- **Model type**: Whether exponential or linear model was selected
- **Per-bin alpha values**: Alpha estimates at each frequency class

Allele Polarization
-------------------

Polarization is the process of determining which allele is **ancestral** (the original state) versus **derived** (a new mutation). This distinction is important for several reasons:

- Calculating derived allele frequency for the asymptotic MK test
- Filtering polymorphisms by frequency (``--min-freq``, ``--no-singletons``)
- Assigning mutations to specific lineages (polarized MK test)

MKado uses outgroup sequences to polarize alleles. The method differs slightly depending on the test type.

Polarization in Standard MK and Asymptotic Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the standard and asymptotic MK tests, MKado uses a **single outgroup** to identify ancestral alleles:

1. **For each polymorphic codon position** in the ingroup, MKado identifies which codons are present and their frequencies
2. **Compare with the outgroup** — find codons that are shared between the ingroup and outgroup
3. **The shared allele is ancestral** — if multiple ingroup alleles match outgroup alleles, the most frequent one is chosen as ancestral
4. **Calculate derived frequency** — derived_freq = 1.0 − frequency(ancestral)

**Example:**

Consider a codon position where the ingroup has:

- Codon AAA at 75% frequency
- Codon AAG at 25% frequency

If the outgroup has codon AAA, then:

- AAA is the **ancestral** allele (shared with outgroup)
- AAG is the **derived** allele (new mutation in ingroup)
- Derived allele frequency = 0.25

**Unpolarizable sites:**

Some polymorphic sites cannot be polarized:

- **No shared allele**: If no ingroup codon matches any outgroup codon, we cannot determine which is ancestral. This can happen when the outgroup is distantly related or multiple mutations have occurred.
- **Ambiguous data**: Codons containing ambiguous bases (N) or gaps (-) are excluded from comparison.

Unpolarizable sites are excluded from frequency-based analyses (asymptotic test, ``--min-freq`` filtering) but are still counted in the standard MK test's Pn/Ps totals.

Polarized MK Test
-----------------

The standard MK test counts fixed differences between ingroup and outgroup but cannot determine **which lineage** the mutation occurred on. Did the ingroup change from the ancestral state, or did the outgroup?

The polarized MK test uses a **second outgroup** (typically a more distantly related species) to answer this question:

.. code-block:: bash

   # Use a second outgroup for polarization
   mkado test alignment.fa -i ingroup -o outgroup1 --polarize-match outgroup2

Polarizing Fixed Differences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each fixed difference between ingroup and outgroup1, MKado applies a **majority rule** using outgroup2:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Scenario
     - Ancestral State
     - Interpretation
   * - outgroup1 = outgroup2
     - The shared outgroup codon
     - Mutation on **ingroup** lineage
   * - ingroup = outgroup2
     - The ingroup codon
     - Mutation on **outgroup1** lineage
   * - All three differ
     - Cannot determine
     - **Unpolarizable** — excluded from lineage-specific counts

**Example:**

Consider a codon position with:

- Ingroup: GCT (Ala)
- Outgroup1: GGT (Gly)
- Outgroup2: GGT (Gly)

Both outgroups share GGT, so GGT is likely ancestral. The GCT→GGT change occurred on the **ingroup lineage**.

Now consider:

- Ingroup: GCT (Ala)
- Outgroup1: GGT (Gly)
- Outgroup2: GCT (Ala)

The ingroup matches the more distant outgroup2, so GCT is likely ancestral. The GCT→GGT change occurred on the **outgroup1 lineage**.

Polarizing Polymorphisms
^^^^^^^^^^^^^^^^^^^^^^^^

Polymorphisms in the ingroup are also polarized to determine if they arose on the ingroup lineage. The logic mirrors fixed difference polarization:

1. **Check for ancestral polymorphism**: If outgroup2 contains ALL the ingroup alleles, the polymorphism existed at the common ancestor and cannot be attributed to the ingroup lineage → excluded
2. **Determine ancestral state**: If outgroup1 and outgroup2 share an allele, that allele is ancestral
3. **Verify derived allele**: If the ingroup has alleles not present in outgroup2, those are derived on the ingroup lineage → counted as ingroup polymorphism
4. **Cannot polarize**: If no clear ancestral state can be determined → excluded

**Example:**

Consider a polymorphic site where:

- Ingroup has: AAA (60%), AAG (40%)
- Outgroup1 has: AAA
- Outgroup2 has: AAA

Outgroup1 and outgroup2 both have AAA, so AAA is ancestral. The AAG allele (40%) is derived on the ingroup lineage. This is an **ingroup polymorphism**.

Now consider:

- Ingroup has: AAA (60%), AAG (40%)
- Outgroup1 has: AAA
- Outgroup2 has: AAA, AAG

Outgroup2 has BOTH ingroup alleles. This polymorphism existed at the common ancestor — it's not derived on the ingroup lineage. This polymorphism is **excluded** from the ingroup count.

Interpreting Polarized Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The polarized MK test reports separate statistics for each lineage:

- **Ingroup lineage**: Dn, Ds, Pn, Ps, alpha — mutations that occurred on the ingroup branch
- **Outgroup lineage**: Dn, Ds — mutations that occurred on the outgroup1 branch
- **Unpolarized**: Dn, Ds, Pn, Ps — changes that could not be assigned to a lineage

Unpolarized polymorphisms include:

- Ancestral polymorphisms (present in outgroup2)
- Polymorphisms where the ancestral state cannot be determined

Only the **ingroup lineage** statistics are used to calculate alpha. This provides a cleaner estimate of selection on the ingroup lineage by excluding ancestral variation.

Why Use Polarization?
^^^^^^^^^^^^^^^^^^^^^

Lineage-specific analysis is useful when:

- You want to test for selection specifically on the ingroup lineage
- The outgroup1 may have experienced different selective pressures
- You need to distinguish ancestral polymorphism from derived changes

**Phylogenetic requirements:**

::

   outgroup2 ─────────┐
                      ├─── root
   outgroup1 ─────┐   │
                  ├───┘
   ingroup ───────┘

Outgroup2 should be more distantly related than outgroup1. The more distant outgroup2 helps establish the ancestral state at the root.

Frequency Filtering Options
---------------------------

MKado provides two different frequency-related options that serve distinct purposes:

``--min-freq`` (Standard MK, Polarized MK, and α_TG)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This option **excludes** rare polymorphisms below a frequency threshold from the Pn/Ps counts.

- **Applies to**: Standard MK test, Polarized MK test, α_TG estimator
- **Purpose**: Filter out singletons or very rare variants that may be sequencing errors or very recent mutations
- **How it works**: Polymorphisms with derived allele frequency < ``min_freq`` are not counted

.. code-block:: bash

   # Exclude polymorphisms below 5% frequency
   mkado test alignment.fa -i dmel -o dsim --min-freq 0.05

   # Also works with alpha-tg
   mkado batch alignments/ -i dmel -o dsim --alpha-tg --min-freq 0.05

.. note::

   ``--min-freq`` cannot be used with ``--asymptotic``. The asymptotic test uses ``--freq-cutoffs`` for frequency filtering.

``--no-singletons`` (Standard MK, Polarized MK, and α_TG)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A convenience option that automatically sets ``--min-freq`` to exclude singletons (variants appearing only once).

- **Applies to**: Standard MK test, Polarized MK test, α_TG estimator
- **How it works**: Calculates ``1/n`` where n is the sample size, and uses that as the minimum frequency threshold
- **Sample size**: Uses ingroup count, or ingroup + outgroup if ``--pool-polymorphisms`` is enabled

.. code-block:: bash

   # Exclude singletons automatically
   mkado test alignment.fa -i dmel -o dsim --no-singletons

   # For batch processing (threshold calculated per gene)
   mkado batch alignments/ -i dmel -o dsim --no-singletons

   # Also works with alpha-tg
   mkado batch alignments/ -i dmel -o dsim --alpha-tg --no-singletons

.. note::

   ``--no-singletons`` cannot be used with ``--min-freq`` (they are mutually exclusive) or ``--asymptotic``.

``--freq-cutoffs`` (Asymptotic MK test)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This option defines the frequency range used for **curve fitting** in the asymptotic test — it does not exclude data.

- **Applies to**: Aggregated asymptotic MK test (``mkado batch -a``)
- **Default**: ``0.1,0.9``
- **Purpose**: Avoid fitting to extreme frequency bins where data may be sparse or noisy
- **How it works**:

  1. All polymorphisms are counted and binned by derived frequency
  2. Only bins within the ``[low, high]`` range are used to fit the exponential/linear model
  3. The curve is still extrapolated to frequency = 1.0

.. code-block:: bash

   # Fit model using only bins between 15% and 85% frequency
   mkado batch alignments/ -i dmel -o dsim -a --freq-cutoffs 0.15,0.85

Key Differences
^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Aspect
     - ``--min-freq``
     - ``--freq-cutoffs``
   * - Effect
     - **Excludes** polymorphisms from counts
     - **Restricts** which bins are used for fitting
   * - Scope
     - Per-polymorphism filtering
     - Curve-fitting range only
   * - Data
     - Polymorphisms below threshold not counted
     - All polymorphisms counted; fitting uses subset
   * - Analysis types
     - Standard MK, Polarized MK, α_TG
     - Asymptotic MK (aggregated batch mode)

Genetic Code Tables
-------------------

By default, MKado uses the standard genetic code. For analyses involving mitochondrial or non-standard nuclear genomes, use ``--code-table`` to specify an alternate code by name or `NCBI table ID <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>`_:

.. code-block:: bash

   # By name
   mkado test alignment.fa -i species1 -o species2 --code-table vertebrate-mito
   mkado batch alignments/ -i species1 -o species2 --code-table invertebrate-mito

   # By NCBI table ID
   mkado test alignment.fa -i species1 -o species2 --code-table 2

   # List all available codes
   mkado codes

The genetic code affects codon translation, synonymous site counting, and the classification of changes as synonymous or non-synonymous.

Available Codes
^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 10 40 30

   * - ID
     - Name
     - ``--code-table`` alias
   * - 1
     - Standard
     - ``standard``
   * - 2
     - Vertebrate Mitochondrial
     - ``vertebrate-mito``
   * - 3
     - Yeast Mitochondrial
     - ``yeast-mito``
   * - 4
     - Mold / Protozoan / Coelenterate Mitochondrial; Mycoplasma; Spiroplasma
     - ``mold-mito``, ``protozoan-mito``, ``mycoplasma``
   * - 5
     - Invertebrate Mitochondrial
     - ``invertebrate-mito``
   * - 6
     - Ciliate, Dasycladacean and Hexamita Nuclear
     - ``ciliate``
   * - 9
     - Echinoderm and Flatworm Mitochondrial
     - ``echinoderm-mito``, ``flatworm-mito``
   * - 10
     - Euplotid Nuclear
     - ``euplotid``
   * - 11
     - Bacterial, Archaeal and Plant Plastid
     - ``bacterial``, ``archaeal``, ``plant-plastid``
   * - 12
     - Alternative Yeast Nuclear
     - ``alt-yeast``
   * - 13
     - Ascidian Mitochondrial
     - ``ascidian-mito``
   * - 14
     - Alternative Flatworm Mitochondrial
     - ``alt-flatworm-mito``
   * - 16
     - Chlorophycean Mitochondrial
     - ``chlorophycean-mito``
   * - 21
     - Trematode Mitochondrial
     - ``trematode-mito``
   * - 22
     - Scenedesmus obliquus Mitochondrial
     - ``scenedesmus-mito``
   * - 23
     - Thraustochytrium Mitochondrial
     - ``thraustochytrium-mito``
   * - 24
     - Rhabdopleuridae Mitochondrial
     - ``rhabdopleuridae-mito``
   * - 25
     - Candidate Division SR1 and Gracilibacteria
     - ``sr1``, ``gracilibacteria``
   * - 26
     - Pachysolen tannophilus Nuclear
     - ``pachysolen``
   * - 27
     - Karyorelictea Nuclear
     - ``karyorelictea``
   * - 29
     - Mesodinium Nuclear
     - ``mesodinium``
   * - 30
     - Peritrich Nuclear
     - ``peritrich``
   * - 31
     - Blastocrithidia Nuclear
     - ``blastocrithidia``
   * - 33
     - Cephalodiscidae Mitochondrial UAA-Tyr
     - ``cephalodiscidae-mito``

See the full `NCBI genetic code reference <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>`_ for codon assignment details.

Output Formats
--------------

MKado supports multiple output formats:

.. code-block:: bash

   # Pretty-printed output (default)
   mkado test alignment.fa -i species1 -o species2

   # Tab-separated values
   mkado test alignment.fa -i species1 -o species2 -f tsv

   # JSON format
   mkado test alignment.fa -i species1 -o species2 -f json

Next Steps
----------

- Learn about :doc:`batch-workflow` for processing multiple genes
- Explore :doc:`asymptotic` for the full asymptotic MK methodology
- Understand :doc:`alpha-tg` for weighted multi-gene estimates
- Review :doc:`file-formats` for input requirements
- Explore the :doc:`api` for programmatic access

References
----------

.. _Haller & Messer 2017: https://doi.org/10.1534/g3.117.039693
.. _McDonald & Kreitman 1991: https://doi.org/10.1038/351652a0
.. _Messer & Petrov 2013: https://doi.org/10.1073/pnas.1220835110
.. _Rand & Kann 1996: https://doi.org/10.1093/oxfordjournals.molbev.a025634
.. _Smith & Eyre-Walker 2002: https://doi.org/10.1038/4151022a
