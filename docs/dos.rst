Direction of Selection (DoS)
============================

MKado reports the Direction of Selection (DoS) statistic from `Stoletzki & Eyre-Walker (2011)`_.

Background
----------

The Neutrality Index (NI) has limitations as a measure of selection:

- **Asymmetric**: NI ranges from 0 to infinity, making interpretation difficult
- **Unbounded**: Values can become very large with small counts
- **Undefined**: Cannot be calculated when Dn=0, Ds=0, or Ps=0

DoS addresses these issues by providing a symmetric, bounded measure of selection.

The DoS Formula
---------------

.. math::

   DoS = \frac{D_n}{D_n + D_s} - \frac{P_n}{P_n + P_s}

Where:

- D\ :sub:`n` = Non-synonymous divergence (fixed differences)
- D\ :sub:`s` = Synonymous divergence (fixed differences)
- P\ :sub:`n` = Non-synonymous polymorphisms
- P\ :sub:`s` = Synonymous polymorphisms

Interpretation
--------------

DoS measures the difference between the proportion of non-synonymous changes among divergent sites versus polymorphic sites:

- **DoS = 0**: Neutral evolution (equal proportions)
- **DoS > 0**: Positive selection (excess adaptive substitutions)
- **DoS < 0**: Slightly deleterious polymorphisms (excess non-synonymous polymorphism)

Advantages over NI
------------------

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Property
     - NI
     - DoS
   * - Range
     - 0 to ∞
     - -1 to +1
   * - Neutral expectation
     - 1
     - 0
   * - Symmetry
     - Asymmetric
     - Symmetric around 0
   * - Behavior with small counts
     - Can be extreme
     - Well-behaved
   * - Interpretability
     - Harder
     - Easier (direct comparison of proportions)

Output
------

DoS is automatically included in all MK test outputs:

**Pretty format:**

.. code-block:: text

   MK Test Results:
     Divergence:    Dn=6, Ds=8
     Polymorphism:  Pn=1, Ps=8
     Sites:         Ln=576.00, Ls=195.00
     Fisher's exact p-value: 0.176
     Neutrality Index (NI):  0.1667
     Alpha (α):              0.8333
     DoS:                    0.3175
     omega (dN/dS):          0.2539

**JSON format:**

.. code-block:: json

   {
     "dn": 6,
     "ds": 8,
     "pn": 1,
     "ps": 8,
     "p_value": 0.176,
     "ni": 0.1667,
     "alpha": 0.8333,
     "dos": 0.3175,
     "ln": 576.0,
     "ls": 195.0,
     "omega": 0.2539
   }

**TSV format:**

.. code-block:: text

   Dn  Ds  Pn  Ps  p_value   NI        alpha     DoS       Ln          Ls          omega
   6   8   1   8   0.175957  0.166667  0.833333  0.317460  576.000000  195.000000  0.253906

Example Calculation
-------------------

Using the classic Kreitman Adh data (Dn=6, Ds=8, Pn=1, Ps=8):

.. math::

   DoS = \frac{6}{6+8} - \frac{1}{1+8} = \frac{6}{14} - \frac{1}{9} = 0.429 - 0.111 = 0.318

This positive value indicates an excess of non-synonymous divergence relative to polymorphism, consistent with positive selection.

Edge Cases
----------

DoS handles edge cases gracefully:

- **Zero divergence (Dn+Ds=0)**: The divergence ratio is treated as 0
- **Zero polymorphism (Pn+Ps=0)**: The polymorphism ratio is treated as 0
- **All zeros**: Returns ``None`` (no data to analyze)

This behavior makes DoS more robust than NI when counts are sparse.

Reference
---------

.. _Stoletzki & Eyre-Walker (2011): https://doi.org/10.1093/molbev/msq249

Stoletzki N, Eyre-Walker A (2011) Estimation of the Neutrality Index. *Molecular Biology and Evolution* 28(1):63-70. https://doi.org/10.1093/molbev/msq249
