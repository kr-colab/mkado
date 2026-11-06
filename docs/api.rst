Python API Reference
====================

MKado provides a Python API for programmatic access to MK test functionality.

Quick Start
-----------

.. code-block:: python

   from mkado import mk_test, asymptotic_mk_test, SequenceSet

   # Run standard MK test
   result = mk_test("ingroup.fa", "outgroup.fa")
   print(f"Alpha: {result.alpha}")
   print(f"P-value: {result.p_value}")

   # Run asymptotic MK test
   result = asymptotic_mk_test("ingroup.fa", "outgroup.fa")
   print(f"Asymptotic Alpha: {result.alpha_asymptotic}")

Core Classes
------------

SequenceSet
^^^^^^^^^^^

.. automodule:: mkado.core.sequences
   :members:
   :undoc-members:
   :show-inheritance:

GeneticCode
^^^^^^^^^^^

.. automodule:: mkado.core.codons
   :members:
   :undoc-members:
   :show-inheritance:

AlignedPair
^^^^^^^^^^^

.. automodule:: mkado.core.alignment
   :members:
   :undoc-members:
   :show-inheritance:

Analysis Functions
------------------

Standard MK Test
^^^^^^^^^^^^^^^^

.. automodule:: mkado.analysis.mk_test
   :members:
   :undoc-members:
   :show-inheritance:

Asymptotic MK Test
^^^^^^^^^^^^^^^^^^

.. automodule:: mkado.analysis.asymptotic
   :members:
   :undoc-members:
   :show-inheritance:

Imputed MK Test
^^^^^^^^^^^^^^^^

.. automodule:: mkado.analysis.imputed
   :members:
   :undoc-members:
   :show-inheritance:

Polarized MK Test
^^^^^^^^^^^^^^^^^

.. automodule:: mkado.analysis.polarized
   :members:
   :undoc-members:
   :show-inheritance:

Usage Examples
--------------

Working with SequenceSet
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from mkado import SequenceSet

   # Load sequences from FASTA
   seqs = SequenceSet.from_fasta("alignment.fa")

   # Filter by name pattern
   ingroup = seqs.filter_by_name("dmel")
   outgroup = seqs.filter_by_name("dsim")

   # Get sequence information
   print(f"Number of sequences: {len(seqs)}")
   print(f"Alignment length: {seqs.alignment_length}")

Running MK Tests
^^^^^^^^^^^^^^^^

.. code-block:: python

   from mkado import mk_test, SequenceSet

   # From file paths
   result = mk_test("ingroup.fa", "outgroup.fa")

   # From SequenceSet objects
   all_seqs = SequenceSet.from_fasta("combined.fa")
   ingroup = all_seqs.filter_by_name("species1")
   outgroup = all_seqs.filter_by_name("species2")
   result = mk_test(ingroup, outgroup)

   # Access results
   print(f"Dn={result.dn}, Ds={result.ds}")
   print(f"Pn={result.pn}, Ps={result.ps}")
   print(f"P-value: {result.p_value}")
   print(f"NI: {result.ni}")
   print(f"Alpha: {result.alpha}")

A ``SequenceSet`` keeps the reading frame it was built with, for example
``SequenceSet.from_fasta("combined.fa", reading_frame=2)``. All groups in one test must use
the same frame. The ``reading_frame`` argument of ``mk_test`` and the other tests sets the
frame only for FASTA paths. If you give it with ``SequenceSet`` objects in a different frame,
the test raises ``ValueError``.

Asymptotic Analysis
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from mkado import asymptotic_mk_test

   result = asymptotic_mk_test("ingroup.fa", "outgroup.fa", num_bins=20)

   print(f"Asymptotic Alpha: {result.alpha_asymptotic}")
   print(f"95% CI: [{result.ci_low}, {result.ci_high}]")

   # Access per-bin data
   for freq, alpha in zip(result.alpha_x_values, result.alpha_by_freq):
       print(f"Frequency {freq:.2f}: Alpha = {alpha:.3f}")

Batch Processing
^^^^^^^^^^^^^^^^

.. code-block:: python

   from pathlib import Path
   from mkado import mk_test, SequenceSet

   # Process multiple files
   results = {}
   for fasta_file in Path("alignments/").glob("*.fa"):
       seqs = SequenceSet.from_fasta(fasta_file)
       ingroup = seqs.filter_by_name("species1")
       outgroup = seqs.filter_by_name("species2")

       if len(ingroup) > 0 and len(outgroup) > 0:
           results[fasta_file.stem] = mk_test(ingroup, outgroup)

   # Summarize results
   for gene, result in results.items():
       print(f"{gene}: alpha={result.alpha:.3f}, p={result.p_value:.3f}")

Generating Volcano Plots
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from pathlib import Path
   from mkado import mk_test, SequenceSet
   from mkado.io.plotting import create_volcano_plot

   # Collect batch results
   results = []
   for fasta_file in Path("alignments/").glob("*.fa"):
       seqs = SequenceSet.from_fasta(fasta_file)
       ingroup = seqs.filter_by_name("species1")
       outgroup = seqs.filter_by_name("species2")

       if len(ingroup) > 0 and len(outgroup) > 0:
           result = mk_test(ingroup, outgroup)
           results.append((fasta_file.stem, result))

   # Generate volcano plot
   create_volcano_plot(results, Path("volcano.png"))

Generating Asymptotic Plots
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from pathlib import Path
   from mkado import asymptotic_mk_test, SequenceSet
   from mkado.io.plotting import create_asymptotic_plot

   # Run asymptotic MK test
   seqs = SequenceSet.from_fasta("alignment.fa")
   ingroup = seqs.filter_by_name("species1")
   outgroup = seqs.filter_by_name("species2")
   result = asymptotic_mk_test(ingroup, outgroup)

   # Generate alpha(x) vs frequency plot
   create_asymptotic_plot(result, Path("asymptotic_fit.png"))

Plotting Functions
------------------

.. automodule:: mkado.io.plotting
   :members:
   :undoc-members:
   :show-inheritance:
