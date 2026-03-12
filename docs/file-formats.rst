File Formats and Input Requirements
====================================

MKado accepts aligned coding sequences in FASTA format. This page describes the input requirements and common pitfalls.

FASTA Format
------------

MKado reads standard FASTA format:

.. code-block:: text

   >sequence_name1 optional description
   ATGCATGCATGC...
   >sequence_name2
   ATGCATGCATGC...

Requirements:

- Standard FASTA format with ``>`` headers
- Sequences can span multiple lines
- Any valid sequence characters (ACGT, gaps, ambiguous codes)

Alignment Requirements
----------------------

Sequences must be pre-aligned:

- All sequences must be the same length
- Gaps should be represented as ``-``
- Alignment should be in-frame (codon-aligned)

Codon Alignment
^^^^^^^^^^^^^^^

MKado analyzes coding sequences at the codon level. Your alignment should:

- Start at the first position of a codon
- Be a multiple of 3 nucleotides (excluding gaps)
- Maintain reading frame throughout

Example of proper codon alignment:

.. code-block:: text

   >seq1
   ATGGCC---TAAACT
   >seq2
   ATGGCCTGATAAACT

If your alignment isn't codon-aligned, use the ``-r`` option to specify the reading frame (1, 2, or 3).

Sequence Naming Conventions
---------------------------

For combined file mode, sequence names should contain identifiable patterns for filtering:

.. code-block:: text

   >speciesA_gene1_sample1
   ATGCATGC...
   >speciesA_gene1_sample2
   ATGCATGC...
   >speciesB_gene1_sample1
   ATGCATGC...

Then filter with:

.. code-block:: bash

   mkado test alignment.fa -i "speciesA" -o "speciesB"

The pattern matching is case-sensitive substring matching.

Ingroup and Outgroup
--------------------

Definitions
^^^^^^^^^^^

- **Ingroup**: The species of primary interest (typically polymorphic population samples)
- **Outgroup**: A closely related species used to identify fixed differences

Selection Guidelines
^^^^^^^^^^^^^^^^^^^^

- Outgroup should be close enough to have reliable alignments
- But distant enough to have accumulated fixed differences
- Multiple outgroup sequences can be used (consensus is taken)

For polarized tests:

- **Second outgroup**: A more distant species to determine mutation direction
- Should be divergent from both ingroup and primary outgroup

Handling Special Cases
----------------------

Gaps
^^^^

- Codons containing gaps (``---``) are excluded from analysis
- Partial gaps within codons are handled conservatively

Stop Codons
^^^^^^^^^^^

- Internal stop codons are flagged but not automatically excluded
- Check your alignments if you see unexpected results

Ambiguous Bases
^^^^^^^^^^^^^^^

- ``N`` and other IUPAC ambiguity codes are supported
- Codons with ambiguous bases may be excluded from certain calculations

Common Problems
---------------

Alignment Not in Frame
^^^^^^^^^^^^^^^^^^^^^^

**Symptom**: Unexpected results, many excluded codons

**Solution**: Use ``mkado info`` to check alignment properties, specify ``-r`` for reading frame

.. code-block:: bash

   mkado info alignment.fa
   mkado test alignment.fa -i sp1 -o sp2 -r 2  # Try reading frame 2

Wrong Species Pattern
^^^^^^^^^^^^^^^^^^^^^

**Symptom**: "No sequences found" error

**Solution**: Check sequence names and pattern

.. code-block:: bash

   mkado info alignment.fa  # Lists all sequence names
   mkado test alignment.fa -i "correct_pattern" -o "outgroup"

Unequal Sequence Lengths
^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom**: Error about alignment length

**Solution**: Re-align sequences or check for truncated sequences

Example Data
------------

The ``examples/`` directory contains properly formatted example data:

.. code-block:: bash

   # Examine example file
   mkado info examples/anopheles_batch/AGAP000078.fa

   # Run test on example
   mkado test examples/anopheles_batch/AGAP000078.fa -i gamb -o afun
