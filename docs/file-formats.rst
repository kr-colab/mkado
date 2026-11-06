File Formats and Input Requirements
====================================

MKado supports two input modes: **FASTA** (aligned coding sequences for ``mkado test`` and ``mkado batch``) and **VCF** (variant calls for ``mkado vcf``). This page describes the input requirements for both.

FASTA Format
------------

MKado reads standard FASTA format:

.. code-block:: text

   >sequence_name1 optional description
   ATGCATGCATGC...
   >sequence_name2
   ATGCATGCATGC...

Requirements:

- Standard FASTA format with ``>`` headers. Each header must have a name, and the file must
  start with a header
- Sequences can span multiple lines. MKado ignores spaces and tabs inside a sequence line, so
  a sequence written in blocks of 10 bases reads correctly
- Bases ``A``, ``C``, ``G`` and ``T``, in uppercase or lowercase. MKado reads any other character, such as ``N``, a gap or an ambiguity code, as missing data (see :ref:`missing-data`)
- No digits in a sequence line. A numbered layout, as in GenBank files, is not FASTA, so MKado
  stops with an error that gives the file and the line number

Alignment Requirements
----------------------

Sequences must be pre-aligned:

- All sequences must be the same length. If they are not, MKado stops with an error
  that gives the most common lengths and the names of some sequences that have them.
  The ingroup, outgroup and second outgroup files must also have the same length as
  each other.
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

The file is one alignment, so MKado checks the length of every sequence in it, also the
sequences that no pattern selects.

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

.. _missing-data:

Missing and Ambiguous Data
^^^^^^^^^^^^^^^^^^^^^^^^^^

MKado uses a codon only if all three bases are ``A``, ``C``, ``G`` or ``T``. Any other character makes the codon missing data. This includes ``N``, a gap (``-``), an IUPAC ambiguity code (``R``, ``Y``, ``K`` and so on) and ``?``. One such base removes the full codon, not only that position.

- MKado removes the sequence from that codon only. The other sequences at that codon are still used, and that sequence is still used at all other codons.
- Allele frequencies at that codon use only the sequences that remain.
- If the ingroup or the outgroup has no clean sequence at a codon, MKado leaves that codon out of divergence (Dn, Ds), polymorphism (Pn, Ps) and the site counts (Ln, Ls). The dropped polymorphisms are counted and reported (see :ref:`Unpolarizable sites <unpolarizable-sites>`). The polarized MK test does not use this rule for polymorphism yet.
- Lowercase bases are read as uppercase.

.. _stop-codons:

Stop Codons
^^^^^^^^^^^

MKado reads a stop codon of the genetic code (``--code-table``) as missing data, the same as a
codon with an ``N``:

- A stop codon adds nothing to divergence (Dn, Ds), polymorphism (Pn, Ps) or the site counts
  (Ln, Ls). The terminal stop codon of a gene is left out in the same way.
- A stop codon in one sequence is missing data for that sequence only. The other sequences
  at that codon are counted, and the frequency of each allele is taken over them.
- A change to or from a stop codon is therefore not counted. A premature stop codon is not
  a nonsynonymous polymorphism.

Each gene with an internal stop codon gets one warning line with the count, for example::

   Warning: AGAP000704: 2 internal stop codon(s) treated as missing data

A stop codon is internal when its own sequence has a later codon with an A, C, G or T base. A
stop codon that only gaps or ``N`` follow is the terminal stop of that sequence, so a sequence
that ends early in a gap-padded alignment does not give a warning.

An internal stop codon in many sequences can mean a pseudogene, a wrong reading frame or an
alignment error, so check such genes. The VCF path uses the same rule (see :doc:`vcf-input`).

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

**Symptom**: ``Sequences have different lengths`` or ``Alignment length mismatch``

**Solution**: Re-align the sequences, or remove the truncated sequences. The first error
gives the most common lengths and the names of some sequences that have them. The second
error compares the ingroup with an outgroup file. In combined file mode, a sequence that no
pattern selects also causes the first error.

VCF Mode Input Files
--------------------

The ``mkado vcf`` command requires three file types in addition to a GFF3 annotation. See :doc:`vcf-input` for the full guide.

VCF Files
^^^^^^^^^

Both the ingroup (``--vcf``) and outgroup (``--outgroup-vcf``) inputs are standard VCF files:

- Bgzipped + tabix-indexed is recommended for performance, but uncompressed VCF also works
- The ingroup VCF should be a multi-sample population VCF with genotype fields
- The outgroup VCF should be a single-sample VCF called against the **same reference genome**
- Multi-allelic sites should be decomposed beforehand with for example ``bcftools norm -m-``

Reference FASTA
^^^^^^^^^^^^^^^

The genome assembly both VCFs were called against (``--ref``):

- Must be indexed with ``samtools faidx`` (creates a ``.fai`` file)
- Both plain FASTA and bgzipped FASTA (with ``.gzi`` index) are supported
- Only **bgzip** compression is supported, not plain gzip (random access requires BGZF block structure)

GFF3 Annotation
^^^^^^^^^^^^^^^

A GFF3 file defining gene models (``--gff``):

- Must contain ``gene``, ``mRNA``/``transcript``, and ``CDS`` features linked by ``Parent`` attributes
- Both plain text and gzip-compressed (``.gff3.gz``) files are supported
- MKado selects the longest transcript per gene automatically
- Genes where the total CDS length is not divisible by 3 are skipped
- **GTF format is not supported** — MKado requires GFF3 (``key=value`` attributes, not GTF's ``key "value"`` style). Convert with ``gffread annotation.gtf -o annotation.gff3`` if needed.

Example Data
------------

The ``examples/`` directory contains properly formatted FASTA example data:

.. code-block:: bash

   # Examine example file
   mkado info examples/anopheles_batch/AGAP000078.fa

   # Run test on example
   mkado test examples/anopheles_batch/AGAP000078.fa -i gamb -o afun
