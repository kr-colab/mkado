VCF Input Mode
==============

MKado can run MK tests directly from VCF files, skipping the need to reconstruct per-gene FASTA alignments. The ``mkado vcf`` command takes an ingroup VCF, an outgroup VCF, a reference genome, and a GFF3 annotation, and produces the same results as the FASTA-based workflows.

This is useful when your starting point is a variant calling pipeline (e.g., GATK, bcftools) rather than pre-aligned coding sequences.

VCF support is included in the standard ``mkado`` install via `cyvcf2 <https://github.com/brentp/cyvcf2>`_ (VCF parsing) and `pysam <https://github.com/pysam-developers/pysam>`_ (indexed FASTA access).

Required Input Files
--------------------

``--vcf`` : Ingroup VCF
^^^^^^^^^^^^^^^^^^^^^^^

A multi-sample population VCF file. Bgzipped and tabix-indexed is recommended for performance, but uncompressed VCF also works.

.. code-block:: bash

   # Prepare your VCF (if not already indexed)
   bgzip population.vcf
   tabix -p vcf population.vcf.gz

``--outgroup-vcf`` : Outgroup VCF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A single-sample VCF of the outgroup species, called against the **same reference genome** as the ingroup VCF. This is used to determine divergence (Dn/Ds) and to polarize polymorphisms.

``--ref`` : Reference FASTA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The genome assembly that both VCFs were called against. Must be indexed with ``samtools faidx``. Both plain FASTA and bgzipped FASTA are supported:

.. code-block:: bash

   # Plain FASTA (creates .fai index)
   samtools faidx reference.fa

   # Bgzipped FASTA (creates .fai and .gzi indices)
   bgzip reference.fa
   samtools faidx reference.fa.gz

.. note::

   Only bgzipped (BGZF) compression is supported for the reference FASTA, not plain gzip. This is because random access requires BGZF block structure. Use ``bgzip`` from htslib, not ``gzip``.

``--gff`` : GFF3 Annotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A GFF3 file defining gene models with CDS features. Both plain text and gzip-compressed files are supported:

.. code-block:: bash

   # Plain GFF3
   mkado vcf --gff annotation.gff3 ...

   # Gzipped GFF3
   mkado vcf --gff annotation.gff3.gz ...

The parser extracts CDS features, groups them by transcript (via the ``Parent`` attribute), and selects the longest transcript per gene. Genes where the total CDS length is not divisible by 3 are skipped.

Basic Usage
-----------

.. code-block:: bash

   # Standard MK test across all genes
   mkado vcf \
       --vcf population.vcf.gz \
       --outgroup-vcf outgroup.vcf.gz \
       --ref reference.fa \
       --gff annotation.gff3

   # Asymptotic MK test (aggregated across genes)
   mkado vcf \
       --vcf population.vcf.gz \
       --outgroup-vcf outgroup.vcf.gz \
       --ref reference.fa \
       --gff annotation.gff3 \
       --asymptotic

   # Single gene
   mkado vcf \
       --vcf population.vcf.gz \
       --outgroup-vcf outgroup.vcf.gz \
       --ref reference.fa \
       --gff annotation.gff3 \
       --gene BRCA1

Gene Selection
--------------

By default, all genes in the GFF3 are processed. You can restrict to specific genes:

.. code-block:: bash

   # Single gene by ID
   mkado vcf ... --gene BRCA1

   # Subset of genes from a file (one ID per line)
   mkado vcf ... --gene-list genes_of_interest.txt

Gene IDs are matched against both the ``ID`` and ``Name`` attributes of gene features in the GFF3.

Analysis Modes
--------------

All analysis types available in ``mkado batch`` are supported:

.. code-block:: bash

   # Standard MK test (default)
   mkado vcf ...

   # Asymptotic MK test (aggregated)
   mkado vcf ... --asymptotic

   # Asymptotic MK test (per-gene)
   mkado vcf ... --asymptotic --per-gene

   # Imputed MK test
   mkado vcf ... --imputed

   # Tarone-Greenland weighted alpha
   mkado vcf ... --alpha-tg

   # With frequency filtering
   mkado vcf ... --min-freq 0.05
   mkado vcf ... --no-singletons

See :doc:`tutorial` for details on each analysis type.

How It Works
------------

For each gene defined in the GFF3, ``mkado vcf`` does the following:

Polymorphism Extraction
^^^^^^^^^^^^^^^^^^^^^^^

1. Query the ingroup VCF for biallelic SNPs overlapping the gene's CDS exons
2. Skip indels and multi-allelic sites
3. For each SNP, map it to a codon using the CDS coordinates
4. Reconstruct the reference and alternative codons from the reference FASTA
5. For minus-strand genes, reverse complement the codons
6. Classify each change as synonymous or nonsynonymous using the genetic code
7. Compute the derived allele frequency from genotype counts
8. If the outgroup carries the ALT allele, flip the polarization (derived frequency = 1 - ALT frequency)

Divergence Extraction
^^^^^^^^^^^^^^^^^^^^^

1. For each codon, check if the outgroup VCF has variant(s) at those positions
2. Only count as divergence if the ingroup is monomorphic for the reference allele
3. Reconstruct the outgroup codon and compare to the reference codon
4. Classify using the shortest mutational path (via ``GeneticCode.get_path()``)

The output is ``PolymorphismData`` (the same intermediate format used by the FASTA-based pipeline), which feeds directly into all existing analysis functions.

Edge Cases
----------

- **Indels**: Skipped (count logged per gene)
- **Multi-allelic sites**: Skipped. Pre-decompose with ``bcftools norm -m-`` if needed
- **Missing genotypes**: Excluded from frequency calculation; only non-missing samples counted
- **Two SNPs in the same codon**: Each SNP classified independently (standard MK convention, avoids phasing issues)
- **Overlapping genes**: Each gene processed independently; the same variant can contribute to multiple genes
- **Stop codons**: Reference stop codons skipped; premature stop-creating variants treated as nonsynonymous
- **Incomplete CDS**: Genes where CDS length % 3 != 0 are skipped

Parallel Processing
-------------------

Like ``mkado batch``, the ``vcf`` command supports parallel gene processing:

.. code-block:: bash

   # Auto-detect worker count
   mkado vcf ... -w 0

   # Use 8 workers
   mkado vcf ... -w 8

   # Sequential (useful for debugging)
   mkado vcf ... -w 1

Output Formats
--------------

.. code-block:: bash

   # TSV (default)
   mkado vcf ... -f tsv > results.tsv

   # Pretty-printed
   mkado vcf ... -f pretty

   # JSON
   mkado vcf ... -f json > results.json

Preparing Your Data
-------------------

A typical workflow starting from raw reads:

.. code-block:: bash

   # 1. Align reads to reference
   bwa mem reference.fa reads_R1.fq.gz reads_R2.fq.gz | samtools sort -o aligned.bam

   # 2. Call variants (ingroup)
   bcftools mpileup -f reference.fa sample1.bam sample2.bam ... | \
       bcftools call -m -v -Oz -o population.vcf.gz
   tabix -p vcf population.vcf.gz

   # 3. Call variants (outgroup)
   bcftools mpileup -f reference.fa outgroup.bam | \
       bcftools call -m -v -Oz -o outgroup.vcf.gz
   tabix -p vcf outgroup.vcf.gz

   # 4. Decompose multi-allelic sites (recommended)
   bcftools norm -m- population.vcf.gz -Oz -o population.norm.vcf.gz
   tabix -p vcf population.norm.vcf.gz

   # 5. Index reference
   samtools faidx reference.fa

   # 6. Run MK tests
   mkado vcf \
       --vcf population.norm.vcf.gz \
       --outgroup-vcf outgroup.vcf.gz \
       --ref reference.fa \
       --gff annotation.gff3.gz

Comparison with FASTA Mode
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Aspect
     - ``mkado batch`` (FASTA)
     - ``mkado vcf``
   * - Input
     - Pre-aligned coding FASTA per gene
     - VCF + reference + GFF3
   * - Alignment
     - User must align sequences
     - Uses reference coordinates
   * - Frequency data
     - Derived from aligned sequences
     - Derived from VCF genotypes
   * - Reading frame
     - User specifies (``-r``)
     - From GFF3 CDS phase
   * - Splicing
     - User must handle
     - Automatic from GFF3 exon structure
   * - Dependencies
     - None beyond base mkado
     - cyvcf2, pysam (included in standard install)
