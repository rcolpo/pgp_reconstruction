# Inputs accepted and how to improve the quality of models

This folder contains examples of genomic information that PGP-Reconstruction accepts as input files. Additionally, it provides an example of a `constraint.tsv` file for manually curating your model, and how to use a reference model.

## Input Files

In this folder you find 3 types of files accepted has input files, all [retrived from NCBI](https://www.ncbi.nlm.nih.gov/genome/167?genome_assembly_id=368427) and below you find their description.
- `Escherichia_coli_strain_2014C-3716.faa`: A FASTA file containing nucleic acid sequences. PGP-Reconstruction uses Prodigal to identify Open Reading Frames (ORFs) and translate them into amino acid sequences, which are then aligned with UniRef90.
- `Escherichia_coli_strain_2014C-3716.fna`: A FASTA file with amino acid sequences. Some sequences contain annotation information that PGP-Reconstruction can read. This input method is preferred over nucleic acid sequences.
- `Escherichia_coli_strain_2014C-3716.gbff`: A GenBank file with both amino acid and nucleic acid sequences. Despite only the amino acid sequences being read, this method is preferred over FASTA files containing amino acid sequences, as the annotation information is written in a more standardized format.

## Constraints File

The `constraint.tsv` file allows you to manually curate the model by including reactions with experimental evidence of existing (or of not existing). This file can include constraints for both metabolites and reactions.

The file has four columns, separated by a TAB: `M_Metabolite/R_Reaction Id`, `Constraint type`, `Score`, and `Consumed/Produced`.

### Metabolite Constraints

- `M_Metabolite`: Accepts IDs from ChEBI, KEGG, BiGG, and SEED, preceded by an `M_` prefix.
- `Constraint type`: Can be either `Soft` or `Hard`. With `Soft`, PGP-Reconstruction attempts to satisfy the constraint. With `Hard`, PGP-Reconstruction returns an error message if it can't satisfy the constraint (infeasible optimization).
- `Score`: Use a positive score for metabolites you wish to be consumed/produced. Use a negative score for metabolites you wish to not be consumed/produced. To avoid infeasible optimization problems potentially arising from a `Hard` constraint type, you may combine a very high (or very low) score with a `Soft` constraint type. 10 (or -10) is an example of a very high (or very low) score.
- `Consumed/Produced`: Only valid for metabolites. Can be either `Media` or `Product`. With `Media`, PGP-Reconstruction attempts to uptake the metabolite. With `Product`, it tries to produce and secrete the metabolite into the external media.

### Reaction Constraints

- `R_Reaction`: Accepts IDs from Rhea, KEGG, BiGG, and SEED, preceded by an `R_` prefix.
- `Constraint type`: Can be either `Soft` or `Hard`. With `Soft`, PGP-Reconstruction tries to satisfy the constraint. With `Hard`, it returns an error message if it can't satisfy the constraint (infeasible optimization).
- `Score`: Use a positive score for reactions that should be present. Use a negative score for reactions that should not be present. To avoid infeasible optimization problems potentially arising from a `Hard` constraint type, you may combine a very high (or very low) score with a `Soft` constraint type. 10 (or -10) is an example of a very high (or very low) score.

## Reference Model

If your organism of interest already has a manually curated model, you can include it as input to PGP-Reconstruction, regardless of its namespace. The model will be treated as a reference model, and all reactions in the model will be translated to the Rhea namespace and receive a slightly negative score (-0.1) if they do not already have a positive score. This indicates to PGP-Reconstruction to use the reactions from the reference model as the first option when "gap-filling" the network. 

Included in this folder is [iEC1344_C.xml](http://bigg.ucsd.edu/models/iEC1344_C), a manually curated model of Escherichia coli, reconstructed using the BiGG namespace.

## Usage

Here's an example of how to use the files in this folder:

```bash
pgprec Escherichia_coli_strain_2014C-3716.gbff --constraints constraints.tsv --reference "iEC1344_C.xml"
```