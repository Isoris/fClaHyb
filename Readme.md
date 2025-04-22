# 🧬 Dual Reference Genomes from F1 Hybrids: Phased Assembly of North African Catfish and Bighead Catfish with Hi-C Data

This repository documents the high-quality, haplotype-resolved genome assemblies of two species — the North African catfish (*Clarias gariepinus*) and the Bighead catfish (*Clarias macrocephalus*) — using a single F1 hybrid individual. By leveraging Hi-C chromatin capture along with long- and short-read sequencing technologies, this project provides a dual reference for aquaculture genomics, hybridization studies, and evolutionary research.

---

## 📊 Project Summary

- 🧬 Haplotype-resolved diploid assembly from an F1 hybrid
- 🧪 PacBio HiFi, ONT, Illumina PE150, and Hi-C data
- 🔬 Phased parental genomes: *C. macrocephalus* (hap1) and *C. gariepinus* (hap2)
- 🧩 Fully scaffolded and manually curated

This resource is especially relevant for aquaculture in **Thailand** and the **Mekong Basin**.

---

## 📂 Data Access and Records

### 🔖 NCBI BioProject  
- 📦 [PRJNA1153495 – Hybrid Genome Project](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1153495)

### 🧬 Genome Assemblies (NCBI Nucleotide)

| Species | Accession | Description |
|---------|-----------|-------------|
| 🐠 Haplotype 1 (*C. macrocephalus*) | [JBLWFY000000000.1](https://www.ncbi.nlm.nih.gov/nuccore/JBLWFY000000000.1) | Bighead catfish haplotype |
| 🐟 Haplotype 2 (*C. gariepinus*) | [JBLWFZ000000000.1](https://www.ncbi.nlm.nih.gov/nuccore/JBLWFZ000000000.1) | North African catfish haplotype |

### 🧪 Raw Sequencing Reads (NCBI SRA)

| Data Type | Accession | Description |
|-----------|-----------|-------------|
| 🔬 Nanopore | [SRR30599638](https://www.ncbi.nlm.nih.gov/sra/SRR30599638) | Long reads with high error rate (~20%) |
| 🧪 HiFi | [SRR30599641](https://www.ncbi.nlm.nih.gov/sra/SRR30599641) | High-fidelity long reads |
| 🖥️ Illumina PE150 | [SRR30599640](https://www.ncbi.nlm.nih.gov/sra/SRR30599640) | Short-read data |
| 🧲 Hi-C PE150 | [SRR30599639](https://www.ncbi.nlm.nih.gov/sra/SRR30599639) | Hi-C chromatin capture |

### 🧫 BioSamples

| Sample | Accession | Description |
|--------|-----------|-------------|
| 🧬 Hybrid F1 | [SAMN43395848](https://www.ncbi.nlm.nih.gov/biosample/SAMN43395848) | F1 hybrid individual (same but combined)|
| 👩‍🔬 Female *C. macrocephalus* | [SAMN42503781](https://www.ncbi.nlm.nih.gov/biosample/SAMN42503781) | Haplotype 1 (C. macrocephalus) (Bio)sample |
| 👨‍🔬 Male *C. gariepinus* | [SAMN43548335](https://www.ncbi.nlm.nih.gov/biosample/SAMN43548335) | Haplotype 2 (C. gariepinus) (Bio)sample |
| 🧬 Additional *C. gariepinus* | [SAMN46907808](https://www.ncbi.nlm.nih.gov/biosample/SAMN46907808), [SAMN46730924](https://www.ncbi.nlm.nih.gov/biosample/SAMN46730924) | Additional references |

---

## 🧪 Final Assemblies

| Version | Description | Output File |
|---------|-------------|-------------|
| 🧬 v19+ | Haplotype 1 – Bighead catfish | `fClaHyb_Mac*.fa` |
| 🧬 v19+ | Haplotype 2 – North African catfish | `fClaHyb_Gar*.fa` |
| 📦 Flye | Collapsed diploid assembly | `CGAF_Flye.assembly.fasta` |
| 🧩 GreenHill | Phased consensus assembly | `out_afterPhase.fa`, `out_ConsensusOutput.fa` |

---

## 🌱 Applications

- Comparative genomics across *Clarias* species  
- Genetic studies of hybridization and adaptation  
- Selective breeding and aquaculture development  
- Support for genome graph models and pangenomics

This dataset contributes to genetic sustainability in aquaculture and biodiversity in Southeast Asia.

---

## 📚 Citation & Reuse

All datasets and scripts are available on Zenodo:

🔗 **Zenodo DOI:** [10.5281/zenodo.14864260](https://doi.org/10.5281/zenodo.14864260)

> Please cite the Zenodo record when reusing the data or referencing this project in your work.

---

 Maintained by Quentin L.S. Andres, PhD. Project under Kasetsart University, THAILAND, Animal Genomics and Bioresources Research Unit.
