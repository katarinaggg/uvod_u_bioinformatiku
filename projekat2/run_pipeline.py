import argparse
import subprocess
from pathlib import Path
import json


def run_step(step_name, command, step_dir):
    step_dir.mkdir(parents=True, exist_ok=True)
    done_file = step_dir / f".{step_name}.done"
    log_file = step_dir / f"{step_name}.log"

    if done_file.exists():
        print(f"[SKIP] {step_name} je već završen.")
        return

    for tmp in step_dir.glob("*.tmp.*"):
        tmp.unlink()
        print(f"[CLEAN] Obrisan nepotpun fajl: {tmp}")

    print(f"[RUN] {step_name}")

    with open(log_file, "w") as log:
        result = subprocess.run(
            command,
            stdout=log,
            stderr=log,
        )

    if result.returncode != 0:
        print(f"[ERROR] {step_name} nije uspeo. Pogledaj log: {log_file}")
        raise RuntimeError(f"Korak {step_name} nije uspeo.")

    done_file.touch()
    print(f"[DONE] {step_name}")


def check_index_files(reference, known_sites):
    ref = Path(reference)
    ks = Path(known_sites)

    potrebni = {
        "reference .fai": Path(str(ref) + ".fai"),
        "reference .dict": ref.with_suffix(".dict"),
        "reference .64.bwt": Path(str(ref) + ".64.bwt"),
        "reference .64.sa": Path(str(ref) + ".64.sa"),
        "reference .64.pac": Path(str(ref) + ".64.pac"),
        "reference .64.amb": Path(str(ref) + ".64.amb"),
        "reference .64.ann": Path(str(ref) + ".64.ann"),
        "known_sites .tbi": Path(str(ks) + ".tbi"),
    }

    for naziv, putanja in potrebni.items():
        if not putanja.exists():
            raise FileNotFoundError(f"Nedostaje indeks fajl ({naziv}): {putanja}")

    print("[OK] Svi indeks fajlovi postoje.")


def main():
    parser = argparse.ArgumentParser(description="Bioinformatics pipeline")
    parser.add_argument("--fastq1", required=True)
    parser.add_argument("--fastq2", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--known-sites", required=True)
    parser.add_argument("--output-dir", required=True)

    args = parser.parse_args()

    fastq1 = Path(args.fastq1)
    fastq2 = Path(args.fastq2)
    reference = Path(args.reference)
    known_sites = Path(args.known_sites)
    output_dir = Path(args.output_dir)

    # provera ulaznih fajlova
    if not fastq1.exists():
        raise FileNotFoundError(f"Ne postoji FASTQ1 fajl: {fastq1}")
    if not fastq2.exists():
        raise FileNotFoundError(f"Ne postoji FASTQ2 fajl: {fastq2}")
    if not reference.exists():
        raise FileNotFoundError(f"Ne postoji referentni genom: {reference}")
    if not known_sites.exists():
        raise FileNotFoundError(f"Ne postoje known sites: {known_sites}")

    # provera indeks fajlova
    check_index_files(reference, known_sites)


    # FastQc

    fastqc_dir = output_dir / "01_fastqc"
    run_step(
        "fastqc",
        [
            "fastqc",
            str(fastq1),
            str(fastq2),
            "-o",
            str(fastqc_dir)
        ],
        fastqc_dir
    )

    # BWA-MEM + sortiranje

    alignment_dir = output_dir / "02_alignment"
    sorted_bam = alignment_dir / "aligned.sorted.bam"
    tmp_bam = alignment_dir / "aligned.sorted.tmp.bam"

    rg = r"@RG\tID:sample\tSM:sample\tPL:ILLUMINA\tLB:lib1"

    run_step(
        "alignment",
        [
            "bash", "-c",
            f"set -o pipefail; bwa mem -R '{rg}' {reference} {fastq1} {fastq2} "
            f"| samtools sort -o {tmp_bam} && mv {tmp_bam} {sorted_bam}"
        ],
        alignment_dir
    )
    # indeksiranje posebno
    run_step(
        "alignment_index",
        [
            "samtools",
            "index",
            str(sorted_bam)
        ],
        alignment_dir
    )   

    # mark duplicates

    markdup_dir = output_dir / "03_mark_duplicates"
    marked_bam = markdup_dir / "marked_duplicates.bam"
    tmp_marked_bam = markdup_dir / "marked_duplicates.tmp.bam"
    metrics_file = markdup_dir / "marked_duplicates_metrics.txt"

    run_step(
        "mark_duplicates",
        [
            "bash", "-c",
            f"gatk MarkDuplicates "
            f"-I {sorted_bam} "
            f"-O {tmp_marked_bam} "
            f"-M {metrics_file} "
            f"&& mv {tmp_marked_bam} {marked_bam}"
        ],
        markdup_dir
    )
    # indeksiranje posebno
    run_step(
        "markdup_index",
        ["samtools", "index", str(marked_bam)],
        markdup_dir
    )


    # BSQR

    bqsr_dir = output_dir / "04_bqsr"
    recal_table = bqsr_dir / "recal_data.table"
    recal_bam = bqsr_dir / "recalibrated.bam"
    tmp_recal_bam = bqsr_dir / "recalibrated.tmp.bam"

    run_step(
        "base_recalibrator",
        [
            "gatk",
            "BaseRecalibrator",
            "-I", str(marked_bam),
            "-R", str(reference),
            "--known-sites", str(known_sites),
            "-O", str(recal_table)
        ],
        bqsr_dir
    )
    # primena recalibration tabele
    run_step(
        "apply_bqsr",
        [
            "bash", "-c",
            f"gatk ApplyBQSR "
            f"-I {marked_bam} "
            f"-R {reference} "
            f"--bqsr-recal-file {recal_table} "
            f"-O {tmp_recal_bam} "
            f"&& mv {tmp_recal_bam} {recal_bam} "
            f"&& rm -f {bqsr_dir / 'recalibrated.tmp.bai'}"
        ],
        bqsr_dir
    )
    #indeksiranje bam-a
    run_step(
        "bqsr_index",
        ["samtools", "index", str(recal_bam)],
        bqsr_dir
    )

    # HaplotypeCaller 

    variants_dir = output_dir / "05_variant_calling"
    raw_vcf = variants_dir / "raw_variants.vcf.gz"
    tmp_raw_vcf = variants_dir / "raw_variants.tmp.vcf.gz"

    run_step(
        "haplotype_caller",
        [
            "bash", "-c",
            f"gatk HaplotypeCaller "
            f"-R {reference} "
            f"-I {recal_bam} "
            f"-O {tmp_raw_vcf} "
            f"&& mv {tmp_raw_vcf} {raw_vcf} "
            f"&& mv {tmp_raw_vcf}.tbi {raw_vcf}.tbi"
        ],
        variants_dir
    )

    # hard filtering

    filter_dir = output_dir / "06_hard_filtering"

    # odvajanje SNP-ova
    snp_vcf = filter_dir / "snps.vcf.gz"
    tmp_snp_vcf = filter_dir / "snps.tmp.vcf.gz"

    run_step(
        "select_snps",
        [
            "bash", "-c",
            f"gatk SelectVariants "
            f"-R {reference} "
            f"-V {raw_vcf} "
            f"--select-type-to-include SNP "
            f"-O {tmp_snp_vcf} "
            f"&& mv {tmp_snp_vcf} {snp_vcf} "
            f"&& mv {tmp_snp_vcf}.tbi {snp_vcf}.tbi"
        ],
        filter_dir
    )

    # odvajanje INDEL-a
    indel_vcf = filter_dir / "indels.vcf.gz"
    tmp_indel_vcf = filter_dir / "indels.tmp.vcf.gz"

    run_step(
        "select_indels",
        [
            "bash", "-c",
            f"gatk SelectVariants "
            f"-R {reference} "
            f"-V {raw_vcf} "
            f"--select-type-to-include INDEL "
            f"-O {tmp_indel_vcf} "
            f"&& mv {tmp_indel_vcf} {indel_vcf} "
            f"&& mv {tmp_indel_vcf}.tbi {indel_vcf}.tbi"
        ],
        filter_dir
    )

    # filtriranje SNP-ova
    filtered_snp_vcf = filter_dir / "snps_filtered.vcf.gz"
    tmp_filtered_snp_vcf = filter_dir / "snps_filtered.tmp.vcf.gz"

    run_step(
        "filter_snps",
        [
            "bash", "-c",
            f"gatk VariantFiltration "
            f"-R {reference} "
            f"-V {snp_vcf} "
            f"--filter-expression 'QD < 2.0' --filter-name 'QD2' "
            f"--filter-expression 'FS > 60.0' --filter-name 'FS60' "
            f"--filter-expression 'MQ < 40.0' --filter-name 'MQ40' "
            f"--filter-expression 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' "
            f"--filter-expression 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' "
            f"-O {tmp_filtered_snp_vcf} "
            f"&& mv {tmp_filtered_snp_vcf} {filtered_snp_vcf} "
            f"&& mv {tmp_filtered_snp_vcf}.tbi {filtered_snp_vcf}.tbi"
        ],
        filter_dir
    )

    # filtriranje INDEL-a
    filtered_indel_vcf = filter_dir / "indels_filtered.vcf.gz"
    tmp_filtered_indel_vcf = filter_dir / "indels_filtered.tmp.vcf.gz"

    run_step(
        "filter_indels",
        [
            "bash", "-c",
            f"gatk VariantFiltration "
            f"-R {reference} "
            f"-V {indel_vcf} "
            f"--filter-expression 'QD < 2.0' --filter-name 'QD2' "
            f"--filter-expression 'FS > 200.0' --filter-name 'FS200' "
            f"--filter-expression 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' "
            f"-O {tmp_filtered_indel_vcf} "
            f"&& mv {tmp_filtered_indel_vcf} {filtered_indel_vcf} "
            f"&& mv {tmp_filtered_indel_vcf}.tbi {filtered_indel_vcf}.tbi"
        ],
        filter_dir
    )

    # spajanje filtriranih SNP-ova i INDEL-a
    merged_vcf = filter_dir / "filtered_variants.vcf.gz"
    tmp_merged_vcf = filter_dir / "filtered_variants.tmp.vcf.gz"

    run_step(
        "merge_filtered",
        [
            "bash", "-c",
            f"gatk MergeVcfs "
            f"-I {filtered_snp_vcf} "
            f"-I {filtered_indel_vcf} "
            f"-O {tmp_merged_vcf} "
            f"&& mv {tmp_merged_vcf} {merged_vcf} "
            f"&& mv {tmp_merged_vcf}.tbi {merged_vcf}.tbi"
        ],
        filter_dir
    )

    # izdvajanje samo PASS varijanti

    pass_dir = output_dir / "07_pass_variants"
    pass_vcf = pass_dir / "pass_variants.vcf.gz"
    tmp_pass_vcf = pass_dir / "pass_variants.tmp.vcf.gz"

    run_step(
        "extract_pass",
        [
            "bash", "-c",
            f"gatk SelectVariants "
            f"-R {reference} "
            f"-V {merged_vcf} "
            f"--exclude-filtered true "
            f"-O {tmp_pass_vcf} "
            f"&& mv {tmp_pass_vcf} {pass_vcf} "
            f"&& mv {tmp_pass_vcf}.tbi {pass_vcf}.tbi"
        ],
        pass_dir
    )

    # json manifest

    manifest_dir = output_dir / "08_manifest"
    manifest_dir.mkdir(parents=True, exist_ok=True)

    fastq1_base = fastq1.name.replace(".fq.gz", "").replace(".fastq.gz", "")
    fastq2_base = fastq2.name.replace(".fq.gz", "").replace(".fastq.gz", "")

    manifest = {
        "fastqc_r1": str(fastqc_dir / f"{fastq1_base}_fastqc.html"),
        "fastqc_r2": str(fastqc_dir / f"{fastq2_base}_fastqc.html"),
        "bam_sorted": str(sorted_bam),
        "bam_markdup": str(marked_bam),
        "bam_recalibrated": str(recal_bam),
        "vcf_raw": str(raw_vcf),
        "vcf_filtered": str(merged_vcf),
        "vcf_pass": str(pass_vcf)
    }

    manifest_path = manifest_dir / "manifest.json"

    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"[DONE] Manifest saved to {manifest_path}")


if __name__ == "__main__":
    main()