# Projekat 2 – BWA-GATK Pipeline
Katarina Gordić 41/2023 RN

## Opis

Ovaj projekat implementira reproducibilan Python pipeline za obradu uparenih FASTQ fajlova prema BWA-GATK protokolu.

Pipeline izvršava sledeće korake:

1. FastQC kontrola kvaliteta
2. Mapiranje read-ova pomoću BWA-MEM
3. Sortiranje i indeksiranje BAM fajla
4. Markiranje duplikata
5. Base Quality Score Recalibration (BQSR)
6. Pozivanje varijanti pomoću GATK HaplotypeCaller
7. Hard filtering SNP i INDEL varijanti
8. Izdvajanje PASS varijanti
9. Generisanje JSON manifesta rezultata

## Korišćeni alati

|   Alat   |   Verzija   |
|----------|-------------|
|  FastQC  | v0.12.1     |
|   BWA    | 0.7.19-r1273|
| samtools | 1.22.1      |
|   GATK   | 4.6.2.0     |
|  Python  | 3.14.4      |


## Struktura projekta

- `run_pipeline.py` – glavna pipeline skripta
- `analysis/` – notebook-ovi za analizu rezultata
- `results/` – rezultatni fajlovi pipeline-a
- `izvestaj.pdf` – izveštaj analize
- `README.md` – dokumentacija projekta

**Napomena:**  ulazni FASTQ fajlovi, referentni genom, known-sites VCF i indeks fajlovi nisu uključeni u repozitorijum u skladu sa zahtevima zadatka.


## Potrebni ulazni fajlovi

Za pokretanje pipeline-a potrebno je obezbediti:

- uparene FASTQ fajlove (`R1` i `R2`)
- referentni genom u FASTA formatu
- VCF fajl sa poznatim varijantama za BQSR (`known-sites`)


## Potrebni indeks fajlovi

Pipeline očekuje da sledeći indeks fajlovi već postoje:

- referentni genom `.fai`
- referentni genom `.dict`
- BWA indeks fajlovi (`.64.amb`, `.64.ann`, `.64.bwt`, `.64.pac`, `.64.sa`)
- `.tbi` indeks za known-sites VCF

Ako neki od potrebnih indeks fajlova nedostaje, izvršavanje se prekida uz odgovarajuću poruku o grešci.


## Pokretanje pipeline-a

Primer pokretanja:

```bash
python run_pipeline.py \
  --fastq1 sample_R1.fastq.gz \
  --fastq2 sample_R2.fastq.gz \
  --reference hg38.fa \
  --known-sites known_sites.vcf.gz \
  --output-dir results
```

## Argumenti

| Argument | Opis |
|-----------|------|
| `--fastq1` | FASTQ fajl sa readovima R1 |
| `--fastq2` | FASTQ fajl sa readovima R2 |
| `--reference` | Referentni genom |
| `--known-sites` | VCF fajl sa poznatim varijantama za BQSR |
| `--output-dir` | Izlazni direktorijum |


## Idempotentnost i nastavak izvršavanja

Pipeline koristi `.done` marker fajlove za evidenciju uspešno završenih koraka.

Ako je korak već uspešno izvršen, pri narednom pokretanju biće automatski preskočen.

Rezultati se prvo upisuju u privremene (`.tmp`) fajlove, a tek nakon uspešnog završetka preimenuju u finalne izlaze.

Privremeni fajlovi nastali tokom neuspelog ili prekinutog izvršavanja uklanjaju se pre ponovnog pokretanja odgovarajućeg koraka.

Sav `stdout` i `stderr` izlaz spoljašnjih alata čuva se u log fajlovima.


## Glavni rezultatni fajlovi

Pipeline generiše:

- FastQC HTML izveštaje
- sortiran i indeksiran BAM fajl
- BAM fajl nakon markiranja duplikata
- BAM fajl nakon BQSR obrade
- sirovi VCF fajl
- filtrirani VCF fajl sa popunjenom FILTER kolonom
- PASS VCF fajl
- JSON manifest sa putanjama do glavnih rezultata


## Analiza rezultata

Kod za izvlačenje metrika i generisanje grafikona nalazi se u direktorijumu `analysis/`:

- `metrics.ipynb`
- `variant_analysis.ipynb`

PDF izveštaj analize nalazi se u fajlu:

- `izvestaj.pdf`