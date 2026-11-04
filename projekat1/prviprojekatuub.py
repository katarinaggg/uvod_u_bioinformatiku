import os
import gc
import glob
import pysam
import pandas as pd
import numpy as np
import math
from scipy.stats import chi2_contingency
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt

# putanje fajlova

bcf_putanja = "podaci/1KG.subset.bcf"

chunk_folder = "rezultati/genotipovi_chunkovi"

finalni_feather_putanja = "rezultati/genotipovi.feather"

status_putanja = "podaci/status_9.tsv"

gwas_feather_putanja = "rezultati/gwas_sve.feather"


def napravi_folder(putanja_foldera):
    os.makedirs(putanja_foldera, exist_ok=True)


def preuzmi_hromozome_iz_bcf(bcf_putanja):
    fajl = pysam.VariantFile(bcf_putanja)
    hromozomi = list(fajl.header.contigs.keys())  #iz headera
    fajl.close()
    return hromozomi


def putanja_za_chunk(chunk_folder, hromozom, broj_chunka):
    naziv = f"{hromozom}_chunk_{broj_chunka:04d}.feather"
    return os.path.join(chunk_folder, naziv)

def sacuvaj_chunk_u_feather(redovi, indeksi, uzorci, izlazna_putanja):
    df = pd.DataFrame(redovi, index=indeksi, columns=uzorci, dtype=np.uint8)
    df = df.reset_index().rename(columns={"index": "varijanta"})
    df.to_feather(izlazna_putanja)

    del df
    gc.collect()  # oslobadjanje memorije jer imamo dosta chunkova


def ispisi_ram():
    try:
        import psutil
        import os as _os
        proces = psutil.Process(_os.getpid())
        memorija_mb = proces.memory_info().rss / 1024 / 1024
        print(f"RAM procesa: {memorija_mb:.2f} MB")
    except Exception:
        print("Nije uspelo ocitavanje RAM-a.")


def ucitaj_hromozom_u_chunkove(
    bcf_putanja,
    hromozom,
    chunk_folder,
    chunk_size=10000,
    korak_ispisa=20000
):

    napravi_folder(chunk_folder)

    fajl = pysam.VariantFile(bcf_putanja)
    uzorci = list(fajl.header.samples)
    broj_uzoraka = len(uzorci)

    # bafer za jedan chunk
    redovi_chunka = []
    indeksi_chunka = []  # chr-pos

    # za provere
    broj_procitanih = 0
    broj_ubacenih = 0
    broj_preskocenih_multialelnih = 0
    broj_preskocenih_nevalidnih = 0

    broj_chunka = 1

    print(f"\n=== Obrada hromozoma {hromozom} ===")
    print(f"Broj uzoraka: {broj_uzoraka}")
    print(f"Chunk size: {chunk_size}")

    # optimizacije
    np_uint8 = np.uint8
    append_red = redovi_chunka.append
    append_indeks = indeksi_chunka.append

    # prolazak kroz sve varijante hromozoma
    for zapis in fajl.fetch(hromozom):
        broj_procitanih += 1

        # preskacemo ako ima vise od 1 alt
        alts = zapis.alts
        if alts is None or len(alts) != 1:
            broj_preskocenih_multialelnih += 1
            continue

        red = np.empty(broj_uzoraka, dtype=np_uint8)
        # ako je neki uzorak nevalidan preskacemo
        preskoci_red = False

        # prolazak kroz sve uzorke
        for i, uzorak in enumerate(uzorci):
            # uzimamo genotip, oblik: (x,y)
            gt = zapis.samples[uzorak]["GT"]

            if gt is None or len(gt) != 2:
                preskoci_red = True
                break

            a1, a2 = gt

            if a1 is None or a2 is None:
                preskoci_red = True
                break

            if a1 < 0 or a2 < 0:
                preskoci_red = True
                break

            if a1 > 1 or a2 > 1:
                preskoci_red = True
                break

            red[i] = (a1 == 1) + (a2 == 1)

        if preskoci_red:
            broj_preskocenih_nevalidnih += 1
            continue

        append_indeks(f"{zapis.chrom}-{zapis.pos}")
        append_red(red)
        broj_ubacenih += 1

        if len(redovi_chunka) >= chunk_size:
            izlazna_putanja = putanja_za_chunk(chunk_folder, hromozom, broj_chunka)

            # ako postoji preskacemo snimanje u fajl
            if not os.path.exists(izlazna_putanja):
                sacuvaj_chunk_u_feather(
                    redovi=redovi_chunka,
                    indeksi=indeksi_chunka,
                    uzorci=uzorci,
                    izlazna_putanja=izlazna_putanja
                )
                print(f"Sacuvan chunk: {izlazna_putanja} | redova: {len(redovi_chunka)}")
                ispisi_ram()
            else:
                print(f"Chunk vec postoji, preskacem snimanje: {izlazna_putanja}")

            redovi_chunka.clear()
            indeksi_chunka.clear()
            gc.collect()
            broj_chunka += 1

        if broj_procitanih % korak_ispisa == 0:
            print(
                f"Hromozom: {hromozom} | "
                f"Procitano: {broj_procitanih} | "
                f"Ubaceno: {broj_ubacenih} | "
                f"Preskoceno multialelnih: {broj_preskocenih_multialelnih} | "
                f"Preskoceno nevalidnih: {broj_preskocenih_nevalidnih} | "
                f"Trenutno u bufferu: {len(redovi_chunka)}"
            )
            ispisi_ram()

    # poslednji chunk
    if len(redovi_chunka) > 0:
        izlazna_putanja = putanja_za_chunk(chunk_folder, hromozom, broj_chunka)

        if not os.path.exists(izlazna_putanja):
            sacuvaj_chunk_u_feather(
                redovi=redovi_chunka,
                indeksi=indeksi_chunka,
                uzorci=uzorci,
                izlazna_putanja=izlazna_putanja
            )
            print(f"Sacuvan poslednji chunk: {izlazna_putanja} | redova: {len(redovi_chunka)}")
            ispisi_ram()
        else:
            print(f"Poslednji chunk vec postoji, preskacem snimanje: {izlazna_putanja}")

        redovi_chunka.clear()
        indeksi_chunka.clear()
        gc.collect()

    fajl.close()

    print(f"\n=== Zavrsen hromozom {hromozom} ===")
    print(f"Ukupno procitano: {broj_procitanih}")
    print(f"Ukupno ubaceno: {broj_ubacenih}")
    print(f"Preskoceno multialelnih: {broj_preskocenih_multialelnih}")
    print(f"Preskoceno nevalidnih: {broj_preskocenih_nevalidnih}")



def ucitaj_i_sacuvaj_sve_hromozome_u_chunkove(
    bcf_putanja,
    chunk_folder,
    chunk_size=10000,
    korak_ispisa=20000
):
    hromozomi = preuzmi_hromozome_iz_bcf(bcf_putanja)
    print("Pronadjeni hromozomi:", hromozomi)

    # prolazimo kroz hromozome i ucitavamo u chunkove
    for hromozom in hromozomi:
        ucitaj_hromozom_u_chunkove(
            bcf_putanja=bcf_putanja,
            hromozom=hromozom,
            chunk_folder=chunk_folder,
            chunk_size=chunk_size,
            korak_ispisa=korak_ispisa
        )

    print("\n=== Svi hromozomi su obradjeni u chunkove ===")



def ucitaj_chunk_u_dataframe(feather_putanja):
    df = pd.read_feather(feather_putanja)
    df = df.set_index("varijanta")
    df = df.astype(np.uint8)
    return df


def preuzmi_sortirane_chunk_putanje(chunk_folder):
    putanje = glob.glob(os.path.join(chunk_folder, "*.feather"))
    putanje.sort()
    return putanje


def spoji_chunkove_u_jedan_dataframe(chunk_folder):
    chunk_putanje = preuzmi_sortirane_chunk_putanje(chunk_folder)

    if len(chunk_putanje) == 0:
        raise ValueError("Nema feather chunk fajlova za spajanje.")

    delovi = []

    for putanja in chunk_putanje:
        print(f"Ucitavam chunk: {putanja}")
        df = ucitaj_chunk_u_dataframe(putanja)
        delovi.append(df)

    spojeni_df = pd.concat(delovi, axis=0)

    print("Spojen DataFrame. Oblik:", spojeni_df.shape)
    return spojeni_df


def sacuvaj_finalni_feather_iz_chunkova(chunk_folder, finalni_feather_putanja):
    df = spoji_chunkove_u_jedan_dataframe(chunk_folder)

    df_za_cuvanje = df.reset_index().rename(columns={"index": "varijanta"})
    df_za_cuvanje.to_feather(finalni_feather_putanja)

    print(f"\nFinalni feather sacuvan na: {finalni_feather_putanja}")
    print(f"Oblik finalnog DataFrame-a: {df.shape}")

    del df
    del df_za_cuvanje
    gc.collect()


# Pomocne funkcije

# dominantni model
# za 1 i 2 alela vraca 1
def primeni_dominantni_model(genotip_red):
    return (genotip_red >= 1).astype(np.uint8)

# recesivni model
# za 2 alela vraca 1
def primeni_recesivni_model(genotip_red):
    return (genotip_red == 2).astype(np.uint8)


def napravi_kontingentnu_tabelu(genotip_model, status_niz):
    nema_zdravi = np.sum((genotip_model == 0) & (status_niz == 0))
    nema_bolesni = np.sum((genotip_model == 0) & (status_niz == 1))
    ima_zdravi = np.sum((genotip_model == 1) & (status_niz == 0))
    ima_bolesni = np.sum((genotip_model == 1) & (status_niz == 1))

    tabela = np.array([
        [nema_zdravi, nema_bolesni],
        [ima_zdravi, ima_bolesni]
    ], dtype=np.int64)

    return tabela


def izracunaj_p_vrednost(tabela):
    tabela_korigovana = tabela + 1   #Laplacova korekcija
    chi2, p_vrednost, dof, expected = chi2_contingency(tabela_korigovana)
    return p_vrednost

def podeli_dataframe_na_blokove(genotip_df, velicina_bloka):
    blokovi = []

    broj_redova = len(genotip_df)

    for pocetak in range(0, broj_redova, velicina_bloka):
        kraj = min(pocetak + velicina_bloka, broj_redova)
        blok = genotip_df.iloc[pocetak:kraj]
        blokovi.append(blok)

    return blokovi

def obradi_blok(args):
    blok_df, status_niz = args

    rezultati = []

    for indeks, red in blok_df.iterrows():
        genotip_red = red.to_numpy(dtype=np.uint8)

        dominantni = primeni_dominantni_model(genotip_red)
        tabela_dom = napravi_kontingentnu_tabelu(dominantni, status_niz)
        p_dom = izracunaj_p_vrednost(tabela_dom)

        recesivni = primeni_recesivni_model(genotip_red)
        tabela_rec = napravi_kontingentnu_tabelu(recesivni, status_niz)
        p_rec = izracunaj_p_vrednost(tabela_rec)

        rezultati.append((indeks, p_dom, p_rec))

    rezultat_bloka = pd.DataFrame(
        rezultati,
        columns=["varijanta", "p_dominantni", "p_recesivni"]
    ).set_index("varijanta")

    return rezultat_bloka

def gwas_analiza(genotip_df, status_serija, bonferroni=True):
    """
    Funkcija vrsi GWAS analizu nad genotipskim podacima i statusom bolesti.

    Implementirani modeli:
    1. Dominantni model - genotipovi 1 i 2 se grupisu zajedno, jer je dovoljno prisustvo BAR JEDNOG alternativnog alela.
    2. Recesivni model - samo genotip 2 se smatra izdvojenom kategorijom, jer su potrebna DVA alternativna alela.

    Ova dva modela su izabrana jer predstavljaju dve osnovne i bioloski intuitivne pretpostavke o uticaju alternativnog alela na fenotip.

    Parametri:
    - genotip_df: pandas DataFrame iz zadatka 1
    - status_serija: pandas Series sa vrednostima 0 (zdrava) i 1 (obolela)
    - bonferroni: da li primeniti Bonferroni filtriranje

    Povratna vrednost:
    - pandas DataFrame sa dve kolone p-vrednosti, po jednom za svaki model, sortiran po p-vrednosti dominantnog modela
    """

    paralelno = True
    broj_procesa = 2
    velicina_bloka = 10000

    status_niz = status_serija.to_numpy(dtype=np.uint8)

    if genotip_df.shape[1] != len(status_niz):
        raise ValueError("Broj kolona genotipskog DataFrame-a mora biti jednak broju status vrednosti.")

    blokovi = podeli_dataframe_na_blokove(genotip_df, velicina_bloka)

    if paralelno:
        if broj_procesa is None:
            broj_procesa = max(1, cpu_count() - 1)

        ulazi = [(blok, status_niz) for blok in blokovi]

        with Pool(processes=broj_procesa) as pool:
            rezultati_blokova = pool.map(obradi_blok, ulazi)
    else:
        rezultati_blokova = []
        for blok in blokovi:
            rezultat_bloka = obradi_blok((blok, status_niz))
            rezultati_blokova.append(rezultat_bloka)

    rezultat_df = pd.concat(rezultati_blokova, axis=0)

    rezultat_df = rezultat_df.sort_values(by="p_dominantni")

    if bonferroni:
        n = len(rezultat_df)
        prag = 0.05 / n
        rezultat_df = rezultat_df[
            (rezultat_df["p_dominantni"] < prag) |
            (rezultat_df["p_recesivni"] < prag)
        ]

    return rezultat_df



def nacrtaj_menhetn_grafikon(rezultat_df, kolona_p, naziv_grafika):
    df_plot = rezultat_df.copy()

    # ako je index umesto varijanta
    if "index" in df_plot.columns and "varijanta" not in df_plot.columns:
        df_plot = df_plot.rename(columns={"index": "varijanta"})

    # parsiranje
    df_plot[["hromozom", "pozicija"]] = df_plot["varijanta"].str.split("-", expand=True)
    df_plot["pozicija"] = df_plot["pozicija"].astype(int)

    df_plot["hrom_broj"] = df_plot["hromozom"].str.replace("chr", "", regex=False).astype(int)

    # sortiranje
    df_plot = df_plot.sort_values(["hrom_broj", "pozicija"]).reset_index(drop=True)

    # logaritmujemo p vrednosti da bi se bolje videle na grafiku
    df_plot["log_p"] = -np.log10(df_plot[kolona_p])

    offset = 0
    offsets = {}
    sredine = {}

    for chrom in sorted(df_plot["hromozom"].unique(), key=lambda x: int(x.replace("chr",""))):
        df_chr = df_plot[df_plot["hromozom"] == chrom]

        offsets[chrom] = offset
        sredine[chrom] = offset + df_chr["pozicija"].median()

        offset += df_chr["pozicija"].max()

    # nova x osa
    df_plot["x"] = df_plot.apply(lambda row: row["pozicija"] + offsets[row["hromozom"]], axis=1)

    # Bonferroni prag
    n = len(df_plot)
    prag = 0.05 / n
    log_prag = -np.log10(prag)

    boje = {
        "chr1": "#C084FC",
        "chr2": "#5EEAD4"
    }

    plt.figure(figsize=(14, 5))

    for chrom, grupa in df_plot.groupby("hromozom", sort=False):
        plt.scatter(
            grupa["x"],
            grupa["log_p"],
            s=6,
            color=boje.get(chrom, "#64748B"),
            alpha=0.8
        )

    # Bonferroni
    plt.axhline(y=log_prag, color="#FB7185", linestyle="--", linewidth=1.8)

    # oznake hromozoma na sredini
    plt.xticks(
        list(sredine.values()),
        list(sredine.keys())
    )

    plt.xlabel("Hromozom")
    plt.ylabel("-log10(p)")
    plt.title(naziv_grafika)

    plt.tight_layout()
    plt.show()


# HMM parametri

# matrica tranzicija
A = np.array([
    [0.9999, 0.0001],
    [0.01,   0.99]
], dtype=np.float64)

# matrica emisija
E = np.array([
    [0.90, 0.07, 0.02, 0.01],
    [0.10, 0.20, 0.30, 0.40]
], dtype=np.float64)

# pocetne vrednosti
pi = np.array([0.9999, 0.0001], dtype=np.float64)

# 4 kategorije signala
def diskretizuj_signal(vrednost):
    if vrednost < 4:
        return 0
    elif vrednost < 6:
        return 1
    elif vrednost < 8:
        return 2
    else:
        return 3

def viterbi(opservacije, A, E, pi):

    opservacije = np.asarray(opservacije, dtype=np.int64)

    broj_pozicija = len(opservacije)
    broj_stanja = A.shape[0]

    # prebacujemo sve u log prostor
    logA = np.log(A)
    logE = np.log(E)
    logPi = np.log(pi)

    # cuvamo najbolje verovatnoce do trenutne pozicije iz svakog stanja
    dp = np.full((broj_pozicija, broj_stanja), -np.inf, dtype=np.float64)

    # odakle smo dosli
    prethodno = np.full((broj_pozicija, broj_stanja), -1, dtype=np.int64)

    # za prvu na osnovu pocetnih verovatnoca
    prva_opservacija = opservacije[0]
    for stanje in range(broj_stanja):
        dp[0, stanje] = logPi[stanje] + logE[stanje, prva_opservacija]


    # kroz sve opservacije
    for t in range(1, broj_pozicija):
        obs = opservacije[t]

        for trenutno_stanje in range(broj_stanja):
            najbolje_prethodno_stanje = -1
            najbolja_vrednost = -np.inf

            for prethodno_stanje in range(broj_stanja):
                kandidat = (dp[t - 1, prethodno_stanje] + logA[prethodno_stanje, trenutno_stanje] + logE[trenutno_stanje, obs])

                if kandidat > najbolja_vrednost:
                    najbolja_vrednost = kandidat
                    najbolje_prethodno_stanje = prethodno_stanje

            dp[t, trenutno_stanje] = najbolja_vrednost
            prethodno[t, trenutno_stanje] = najbolje_prethodno_stanje


    najbolja_stanja = np.full(broj_pozicija, -1, dtype=np.int64)
    najbolja_stanja[-1] = np.argmax(dp[-1])

    for t in range(broj_pozicija - 2, -1, -1):
        najbolja_stanja[t] = prethodno[t + 1, najbolja_stanja[t + 1]]

    return najbolja_stanja

def izdvoji_regione_iz_stanja(df_hromozom, stanja, kolona_p):

    df_hromozom = df_hromozom.reset_index(drop=True).copy()

    if len(df_hromozom) != len(stanja):
        raise ValueError("Dužina DataFrame-a i niza stanja mora biti ista.")

    regioni = []
    u_regionu = False
    pocetak = None

    for i, stanje in enumerate(stanja):

        if stanje == 1 and not u_regionu:
            u_regionu = True
            pocetak = i

        elif stanje == 0 and u_regionu:
            kraj = i - 1
            region_df = df_hromozom.iloc[pocetak:kraj + 1]

            # lead snp (najmanja p)
            indeks_lead = region_df[kolona_p].idxmin()
            lead_red = df_hromozom.loc[indeks_lead]

            regioni.append({
                "chrom": region_df["chrom"].iloc[0],
                "start": int(region_df["pos"].min()),
                "end": int(region_df["pos"].max()),
                "n_snps": int(len(region_df)),
                "lead_snp": lead_red["varijanta"],
                "lead_p": float(lead_red[kolona_p]),
                "snps": region_df["varijanta"].tolist()
            })

            u_regionu = False
            pocetak = None

    # ako se region zavrsava na poslednjoj poziciji
    if u_regionu:
        kraj = len(stanja) - 1
        region_df = df_hromozom.iloc[pocetak:kraj + 1]

        indeks_lead = region_df[kolona_p].idxmin()
        lead_red = df_hromozom.loc[indeks_lead]

        regioni.append({
            "chrom": region_df["chrom"].iloc[0],
            "start": int(region_df["pos"].min()),
            "end": int(region_df["pos"].max()),
            "n_snps": int(len(region_df)),
            "lead_snp": lead_red["varijanta"],
            "lead_p": float(lead_red[kolona_p]),
            "snps": region_df["varijanta"].tolist()
        })

    regioni_df = pd.DataFrame(regioni)

    if not regioni_df.empty:
        regioni_df = regioni_df.sort_values("lead_p").reset_index(drop=True)

    return regioni_df

def detektuj_regione_hmm(df, kolona_modela, A, E, pi):
    rezultat_df = df.copy()

    if "index" in rezultat_df.columns and "varijanta" not in rezultat_df.columns:
        rezultat_df = rezultat_df.rename(columns={"index": "varijanta"})

    # parsiranje kolone
    rezultat_df[["chrom", "pos"]] = rezultat_df["varijanta"].str.split("-", n=1, expand=True)
    rezultat_df["pos"] = rezultat_df["pos"].astype(int)

    # logaritmovanje p vrednosti
    eps = 1e-300
    rezultat_df[kolona_modela] = rezultat_df[kolona_modela].clip(lower=eps)
    # dodajemo u novu kolonu signal
    rezultat_df["signal"] = -np.log10(rezultat_df[kolona_modela])

    # lista opservacija (po kategorijama)
    rezultat_df["opservacija"] = rezultat_df["signal"].apply(diskretizuj_signal).astype(np.int64)

    # sortiramo po hromozomima i onda po pozicijama unutar njih
    rezultat_df = rezultat_df.sort_values(["chrom", "pos"]).reset_index(drop=True)

    # za svaki hromozom trazimo regione
    svi_regioni = []

    for chrom in rezultat_df["chrom"].unique():
        df_hromozom = rezultat_df[rezultat_df["chrom"] == chrom].copy().reset_index(drop=True)

        opservacije = df_hromozom["opservacija"].values
        stanja = viterbi(opservacije, A, E, pi)

        regioni_hromozoma = izdvoji_regione_iz_stanja(df_hromozom, stanja, kolona_modela)

        if not regioni_hromozoma.empty:
            svi_regioni.append(regioni_hromozoma)

    # spajanje svih regiona u tabelu
    if len(svi_regioni) == 0:
        return pd.DataFrame(columns=["chrom", "start", "end", "n_snps", "lead_snp", "lead_p", "snps"])

    svi_regioni_df = pd.concat(svi_regioni, ignore_index=True)
    # sortitamo po verovatnoci
    svi_regioni_df = svi_regioni_df.sort_values("lead_p").reset_index(drop=True)

    return svi_regioni_df


def nacrtaj_signal_i_regione(rezultat_df, kolona_p, regioni_df=None):
    df = rezultat_df.copy()

    if "index" in df.columns and "varijanta" not in df.columns:
        df = df.rename(columns={"index": "varijanta"})

    df[["chrom", "pos"]] = df["varijanta"].str.split("-", n=1, expand=True)
    df["pos"] = df["pos"].astype(int)

    eps = 1e-300
    df[kolona_p] = df[kolona_p].clip(lower=eps)
    df["signal"] = -np.log10(df[kolona_p])

    def chrom_sort_key(c):
        c = str(c).replace("chr", "")
        return (0, int(c)) if c.isdigit() else (1, c)

    hromozomi = sorted(df["chrom"].unique(), key=chrom_sort_key)
    df = df.sort_values(["chrom", "pos"]).reset_index(drop=True)

    ymax = df["signal"].max() + 1

    BOJA_POZADINA = "#aaaaaa"
    BOJA_SIGNAL   = "#e74c3c"
    BOJA_REGION   = "#f39c12"
    BOJA_LEAD     = "#8e44ad"
    GRANICE       = [4, 6, 8]

    for chrom in hromozomi:
        df_chr = df[df["chrom"] == chrom].copy()
        regioni_chr = regioni_df[regioni_df["chrom"] == chrom] if regioni_df is not None else None

        # HMM stanje svake tacke (S/B)
        df_chr["u_regionu"] = False
        if regioni_chr is not None and len(regioni_chr) > 0:
            for _, r in regioni_chr.iterrows():
                maska = (df_chr["pos"] >= r["start"]) & (df_chr["pos"] <= r["end"])
                df_chr.loc[maska, "u_regionu"] = True

        fig, ax = plt.subplots(figsize=(14, 5))

        for g in GRANICE:
            ax.axhline(g, color="steelblue", linewidth=0.8, linestyle="--", alpha=0.5, zorder=1)

        if regioni_chr is not None:
            for _, r in regioni_chr.iterrows():
                ax.axvspan(r["start"], r["end"], alpha=0.15, color=BOJA_REGION, zorder=2)

        # B tacke
        pozadina = df_chr[~df_chr["u_regionu"]]
        ax.scatter(pozadina["pos"], pozadina["signal"],
                   s=5, alpha=0.5, color=BOJA_POZADINA, label="Pozadina (B)", zorder=3)

        # S tacke (SIGNAL)
        signal = df_chr[df_chr["u_regionu"]]
        if len(signal) > 0:
            ax.scatter(signal["pos"], signal["signal"], s=8, alpha=0.8, color=BOJA_SIGNAL, label="Signal (S)", zorder=4)

        # lead SNP
        if regioni_chr is not None:
            for _, r in regioni_chr.iterrows():
                lead_pos_str = str(r["lead_snp"]).split("-")[-1]
                try:
                    lead_pos = int(lead_pos_str)
                    lead_tacka = df_chr[df_chr["pos"] == lead_pos]
                    if len(lead_tacka) > 0:
                        ly = lead_tacka["signal"].values[0]
                        ax.scatter(lead_pos, ly, s=30, color=BOJA_LEAD, zorder=5, label="Lead SNP" if _ == regioni_chr.index[0] else "")
                        ax.annotate(r["lead_snp"], (lead_pos, ly),
                                    textcoords="offset points", xytext=(0, 6),
                                    fontsize=7, ha="center", color=BOJA_LEAD)
                except ValueError:
                    pass

        ax.set_ylim(0, ymax)
        ax.set_xlabel("Pozicija")
        ax.set_ylabel("−log₁₀(p)")
        ax.set_title(f"Signal duž {chrom}")
        ax.legend(loc="upper right", markerscale=2, fontsize=8)
        plt.tight_layout()
        plt.show()

    # zoom grafikoni za svaki region
    if regioni_df is not None:
        for i, r in regioni_df.iterrows():
            prosirenje = 50_000

            podskup = df[
                (df["chrom"] == r["chrom"]) &
                (df["pos"] >= r["start"] - prosirenje) &
                (df["pos"] <= r["end"] + prosirenje)
            ].copy()

            podskup["u_regionu"] = ((podskup["pos"] >= r["start"]) & (podskup["pos"] <= r["end"]))

            fig, ax = plt.subplots(figsize=(10, 4))

            for g in GRANICE:
                ax.axhline(g, color="steelblue", linewidth=0.8, linestyle="--", alpha=0.5)

            ax.axvspan(r["start"], r["end"], alpha=0.15, color=BOJA_REGION)

            ax.scatter(podskup[~podskup["u_regionu"]]["pos"], podskup[~podskup["u_regionu"]]["signal"], s=10, color=BOJA_POZADINA, alpha=0.6, label="Pozadina (B)")

            ax.scatter(podskup[podskup["u_regionu"]]["pos"], podskup[podskup["u_regionu"]]["signal"], s=14, color=BOJA_SIGNAL, alpha=0.9, label="Signal (S)")

            # lead SNP
            lead_pos_str = str(r["lead_snp"]).split("-")[-1]
            try:
                lead_pos = int(lead_pos_str)
                lead_tacka = podskup[podskup["pos"] == lead_pos]
                if len(lead_tacka) > 0:
                    ly = lead_tacka["signal"].values[0]
                    ax.scatter(lead_pos, ly, s=40, color=BOJA_LEAD, zorder=5, label="Lead SNP")
                    ax.annotate(r["lead_snp"], (lead_pos, ly),
                                textcoords="offset points", xytext=(0, 7),
                                fontsize=8, ha="center", color=BOJA_LEAD)
            except ValueError:
                pass

            ax.set_title(f"Zoom — Region {i+1}: {r['chrom']} ({r['start']:,} – {r['end']:,})"
                         f"\nLead SNP: {r['lead_snp']}  (p = {r['lead_p']:.2e})")
            ax.set_xlabel("Pozicija")
            ax.set_ylabel("−log₁₀(p)")
            ax.legend(loc="upper right", fontsize=8)
            plt.tight_layout()
            plt.show()


def main():

    # 1. ZADATAK
    ucitaj_i_sacuvaj_sve_hromozome_u_chunkove(
        bcf_putanja=bcf_putanja,
        chunk_folder=chunk_folder,
        chunk_size=10000,
        korak_ispisa=20000
    )

    df = spoji_chunkove_u_jedan_dataframe(chunk_folder)

    print("Oblik genotipskog DataFrame-a:", df.shape)

    # 2. ZADATAK
    status_df = pd.read_csv(status_putanja, sep='\t')
    status_serija = status_df["disease"]

    rezultat_df = gwas_analiza(df, status_serija, bonferroni=False)

    # provere
    print("Oblik GWAS rezultata:", rezultat_df.shape)
    print("Kolone:", rezultat_df.columns.tolist())
    print("Jedinstven broj varijanti:", rezultat_df.index.nunique())
    print("Broj NaN vrednosti:")
    print(rezultat_df.isna().sum())
    print("p_dominantni u [0,1]:",
          (rezultat_df["p_dominantni"].between(0, 1)).all())
    print("p_recesivni u [0,1]:",
          (rezultat_df["p_recesivni"].between(0, 1)).all())
    
    nacrtaj_menhetn_grafikon(rezultat_df.reset_index(), "p_dominantni", "Menhetn grafikon - dominantni model")
    nacrtaj_menhetn_grafikon(rezultat_df.reset_index(), "p_recesivni", "Menhetn grafikon - recesivni model")

    # 3. ZADATAK
    kolona_modela = "p_recesivni"

    regioni_df = detektuj_regione_hmm(rezultat_df.reset_index(), kolona_modela, A, E, pi)

    print(regioni_df)
    print("Ukupan broj regiona:", len(regioni_df))

    nacrtaj_signal_i_regione(rezultat_df.reset_index(), kolona_modela, regioni_df)


if __name__ == "__main__":
    main()