HOW TO CALCULATE OVERLAP IN DADA

When working with paired-end sequencing, the forward ($R1$) and reverse ($R2$) reads start from opposite ends of the exact same DNA fragment (amplicon) and read toward each other. To stitch them back into a single continuous sequence, DADA2 looks for a identical sequence region where the two reads meet: the **overlap**.

---

### The Fundamental Equation

To calculate how much overlap you have left after trimming your reads, use this simple formula:

$$\text{Overlap} = (\text{Truncated } R1 + \text{Truncated } R2) - \text{Amplicon Length}$$

Where:

* **Truncated $R1$**: `--p-trunc-len-f` (minus any `--p-trim-left-f`)
* **Truncated $R2$**: `--p-trunc-len-r` (minus any `--p-trim-left-r`)
* **Amplicon Length**: The total physical length of the targeted target genomic region (excluding primers, if you already stripped them with Cutadapt).

---

### Why the Math Matters (3 Scenarios)

**1. Adequate Overlap ($\ge 20\text{ bp}$)**

* **Math:** $250\text{ bp } (R1) + 230\text{ bp } (R2) = 480\text{ bp}$.
* **Overlap:** $480 - 460 = 20\text{ bp}$.
* **Outcome:** DADA2 finds the identical $20\text{ bp}$ sequence shared between $R1$ and $R2$, merges them successfully, and preserves your sample throughput.

**2. Zero or Negative Overlap ($\le 0\text{ bp}$)**

* **Math:** $240\text{ bp } (R1) + 210\text{ bp } (R2) = 450\text{ bp}$.
* **Overlap:** $450 - 460 = -10\text{ bp}$.
* **Outcome:** The reads physically do not reach each other. There is a $10\text{ bp}$ gap of unsequenced DNA in the middle. DADA2 cannot stitch them together and **drops 100% of the reads** during the merge step, leaving you with an empty feature table.

**3. Insufficient Overlap ($< 12\text{ bp}$ default)**

* **Math:** $240\text{ bp } (R1) + 230\text{ bp } (R2) = 470\text{ bp}$.
* **Overlap:** $470 - 460 = 10\text{ bp}$.
* **Outcome:** Although they overlap slightly, DADA2 requires a default minimum overlap of $12\text{ bp}$ (or your custom target of $20\text{ bp}$ for safety) to prevent false-positive merges. Most reads will be discarded.

---

### Rule of Thumb across Regions

| Hypervariable Region | Avg. Amplicon Length | Minimum Combined ($R1 + R2$) Needed for $\approx 20\text{ bp}$ Overlap |
| --- | --- | --- |
| **V4** | $\approx 250\text{ bp}$ | $\ge 270\text{ bp}$ (Easily achieved with $2 \times 150\text{ bp}$ or $2 \times 250\text{ bp}$) |
| **V3-V4** | $\approx 460\text{ bp}$ | $\ge 480\text{ bp}$ (Requires aggressive retention: e.g., $250 + 230\text{ bp}$) |
| **V1-V3** | $\approx 490\text{ bp}$ | $\ge 510\text{ bp}$ (Extremely tight on standard $2 \times 300\text{ bp}$ Illumina runs) |

---

### The Quality vs. Overlap Trade-off

Because reverse reads ($R2$) typically drop in quality faster than forward reads ($R1$), it is tempting to truncate $R2$ aggressively (e.g., down to $150\text{ bp}$).

* **For short amplicons (V4):** You can trunc $R2$ heavily because $200 + 150 = 350\text{ bp}$, which is way higher than the required $270\text{ bp}$.
* **For long amplicons (V3-V4):** You **cannot** truncate $R2$ aggressively. You must accept a slightly lower quality score at the 3' end of $R2$ to preserve enough physical length so the math holds true ($\ge 480\text{ bp}$).



### ----------------------------------------------------------------------------
### ----------------------------------------------------------------------------
### ----------------------------------------------------------------------------
## CALCULUS FOR HAIMOV
### ----------------------------------------------------------------------------
### ----------------------------------------------------------------------------

Looking at your quality profiles for **Haimov 2022**, the reads drop off sharply in quality much earlier than standard Illumina runs.

Here is the breakdown of the quality scores, the overlap math, and the recommended truncation parameters depending on which **16S hypervariable region** this dataset targeted.

---

### **1. Quality Score Diagnostics**

* **Forward Reads ($R1$ - Left Plot):**
* Quality remains high ($Q > 30$) until roughly **$\approx 200\text{ bp}$**.
* After $200\text{ bp}$, quality drops into the $Q20-Q25$ range, and precipitously degrades past **$240\text{ bp}$**.


* **Reverse Reads ($R2$ - Right Plot):**
* Quality drops significantly much earlier.
* Beyond **$\approx 150\text{ bp}$**, the 25th percentile dips down to $Q15$, and past **$220\text{ bp}$**, it hits $Q12-Q15$ across the board.



---

### **2. Scenario Analysis & Overlap Math**

To make sure DADA2 doesn't drop 100% of your reads during the merge step, your truncation parameters must be calculated against the target amplicon size:

$$\text{Overlap} = (\text{Truncated } R1 + \text{Truncated } R2) - \text{Amplicon Length} \ge 20\text{ bp}$$

#### **Scenario A: Targeted Region is V4 ($\approx 252\text{ bp}$ amplicon)**

If Haimov targeted the V4 region (e.g., $515\text{F} / 806\text{R}$):

* **Recommended Truncation:**
* `--p-trunc-len-f`: **200**
* `--p-trunc-len-r`: **150**


* **Math Check:**

$$\text{Combined Length} = 200 + 150 = 350\text{ bp}$$


$$\text{Overlap} = 350\text{ bp} - 252\text{ bp} = \mathbf{98\text{ bp}}$$


* **Verdict:** **PERFECT.** You have massive overlap ($98\text{ bp} \gg 20\text{ bp}$). Truncating $R2$ aggressively at $150\text{ bp}$ cuts out all the low-quality noise while guaranteeing effortless merging.

---

#### **Scenario B: Targeted Region is V3-V4 ($\approx 460\text{ bp}$ amplicon)**

If Haimov targeted the V3-V4 region (e.g., $341\text{F} / 805\text{R}$):

* **Required Combined Length:** $460\text{ bp} + 20\text{ bp} = \mathbf{480\text{ bp}}$.
* **The Problem:** To reach $480\text{ bp}$, you would need parameters like `--p-trunc-len-f 240` and `--p-trunc-len-r 240`. But looking at the $R2$ plot, the quality at $240\text{ bp}$ is atrocious ($Q \approx 13$). DADA2 will throw out those low-quality reads during error modeling.
* **Recommended Truncation (Compromise):**
* `--p-trunc-len-f`: **240**
* `--p-trunc-len-r`: **240**


* **Math Check:**

$$\text{Combined Length} = 240 + 240 = 480\text{ bp}$$


$$\text{Overlap} = 480\text{ bp} - 460\text{ bp} = \mathbf{20\text{ bp}}$$


* **Verdict:** **RISKY.** If Haimov is V3-V4, this dataset will have low retention after DADA2 due to the poor $R2$ quality tail. *(Note: If read retention is under $10\%$, you may need to process Haimov as single-end using only $R1$ truncated at $230\text{ bp}$).*

---

### **3. Config TSV Recommendation**

Assuming **Haimov 2022** targets the **V4 region** (which matches your initial array values `TRUNC_LEN_F=200` and `TRUNC_LEN_R=150` in `05_dada2_denoise_datasets.sh`), your line in `0_preprocessing/config/datasets_config.tsv` should be:

| study_id | fwd_primer | rev_primer | trunc_len_f | trunc_len_r | trim_left_f | trim_left_r | region |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **haimov_2022** | *(your_fwd)* | *(your_rev)* | **200** | **150** | **0** | **0** | **V4** |




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
####        LIUB_2019
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

La calidad general de **`liub_2019`** es alta (la mediana de Phred se mantiene holgadamente sobre Q35-Q40 a lo largo de casi toda la lectura).

Sin embargo, hay dos aspectos clave a considerar para configurar DADA2:

1. **Caída puntual alrededor de la base ~260-280:** Hay una ligera bajada en el percentil inferior, pero la mediana (la línea negra central) sigue estando muy por encima de Q30.
2. **Longitud mínima variable:** El texto indica que la longitud mínima encontrada en el submuestreo es de **221 bases**. Si aplicas un truncamiento rígido por encima de 221 pb (ej. `--p-trunc-len 250`), DADA2 descartará por completo todas las lecturas cuya longitud sea menor a ese límite.

---

### **Parámetros Recomendados para DADA2 Single-End**

Usa `qiime dada2 denoise-single` con la siguiente configuración:

* **`--p-trim-left 0`**: Los primers ya fueron removidos eficientemente por Cutadapt.
* **`--p-trunc-len 220`**: Cortar a 220 pb retiene más del 99% de las lecturas evitando que se descarten secuencias cortas, al tiempo que elimina la caída de calidad del extremo 3' (bases 260-280).

```bash
qiime dada2 denoise-single \
  --i-demultiplexed-seqs liub_2019/demux_trimmed.qza \
  --p-trim-left 0 \
  --p-trunc-len 220 \
  --p-n-threads 4 \
  --o-table liub_2019/table.qza \
  --o-representative-sequences liub_2019/rep-seqs.qza \
  --o-denoising-stats liub_2019/stats.qza

```


 
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ###
# ###
# ###
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

