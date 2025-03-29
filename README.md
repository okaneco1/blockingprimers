# Blocking Primer Effectiveness for Dietary DNA Metabarcoding of Sea Lamprey

This repository accompanies the manuscript (currently in review):

**"Development of PCR Blocking Primers Enabling DNA Metabarcoding Analysis of Dietary Composition in Hematophagous Sea Lamprey"**  
*Conor O’Kane, Nicholas S. Johnson, Kim T. Scribner, Jeannette Kanefsky, Weiming Li, John D. Robinson*

---

## Overview

This study employs molecular techniques to detect trace amounts of host DNA in sea lamprey gut contents using universal primers that target the mitochondrial 12S rRNA gene. However, a major technical challenge with universal primers is the overwhelming amplification of sea lamprey DNA, which can overshadow host sequences in sequencing outputs. To address this, we developed and tested **eight custom PCR blocking primers** specifically designed to suppress sea lamprey 12S amplification during PCR, while still allowing amplification of prey fish DNA.

These primers vary in sequence length (34bp vs 36bp), end modification (C3 spacer vs inverted dT), and purification method (standard vs HPLC). Their effectiveness was evaluated using quantitative PCR (qPCR), gel electrophoresis, and high-throughput sequencing on both mock DNA communities and wild-caught sea lamprey digestive samples.

This repository contains the data, analysis scripts, and documentation used to:

- Quantify the effectiveness of each blocking primer
- Compare 25 vs 40 PCR cycle outputs
- Extract and examine OTU composition per sample
- Visualize sequence similarity between host species for primer design

These tools are intended to facilitate future studies applying blocking primers in dietary DNA metabarcoding, particularly when developing novel blocking primers and assessing the balance between non-target DNA amplificaion and target DNA suppression.

