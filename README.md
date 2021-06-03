# MetaboMgITC

## Metabolite binding to divalent Mg cations at biologically relevant pH and ionic strength

### Authors: Jacob P. Sieg, Ryota Yamagami, and Philip C. Bevilacqua

## 1. What is MetaboMgITC?

#### RNA regulates myriad cellular events such as transcription, translation, and splicing. To perform these essential functions, RNA folds into complex tertiary structures where a negatively charged ribose-phosphate backbone interacts with metal ions. Magnesium, the most abundant divalent metal ion in cells, neutralizes the backbone thereby playing essential roles in RNA folding and function. In the cell, some Mg2+ ions are chelated by metabolites. Recently, we have demonstrated that this metabolite chelated magnesium (MCM) pool can effect RNA folding and function, depending on the strength of the interaction between Mg2+ and metabolites. The binding of Mg2+ to metabolites is highly dependent on pH and ionic strength, which changes the population of Mg2+ binding competent protonation states for a metabolite ligand and changes the activity o ions in aqueose solution respectivly, further complicating the actual Mg2+ status in cells. in an effort towards gaining insight into the role of biologically chelated magnesium ions in RNA chemistry and biology, we have:

#### (1) Currated a database of 292 absolute metabolite concentration from E. coli, Yeast, and mouse iMBK cells.
#### (2) Currated a database of pKa's, Mg2+ binding constants, ionic strengths, and temperatures for 120 metabolite ligands. This data, 582 constants in total, is used to calculate apparant metabolite/Mg2+ dissacociation constants at specific pH's and ionic strengths.
#### (3) R functions to provide a user with a relatively facile way to acess and interperate the data. Most importantly, we provide "Kd.app.calc", which enables the user to calculate apparent metabolite/Mg2+ dissacociation constants at specific pH's and ionic strengths.
#### (4) Isothermal titration calorimetry (ITC) data analysis tools to experimentally determine apparent metabolite/Mg2+ dissacociation constants at specific pH's and ionic strengths.

#### An example implementation of MetaboMgITC is dementrated below, where MCM levels are projected for E. coli, Yeast, and mouse iMBK cells at pH 7.5 and ionic strength = 0.15 M.

![Final](https://user-images.githubusercontent.com/63312483/120546207-96225600-c3bd-11eb-941d-49a7ace89fa9.png)

# 2. Using MetaboMgITC on your own console

## Dependencies

#### R version version 4.1.0
#### devtools (R package)
#### pracma (R package)
#### dplyr (R package)

## Recomended packages and programs

#### R studio
#### tidyverse (R package)

## Video Tutorials

#### 1.) Reconstituting MetaboMgITC on your own console
#### 2.) Using Kd.app.calc
#### 3.) ITC data analysis

# 3. Methods

## 3.1 Apparant dissacociation constant approximation

#### The apparent disassociation constant (KD') for a metal ion  binding to a metabolite at a pH and ionic strength is:



#### Where [M] is the concentration of metal ions. [L] and [ML] are the sum of the concentration of all protonation states for the metabolite and the metabolite magnesium complex respectively. The KD' for each ligand at a given ionic strength and pH was calculated using absolute metal ion binding constants for a given ligand protonation state and protonation constants (pKa’s) from the literature using equation 2.

## 3.2 Calculating MCM levels in cells

## 3.3 ITC data collection and analysis

#### Isothermal Titration Calorimetry. Buffers were prepared by dissolving high purity (>99.99%) salts in purified (18 mΩ) water in volumetric flasks to the final concentrations described in Supplementary Table 1. A standard 1 M (±0.01) magnesium chloride solution from Millipore Sigma was diluted to ensure accurate magnesium chloride concentrations. Samples were degassed at 25°C using a ThermoVac (MicroCal, LLC) before loading into a VP-ITC MicroCalorimeter (MicroCal, LLC) according to the manufactures recommendations with pure (18 mΩ) water in the reference cell. Titration was performed at 25°C with a 10 µcal/sec reference power and a stirring speed of 310 rpm. 29 total injections (one 2 μL injection followed by twenty eight 10 μL injections) were performed at an injection rate of 0.5 μL/sec with a 200 sec space following each injection. The first injection was not included in subsequent analysis.
#### Data were analyzed in R (version 4.1.0). Briefly, raw ITC file were parsed using a custom function. Heats of injections were determined by integrating the raw differential power curves with the polyarea function in pracma (https://CRAN.R-project.org/package=pracma). Curves were fit to a 1:1 Wisman Isotherm1 using the nls non-linear regression function in base R 4.1.0 to determine apparent association and disassociation constants. Data, standalone functions, and analysis are available on Github with help documentation for reconstitution on any console (https://github.com/JPSieg/MetaboMgITC).

# 4. References



