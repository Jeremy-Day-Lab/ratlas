# The Ratlas datasets

The Ratlas single-nucleus RNA-seq (snRNA-seq) datasets are currently composed of the following:

__1) The rat nucleus accumbens (NAc), with samples taken from four experimental groups:__

* Male adult rats, 1hr following saline injection
* Male adult rats, 1hr following cocaine injection (20mg/kg, I.P.),
* Female adult rats, 1hr following first saline injection
* Female adult rats, 1hr following first cocaine injection (20mg/kg, I.P.)
<br>

<img src="Adult_NAc_snRNA-seq.jpg" align="center" width="850px" />

<br>
<br>

__2) The rat nucleus accumbens (NAc), with samples taken from two experiments and eight experimental groups:__

&nbsp;&nbsp;&nbsp;&nbsp; **Acute** - Rats received a single dose of saline or cocaine

* Male adult rats, 1hr following first saline injection
* Male adult rats, 1hr following first cocaine injection (20mg/kg, I.P.)
* Female adult rats, 1hr following first saline injection
* Female adult rats, 1hr following first cocaine injection (20mg/kg, I.P.)

&nbsp;&nbsp;&nbsp;&nbsp; **Repeated** - Rats received seven consecutive doses of saline or cocaine

* Male adult rats, 1hr following seventh saline injection
* Male adult rats, 1hr following seventh cocaine injection (20mg/kg, I.P.)
* Female adult rats, 1hr following seventh saline injection
* Female adult rats, 1hr following seventh cocaine injection (20mg/kg, I.P.)
<br>

<img src="Acute_Repeated_Adult_NAc_snRNA-seq.jpg" align="center" width="850px" />

<br>
<br>

__3) Primary striatal neurons (mixed from male and female E18 rat brains and cultured to DIV11) from four experimental groups:__

* Vehicle (media alone, 1hr)
* Dopamine (50µM, 1hr)
* SKF-38393 (1µM, 1hr)
* Potassium chloride (25mM, 1hr)
<br>

<img src="Primary_culture_snRNA-seq.jpg" align="center" width="850px" />

<br>
<br>

__4) Adult rat ventral tegmental area (VTA) from naive male and female rats.__

<br>

<img src="Adult_VTA_snRNA-seq.jpg" align="center" width="850px" />

<br>

----
### Citation:

Savell, K.E.\*, Tuscher, J.J.\*, Zipperly, M.E\*, Duke, C.G.\*, Phillips III, R.A.\*, Bauman, A.J., Thukral, S., Sultan, F.A., Goska, N.A, Ianov, L. & Day, J.J. (Science Advances, June, 2020). [_A dopamine-induced gene expression signature regulates neuronal function and cocaine response_](https://advances.sciencemag.org/content/6/26/eaba4221)  DOI: 10.1126/sciadv.aba4221

Phillips III, R.A.\* Tuscher, J.J.\*, Black, S.L., Andraka E., Fitzgerald, N.D., Ianov, L., & Day, J.J. (Cell Reports, April, 2022). [_An atlas of transcriptionally defined cell populations in the rat ventral tegmental area._](https://www.cell.com/cell-reports/fulltext/S2211-1247%2822%2900364-3)  DOI: <https://doi.org/10.1016/j.celrep.2022.110616>

Phillips III, R.A.\* Tuscher, J.J., Wan E., Fitzgerald, N.D., Zipperly, M.E, Duke, C.G., Ianov, L., & Day, J.J. (Molecular and Cellular Neuroscience, June, 2023). [_Distinct subpopulations of D1 medium spiny neurons exhibit unique transcriptional responsiveness to cocaine_](https://doi.org/10.1016/j.mcn.2023.103849)  DOI: <https://doi.org/10.1016/j.mcn.2023.103849>

### Links:

Direct link to analytical code for the work cited above (including VTA Dockerfiles) can be found at <https://github.com/Jeremy-Day-Lab>

All Day lab resources may be found at the [Day Lab website](http://day-lab.org/resources)

### FAQ:

1. Why can't I find my gene?

	Users of the app should be aware of the feature to search by gene IDs for cases when a gene name is not assigned during the generation of the cell matrices (gene ids are typically associated with more novel genes in the rat assembly. You will normally find either a gene name or ID, not both). Please see the "Gene name search" help icon at the upper right corner of the `Choose a gene` option for more information including reference genome source and versions. We add special emphasis on this section since the addition of the rn7 mapping.
	
2. What is the difference between `rn6` and `rn7` tabs?

	The difference between `rn6` vs `rn7` in the first three datasets published in Ratlas (adult acute NAc, primary striatal neurons and VTA) is the Rat assembly. The rat `rn6` assembly was implemented in the original published data across these three datasets, and later, we have re-mapped the data and updated the objects to contain the updated rat `rn7` assembly. [This Twitter thread, written by Dr. Jeremy Day](https://twitter.com/DayLabUAB/status/1542635405542957058), nicely summarizes some of the key differences we have found between the assemblies.
	
	Notably, newer datasets published since then (adult acute and repeated NAc), have only been mapped with the newer rn7 assembly.

### About the app:

This app was developed and is actively maintained by Lara Ianov, Ph.D., Managing Director of the [UAB Biological Data Science Core](https://www.uab.edu/cores/ircp/bds), and bioinformatics specialist for the Civitan International Research Center, at the University of Alabama at Birmingham. For questions about the app, please send an email to lianov@uab.edu
