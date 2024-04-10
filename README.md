# Ratlas

`Ratlas` is an R Shiny web application which hosts the Day lab published single-nuclei datasets. **You may use the `Ratlas` app by visiting the following URL: <https://day-lab.shinyapps.io/ratlas/>.** All information regarding the currently hosted datasets can be found at the app's homepage.

## FAQs:

<details> 

<summary> 1. Why can't I find my gene? </summary>

Users of the app should be aware of the feature to search by gene IDs for cases when a gene name is not assigned during the generation of the cell matrices (gene ids are typically associated with more novel genes in the rat assembly. You will normally find either a gene name or ID, not both). Please see the "Gene name search" help icon at the upper right corner of the `Choose a gene` option for more information including reference genome source and versions. We add special emphasis on this section since the addition of the rn7 mapping.
</details> 

<details> 
<summary> 2. What is the difference between `rn6` and `rn7` tabs?  </summary>

The difference between `rn6` vs `rn7` in the first three datasets published in Ratlas (adult acute NAc, primary striatal neurons and VTA) is the Rat assembly. The rat `rn6` assembly was implemented in the original published data across these three datasets, and later, we have re-mapped the data and updated the objects to contain the updated rat `rn7` assembly. [This Twitter thread, written by Dr. Jeremy Day](https://twitter.com/DayLabUAB/status/1542635405542957058), nicely summarizes some of the key differences we have found between the assemblies.

Notably, newer datasets published since then (adult acute and repeated NAc), have only been mapped with the newer rn7 assembly.

</details> 

## Credit to data sources

The data to the following publication are hosted in `Ratlas`:

Savell, K.E.\*, Tuscher, J.J.\*, Zipperly, M.E\*, Duke, C.G.\*, Phillips III, R.A.\*, Bauman, A.J., Thukral, S., Sultan, F.A., Goska, N.A, Ianov, L. & Day, J.J. (Science Advances, June, 2020). [_A dopamine-induced gene expression signature regulates neuronal function and cocaine response_](https://advances.sciencemag.org/content/6/26/eaba4221)  DOI: 10.1126/sciadv.aba4221

Phillips III, R.A.\* Tuscher, J.J.\*, Black, S.L., Andraka E., Fitzgerald, N.D., Ianov, L., & Day, J.J. (Cell Reports, April, 2022). [_An atlas of transcriptionally defined cell populations in the rat ventral tegmental area._](https://www.cell.com/cell-reports/fulltext/S2211-1247%2822%2900364-3)  DOI: <https://doi.org/10.1016/j.celrep.2022.110616>

Phillips III, R.A.\* Tuscher, J.J., Wan E., Fitzgerald, N.D., Zipperly, M.E, Duke, C.G., Ianov, L., & Day, J.J. (Molecular and Cellular Neuroscience, June, 2023). [_Distinct subpopulations of D1 medium spiny neurons exhibit unique transcriptional responsiveness to cocaine_](https://doi.org/10.1016/j.mcn.2023.103849)  DOI: <https://doi.org/10.1016/j.mcn.2023.103849>

## Citation to Ratlas

If this application benefits your work, we kindly ask to acknowledge the app by including the following DOI: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10957420.svg)](https://doi.org/10.5281/zenodo.10957420)