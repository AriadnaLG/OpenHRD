# OpenHRD
OpenHRD is an R scrip that at incorporates open-source packages combined with a custom-made script to calculate genomic instability scars from SNP microarrays data. Using Affimetrix OncoScan “CEL” files as unique input.
The packages used included EaCoN (Easy Copy Number) (Job, 2021) (for normalization and segmentation, ASCAT (Allele-Specific Copy Number Analysis of Tumors) (Ross et al., 2021; Van Loo et al., 2010) for window segmentation, ploidy, and allele-specific copy number profiles; and scarHRD (Sztupinszki et al., 2018) for the determination of genomic instability scars.

### Genomic instability scars
+ **Telomeric allelic imbalance - TAI:** Number of regions with allelic imbalance that extend to the telomeric region (Birkbak et al., 2012).
+ **Large scale state transitions - LST:** Number of chromosomal breaks with a minimum size of 10 Mb and a distance between these breaks of 3 Mb (Popova et al., 2012).
+ **Loss of heterozigosity - LOH:** Number of regions with any type of loss of heterozygosity larger than 15 Mb (Abkevich et al., 2012).

The use of SNP arrays also allows for the assessment of copy number neutrality. The sum of the scars constituted the final GIS score.

### *For scientific use only 


# References
Abkevich, V., Timms, K. M., Hennessy, B. T., Potter, J., Carey, M. S., Meyer, L. A., Smith-McCune, K., Broaddus, R., Lu, K. H., Chen, J., Tran, T. V., Williams, D., Iliev, D., Jammulapati, S., FitzGerald, L. M., Krivak, T., DeLoia, J. A., Gutin, A., Mills, G. B., & Lanchbury, J. S. (2012). Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. British Journal of Cancer, 107(10), 1776–1782. https://doi.org/10.1038/bjc.2012.451 

Birkbak, N. J., Wang, Z. C., Kim, J.-Y., Eklund, A. C., Li, Q., Tian, R., Bowman-Colin, C., Li, Y., Greene-Colozzi, A., Iglehart, J. D., Tung, N., Ryan, P. D., Garber, J. E., Silver, D. P., Szallasi, Z., & Richardson, A. L. (2012). Telomeric Allelic Imbalance Indicates Defective DNA Repair and Sensitivity to DNA-Damaging Agents. Cancer Discovery, 2(4), 366–375. https://doi.org/10.1158/2159-8290.CD-11-0206

Job, B. (2021). EaCoN : Easy Copy Number ! (Version 0.3.6-2) [Computer software]. https://github.com/gustaveroussy/EaCoN

Popova, T., Manié, E., Rieunier, G., Caux-Moncoutier, V., Tirapo, C., Dubois, T., Delattre, O., Sigal-Zafrani, B., Bollet, M., Longy, M., Houdayer, C., Sastre-Garau, X., Vincent-Salomon, A., Stoppa-Lyonnet, D., & Stern, M.-H. (2012). Ploidy and Large-Scale Genomic Instability Consistently Identify Basal-like Breast Carcinomas with BRCA1/2 Inactivation. Cancer Research, 72(21), 5454–5462. https://doi.org/10.1158/0008-5472.CAN-12-1470

Ross, E. M., Haase, K., Van Loo, P., & Markowetz, F. (2021). Allele-specific multi-sample copy number segmentation in ASCAT. Bioinformatics, 37(13), 1909–1911. https://doi.org/10.1093/bioinformatics/btaa538

Sztupinszki, Z., Diossy, M., Krzystanek, M., Reiniger, L., Csabai, I., Favero, F., Birkbak, N. J., Eklund, A. C., Syed, A., & Szallasi, Z. (2018). Migrating the SNP array-based homologous recombination deficiency measures to next generation sequencing data of breast cancer. Npj Breast Cancer, 4(1), 16. https://doi.org/10.1038/s41523-018-0066-6

Van Loo, P., Nordgard, S. H., Lingjærde, O. C., Russnes, H. G., Rye, I. H., Sun, W., Weigman, V. J., Marynen, P., Zetterberg, A., Naume, B., Perou, C. M., Børresen-Dale, A.-L., & Kristensen, V. N. (2010). Allele-specific copy number analysis of tumors. Proceedings of the National Academy of Sciences of the United States of America, 107(39), 16910–16915. https://doi.org/10.1073/pnas.1009843107
