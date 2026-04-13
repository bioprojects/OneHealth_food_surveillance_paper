
# Source data and custom codes for the One Health genomic surveillance paper 
This page includes input files to recreate some figures (not included in the Source Data file), and explains custom codes used in the paper "One Health genomic surveillance identifies potential genetic links of antimicrobial-resistant food-borne pathogens across sectors" by Norikazu Kitamura#, Hiroto Shinomiya#, Michiko Kawanishi#, Noriko Konishi#, Kenji Ohya#, Haruyoshi Tomita#, Yukiko Asano, Chizuko Kimura, Hideto Sekiguchi, Manao Ozawa, Kanako Ishihara, Jun Kurushima, Koichi Tanimoto, Takahiro Nomura, Yusuke Hashimoto, Foodborne AMR Research Group of Public Health Institutes, Tetsuo Asai, Haruo Watanabe, Yo Sugawara, Koji Yahara#, and Motoyuki Sugai.

The following figures and tables were created using the input files and the custom codes as explained below.
- Fig. 1: Phylogenetic tree constructed from core genomes of 1,613 non-typhoidal *Salmonella enterica* isolates.

- Fig. 3: Core genome trees highlighting closely related isolate pairs detected in specific serovars, Blockley and Agona.
  
- Fig. 4: Phylogenetic trees constructed from core genomes of *Campylobacter* species isolates from human, food, and animal.
  
- Supplementary Fig. 3: Core genome phylogenetic tree of *S*. Heldeiberg (N=248), including publicly available genomes lacking detailed sample origin information.
  
- Supplementary Fig. 4: Core genome phylogenetic tree of *S*. Schwarzengrund (N=853), including publicly available genomes lacking detailed sample origin information
  
- Supplementary Fig. 5: Core genome phylogenetic tree of *S*. Infantis (N=1061), including publicly available genomes lacking detailed sample origin information.
  
- Supplementary Fig. 8: Expanded views (right) of parts of the core genome tree (left, same as fig. 4), where food-derived isolates (collected in 2023 or 2024) and human-derived isolates (collected in 2022 or 2023) exhibit lower divergence.

- Supplementary Fig. 9: Core genome tree of 38 *Enterococcus* spp. isolates derived from humans and food.

- Table 2:  *Salmonella* pairwise SNPs within each serovar

- Table 3: *Campylobacter* pairwise SNPs within each clonal complex or sequence type

  

This page explains the custome codes according to the guildeline "For manuscripts utilizing custome algorithms or software that are central to the research but not yet described in published literature" of Nature Research.

This page is written based on "Code and Software Submission Checklist" of Nature Research.

## System requirements 
- R and Rstudio

  - packages
    - None
  





#### Any required non-standard hardware

None for the custom codes.



## Instruction for use

Please go to the directory of each figure and table, and upload the files or execute the code.


------
##### Fig. 1: Phylogenetic tree constructed from core genomes of 1,613 non-typhoidal *Salmonella enterica* isolates.

- The figure was created using Phandango by uploading nwk and csv files  (https://jameshadfield.github.io/phandango/#/).
  - (expected run time: within a few seconds)


##### Fig. 3: Core genome trees highlighting closely related isolate pairs detected in specific serovars, Blockley and Agona.

- The figure was created using Phandango by uploading nwk and csv files (https://jameshadfield.github.io/phandango/#/).
  - (expected run time: within a few seconds)

##### Fig. 4: Phylogenetic trees constructed from core genomes of *Campylobacter* species isolates from human, food, and animal.

- The figure was created using Phandango by uploading nwk and csv files (https://jameshadfield.github.io/phandango/#/).
  - (expected run time: within a few seconds)

##### Supplementary Fig. 3: Core genome phylogenetic tree of *S*. Heidelberg (N=248), including publicly available genomes lacking detailed sample origin information.

- The figure was created using Phandango by uploading nwk and csv files (https://jameshadfield.github.io/phandango/#/).
  - (expected run time: within a few seconds)

##### Supplementary Fig. 4: Core genome phylogenetic tree of *S*. Schwarzengrund (N=853), including publicly available genomes lacking detailed sample origin information

- The figure was created using Phandango by uploading nwk and csv files (https://jameshadfield.github.io/phandango/#/).
  - (expected run time: within a few seconds)

##### Supplementary Fig. 5: Core genome phylogenetic tree of *S*. Infantis (N=1061), including publicly available genomes lacking detailed sample origin information.

- The figure was created using Phandango by uploading nwk and csv files (https://jameshadfield.github.io/phandango/#/).
  - (expected run time: within a few seconds)


##### Supplementary Fig. 8: Expanded views (right) of parts of the core genome tree (left, same as fig. 4), where food-derived isolates (collected in 2023 or 2024) and human-derived isolates (collected in 2022 or 2023) exhibit lower divergence. 

- The figure was created using Phandango by uploading nwk and csv files (https://jameshadfield.github.io/phandango/#/).
  - (expected run time: within a few seconds)

##### Supplementary Fig. 9: Core genome tree of 38 *Enterococcus* spp. isolates derived from humans and food.

- The figure was created using Phandango by uploading nwk and csv files (https://jameshadfield.github.io/phandango/#/).
  - (expected run time: within a few seconds)

##### Table 2: *Salmonella* pairwise SNPs within each serovar

- Table2.R
  - (expected run time: within a few seconds)

##### Table 3: *Campylobacter* pairwise SNPs within each clonal complex or sequence type

- Table3.R
  - (expected run time: within a few seconds)
