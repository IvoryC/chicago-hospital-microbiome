# Chicago Hospital Dataset

## Reproducing this pipeline

**Note** Before you begin, you will need to download data manually, see instructions in `01_download_data/`.

Each script is meant to run with the working directory set to the script's location.

Each R markdown (`.Rmd`) can be run from Rstudio using the `knit` button (preferred), or run from the command line:

> Rscript -e 'library(rmarkdown); rmarkdown::render("/path/to/test.Rmd", "html_document")'

Either way, you need to have the `rmarkdown` package.



## Publication

Lax et al 2017 
["Bacterial colonization and succession in a newly opened hospital"](https://www.science.org/doi/10.1126/scitranslmed.aah6500)

**Abstract**

Patients share their microbiota with their rooms and with nursing staff, and this shapes the microbial ecology of the hospital environment.  A new hospital teems with life
Lax et al. conducted a yearlong survey of the bacterial diversity associated with the patients, staff, and built surfaces in a newly opened hospital. 
They found that the bacterial communities on patient skin strongly resembled those found in their rooms. The authors demonstrated that the patient skin 
microbial communities were shaped by a diversity of clinical and environmental factors during hospitalization. They found little effect of intravenous or 
oral antibiotic treatment on the skin microbiota of patients.  
The microorganisms that inhabit hospitals may influence patient recovery and outcome, although the complexity and diversity of these bacterial communities 
can confound our ability to focus on potential pathogens in isolation. To develop a community-level understanding of how microorganisms colonize and move through 
the hospital environment, we characterized the bacterial dynamics among hospital surfaces, patients, and staff over the course of 1 year as a new hospital became operational. 
The bacteria in patient rooms, particularly on bedrails, consistently resembled the skin microbiota of the patient occupying the room. 
Bacterial communities on patients and room surfaces became increasingly similar over the course of a patient’s stay. 
Temporal correlations in community structure demonstrated that patients initially acquired room-associated taxa that predated their stay but that their own microbial 
signatures began to influence the room community structure over time. The α- and β-diversity of patient skin samples were only weakly or nonsignificantly associated with 
clinical factors such as chemotherapy, antibiotic usage, and surgical recovery, and no factor except for ambulatory status affected microbial similarity between the 
microbiotas of a patient and their room. Metagenomic analyses revealed that genes conferring antimicrobial resistance were consistently more abundant on room surfaces 
than on the skin of the patients inhabiting those rooms. In addition, persistent unique genotypes of Staphylococcus and Propionibacterium were identified. 
Dynamic Bayesian network analysis suggested that hospital staff were more likely to be a source of bacteria on the skin of patients than the reverse but that there 
were no universal patterns of transmission across patient rooms.


**Data**

The sequences for this dataset are in multiple SRA projects.  
The paper describes >6,000 samples.  The SRA project PRJEB14474 has 3,079. 
Project PRJEB13117 has 55 air samples.
