#!/bin/bash
## This script should be run in the DM0 folder on Turing (or your friendly local high-performance computing cluster).
## This script uses the slurm job scheduler. Edit throughout if your HPC uses something other than slurm to schedule jobs on the cluster.

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/DM0-evolved-reruns/CZB151_R1.fastq.gz sequence-data/DM0-evolved-reruns/CZB151_R2.fastq.gz -o consensus/CZB151"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/DM0-evolved-reruns/CZB151_R1.fastq.gz sequence-data/DM0-evolved-reruns/CZB151_R2.fastq.gz -o polymorphism/CZB151"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/DM0-evolved-reruns/CZB152_R1.fastq.gz sequence-data/DM0-evolved-reruns/CZB152_R2.fastq.gz sequence-data/CZB152_combined_reads.fastq -o consensus/CZB152"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/DM0-evolved-reruns/CZB152_R1.fastq.gz sequence-data/DM0-evolved-reruns/CZB152_R2.fastq.gz sequence-data/CZB152_combined_reads.fastq -o polymorphism/CZB152"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/DM0-evolved-reruns/CZB154_R1.fastq.gz sequence-data/DM0-evolved-reruns/CZB154_R2.fastq.gz sequence-data/CZB154_reads.fastq -o consensus/CZB154"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/DM0-evolved-reruns/CZB154_R1.fastq.gz sequence-data/DM0-evolved-reruns/CZB154_R2.fastq.gz sequence-data/CZB154_reads.fastq -o polymorphism/CZB154"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDB67_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDB67_R2.fastq.gz -o consensus/ZDB67"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDB67_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDB67_R2.fastq.gz -o polymorphism/ZDB67"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDB68_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDB68_R2.fastq.gz -o consensus/ZDB68"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDB68_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDB68_R2.fastq.gz -o polymorphism/ZDB68"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDB69_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDB69_R2.fastq.gz -o consensus/ZDB69"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDB69_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDB69_R2.fastq.gz -o polymorphism/ZDB69"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp871_combined_reads.fastq -o consensus/ZDBp871"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp871_combined_reads.fastq -o polymorphism/ZDBp871"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp874_combined_reads.fastq -o consensus/ZDBp874"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp874_combined_reads.fastq -o polymorphism/ZDBp874"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDBp875_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDBp875_R2.fastq.gz -o consensus/ZDBp875"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDBp875_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDBp875_R2.fastq.gz -o polymorphism/ZDBp875"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDBp877_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDBp877_R2.fastq.gz sequence-data/ZDBp877.solexa.fastq -o consensus/ZDBp877"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDBp877_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDBp877_R2.fastq.gz sequence-data/ZDBp877.solexa.fastq -o polymorphism/ZDBp877"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp880_combined_reads.fastq -o consensus/ZDBp880"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp880_combined_reads.fastq -o polymorphism/ZDBp880"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDBp883_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDBp883_R2.fastq.gz sequence-data/ZDBp883_combined_reads.fastq -o consensus/ZDBp883"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDBp883_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDBp883_R2.fastq.gz sequence-data/ZDBp883_combined_reads.fastq -o polymorphism/ZDBp883"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp886_combined_reads.fastq -o consensus/ZDBp886"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp886_combined_reads.fastq -o polymorphism/ZDBp886"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDBp889_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDBp889_R2.fastq.gz sequence-data/ZDBp889_reads.fastq -o consensus/ZDBp889"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDBp889_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDBp889_R2.fastq.gz sequence-data/ZDBp889_reads.fastq -o polymorphism/ZDBp889"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp892_combined_reads.fastq -o consensus/ZDBp892"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp892_combined_reads.fastq -o polymorphism/ZDBp892"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDBp895_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDBp895_R2.fastq.gz sequence-data/ZDBp895_reads.fastq  -o consensus/ZDBp895"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/DM0-evolved-reruns/ZDBp895_R1.fastq.gz sequence-data/DM0-evolved-reruns/ZDBp895_R2.fastq.gz sequence-data/ZDBp895_reads.fastq  -o polymorphism/ZDBp895"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp898_combined_reads.fastq -o consensus/ZDBp898"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp898_combined_reads.fastq -o polymorphism/ZDBp898"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp901_combined_reads.fastq -o consensus/ZDBp901"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp901_combined_reads.fastq -o polymorphism/ZDBp901"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp904_combined_reads.fastq -o consensus/ZDBp904"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp904_combined_reads.fastq -o polymorphism/ZDBp904"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp910_combined_reads.fastq -o consensus/ZDBp910"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp910_combined_reads.fastq -o polymorphism/ZDBp910"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp911_combined_reads.fastq -o consensus/ZDBp911"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp911_combined_reads.fastq -o polymorphism/ZDBp911"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp912_combined_reads.fastq -o consensus/ZDBp912"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp912_combined_reads.fastq -o polymorphism/ZDBp912"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp913_combined_reads.fastq -o consensus/ZDBp913"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp913_combined_reads.fastq -o polymorphism/ZDBp913"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp914_combined_reads.fastq -o consensus/ZDBp914"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp914_combined_reads.fastq -o polymorphism/ZDBp914"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp915_combined_reads.fastq -o consensus/ZDBp915"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp915_combined_reads.fastq -o polymorphism/ZDBp915"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp916_combined_reads.fastq -o consensus/ZDBp916"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp916_combined_reads.fastq -o polymorphism/ZDBp916"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp917_combined_reads.fastq -o consensus/ZDBp917"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp917_combined_reads.fastq -o polymorphism/ZDBp917"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp918_combined_reads.fastq -o consensus/ZDBp918"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp918_combined_reads.fastq -o polymorphism/ZDBp918"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp919_combined_reads.fastq -o consensus/ZDBp919"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp919_combined_reads.fastq -o polymorphism/ZDBp919"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp920_combined_reads.fastq -o consensus/ZDBp920"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp920_combined_reads.fastq -o polymorphism/ZDBp920"

sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -r LCA.gff3 sequence-data/ZDBp921_combined_reads.fastq -o consensus/ZDBp921"
sbatch -p main -t 24:00:00 --mem=12G --wrap="breseq -p -r LCA.gff3 sequence-data/ZDBp921_combined_reads.fastq -o polymorphism/ZDBp921"
