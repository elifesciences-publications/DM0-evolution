## these shell command make the html files for the S4 Supplementary Data File.
## I set up the directory structure by hand, and copied the input DM0-gd and DM25-gd directories by hand from
## results/genome-analysis/environment-comparison .

gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-genome-summary.html DM0-gd/*.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-genome-summary.html DM25-gd/*.gd

gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-html/ZDBp871_minus_CZB151.html DM0-gd/ZDBp871_minus_CZB151.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-html/ZDBp875_minus_CZB151.html DM0-gd/ZDBp875_minus_CZB151.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-html/ZDBp877_minus_CZB152.html DM0-gd/ZDBp877_minus_CZB152.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-html/ZDBp880_minus_CZB152.html DM0-gd/ZDBp880_minus_CZB152.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-html/ZDBp883_minus_CZB154.html DM0-gd/ZDBp883_minus_CZB154.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-html/ZDBp886_minus_CZB154.html DM0-gd/ZDBp886_minus_CZB154.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-html/ZDBp889_minus_ZDB67.html DM0-gd/ZDBp889_minus_ZDB67.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-html/ZDBp892_minus_ZDB67.html DM0-gd/ZDBp892_minus_ZDB67.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-html/ZDBp895_minus_ZDB68.html DM0-gd/ZDBp895_minus_ZDB68.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-html/ZDBp898_minus_ZDB68.html DM0-gd/ZDBp898_minus_ZDB68.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-html/ZDBp901_minus_ZDB69.html DM0-gd/ZDBp901_minus_ZDB69.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM0-html/ZDBp904_minus_ZDB69.html DM0-gd/ZDBp904_minus_ZDB69.gd

gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp910_minus_CZB151.html DM25-gd/ZDBp910_minus_CZB151.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp910_minus_CZB151.html DM25-gd/ZDBp910_minus_CZB151.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp911_minus_CZB151.html DM25-gd/ZDBp911_minus_CZB151.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp912_minus_CZB152.html DM25-gd/ZDBp912_minus_CZB152.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp913_minus_CZB152.html DM25-gd/ZDBp913_minus_CZB152.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp914_minus_CZB154.html DM25-gd/ZDBp914_minus_CZB154.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp915_minus_CZB154.html DM25-gd/ZDBp915_minus_CZB154.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp916_minus_ZDB67.html DM25-gd/ZDBp916_minus_ZDB67.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp917_minus_ZDB67.html DM25-gd/ZDBp917_minus_ZDB67.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp918_minus_ZDB68.html DM25-gd/ZDBp918_minus_ZDB68.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp919_minus_ZDB68.html DM25-gd/ZDBp919_minus_ZDB68.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp920_minus_ZDB69.html DM25-gd/ZDBp920_minus_ZDB69.gd
gdtools ANNOTATE -f HTML -r ../../genomes/curated-diffs/LCA.gbk -o DM25-html/ZDBp921_minus_ZDB69.html DM25-gd/ZDBp921_minus_ZDB69.gd