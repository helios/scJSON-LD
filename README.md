scJSON-LD is inteded for Single Cell analysis in Python and R. For Python it's using [Scanpy](https://github.com/theislab/scanpy) while for R it's using [Seurat](https://github.com/satijalab/seurat).

During the [BioHackathon 2017](http://2017.biohackathon.org/) scJSON-LD was an idea to keep track of the analysis workflow and data generated along the analyses. From the begining was clear that describing raw data measures was quite complicated, JSON object were huge and required a definition of new data object which was not really easy to design. At that time Single Cells libraries both Python and R , respectivaly [Scanpy](https://github.com/theislab/scanpy) and  [Seurat](https://github.com/satijalab/seurat) were changing quickly and were not providing a stable development environement.

More recentrly bot Seurat and ScanPy stabilized so the time is mature to smoothly inject JSON-LD concepts into the Single Cell workflows.

