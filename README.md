
### Seed banks in Clarkia xantiana

Siegmund, G.-F., and M. A. Geber.

-----

### Abstract

….

-----

### Repository Directory

1.  `library` contains the \[…\] datasets associated with the paper. I
    have included the metadata for each dataset below.

2.  `reshapeDataScripts` contains the R scripts used to download and
    reshape the datasets.

3.  `modelBuild` contains the R scripts to fit models.

4.  `modelCodeJAGS` contains the JAGS scripts.

5.  `analysisCode` contains scripts to analyze model output.

6.  `products` contains a folder with the R scripts used to make all
    figures and tables, and a folder with the manuscript.

-----

### Library Metadata

The table below summarizes the data generated in this study (also
available [as a pdf at this
link](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/products/tables/data-summary.pdf))

|                                        |                  |               |           |
| :------------------------------------- | :--------------: | :-----------: | :-------: |
| **Seed vital rates**                   |        —         |       —       |     —     |
| Seed survival and germination          | Seed bag burial  | Y<sub>1</sub> | 2006-2009 |
| Seed viability                         | Viability trials | Y<sub>2</sub> | 2006-2009 |
| Seed survival and germination          |    Seed pots     | Y<sub>3</sub> | 2013-2019 |
| **Seedling survival** -                |        –         |       —       |     —     |
| Seedling survival to fruiting          |  Field surveys   | Y<sub>4</sub> | 2006-2019 |
| **Fruits per plant**                   |        —         |       —       |     —     |
| Total fruit equivalents per plant      |  Field surveys   | Y<sub>5</sub> | 2006-2012 |
| Undamaged and damaged fruits per plant |  Field surveys   | Y<sub>6</sub> | 2013-2019 |
| Total fruit equivalents per plant      |   Extra plots    | Y<sub>7</sub> | 2006-2012 |
| Undamaged and damaged fruits per plant |   Extra plots    | Y<sub>8</sub> | 2013-2019 |
| **Seeds per fruit**                    |        —         |       —       |     —     |
| Seeds per undamaged fruit              |    Lab counts    | Y<sub>9</sub> | 2006-2019 |
| Seeds per damaged fruit                |    Lab counts    | Y<sub>0</sub> | 2013-2019 |

Links to all the datasets follow:

  - [Seed bag
    burial](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/library/dataProcessingWorkflow.md#seed-bag-data)
  - [Viability
    trials](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/library/dataProcessingWorkflow.md#viability-trial-data)
  - [Seed
    pots](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/library/dataProcessingWorkflow.md#seed-pot-data)
    - needs to be updated
  - [Field surveys for seedling survival to
    fruiting](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/library/dataProcessingWorkflow.md#seedlings-and-fruiting-plant-data)
  - [Field surveys for total fruit equivalents per plant, and
    undamaged/damaged fruits per
    plant](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/library/dataProcessingWorkflow.md#fruits-per-plant-data-extra-plots)
  - [Lab counts for seeds per fruit from undamaged and damaged
    fruits](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/library/dataProcessingWorkflow.md#seeds-per-fruit-data)

### Summary

The table below summarizes the parameters and years for which I have
estimates.

| Parameter | 2006 | 2007 | 2008 | 2009 | 2010 | 2011 | 2012 | 2013 | 2014 | 2015 | 2016 | 2017 | 2018 | 2019 |
| :-------- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| s0        | N    | N    | N    | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   |
| s1        | Y    | Y    | Y    | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   |
| g1        | Y    | Y    | Y    | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   |
| s2        | Y    | Y    | Y    | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   | \-   |
| sigma     | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    |
| Fec       | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | N    | N    | N    | N    | N    | N    |
| phi       | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | Y    | N    |

#### sigmaDF

Metadata for the seedling survival to fruiting dataset described in the
Methods and Appendices XXXX.

| field            | description                                                      |
| :--------------- | :--------------------------------------------------------------- |
| site             | name of site, see methods and Appendix XX                        |
| year             | year of …                                                        |
| transect         | transect identifier, see Methods and Appendix XX                 |
| position         | position within transect identifier, see Methods and Appendix XX |
| noSeedlings      | number of seedlings at plot                                      |
| noFruitingPlants | number of fruiting plants at plot                                |

<br>
