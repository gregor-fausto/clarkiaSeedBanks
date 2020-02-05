
### Data processing workflow

I’m trying to document the data processing workflow. I am using the
packages `tidyverse` and `readxl` (documentation:
<https://readxl.tidyverse.org/>). Link to the Issue documentation for
the data processing workflow:
<https://github.com/gregor-fausto/clarkiaSeedBanks/issues/6>.

Load the libraries for data processing (see
<https://github.com/r-lib/rlang/issues/669> for the overwrite message I
am suppressing)

``` r
library(tidyverse)
library(readxl)
library(knitr)
library(lubridate)
library(tidyxl)
```

Scripts for the original processing I did is located in the following
folder. I want to document when I last updated these files before I
start updating the workflow. None of the original scripts are going to
be changed.

``` r
dir = c("/Users/Gregor/Dropbox/projects/clarkiaScripts/code/reshapeDataScripts/")
files = list.files(path=dir)
dt <- data.frame(files)
fileModified <- function(x){file.info(paste0(dir,x))$mtime}

for(i in 1:dim(dt)[1]){
  dt[i,2] <- ymd_hms(fileModified(files[i])) %>%
    data.frame() %>% rename(dateTime='.') 
}

kable(dt)
```

| files                                 | dateTime            |
| :------------------------------------ | :------------------ |
| Estimation of Clarkia vital rates.doc | 2018-08-05 10:03:31 |
| germ\_and\_viability.xls              | 2016-10-24 09:59:32 |
| germ\&viability.xls                   | 2016-10-24 09:59:32 |
| reshapeAllClimate.R                   | 2017-11-07 10:39:30 |
| reshapeFruitAndSeedData.R             | 2017-12-18 14:09:52 |
| reshapeFruits.R                       | 2016-12-23 19:53:00 |
| reshapeFruitsPerPlantFromPlots.R      | 2016-10-18 08:55:23 |
| reshapeMasterFile.R                   | 2016-02-02 17:33:17 |
| reshapeSeedBags-VersionThree.R        | 2018-08-01 19:01:20 |
| reshapeSeedBags-VersionTwo.R          | 2017-02-14 13:57:21 |
| reshapeSeedBags.R                     | 2017-02-03 09:50:15 |
| reshapeSeedBagsToSeedPots.R           | 2017-11-07 11:31:11 |
| reshapeSeedPots.R                     | 2017-11-03 14:44:45 |
| reshapeSeeds.R                        | 2016-12-23 19:52:24 |
| reshapeSiteVariables.R                | 2017-12-08 14:23:44 |
| reshapeSurvivalFecundity.R            | 2016-10-18 08:56:25 |

Start with the seed data. This is file `.../reshapeSeeds.R` file in the
list above

``` r
setwd("/Users/Gregor/Dropbox/Clarkia-LTREB/20_demography_sites")

filenames<-list.files(pattern=paste("^seeds_"), recursive=TRUE)
worksheet_2006.2010 <- filenames[1:5]
worksheet_2011.2012 <- filenames[6:7]
worksheet_2013.2015 <- filenames[8:10]
print(c(worksheet_2006.2010,worksheet_2011.2012,worksheet_2013.2015))
```

    ##  [1] "seeds_2006.xls" "seeds_2007.xls" "seeds_2008.xls" "seeds_2009.xls"
    ##  [5] "seeds_2010.xls" "seeds_2011.xls" "seeds_2012.xls" "seeds_2013.xls"
    ##  [9] "seeds_2014.xls" "seeds_2015.xls"

read excel data

``` r
seedFiles<-c(worksheet_2006.2010,worksheet_2011.2012,worksheet_2013.2015)
fileList<-paste0("/Users/Gregor/Dropbox/Clarkia-LTREB/20_demography_sites/",seedFiles)
```

Check and count which spreadsheets have multiple sheets:

``` r
nSheets<-lapply(lapply(fileList,excel_sheets),length)
```

Check the number of columns in each spreadshseet:

``` r
lapply(lapply(fileList[1:7],read_excel),ncol)
```

    ## [[1]]
    ## [1] 2
    ## 
    ## [[2]]
    ## [1] 2
    ## 
    ## [[3]]
    ## [1] 2
    ## 
    ## [[4]]
    ## [1] 2
    ## 
    ## [[5]]
    ## [1] 2
    ## 
    ## [[6]]
    ## [1] 4
    ## 
    ## [[7]]
    ## [1] 4

``` r
lapply(lapply(fileList[8:10],read_excel,sheet=1),ncol)
```

    ## New names:
    ## * `` -> ...4

    ## New names:
    ## * `` -> ...4
    ## * `` -> ...5

    ## [[1]]
    ## [1] 4
    ## 
    ## [[2]]
    ## [1] 5
    ## 
    ## [[3]]
    ## [1] 6

``` r
lapply(lapply(fileList[8:10],read_excel,sheet=2),ncol)
```

    ## New names:
    ## * `` -> ...3
    ## * `` -> ...4

    ## New names:
    ## * `` -> ...3
    ## * `` -> ...4
    ## * `` -> ...5

    ## New names:
    ## * `` -> ...3
    ## * `` -> ...4

    ## [[1]]
    ## [1] 5
    ## 
    ## [[2]]
    ## [1] 6
    ## 
    ## [[3]]
    ## [1] 5

From 2006-2012, we collected data on the number of seeds per undamaged
fruit. At each site, undamaged fruits were collected from plants spread
across the entire site. The data files thus include 1 sheet with the
following data: 1 column for the site at which the undamaged fruit was
collected, and 1 column for the number of seeds in the undamaged fruit.

The data files for 2011-2012 have three data columns: 1 column for the
site at which the undamaged fruit was collected, 1 column for the number
of seeds in the undamaged fruit, and 1 column indicating whether the
fruit was collected from a permanent plot or from across the site. A
fourth column has notes: some seed counts were randomly resampled from
previous years’ counts and some sites did not have any undamaged fruits
in the given year.

Start by processing these files; to do is including information on
permanent plot as NA, removing the notes, talk to Monica about putting
any files with additional formatting into xlsx?

From 2013-2015, we collected data on the number of seeds per undamaged
fruit and the number of seeds per damaged fruit. The data files thus
include two sheets, 1 for the number of seeds per undamaged fruit and 1
for the numbe of seeds per damaged fruit. At each site, undamaged fruits
were collected from plants spread across the entire site. The data files
thus include 1 column for the site at which the undamaged fruit was
collected, and 1 column for the number of seeds in the fruit.

The sheets for undamaged fruits have the following columns: 1 column for
the site at which the undamaged fruit was collected, 1 column for the
number of seeds in the undamaged fruit, and 1 column indicating whether
the fruit was collected from a permanent plot or from across the site. A
fourth column has notes: some seed counts were randomly resampled from
previous years’ counts and some sites did not have any undamaged fruits
in the given year. Some of the data files have empty columns between the
data and the notes.

The sheets for damaged fruits have the following columns: 1 column for
the site at which the damaged fruit was collected, and 1 column for the
number of seeds in the damaged fruit. Another column has notes: some
sites did not have any damaged fruits in the given year. Some of the
data files have empty columns between the data and the notes.

One issue is that formatting is used in these data files. The tidyxl
packages retains cell information.

``` r
library(tidyxl)
#lapply(fileList[8:10],xlsx_cells,sheet=2)
```
