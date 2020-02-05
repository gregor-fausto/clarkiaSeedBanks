
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
