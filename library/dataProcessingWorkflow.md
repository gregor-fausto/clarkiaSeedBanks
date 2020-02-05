
### Data processing workflow

I’m trying to document the data processing workflow. I am using the
packages `tidyverse` and `readxl` (documentation:
<https://readxl.tidyverse.org/>). Link to the Issue documentation for
the data processing workflow:
<https://github.com/gregor-fausto/clarkiaSeedBanks/issues/6>.

  - [File directories](#file-directory-link)
  - [Processing data](#processing-data-link)
  - [Seeds per fruit data](#seeds-per-fruit-data-link)

### File directories

The files for data analysis are originally found in the shared Dropbox
folder `.../Dropbox/Clarkia-LTREB/20_demography_sites/`. I created a
folder on my local Dropbox `/Users/Gregor/Dropbox/dataLibrary/` to hold
a copy of these files. To copy the files, I run the following in my
Terminal:

`cp /Users/Gregor/Dropbox/Clarkia-LTREB/20_demography_sites/*
/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/`

I last ran this on 02/05/20 at 3:18 PM. The files in this directory are
now a copy of those in the shared Dropbox folder. For the time being,
the copy of the files remains outside of the folder that gets updated on
Git.

### Processing data

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

### Seeds per fruit data

Start with the seed data. This is file `.../reshapeSeeds.R` file in the
list above

``` r
setwd("/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites")

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
fileList<-paste0("/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/",seedFiles)
```

Check and count which spreadsheets have multiple sheets:

``` r
nSheets<-lapply(lapply(fileList,excel_sheets),length)
unlist(nSheets)
```

    ##  [1] 1 1 1 1 1 1 1 2 2 2

Check the number of columns in each spreadshseet:

``` r
# files for 2006-2012
unlist(lapply(lapply(fileList[1:7],read_excel),ncol))
```

    ## [1] 2 2 2 2 2 4 4

``` r
# files for 2013-2015, sheet 1
unlist(lapply(lapply(fileList[8:10],read_excel,sheet=1),ncol))
```

    ## [1] 4 5 6

``` r
# files for 2013-2015, sheet 2
unlist(lapply(lapply(fileList[8:10],read_excel,sheet=2),ncol))
```

    ## [1] 5 6 5

The data files for 2006-2011 have two columns. The data files include 1
sheet with the following data: 1 column for the site at which the
undamaged fruit was collected, and 1 column for the number of seeds in
the undamaged fruit.

``` r
yearExtract<-function(x){
  tmp<-as.numeric(sapply(strsplit(x, c("[seeds_.xls]")), "["))
  tmp<-tmp[!is.na(tmp)]
  tmp
}

years <- unlist(lapply(fileList[1:5],yearExtract))
```

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion
    
    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion
    
    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion
    
    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion
    
    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

``` r
reshapeFun <- function(x){
x %>%
  dplyr::rename(site='Site') %>%
  dplyr::rename(seedCount = `Seed no per fruit`) %>%
  dplyr::mutate(damageStatusBinary = 0, permanentPlotBinary = NA)
}

fileDfs<-lapply(fileList[1:5],read_excel)
fileDfs<-lapply(fileDfs,reshapeFun)

for(i in 1:length(years)){
  fileDfs[[i]] <- fileDfs[[i]] %>%
      dplyr::mutate(year=years[i])
}

seed20062010 <- fileDfs %>% 
  purrr::reduce(full_join) 
```

    ## Joining, by = c("site", "seedCount", "damageStatusBinary", "permanentPlotBinary", "year")
    ## Joining, by = c("site", "seedCount", "damageStatusBinary", "permanentPlotBinary", "year")
    ## Joining, by = c("site", "seedCount", "damageStatusBinary", "permanentPlotBinary", "year")
    ## Joining, by = c("site", "seedCount", "damageStatusBinary", "permanentPlotBinary", "year")

The data files for 2011-2012 have four columns: (1) a column for the
site at which the undamaged fruit was collected, (2) a column for the
number of seeds in the undamaged fruit, and (3) a column indicating
whether the fruit was collected from a permanent plot or from across the
site. (4) An additional column has notes: some seed counts were randomly
resampled from previous years’ counts and some sites did not have any
undamaged fruits in the given year.

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

Write a function to take data frames for 2006-2010 and rename the two
columns to lowercase and camelCase, add a binary variable for whether or
not the seeds were damaged (all of these fruits are undamaged), and add
a binary variable for whether or not the fruits were collected in a
permanent plot.
