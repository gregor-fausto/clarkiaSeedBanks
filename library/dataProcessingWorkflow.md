
### Data processing workflow

I’m trying to document the data processing workflow. I am using the
packages `tidyverse` and `readxl` (documentation:
<https://readxl.tidyverse.org/>). Link to the Issue documentation for
the data processing workflow:
<https://github.com/gregor-fausto/clarkiaSeedBanks/issues/6>.

  - [File directories](#file-directories)
  - [Processing data](#processing-data)
  - [Seeds per fruit data](#seeds-per-fruit-data)
  - [Fruits per plant data](#fruits-per-plant-data)
  - [Seedlings and fruiting plant
    data](#seedlings-and-fruiting-plant-data)

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

Additional changes:

  - saved
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2011.xls`
    as
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2011.xlsx`

  - saved
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2012.xls`
    as
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2012.xlsx`

  - deleted empty rows at the bottom of the file
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2012.xlsx`;
    these are filled with some formatting and otherwise may give
    problems with data importing

  - saved
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2013.xls`
    as
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2013.xlsx`

  - deleted empty rows at the bottom of the file
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2013.xlsx`;
    these are filled with some formatting and otherwise may give
    problems with data importing

  - saved
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2014.xls`
    as
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2014.xlsx`

  - deleted empty rows at the bottom of the file
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2014.xlsx`,
    sheet 2; these are filled with some formatting and otherwise may
    give problems with data importing

  - saved
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2015.xls`
    as
    `/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/seeds_2015.xlsx`

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
list above.

``` r
setwd("/Users/Gregor/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites")

filenames<-list.files(pattern=paste("^seeds_"), recursive=TRUE)
worksheet_2006.2010 <- filenames[1:5]
worksheet_2011.2012 <- filenames[c(7,9)]
worksheet_2013.2015 <- filenames[c(11,13,15)]
print(c(worksheet_2006.2010,worksheet_2011.2012,worksheet_2013.2015))
```

    ##  [1] "seeds_2006.xls"  "seeds_2007.xls"  "seeds_2008.xls"  "seeds_2009.xls" 
    ##  [5] "seeds_2010.xls"  "seeds_2011.xlsx" "seeds_2012.xlsx" "seeds_2013.xlsx"
    ##  [9] "seeds_2014.xlsx" "seeds_2015.xlsx"

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

The data files for 2006-2010 have two columns. The data files include 1
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

``` r
yearExtract<-function(x){
  tmp<-as.numeric(sapply(strsplit(x, c("[seeds_.xls]")), "["))
  tmp<-tmp[!is.na(tmp)]
  tmp
}

years <- unlist(lapply(fileList[6:7],yearExtract))
```

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion
    
    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

``` r
reshapeFun <- function(x){
x %>%
  dplyr::rename(site='Site') %>%
  dplyr::rename(seedCount = `Seed no per fruit`) %>%
  dplyr::mutate(damageStatusBinary = 0)
}

appendFontColor <- function(x){

# Step 1: import the table taking only cell values and ignoring the formatting
tidyExcel <- read_excel(x,na="NA") 
tidyExcel <- reshapeFun(tidyExcel)

# Step 2: import one column of the table, taking only the formatting and not the
# cell values

# `formats` is a pallette of fill colours that can be indexed by the
# `local_format_id` of a given cell to get the fill colour of that cell
font_colours <- xlsx_formats(x)$local$font$color$rgb

# Import all the cells, filter out the header row, filter for the first column,
# and create a new column `font_colour` of the font colours, by looking up the
# local_format_id of each cell in the `font_colours` pallette.
fonts <- xlsx_cells(x, sheet = 1) %>%
  filter(row >= 2, col == 1) %>% # Omit the header row
  mutate(font_colour = font_colours[local_format_id]) %>%
  select(font_colour)

# Step 3: append the `font` column to the rest of the data
tmp <- bind_cols(tidyExcel, fonts)

out <- tmp %>% 
  dplyr::select(site,seedCount,permanentPlot,damageStatusBinary,font_colour) %>%
  dplyr::mutate(fieldData = ifelse(font_colour=="FF000000",1,0)) %>%
  dplyr::select(-font_colour)

return(out)

}

fileDfs<-lapply(fileList[6:7],appendFontColor)

for(i in 1:length(years)){
  fileDfs[[i]] <- fileDfs[[i]] %>%
      dplyr::mutate(year=years[i])
}

seed20112012 <- fileDfs %>% 
  purrr::reduce(full_join) 
```

    ## Joining, by = c("site", "seedCount", "permanentPlot", "damageStatusBinary", "fieldData", "year")

``` r
seed20112012 %>% 
  dplyr::filter(permanentPlot!=fieldData)
```

    ## # A tibble: 1 x 6
    ##   site  seedCount permanentPlot damageStatusBinary fieldData  year
    ##   <chr>     <dbl>         <dbl>              <dbl>     <dbl> <dbl>
    ## 1 LCE          34             1                  0         0  2011

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

``` r
yearExtract<-function(x){
  tmp<-as.numeric(sapply(strsplit(x, c("[seeds_.xls]")), "["))
  tmp<-tmp[!is.na(tmp)]
  tmp
}

years <- unlist(lapply(fileList[8:10],yearExtract))
```

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion
    
    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion
    
    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

``` r
reshapeFun <- function(x){
x %>%
  dplyr::rename(site='Site') %>%
  dplyr::rename(seedCount = `Seed no per undamaged fruit`) %>%
  dplyr::mutate(damageStatusBinary = 0)
}

appendFontColor <- function(x){

# Step 1: import the table taking only cell values and ignoring the formatting
tidyExcel <- read_excel(x,na="NA",sheet=1) 
tidyExcel <- reshapeFun(tidyExcel)

# Step 2: import one column of the table, taking only the formatting and not the
# cell values

# `formats` is a pallette of fill colours that can be indexed by the
# `local_format_id` of a given cell to get the fill colour of that cell
font_colours <- xlsx_formats(x)$local$font$color$rgb

# Import all the cells, filter out the header row, filter for the first column,
# and create a new column `font_colour` of the font colours, by looking up the
# local_format_id of each cell in the `font_colours` pallette.
fonts <- xlsx_cells(x, sheet = 1) %>%
  filter(row >= 2, col == 1) %>% # Omit the header row
  mutate(font_colour = font_colours[local_format_id]) %>%
  select(font_colour)

# Step 3: append the `font` column to the rest of the data
tmp <- bind_cols(tidyExcel, fonts)

out <- tmp %>% 
  dplyr::select(site,seedCount,permanentPlot,damageStatusBinary,font_colour) %>%
  dplyr::mutate(fieldData = ifelse(font_colour=="FF000000"|font_colour=="FF008000",1,0)) %>%
  dplyr::select(-font_colour)

return(out)

}

fileDfs<-lapply(fileList[8:10],appendFontColor)
```

    ## New names:
    ## * `` -> ...4

    ## New names:
    ## * `` -> ...4
    ## * `` -> ...5

``` r
for(i in 1:length(years)){
  fileDfs[[i]] <- fileDfs[[i]] %>%
      dplyr::mutate(year=years[i])
}

seedUndamaged20132015 <- fileDfs %>% 
  purrr::reduce(full_join) 
```

    ## Joining, by = c("site", "seedCount", "permanentPlot", "damageStatusBinary", "fieldData", "year")

    ## Joining, by = c("site", "seedCount", "permanentPlot", "damageStatusBinary", "fieldData", "year")

``` r
seedUndamaged20132015 %>% 
  dplyr::filter(permanentPlot!=fieldData)
```

    ## # A tibble: 0 x 6
    ## # … with 6 variables: site <chr>, seedCount <dbl>, permanentPlot <dbl>,
    ## #   damageStatusBinary <dbl>, fieldData <dbl>, year <dbl>

The sheets for damaged fruits have the following columns: 1 column for
the site at which the damaged fruit was collected, and 1 column for the
number of seeds in the damaged fruit. Another column has notes: some
sites did not have any damaged fruits in the given year. Some of the
data files have empty columns between the data and the notes.

``` r
yearExtract<-function(x){
  tmp<-as.numeric(sapply(strsplit(x, c("[seeds_.xls]")), "["))
  tmp<-tmp[!is.na(tmp)]
  tmp
}

years <- unlist(lapply(fileList[8:10],yearExtract))
```

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion
    
    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion
    
    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

``` r
reshapeFun <- function(x){
x %>%
  dplyr::rename(site='Site') %>%
  dplyr::rename(seedCount = `Seed no per damaged fruit`) %>%
  dplyr::mutate(damageStatusBinary = 1, permanentPlot = NA)
}

appendFontColor <- function(x){

# Step 1: import the table taking only cell values and ignoring the formatting
tidyExcel <- read_excel(x,na=c("NA","MA"),sheet=2) 
tidyExcel <- reshapeFun(tidyExcel)

# Step 2: import one column of the table, taking only the formatting and not the
# cell values

# `formats` is a pallette of fill colours that can be indexed by the
# `local_format_id` of a given cell to get the fill colour of that cell
font_colours <- xlsx_formats(x)$local$font$color$rgb

# Import all the cells, filter out the header row, filter for the first column,
# and create a new column `font_colour` of the font colours, by looking up the
# local_format_id of each cell in the `font_colours` pallette.
fonts <- xlsx_cells(x, sheet = 2) %>%
  filter(row >= 2, col == 1) %>% # Omit the header row
  mutate(font_colour = font_colours[local_format_id]) %>%
  select(font_colour)

# Step 3: append the `font` column to the rest of the data
tmp <- bind_cols(tidyExcel, fonts)

out <- tmp %>% 
  dplyr::select(site,seedCount,permanentPlot,damageStatusBinary,font_colour) %>%
  dplyr::mutate(fieldData = ifelse(font_colour=="FF000000"|font_colour=="FF008000",1,0)) %>%
  dplyr::select(-font_colour)

return(out)

}

fileDfs<-lapply(fileList[8:10],appendFontColor)
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

``` r
for(i in 1:length(years)){
  fileDfs[[i]] <- fileDfs[[i]] %>%
      dplyr::mutate(year=years[i])
}

seedDamaged20132015 <- fileDfs %>% 
  purrr::reduce(full_join) 
```

    ## Joining, by = c("site", "seedCount", "permanentPlot", "damageStatusBinary", "fieldData", "year")

    ## Joining, by = c("site", "seedCount", "permanentPlot", "damageStatusBinary", "fieldData", "year")

``` r
seedDamaged20132015 %>% 
  dplyr::filter(permanentPlot!=fieldData)
```

    ## # A tibble: 0 x 6
    ## # … with 6 variables: site <chr>, seedCount <dbl>, permanentPlot <lgl>,
    ## #   damageStatusBinary <dbl>, fieldData <dbl>, year <dbl>

It might be useful to summarize the datasets:

``` r
summary20062010<-seed20062010 %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(count = sum(!is.na(seedCount))) %>%
  tidyr::spread(key="year",value="count")

summary20112012<-seed20112012 %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(count = sum(!is.na(seedCount))) %>%
  tidyr::spread(key="year",value="count")

summaryUndamaged20132012<-seedUndamaged20132015  %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(count = sum(!is.na(seedCount))) %>%
  tidyr::spread(key="year",value="count")

summaryDamaged20132012<-seedDamaged20132015  %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(count = sum(!is.na(seedCount))) %>%
  tidyr::spread(key="year",value="count")

summaryTable <- summary20062010 %>%
  dplyr::full_join(summary20112012) %>%
  dplyr::full_join(summaryUndamaged20132012)
```

    ## Joining, by = "site"
    ## Joining, by = "site"

``` r
kable(summaryTable, caption="Summary table of the number of counts of seeds from undamaged fruits")
```

| site | 2006 | 2007 | 2008 | 2009 | 2010 | 2011 | 2012 | 2013 | 2014 | 2015 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BG   |   21 |   21 |   41 |   30 |   30 |   28 |   30 |   30 |   30 |   29 |
| BR   |   20 |   20 |   32 |   30 |   29 |   30 |   30 |   30 |   31 |   31 |
| CF   |   20 |   20 |   30 |   29 |   34 |   30 |   30 |   30 |   30 |   26 |
| CP3  |   20 |   20 |   41 |   30 |   30 |   29 |   30 |   30 |   30 |   21 |
| DEM  |   20 |   20 |   29 |   30 |   32 |   27 |   30 |   30 |   24 |   28 |
| DLW  |   20 |   20 |   22 |   30 |   31 |   28 |   30 |   30 |    1 |   30 |
| EC   |   20 |   20 |   29 |   30 |   31 |   26 |   30 |   30 |   31 |   31 |
| FR   |   20 |   20 |   31 |   30 |   31 |   31 |   30 |   30 |   46 |   30 |
| GCN  |   20 |   20 |   29 |   30 |   32 |   30 |   30 |   30 |   28 |   29 |
| KYE  |   20 |   20 |   30 |   30 |   30 |   30 |   30 |   30 |   30 |   29 |
| LCE  |   20 |   20 |   30 |   30 |   32 |   30 |    0 |    0 |   30 |   38 |
| LCW  |   20 |   20 |   28 |   30 |   35 |   32 |   30 |   30 |    1 |    0 |
| LO   |   32 |   32 |   30 |   30 |   37 |   32 |   30 |   30 |   30 |    0 |
| MC   |   20 |   20 |   29 |   30 |   35 |   30 |   30 |   30 |   46 |   35 |
| OKRE |   20 |   20 |   26 |   30 |   30 |   28 |   30 |   30 |   18 |   24 |
| OKRW |   20 |   20 |   33 |   30 |   34 |   28 |   30 |   30 |   30 |    0 |
| OSR  |   20 |   20 |   32 |   30 |   30 |   28 |   30 |   30 |   30 |   37 |
| S22  |   20 |   20 |   33 |   30 |   28 |   23 |   30 |   30 |   23 |   30 |
| SM   |   20 |   20 |   31 |   29 |   32 |   30 |   30 |   30 |   30 |    8 |
| URS  |   18 |   30 |   25 |   30 |   30 |   27 |   30 |   30 |    1 |    0 |

Summary table of the number of counts of seeds from undamaged fruits

``` r
kable(summaryDamaged20132012, caption="Summary table of the number of counts of seeds from damaged fruits")
```

| site | 2013 | 2014 | 2015 |
| :--- | ---: | ---: | ---: |
| BG   |   17 |   20 |   11 |
| BR   |   24 |   25 |   23 |
| CF   |   22 |   29 |   27 |
| CP3  |   23 |   11 |    9 |
| DEM  |    5 |   14 |   25 |
| DLW  |    8 |    0 |   30 |
| EC   |   12 |   22 |    8 |
| FR   |    2 |   25 |   15 |
| GCN  |    0 |    0 |    3 |
| KYE  |   23 |   34 |   15 |
| LCE  |    0 |   11 |   15 |
| LCW  |    0 |    0 |    0 |
| LO   |    4 |   14 |    0 |
| MC   |    4 |   15 |   15 |
| OKRE |   13 |    8 |    9 |
| OKRW |    0 |    4 |    0 |
| OSR  |    1 |   19 |   26 |
| S22  |    1 |    3 |    2 |
| SM   |    1 |    3 |    0 |
| URS  |    0 |    0 |    0 |

Summary table of the number of counts of seeds from damaged fruits

### Fruits per plant data

### Seedlings and fruiting plant data
