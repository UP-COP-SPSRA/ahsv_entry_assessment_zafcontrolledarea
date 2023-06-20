# An entry risk assessment of African horse sickness virus into the controlled area of South Africa through the legal movement of equids: Datasets and analysis code

## Introduction

This repository contains datasets and code used for the manuscript as titled above. All details regarding the data source and considerations that should be noted are discussed in the manuscript - these are important to understand before embarking on any analysis of these data.

## Referencing this dataset

John D. Grewar, Johann L. Kotze, Beverly J. Parker, Lesley S. van Helden, Gary Buhrmann, Camilla T. Weyer (2021) An entry risk assessment of African horse sickness virus into the controlled area of South Africa through the legal movement of equids [Dataset]. University of Pretoria - Community of Practice in Sanitary and Phytosanitary Risk Assessment. <https://github.com/UP-COP-SPSRA/ahsv_entry_assessment_zafcontrolledarea>

## R Code

The descriptive analysis for the manuscript was performed using R. The full code is available in the *manuscriptcode.R* file. Instructions to download this repository for the required datasets used are included in the code as well as the required libraries.

## Datasets

### CSV files

Associated spatial data relating to this repository is stored in the */datafiles* folder.

#### Outbreak case and herd level data for outbreaks in the AHS controlled area of South Africa

*data_outbreaks_ca.csv* contains actual numeric data for outbreaks  
*references_data_outbreaks_ca.csv* are the relevant references for said data, doi's have been included where available

#### Movement data for analysis of 2019 movement patterns within South Africa

*data_2019movements.csv* contains all movements associated with the AHS controlled area (CA) and Infected zone (IZ) where the origin of movement occurs in the IZ and for the 2019 calendar year.

-   *mid* refers to the unique identifier for that movement
-   *lmgid* refers to the unique identifier of the local municipality that the movement occurred from
-   *movementdate* refers to the date of movement
-   *originholdingid* refers to the unique identifier of the holding of origin that the movement occurred from
-   *destinationholdingid* refers to the unique identifier of the holding of destination that the movement occurred to
-   *totalmoved* refers to the total equids that moved with *movementspecies* the species of equid associated with the movement
-   *movementtype* relates to the type of movement that occurs, each relating to varying levels of risk of AHS introduction into the AHS CA
    -   *standard* movements are directly between the AHS IZ and the AHS CA
    -   *soq* movements are movements where horses stand at a stop-over quarantine facility in a low risk area for at least 14 days prior to an AHS PCR test and onward movement
    -   *zebq* movements are specifically associated with zebra movements where control measures are in place to facilitate safe movement, including quarantine and testing
    -   *vpqo* movements are movements where the origin of movement is in a quarantine facility that is vector protected but in the AHS infected zone.
    -   *vpsoq* movements are stop over movements where the quarantine facility is both vector protected and situated in the AHS protection zone (within the controlled area).
-   All information regarding control of equine movements in South Africa can be obtained at [here](https://www.myhorse.org.za/ahsvpn/)

#### AHS case data for 2019 based on the ECOD disease reporting system in South Africa

*data_2019cases_ecod.csv* contains monthly aggregated case totals from 2019 for AHS cases as reported on the Equine Cause of Disease (ECOD) system which is a private veterinary association reporting system. Data is primarily based on laboratory reports and a review of the AHS case reporting for 2019 can be found [here](http://jdata.co.za/myhorse/documents/infographics/Reports/2019%20General%20AHS%20surveillance%20and%20testing%20report.pdf), while the ECOD system is available through the disease reporting links on the SAEVA [website](www.saeva.co.za).  
The dataset includes the unique identifiers for the local municiaplity (*lmgid*) within which the cases occurred, allowing linking to the *data_2019movements.csv* and spatial datasets associated with this manuscript.

#### Equine census data

The *data_census_rsa.csv* dataset includes the *totalhorses* and the associated local municipality (*lmgid*) unique identifier which can be linked to other datasets in this repository. Please take note to read the materials and methods section in the manuscript before using these data for other projects. The *data_horsepopulation_factored.tif* file is the raster underlying the local municipality aggregation. Please see the ***Population data*** section of the manuscript for details on the derivation of this raster.

### Associated spatial data

Associated spatial data relating to this repository is stored in the */gis* folder. Unless otherwise indicated all files are shapefiles with a CRS of 4326. For access to maintained administrative boundary datasets in South Africa please go to the South African Demarcation Board's website which is [here](http://www.demarcation.org.za/) and observe required usage rights.

-   *lm* Local municipalities (2011) of South Africa. Movements and AHS case data are aggregated at this level.
