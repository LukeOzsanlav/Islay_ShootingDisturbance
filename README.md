# Contrasting effects of shooting disturbance on the movement and behaviour of sympatric wildfowl species 🔫 🦆
The impact of shooting on time-activity budgets, movement and habitat selection in Barnacle Geese *Branta leucopsis* and White-fronted Geese *Anser albifrons flavirostris* on Islay, Scotland. 

## Authors 🖊️
- Luke Ozsanlav-Harris <a itemprop="sameAs" content="https://orcid.org/0000-0003-3889-6722" href="https://orcid.org/0000-0003-3889-6722" target="orcid.widget" rel="noopener" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" target="_blank" style="width:1em;margin-right:.5em;"/></a>
- Aimée L. S. Mcintosh <a itemprop="sameAs" content="https://orcid.org/0000-0002-4975-3682" href="https://orcid.org/0000-0002-4975-3682" target="orcid.widget" rel="noopener" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" target="_blank" style="width:1em;margin-right:.5em;"/></a>
- Larry R. Griffin
- Geoff M. Hilton <a itemprop="sameAs" content="https://orcid.org/0000-0001-9062-3030" href="https://orcid.org/0000-0001-9062-3030" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
- Jessica M. Shaw <a itemprop="sameAs" content="https://orcid.org/0000-0003-0862-9260" href="https://orcid.org/0000-0003-0862-9260" target="orcid.widget" rel="noopener" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" target="_blank" style="width:1em;margin-right:.5em;"/></a>
- Stuart Bearhop <a itemprop="sameAs" content="https://orcid.org/0000-0002-5864-0129" href="https://orcid.org/0000-0002-5864-0129" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>

## Code description 👨‍💻
- `R code/1- Create maps for publication.R`: Create maps of Islay, including habitat and roads
- `R code/2- Displacement Model.R`: Objective 1, modelling the immediate displacement in response to shooting
- `R code/3- Daily Behaviour Model.R`: Objective 4, modelling behavioral changes in response to shooting disturbance (GWfG only)
- `R code/4- Daily ODBA Model.R`: Objective 5, modelling change in daily ODBA due to shooting disturbance
- `R code/5- Daily Foraging Distance.R`: Objective 3, modelling changed in total daily movement in response to shooting disturbance 

## Data description 📊
- `Derived data/Foraging_Distance_Model_Data.csv`
- `Derived data/GBG_INDIV_SHOOT_DATE_797m.csv`
- `Derived data/GWFG_INDIV_SHOOT_DATE_797m.csv`
- `Biologging Data/Script2_BiologgingData.RDS`: contains the general GPS locations for both species with the step lengths (measured in km) calculated from the precise locations. If the bird was within 4km of a shooting event when it occurred then the information for that shooting event and the distance of the bird to the shooting event is appended onto the GPS fix immediately prior to the shooting event.
- `Biologging Data/Script3_BiologgingData.RDS`: contains accelerometer data classified into behaviours of GWfG on Islay. 
- `Biologging Data/Script4_BiologgingData.RDS`: contains ODBA values calculated from accelerometers for GBG and GWfG on Islay.
- `Biologging Data/Script3_ShootingProximity.RDS`: contains all the instances that GWfG were close spatially and temporally to a shooting event while the device was collecting accelerometer data. This can be used to work out which days birds were disturbed by shooting on
- `Biologging Data/Script4_ShootingProximity.RDS`: contains all the instances that GWfG and GBG were close spatially and temporally to a shooting event while the device was collecting accelerometer data. This can be used to work out which days birds were disturbed by shooting on
- `Landcover Data/Islay landcover data`: habitat data from Islay from the years 2015, 2017, 2018, 2019 amd 2020 as rasters. This was cropped from larger UK wide raster. This data is from UKCEH and wider datsets can be accessed [here](https://www.ceh.ac.uk/data/ukceh-land-cover-maps).
- `Landcover Data/High-Res Coastline`: high resolution coast outline for the UK and Islay as well as bounding boxes for both of these regions.
- `Landcover Data/88090_ISLAY_GMS_FIELD_BOUNDARY`: shapefile of all the agricultural field boundaries on Islay.
- `Landcover Data/Ramsar Outline`: shapefiles of the two RAMSAR areas on Islay. One at Gruinart and one at the Oa.
- `Landcover Data/Roads`: shapefile of all classified public roads in Scotland.
- `MetaData/Tagged bird summary data new.csv` metadata file that contains the sex, ringing location and deployment dates of all tagged GWfG
- `Shooting logs/All_logs_cleaned.csv`: All shooting logs from Islay containing the spatiotemporal information on each shooting event


