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
- `R code/0- Create maps for publication.R`: Create maps of Islay, including habitat and roads
- `R code/1- Combine tracks + Add landcover.R`: Prepare GWfG tracking data for analysis and add landcover associated with each GPS locations
- `R code/2- Displacement Model.R`: Objective 1, modelling the immediate displacement in response to shooting
- `R code/3- Daily Behaviour Model.R`: Objective 4, modelling behavioral changes in response to shooting disturbance (GWfG only)
- `R code/4- Daily ODBA Model.R`: Objective 5, modelling change in daily ODBA due to shooting disturbance
- `R code/5- Daily Foraging Distance.R`: Objective 3, modelling changed in total daily movement in response to shooting disturbance 
- `R code/Useful Functions` folder containing useful functions for this project, it is not currently used in the main work flow

## Data description 📊
- `Derived data/All_winter_GPS_with_habitat.RDS` R data file containing all of the winter Islay GPS fixes for GWfG with the habitats appended from `Landcover Data/Islay landcover data`
- `Derived data/Foraging_Distance_Model_Data.csv`
- `Derived data/GBG_INDIV_SHOOT_DATE_797m.csv`
- `Derived data/GWFG_INDIV_SHOOT_DATE_797m.csv`
- `Landcover Data/Islay landcover data`: habitat data from Islay from the years 2015, 2017, 2018, 2019 amd 2020 as rasters. This was cropped from larger UK wide raster. This data is from UKCEH and wider datsets can be accessed [here](https://www.ceh.ac.uk/data/ukceh-land-cover-maps).
- `Landcover Data/High-Res Coastline`: high resolution coast outline for the UK and Islay as well as bounding boxes for both of these regions.
- `Landcover Data/88090_ISLAY_GMS_FIELD_BOUNDARY`: shapefile of all the agricultural field boundaries on Islay.
- `Landcover Data/Ramsar Outline`: shapefiles of the two RAMSAR areas on Islay. One at Gruinart and one at the Oa.
- `Landcover Data/Roads`: shapefile of all classified public roads in Scotland.
- `MetaData/Tagged bird summary data new.csv` metadata file that contains the sex, ringing location and deployment dates of all tagged GWfG
- `Shooting logs/All_logs_cleaned.csv`: All shooting logs from Islay containing the spatiotemporal information on each shooting event


