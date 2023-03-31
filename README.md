# Contrasting effects of shooting disturbance on the movement and behaviour of sympatric wildfowl species 🔫 🦆
The impact of shooting on time-activity budgets, movement and habitat selection in Barnacle Geese *Branta leucopsis* and White-fronted Geese *Anser albifrons flavirostris* on Islay, Scotland. 

## Authors 🖊️
- Luke Ozsanlav-Harris <a itemprop="sameAs" content="https://orcid.org/0000-0003-3889-6722" href="https://orcid.org/0000-0003-3889-6722" target="orcid.widget" rel="noopener" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" target="_blank" style="width:1em;margin-right:.5em;"/></a>
- Aimee L. S. Mcintosh <a itemprop="sameAs" content="https://orcid.org/0000-0002-4975-3682" href="https://orcid.org/0000-0002-4975-3682" target="orcid.widget" rel="noopener" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" target="_blank" style="width:1em;margin-right:.5em;"/></a>
- Larry R. Griffin
- Geoff M. Hilton <a itemprop="sameAs" content="https://orcid.org/0000-0001-9062-3030" href="https://orcid.org/0000-0001-9062-3030" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
- Jessica M. Shaw <a itemprop="sameAs" content="https://orcid.org/0000-0003-0862-9260" href="https://orcid.org/0000-0003-0862-9260" target="orcid.widget" rel="noopener" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" target="_blank" style="width:1em;margin-right:.5em;"/></a>
- Stuart Bearhop <a itemprop="sameAs" content="https://orcid.org/0000-0002-5864-0129" href="https://orcid.org/0000-0002-5864-0129" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>

## Code description 👨‍💻
- `R code/1_Combine tracking data + label with landcover.R`: Prepare GWfG tracking data for analysis
- `R code/2_Shooting proximity effects on Distance travelled.R`: Modelling the immediate movement in reponse to shooting
- `R code/4_Daily behaviour on shooting vs non-shooting days.R`: Modelling behavioral changes in response to shooting disturbance (GWfG only)
- `R code/5_Daily ODBA on shooting vs non-shooting days.R`; Modelling change in daily ODBA due to shooting disturbance
- `R code/X_Create maps for publication.R`: Create maps of Islay, including habitat and roads
- `R code/Useful Functions` folder containing useful code for this project, it is not used in the main work flow

## Data description 📊
- `Landcover Data/Islay landcover data`: habitat data from Islay from the years 2015, 2017, 2018, 2019 amd 2020 as rasters. This was cropped from larger UK raster using `1_Crop RAW landcover data.R`
- `Landcover Data/High-Res Coastline`: high resultion coast outline for the UK and Islay and bounding boxes for both regions
- `Landcover Data/88090_ISLAY_GMS_FIELD_BOUNDARY` shapefile of all the agricultural field boundaries on Islay
- `Landcover Data/Ramsar Outline`: Shapefiles of the two RAMSAR areas on Islay. One at Gruinart and one at the Oa
- `Landcover Data/Roads`: shapefile of all classified roads in Scotland
- `MetaData/Tagged bird summary data new.csv` metadata file that contains the sex, ringing location and deployment dates of all tagged GWfG
- `Shooting logs/All_logs_cleaned.csv`: All cleaned shooting logs from Islay containing the information on each shooting event
- `Derived data/All_winter_GPS_with_habitat.RDS` R data file containing all of the winter Islay GPS fixes for GWfG with the habitats appended from `Landcover Data/Islay landcover data`

