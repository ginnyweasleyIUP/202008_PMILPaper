SISAL iHadCM3 model-data comparison last millennium (810CE-1850CE)

All files contain all entities from the SISALv2 database (https://researchdata.reading.ac.uk/256/) that have at least 2 radiometric or laminated dates and ten oxygen isotope measurements within the last millennium. 

This archive includes the following files:
(1) SISAL_HadCM3_ds.csv
(2) SISAL_HadCM3_xnapa_PMIL_yearly.csv
(3) SISAL_HadCM3_xnapb_PMIL_yearly.csv
(4) SISAL_HadCM3_xnapc_PMIL_yearly.csv

File (1) contains the data that is down-sampled to record resolution of all three ensemble members a,b, and c. 

site_id		ID of the cave site
entity_id 	ID of the individual speleothem entity
interp_age	in years BP, original SISALv2 chronology
d18O_measurement	in permil, original SISALv2 measurements
d18O_dw_eq	in permil, drip-water equivalent conversions from d18O_measurements
TEMP_*		in °C, simulated down-sampled temperature extracted at cave site
PREC_*		in kg/ms^2, simulated down-sampled precipitation extracted at cave site
ISOT_*		in permil, simulated down-sampled oxygen isotope ratio of precipitation extracted at cave site
ITPC_*		in permil, simulated down-sampled precipitation weighted oxygen isotope ratio of precipitation extracted at cave site

File (2), (3) and (4) only differ in the ensemble member in the table. The explanation of (2) holds for (3) and (4) respectively for the other ensemble members. Again 

site_id		ID of the cave site
entity_id 	ID of the individual speleothem entity
year_BP		time in years BP
TEMP		in °C, simulated down-sampled temperature extracted at cave site
PREC		in kg/ms^2, simulated down-sampled precipitation extracted at cave site
ISOT		in permil, simulated down-sampled oxygen isotope ratio of precipitation extracted at cave site
ITPC		in permil, simulated down-sampled precipitation weighted oxygen isotope ratio of precipitation extracted at cave site

