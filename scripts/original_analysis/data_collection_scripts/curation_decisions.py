"""
curation_decisions.py
Manual eligibility decisions (reviewer: Claude, AI assistant; review date 2026-10-31)
for INSDC candidate studies with >=1 public WGS/METAGENOMIC run in ENA, and for
MG-RAST candidate projects. Decisions are based on the study title, ENA run
records (library strategy/source, country) and, for core datasets, the
publications. They are first-pass decisions and REQUIRE AUTHOR CONFIRMATION.

Format: accession: (decision, Sample_Category, reason)
decision in {Included, Excluded, Pending verification, Separate contextual dataset}
Any study not listed below but with WGS/METAGENOMIC runs defaults to
Excluded / "Off-target on title review (host-associated, agricultural, marine or
remote environment, engineered process, or method/benchmark study)".
"""
I, E, P, S = "Included", "Excluded", "Pending verification", "Separate contextual dataset"
MW, HW, CW = "Municipal wastewater", "Hospital wastewater", "Community wastewater"
TS, BE, SO, SW, AIR, FOOD, OTH = "Transit surface", "Built environment", "Urban soil", "Surface water", "Air", "Food-related environment", "Other"

CORE = {
 "PRJEB13831": (I, MW, "Core: global sewage surveillance (Munk et al. 2022 'pilot 2016 reproduce' runs; Hendriksen et al. 2019 project). Secondary accession ERP015409."),
 "ERP015409": (I, MW, "Core: secondary accession of PRJEB13831 (merged with PRJEB13831)."),
 "PRJEB27054": (I, MW, "Core: global sewage 2016 pilot round; 227 of 235 runs listed in Munk et al. 2022 Suppl. Data 1; underlies MGnify MGYS00005846."),
 "PRJEB27621": (I, MW, "Core: global sewage June 2017 round (Munk et al. 2022)."),
 "PRJEB40798": (I, MW, "Core: global sewage Nov 2018 round (Munk et al. 2022)."),
 "PRJEB40815": (I, MW, "Core: global sewage Nov 2017 round (Munk et al. 2022)."),
 "PRJEB40816": (I, MW, "Core: global sewage June 2018 round (Munk et al. 2022)."),
 "PRJEB51229": (I, MW, "Core: global sewage longitudinal 2017-2018 and Tanzania series (Munk et al. 2022)."),
 "PRJNA1415085": (I, CW, "Core: Indian metropolitan open-drain wastewater, 4 cities (Singh et al. 2026); processed data on Zenodo 17590646."),
 "PRJNA1399457": (I, HW, "Core: Dhaka hospital and community wastewater (BCSIR); 4 metagenomes included, 3 isolate genomes excluded."),
 "PRJNA732392": (I, TS, "Core: MetaSUB City Sampling Day 2016-2017 (Danko et al. 2021); SRA subset of published samples."),
 "PRJNA1522440": (S, FOOD, "Food-chain context only: FDA retail lamb; metagenomic and enrichment-based data not distinguishable in run metadata; no urban location."),
}
INCLUDED = {
 # wastewater / sewage
 "PRJDB15925": (MW, "Kolkata sewage project (India)."), "PRJDB20415": (HW, "Japanese hospital sewage."),
 "PRJDB9036": (HW, "Japanese hospital sewage (companion project)."), "PRJDB20715": (HW, "Hospital wastewater sources, Japan."),
 "PRJDB37165": (HW, "Hospital wastewater, Chongqing."), "PRJDB37674": (HW, "Hospital wastewater, Metro Manila."),
 "PRJEB109763": (HW, "Hospital wastewater, Turkey."), "PRJEB13832": (MW, "Local sewage surveillance, Denmark (DTU)."),
 "PRJEB13833": (CW, "Sewage surveillance in Kibera, Nairobi (informal settlement)."), "PRJEB14051": (MW, "Three Swedish sewage treatment plants."),
 "PRJEB60913": (MW, "Paired hospital/community urban wastewater, Scotland (mixed hospital and municipal)."),
 "PRJEB70907": (MW, "Sewage AMR surveillance, South Africa."), "PRJEB79273": (MW, "Urban sewage, Latvian cities."),
 "PRJEB85938": (MW, "AMR in the urban water cycle, Europe."), "PRJNA1077107": (MW, "Urban sewage surveillance, China."),
 "PRJNA1082221": (HW, "Hospital wastewater, Poland."), "PRJNA1091428": (MW, "Kampala wastewater resistome."),
 "PRJNA1120823": (MW, "Sewage resistome, Belem, Brazil."), "PRJNA1168082": (MW, "City-scale sewage network AMR mapping, Spain."),
 "PRJNA1223851": (CW, "Community-scale wastewater resistome, Spain."), "PRJNA1266479": (HW, "Egyptian hospital wastewater (and tap water)."),
 "PRJNA1282049": (MW, "Municipal WWTP functional resistomes, Hungary."), "PRJNA1302737": (MW, "Wastewater influent, Iasi, Romania."),
 "PRJNA1328447": (MW, "Zimbabwe wastewater 2021."), "PRJNA1390513": (MW, "12-month wastewater profiling, Brazil."),
 "PRJNA1475805": (HW, "Hospital and urban wastewater, Dhaka (South Asian; distinct from PRJNA1399457)."),
 "PRJNA1480870": (MW, "Urban wastewater, coastal metropolis."), "PRJNA1494178": (HW, "Hospital wastewater, Paraguay."),
 "PRJNA1511879": (HW, "Hospital wastewater metagenome."), "PRJNA1534091": (HW, "Hospital wastewater ARGs and pathogens."),
 "PRJNA358310": (MW, "Wastewater, Burkina Faso/Cameroon."), "PRJNA400857": (SW, "Urban water samples AMR, Singapore."),
 "PRJNA438174": (MW, "Wastewater AMR, Singapore."), "PRJNA750571": (MW, "Sewage surveillance resistome, Singapore."),
 "PRJNA482680": (MW, "Wastewater metagenome, Canada."), "PRJNA506137": (MW, "Wastewater metagenome, South Korea."),
 "PRJNA532515": (MW, "Urban WWTP, Portugal."), "PRJNA793838": (MW, "Wastewater metagenomics, Italy."),
 "PRJNA904380": (MW, "Wastewater metagenomics, New Zealand."), "PRJNA999488": (HW, "Hospital wastewater disinfection, Poland."),
 "PRJNA413894": (MW, "WWTP effluent resistome, Italy (treated effluent)."), "PRJNA881852": (MW, "Wastewater disinfection resistome, Italy."),
 "PRJNA1329243": (MW, "Full-scale WWTPs, sentinel ARGs, Spain."), "PRJNA1445238": (MW, "Municipal WWTP treatment stages, Israel."),
 "PRJNA718910": (MW, "WWTPs, China."), "PRJNA1031342": (MW, "Sewage metagenome, China."),
 "PRJNA1309100": (MW, "Wastewater and river resistomes, South Africa."), "PRJEB34690": (MW, "Resistome risks of wastewaters/aquatic environments, South Africa."),
 "PRJNA1203039": (CW, "Urban wastewater discharge to rivers, Indian twin cities."), "PRJNA1008397": (CW, "Urban wastewater and Kathajodi river, Cuttack, India."),
 "PRJNA1187900": (CW, "Wastewater, River Mahananda, Siliguri, India."), "PRJEB90074": (SW, "River water and wastewater, UK."),
 # transit / built environment
 "PRJDB14136": (BE, "Urban microbiomes in Japan (matrix to confirm)."), "PRJEB101092": (BE, "Public-place surfaces, Xi'an."),
 "PRJEB27640": (BE, "Controlled vs uncontrolled built environments."), "PRJNA1009672": (TS, "Subway environment, Thailand."),
 "PRJNA1201944": (TS, "Chennai Metro station and train surfaces (India)."), "PRJNA1283065": (TS, "Public transit microbiome, India."),
 "PRJNA301589": (TS, "Urban mass-transit system, USA."), "PRJNA561080": (TS, "gCSD17 subway and urban air (MetaSUB-related; air subset classified per sample)."),
 "PRJNA671748": (BE, "Hong Kong outdoor built-environment surfaces."), "PRJNA1393330": (BE, "Hospital surface microbiome, China."),
 "PRJNA1424689": (BE, "Citywide Liverpool outdoor urban environment."), "PRJNA1053252": (BE, "New Delhi house microbiome (indoor; India)."),
 "PRJNA299404": (BE, "Hospital shower-hose biofilms."), "PRJNA1238272": (BE, "Environmental samples, urban informal settlement (matrix to confirm)."),
 # air
 "PRJEB93969": (AIR, "Airborne ARGs in Nordic urban homes (indoor air)."), "PRJNA1262880": (AIR, "Public air samples, China."),
 "PRJNA1531471": (AIR, "Urban airborne microbiota, Changsha."), "PRJNA527324": (AIR, "Semi-urban/rural outdoor aerosol, Norway."),
 # urban soil
 "PRJEB72286": (SO, "ARGs in soil of a European capital."), "PRJNA1158962": (SO, "Historical and recent urban park soils, Turin."),
 "PRJNA1192737": (SO, "Xiamen urban park."), "PRJNA1335812": (SO, "Urban forest remnant patches, Guiyang."),
 "PRJNA1401303": (SO, "Urban forest soil."), "PRJNA633324": (SO, "Urban vs deep-forest resistome, China."),
 # surface water
 "PRJDB14853": (SW, "River resistome, dry/wet weather, Japan."), "PRJDB41754": (SW, "Five urban lakes, Nanchang."),
 "PRJDB6153": (SW, "Tama River (urban Tokyo)."), "PRJEB22915": (SW, "Anthropic effects on aquatic resistomes, Colombia."),
 "PRJEB37889": (SW, "Polluted urban lake, Bolivia."), "PRJEB41039": (SW, "Contaminated urban lake, Bolivia."),
 "PRJEB57859": (SW, "Tijuana River watershed (US portion)."), "PRJNA1021530": (SW, "Urban mangrove habitats, India."),
 "PRJNA1026655": (SW, "Urban-rural ponds, Bangladesh."), "PRJNA1026961": (SW, "Urban stormwater runoff, USA."),
 "PRJNA1071153": (SW, "Urban stormwater, long-read metagenomes (different read technology)."), "PRJNA1071831": (SW, "WWTP effluent impact on River Iskar, Bulgaria."),
 "PRJNA1091838": (SW, "Tropical urban environment resistome, India (matrix to confirm)."), "PRJNA1095672": (SW, "Urban aquatic environment AMR, Australia."),
 "PRJNA1153130": (SW, "Chinese urban rivers."), "PRJNA1240286": (SW, "Urban river water and sediment, Costa Rica."),
 "PRJNA1364797": (SW, "Three urban wetlands, China."), "PRJNA1403057": (SW, "Urban freshwater ecosystem, Manitoba."),
 "PRJNA1423341": (SW, "Tijuana river ARGs."), "PRJNA1480072": (SW, "Urban rivers, Changzhou."),
 "PRJNA1512102": (SW, "Urban ponds, Shenzhen."), "PRJNA1527017": (SW, "Urban river receiving reclaimed water, Ningbo."),
 "PRJNA267173": (SW, "Urban waterway sediments, Singapore."), "PRJNA531081": (SW, "Yamuna sediment, Wazirabad, New Delhi."),
 "PRJNA531097": (SW, "Yamuna sediment, Faizupur Khadar, New Delhi."), "PRJNA551779": (SW, "Urban river sediments, Zhengzhou."),
 "PRJNA562643": (SW, "Urban watershed AMR, USA."), "PRJNA864791": (SW, "Cunha canal, Rio de Janeiro."),
 "PRJNA940122": (SW, "Gowanus Canal Superfund site, Brooklyn."), "PRJNA976879": (SW, "Kanda River, urban Tokyo."),
}
SEPARATE = {
 "PRJEB63604": (FOOD, "Food-associated resistome in processing and production environments, Ireland."),
 "PRJNA941197": (FOOD, "Food processing environment metagenomes."), "PRJNA340876": (FOOD, "Food contamination metagenome."),
 "PRJNA1140288": (FOOD, "Farm-to-fork poultry microbiomes/resistomes, Maputo City."),
 "PRJNA1417753": (FOOD, "Italian slaughterhouse AMR and ESKAPE pathogens."), "PRJEB2763": (FOOD, "Marinated vs unmarinated broiler meat."),
 "PRJNA298260": (FOOD, "Beef and dairy metagenomes."), "PRJNA1068433": (FOOD, "ARGs in Chinese traditional fermented foods."),
 "PRJNA1274045": (FOOD, "Radish sprouts irrigated with treated wastewater."), "PRJNA593573": (FOOD, "Leafy vegetables from an urban farm, Singapore."),
 "PRJEB85115": (OTH, "Food waste metagenome, India."), "PRJNA1444737": (OTH, "Urban compost, New York City farms."),
 "PRJNA1522458": (OTH, "Urban compost piles."), "PRJNA988310": (OTH, "Aerobic home compost."),
 "PRJNA1484761": (OTH, "Chicken manure and associated wastewater, Pakistan (agricultural, not urban)."),
}
PENDING = {
 "PRJEB34634": "Mixed One Health design (farm animals, human/animal sewage, sewage treatment works); sewage subset needs identification.",
 "PRJEB90275": "Microplastics from hospital wastewater to sea; plastisphere matrix.", "PRJEB11494": "Single UK wastewater run; minimal metadata.",
 "PRJNA1346432": "Activated sludge pathogen/ARG reduction (process study).", "PRJNA1099934": "Pipeline biofilm and wastewater, China.",
 "PRJDB18077": "Extracellular resistome in a WWTP (method/process study).", "PRJDB19716": "Intra/extracellular DNA in a WWTP (method study).",
 "PRJNA1514836": "Cross-sectoral AMR surveillance, Australia (mixed sectors).", "PRJNA1064560": "Citizen-science soil antibiotics study (urban status unclear).",
 "PRJNA1236264": "'Urbanization metagenomic analysis' - matrix unstated.", "PRJEB72287": "ARG soils across countries (urban status unclear).",
 "PRJEB87523": "Riverine metagenome, Spain (urban status unclear).", "PRJNA1066593": "Anthropogenic pollution and aquatic resistome, Italy.",
 "PRJNA1067305": "Urbanized-to-marine protected area transect (coastal).", "PRJNA1255235": "Urban and island beaches, Belem.",
 "PRJNA1330252": "Hong Kong urban rivers incl. microplastic and natural particles.", "PRJNA1432522": "Eutrophication and phthalates in urban water (possible experiment).",
 "PRJNA892917": "Lockwitzbach river (Dresden) microbiome/resistome; urban status unclear.", "PRJNA820475": "Binan River, Philippines.",
 "PRJNA823379": "Laguna Lake tributaries, Philippines.", "PRJNA910298": "Eutrophic South Australian surface waters.",
 "PRJNA1517594": "Waterborne pathogens, coastal Bangladesh (urban status unclear).", "PRJEB6102": "Industrially polluted lake, India (industrial effluent).",
 "PRJNA1399168": "Foraging grounds of urban bats (environmental, matrix mixed).",
}
EXCLUDED_SPECIFIC = {
 "PRJEB41153": "Experimental spike/release design.", "PRJEB41171": "Plasmid-enriched (plasmidome) libraries; different data type.",
 "PRJEB111923": "Phage-targeted study.", "PRJEB113830": "WWTP process-engineering study.",
 "PRJNA593593": "Biofuel bioreactor series (JGI).", "PRJNA593594": "Biofuel bioreactor series (JGI).",
 "PRJNA1447082": "Virus-enriched metagenomes.", "PRJNA681008": "Controlled cleanroom built environment (not urban public).",
 "PRJNA1016375": "Experimental indoor green-wall growth media.", "PRJNA975673": "Farm-exposure house dust (rural focus).",
 "PRJNA1516236": "Industrial smelting-area PM2.5.", "PRJNA1228129": "Aircraft filters and face masks.",
 "PRJNA1298337": "Continental background aerosol; mixed metatranscriptomic design.", "PRJNA1223855": "Culture plate-wash (enrichment) of soil.",
 "PRJNA1028102": "Plastisphere.", "PRJNA1219086": "Plastisphere.", "PRJNA1215273": "Microplastic biofilms.",
 "PRJNA254745": "Viral metagenomics of meat.", "PRJNA564982": "Enrichment-culture metagenomes.",
 "PRJEB27316": "Proficiency/ring trial.", "PRJEB29454": "Proficiency/ring trial.", "PRJEB52722": "Sequencing-method benchmark.",
 "PRJNA746311": "Capture-bait method study.", "PRJNA484080": "Database/platform project.",
 "PRJEB12831": "Archaeological samples.", "PRJEB30834": "Human faecal metagenomes.", "PRJNA349463": "Human microbiome urban/rural comparison.",
 "PRJNA1356541": "Child faecal strain tracking.", "PRJNA300541": "Human/animal faecal microbiota.", "PRJEB9661": "Wild rat microbiota.",
 "PRJNA1365129": "Fly vectors.", "PRJNA385554": "Blowfly/housefly microbiomes.", "PRJNA377922": "Banknote microbiota (fomite; not an environmental matrix in scope).",
}
MGRAST = {
 "mgp16346": (I, SO, "GLUSEEN global urban soil network."),
 "mgp89046": (I, MW, "Swiss wastewater and river bacteriomes (NRP72 river resistome); per-sample matrix to separate."),
 "mgp90475": (I, MW, "Mobile antibiotic resistome in WWTPs, Tychy 2018."), "mgp90740": (I, MW, "Mobile antibiotic resistome in WWTPs, Olsztyn 2018."),
 "mgp92067": (I, MW, "Mobile antibiotic resistome in WWTPs, Tychy 2019."), "mgp92140": (I, MW, "Mobile antibiotic resistome in WWTPs, Olsztyn 2018-19."),
 "mgp93181": (I, MW, "Mobile antibiotic resistome in WWTPs, Mix20."), "mgp87891": (P, BE, "Showerhead biofilms (built environment); urban status unclear."),
 "mgp6602": (E, "", "Not accessible: MG-RAST API returned HTTP 401 (not public)."),
}
DEFAULT_EXCLUDED_REASON = ("Off-target on title review (host-associated, agricultural, marine or remote environment, "
                           "engineered process, or method/benchmark study).")

def decide(acc):
    if acc in CORE: return CORE[acc]
    if acc in INCLUDED: return (I,) + INCLUDED[acc]
    if acc in SEPARATE: return (S,) + SEPARATE[acc]
    if acc in PENDING: return (P, "Unclassified", PENDING[acc])
    if acc in EXCLUDED_SPECIFIC: return (E, "", EXCLUDED_SPECIFIC[acc])
    return (E, "", DEFAULT_EXCLUDED_REASON)
