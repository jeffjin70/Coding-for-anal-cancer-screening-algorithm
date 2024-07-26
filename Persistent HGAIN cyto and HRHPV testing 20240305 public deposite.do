**code written for Stata

version 13
clear all
macro drop _all
set linesize 120
set more off
capture log close

**data file available upon request with review
sort idno nround
drop cohrhpv

**merge with biomarker results
merge 1:1 idno round using "p:\hepp\4. Completed Projects\spanc\spanc data\hpv results\biomarker"
drop _merge

merge 1:1 idno round using "p:\hepp\4. Completed Projects\spanc\spanc data\cyto results\dualstain"
drop _merge

merge 1:1 idno nround using "p:\hepp\4. Completed Projects\spanc\spanc data\hpv results\mrnav3", update
drop _merge

merge 1:1 idno nround using "p:\hepp\4. Completed Projects\spanc\spanc data\hpv results\viral load v3", update
drop _merge

merge 1:1 idno round using "p:\hepp\4. Completed Projects\spanc\spanc data\EuroAssay results\EuroAssay Results"
drop _merge

drop if cyto==. & idno~=2042 //delete participants who did not undergo cyto and histo

egen ascc=max(hist), by(idno)
/*drop if ascc==6
drop if idno==1134 | idno==1138 | idno==1462 //refered for treatment of HSIL*/

gen exc=1 if ascc==6
replace exc=1 if exc==. & (idno==1134 | idno==1462 | idno==1138)

sort idno dov

drop if nround==2
sort idno nround
egen elig=count(hgain), by(idno)
recode elig (1=0) (2/4=1)

stset dov, id(idno) scale(365.25) failure(hgain==0) origin(round==1) enter(hgain==1)
egen clhgain=max(_d), by(idno)

by idno: gen perhgain=1 if (elig==1 & hgain==1 & hgain[_n+1]==1) & round==1
replace perhgain=0 if elig==1 & perhgain==. & round==1

gen plsil=(ncyto>=2) if !missing(ncyto)
sort idno nround
by idno: gen p2lsil=plsil[_n+1] if round==1
gen lsil=(ncyto>=3) if !missing(ncyto)
gen phsil=(ncyto>=4) if !missing(ncyto)
by idno: gen p2hsil=phsil[_n+1] if round==1

by idno: gen perplsil=1 if plsil==1 & p2lsil==1 & round==1
replace perplsil=0 if perplsil==. & !missing(plsil) & !missing(p2lsil)

by idno: gen perphsil=1 if phsil==1 & p2hsil==1 & round==1
replace perphsil=0 if perphsil==. & !missing(phsil) & !missing(p2hsil)

egen cyhpv16=group(plsil hpv16)

recode cyhpv16 (1/3=0) (4=1), g(tricyhpv16)
recode cyhpv16 (1=0) (2/4=1), g(cocyhpv16)

diagt perhgain hpv16 if !missing(ncyto), t
diagt perhgain tricyhpv16, t
diagt perhgain cocyhpv16, t

diagt perhgain hpv16 if !missing(ncyto) & rpthiv==0, t
diagt perhgain hpv16 if !missing(ncyto) & rpthiv==1, t
diagt perhgain tricyhpv16 if rpthiv==0, t
diagt perhgain tricyhpv16 if rpthiv==1, t
diagt perhgain cocyhpv16 if rpthiv==0, t
diagt perhgain cocyhpv16 if rpthiv==1, t

local ohrlist "18 31 33 35 39 45 51 52 56 58 59 68"
sort idno round
gen perohrhpv=0 if !missing(hpv16) & round==1
foreach var in `ohrlist' {
    by idno: replace perohrhpv=1 if hpv`var'==1 & hpv`var'[_n+1]==1 & round==1
}

local ohrlist "18 31 33 35 39 45 51 52 56 58 59 68"
sort idno round
foreach var in `ohrlist' {
	gen perhpv`var'=0 if !missing(hpv16) & round==1
    by idno: replace perhpv`var'=1 if hpv`var'==1 & hpv`var'[_n+1]==1 & round==1
}

replace perohrhpv=0 if perohrhpv==1 & hpv16==1
egen cyohrhpv=group(plsil perohrhpv)
recode cyohrhpv (1/3=0) (4=1), g(tricyohrhpv)
recode cyohrhpv (1=0) (2/4=1), g(cocyohrhpv)

local ohrmrnalist "18 31 33 45"
sort idno round
gen perohrhpvmrna=0 if !missing(hpv16mrna) & round==1
foreach var in `ohrmrnalist' {
    by idno: replace perohrhpvmrna=1 if hpv`var'mrna==1 & hpv`var'mrna[_n+1]==1 & round==1
}

replace perohrhpvmrna=0 if perohrhpvmrna==1 & hpv16mrna==1


diagt perhgain perohrhpv if !missing(ncyto) & hpv16==0, t
diagt perhgain tricyohrhpv if hpv16==0, t
diagt perhgain cocyohrhpv if hpv16==0, t

diagt perhgain perohrhpv if !missing(ncyto) & hpv16==0 & rpthiv==0, t
diagt perhgain perohrhpv if !missing(ncyto) & hpv16==0 & rpthiv==1, t
diagt perhgain tricyohrhpv if hpv16==0 & rpthiv==0, t
diagt perhgain tricyohrhpv if hpv16==0 & rpthiv==1, t
diagt perhgain cocyohrhpv if hpv16==0 & rpthiv==0, t
diagt perhgain cocyohrhpv if hpv16==0 & rpthiv==1, t

gen sigscr=1 if (hpv16==1 | perohrhpv==1) & (round==1 & !missing(perhgain) & !missing(ncyto)) //Algorithm #4
replace sigscr=0 if sigscr==. & !missing(perhgain) & !missing(ncyto) & !missing(hpv16)
egen sigcom=group(sigscr plsil)
recode sigcom (1/3=0) (4=1), g(sigtri)
recode sigcom (1=0) (2/4=1), g(sigco)

egen sig2com=group(sigscr perplsil)
recode sig2com (1/3=0) (4=1), g(sig2tri)

replace perhgain=. if exc==1
replace sigscr=. if exc==1

gen inc=(!missing(perhgain) & !missing(perohrhpv) & !missing(perplsil)) if round==1

**gen inc=!missing(sigscr) if round==1

gen bhrhpv=hrhpv if round==1
replace bhrhpv=2 if hpv16==1

label define bhrhpv 0 "HRHPV Neg" 1 "non16HRHPV Pos" 2 "HPV16 Pos"
label value bhrhpv bhrhpv

diagt perhgain sigscr, t
diagt perhgain sigtri, t
diagt perhgain sigco, t

recode cignow 	(.=0 "Never") ///
				(0=1 "Past") (1=2 "Current"), ///
				g(cigstatus)
				
label variable cigstatus "Have you ever been a smoker"

replace smenlif=. if smenlif==97

recode smenlif 	(1/2=0 "less than 10 men") ///
				(4=1 "11-50 men") ///
				(5=3 "51 - 200 men") ///
				(6=4 "201 - 500 men") ///
				(7/8=5 ">500 men") ///
				(else=.), ///
				g(smenliff)

label var smenliff "Grouped: lifetime male partners"

recode smen6m 	(0=0 "None") ///
				(1=1 "only 1 man") ///
				(2=2 "2-5 men") ///
				(3=3 "6-10 men") ///
				(4/max=4 "Over 10 men"), ///
				g(smen6mf)

label var smen6mf "Grouped: Number of partners in the past 6 months"

recode smen6mf (0/1=1) (2=2) (3/4=3), gen(sm6)

**Alternative screening approach HPV16 +/- cyto PLUS persistent non-16HRHPV + cyto
gen alsigtri=sigtri // Algorithm #2
replace alsigtri=1 if alsigtri==0 & hpv16==1 & plsil==0

gen alsig2tri=sig2tri
replace alsig2tri=1 if alsig2tri==0 & hpv16==1 & perplsil==0

**Alternative screening approach HPV16+/- cyto PLUS non16HRHPV baseline only + cyto
gen al2sigtri=1 if hrhpv==1 & plsil==1 & round==1
replace al2sigtri=1 if al2sigtri==. & hpv16==1 & round==1
replace al2sigtri=0 if al2sigtri==. & round==1 & !missing(perhgain) & !missing(ncyto)
replace al2sigtri=. if inc==0

gen al2sig2tri=1 if hrhpv==1 & perplsil==1 & round==1
replace al2sig2tri=1 if al2sig2tri==. & hpv16==1 & round==1
replace al2sig2tri=0 if al2sig2tri==. & round==1 & !missing(perhgain) & !missing(perplsil)
replace al2sig2tri=. if alsig2tri==.

**Alternative screening approach HPV16+LSIL PLUS persistent non16HRHPV+PLSIL
gen al3sigtri=sigtri 
replace al3sigtri=0 if al3sigtri==1 & (hpv16==1 & ncyto==2)


**Alternative sreening approach PHSIL for cytology cut-off
egen alsigcom=group(sigscr phsil)
recode alsigcom (1/3=0) (4=1), g(al4sigtri)

egen alsig2com=group(sigscr perphsil)
recode alsig2com (1/3=0) (4=1), g(al4sig2tri) //Algorithm #5

**Alternative sreening approach PHSIL for cytology cut-off HPV16 +/- cytology
gen al5sigtri=al4sigtri //Algorithm #3
replace al5sigtri=1 if al5sigtri==0 & (hpv16==1 & phsil==0)

gen al5sig2tri=al4sig2tri
replace al5sig2tri=1 if al5sig2tri==0 & (hpv16==1 & perphsil==0)

**Alternatvie screening approch HPV16+/- cyto PLUS non16HRHPV baesline only+PHSIL
gen ohrhpv=1 if (hrhpv==1 & hpv16==0) 
replace ohrhpv=0 if missing(ohrhpv) & !missing(hpv16)
gen al6sigtri=1 if hpv16==1 | (ohrhpv==1 & phsil==1)
replace al6sigtri=0 if missing(al6sigtri) & !missing(sigtri)
replace al6sigtri=. if missing(sigtri)

**Alternative screening apporach HPV16+ /- cytyo OR non16HRHPV baseline only+PHSIL OR persistent non16HPHPV+PLSIL
gen al7sigtri=al6sigtri //Algorithm #1
replace al7sigtri=1 if al7sigtri==0 & perohrhpv==1 & plsil==1



gen hpvscr=1 if hpv16==1 & inc==1
replace hpvscr=2 if perohrhpv==1 & inc==1
replace hpvscr=0 if hpvscr==. & inc==1

label define hpvscr 0 "HPV screen neg" 1 "HPV16 Pos" 2 "Persist oHRHPV"
label value hpvscr hpvscr

sort idno nround

preserve

drop if nround==3

by idno: gen per2hgain=1 if (elig==1 & hgain==1 & hgain[_n+1]==1) & round==1

by idno: replace per2hgain=0 if elig==1 & per2hgain==. & round==1 & !missing(hgain[_n+1])

keep idno per2hgain nround

sort idno nround
 
save temp, replace

restore

merge 1:1 idno nround using temp
drop _merge

gen per2yrhgain=(perhgain==1 & per2hgain==1) if !missing(per2hgain)



