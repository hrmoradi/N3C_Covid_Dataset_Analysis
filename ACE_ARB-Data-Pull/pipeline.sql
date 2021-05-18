

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7"),
    Covid_positive_persons_LDS=Input(rid="ri.foundry.main.dataset.a917b9e7-b891-4ac5-93dc-de1a3b93c7ad"),
    compiled_gfr=Input(rid="ri.foundry.main.dataset.76d20b33-17ef-4566-9c80-2d19872d35cc"),
    complete_patient_table_with_derived_scores=Input(rid="ri.foundry.main.dataset.d467585c-53b3-4f10-b09e-8e86a4796a1a"),
    covid_medications=Input(rid="ri.foundry.main.dataset.3116094b-57e5-44ba-908c-4494c1883c56"),
    first_diag_table_pos=Input(rid="ri.foundry.main.dataset.d3b1dfd4-9a7c-4fa1-b29a-de58ed191d1b"),
    selected_critical_visits=Input(rid="ri.foundry.main.dataset.a8a5404a-760c-4f8e-bee2-76ccc7b5b781"),
    summary_Table=Input(rid="ri.foundry.main.dataset.5dd50df8-79ce-461c-9397-0a429352f52a")
)
SELECT     c. person_id, 
           c.visit_concept_name,
           c.Severity_Type,
           c.age_at_visit_start_in_years_int AS age_pat,
           CASE
                      WHEN cv.ACE = 1 THEN 1
                      WHEN cv.ARB = 1 THEN 1
                      ELSE 0
           END AS ACE_ARB,
           CASE
                      WHEN c.Race = 'Black or African American' THEN 1
                      ELSE 0
           END AS race_AA,
           c.Ethnicity,
           c.gender_concept_name AS gender_pat,
           CASE 
                      WHEN c.gender_concept_name = 'Female' THEN 0
                      WHEN c.gender_concept_name = 'Male' THEN 1
                      WHEN c.gender_concept_name = 'Other' THEN 2
           END AS gender,
           CASE
                      WHEN c.Severity_Type = 'Dead_w_COVID' THEN 1
                      WHEN c.in_death_table = 1 THEN 1
                      ELSE 0
           END AS mort_hospice,
           CASE
                      WHEN st.diabetes = 1 THEN 1
                      WHEN st.dmcx = 1 THEN 1
                      ELSE 0
           END AS hx_DM,
           st.CHF AS hx_CHF,
           st.hypertension AS hx_HTN,
           cl.GFR, ---FLOOR(cl.GFR/10.00) * 10 as GFR
           MONTHS_BETWEEN( trunc(cp.date_of_earliest_covid_diagnosis, 'MM'), trunc(fdp.first_diag , 'MM'  ) ) as month_diag -- trunc(cp.date_of_earliest_covid_diagnosis, 'MM')
FROM       complete_patient_table_with_derived_scores c
INNER JOIN Covid_positive_persons_LDS cp
ON         cp.person_id = c.person_id
INNER JOIN selected_critical_visits s
ON         s.person_id = c.person_id
INNER JOIN summary_Table st
ON         st.person_id = c.person_id
LEFT JOIN  covid_medications cv
ON         c.person_id = cv.person_id
LEFT JOIN  compiled_gfr cl ---ecreatinine_Labs cl --- using eGFR table with 366k rows instead of GFR 255k
ON         c.person_id = cl.person_id
cross join first_diag_table_pos as fdp
WHERE      s.visit_concept_name = 'Inpatient Visit'

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4"),
    compiled_gfr_neg=Input(rid="ri.foundry.main.dataset.de0b8ba5-9bca-47f1-9e40-8e233ed3c6ad"),
    complete_patient_table_with_derived_scores=Input(rid="ri.foundry.main.dataset.d467585c-53b3-4f10-b09e-8e86a4796a1a"),
    cov_med_365_neg=Input(rid="ri.foundry.main.dataset.cd6a4e5d-4cdc-4438-be6c-33aa604a600b"),
    first_diag_table_neg=Input(rid="ri.foundry.main.dataset.cef5e6f6-2e13-4ab2-adb0-38439afa5b32"),
    neg_earliest=Input(rid="ri.foundry.main.dataset.19034815-c3a7-43e6-85f3-95725a63b850"),
    selected_critical_visits=Input(rid="ri.foundry.main.dataset.a8a5404a-760c-4f8e-bee2-76ccc7b5b781"),
    summary_Table_neg=Input(rid="ri.foundry.main.dataset.7487019f-58ba-4147-a1b3-952fcca48c10")
)
SELECT     c. person_id, 
           c.visit_concept_name,
           c.Severity_Type,
           c.age_at_visit_start_in_years_int AS age_pat,
           CASE
                      WHEN cv.ACE = 1 THEN 1
                      WHEN cv.ARB = 1 THEN 1
                      ELSE 0
           END AS ACE_ARB,
           CASE
                      WHEN c.Race = 'Black or African American' THEN 1
                      ELSE 0
           END AS race_AA,
           c.Ethnicity,
           c.gender_concept_name AS gender_pat,
           CASE 
                      WHEN c.gender_concept_name = 'Female' THEN 0
                      WHEN c.gender_concept_name = 'Male' THEN 1
                      WHEN c.gender_concept_name = 'Other' THEN 2
           END AS gender,
           CASE
                      WHEN c.Severity_Type = 'Dead_w_COVID' THEN 1
                      WHEN c.in_death_table = 1 THEN 1
                      ELSE 0
           END AS mort_hospice,
           CASE
                      WHEN st.diabetes = 1 THEN 1
                      WHEN st.dmcx = 1 THEN 1
                      ELSE 0
           END AS hx_DM,
           st.CHF AS hx_CHF,
           st.hypertension AS hx_HTN,
           cl.GFR, ---FLOOR(cl.GFR/10.00) * 10 as GFR
           MONTHS_BETWEEN( trunc(cp.date_of_earliest_covid_diagnosis, 'MM'), trunc(fdn.first_diag , 'MM'  ) ) as month_diag -- trunc(cp.date_of_earliest_covid_diagnosis, 'MM')
FROM       complete_patient_table_with_derived_scores c
INNER JOIN neg_earliest cn
ON         cn.person_id = c.person_id
INNER JOIN selected_critical_visits s
ON         s.person_id = c.person_id
INNER JOIN summary_Table_neg st
ON         st.person_id = c.person_id
LEFT JOIN  cov_med_365_neg cv
ON         c.person_id = cv.person_id
LEFT JOIN  compiled_gfr_neg cl ---ecreatinine_Labs cl --- using eGFR table with 366k rows instead of GFR 255k
ON         c.person_id = cl.person_id
cross join first_diag_table_neg fdn
WHERE       s.visit_concept_name = 'Inpatient Visit' --  cl.GFR is not null -- and

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.39ec90ea-6062-46e1-bd5a-6a22d7d7aca1"),
    cs=Input(rid="ri.foundry.main.dataset.e670c5ad-42ca-46a2-ae55-e917e3e161b6"),
    selected_Medication_Codesets=Input(rid="ri.foundry.main.dataset.9d244035-667e-4f81-ab43-c257ea44f635")
)
SELECT
    s.concept_Set_Alias AS Alias,
    c.concept_name,
    c.concept_id
FROM cs c
INNER JOIN selected_Medication_Codesets s ON s.codeset_ID = c.codeset_id
WHERE is_most_recent_version = true 

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.7c91b495-099c-4f77-8174-68524facec95"),
    COVID19_Med_Concepts=Input(rid="ri.foundry.main.dataset.39ec90ea-6062-46e1-bd5a-6a22d7d7aca1"),
    drug_exposure=Input(rid="ri.foundry.main.dataset.ec252b05-8f82-4f7f-a227-b3bb9bc578ef")
)
SELECT
    a.person_id,
    a.drug_exposure_start_date,
    COALESCE(a.drug_exposure_end_date, a.drug_exposure_start_date) drug_exposure_end_date,
    a.quantity,
    a.days_supply,
    a.visit_occurrence_id,
    b.Alias,
    b.concept_name
FROM drug_exposure a
inner join (
    select 
        Alias, concept_id, concept_name 
    from COVID19_Med_Concepts
) b on a.drug_concept_id = b.concept_id
where a.drug_exposure_start_date is not null

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.c7a605e0-17bc-4765-8684-9876c0e382bb"),
    COVID19_Medications_both=Input(rid="ri.foundry.main.dataset.7c91b495-099c-4f77-8174-68524facec95"),
    neg_earliest=Input(rid="ri.foundry.main.dataset.19034815-c3a7-43e6-85f3-95725a63b850")
)
select distinct
    person_id,
    Alias, 
    days_covid_to_med
from (
    SELECT
        a.person_id,
        a.Alias,
        min(datediff(a.drug_exposure_start_date, b.date_of_earliest_covid_diagnosis)) days_covid_to_med
    FROM COVID19_Medications_both a
    inner join neg_earliest b on a.person_id = b.person_id
    where a.drug_exposure_end_date >= b.date_of_earliest_covid_diagnosis
    group by a.person_id, b.person_id, a.Alias
) x
where days_covid_to_med <= 365

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.cd6a4e5d-4cdc-4438-be6c-33aa604a600b"),
    cov_med_365_both=Input(rid="ri.foundry.main.dataset.c7a605e0-17bc-4765-8684-9876c0e382bb")
)
SELECT
    person_id,
    sum(case when Alias = 'Amiodarone gtt' then 1 else 0 end) Amiodarone_gtt,
    sum(case when Alias = 'Anakinra' then 1 else 0 end) Anakinra,
    sum(case when Alias = 'Azithromycin' then 1 else 0 end) Azithromycin,
    sum(case when Alias = 'Chloroquine' then 1 else 0 end) Chloroquine,
    sum(case when Alias = 'Dexamethasone' then 1 else 0 end) Dexamethasone,
    sum(case when Alias = 'dialysis CRRT/HD' then 1 else 0 end) dialysis_CRRT_HD,
    sum(case when Alias = 'Dobutamine gtt' then 1 else 0 end) Dobutamine_gtt,
    sum(case when Alias = 'Dopamine gtt' then 1 else 0 end) Dopamine_gtt,
    sum(case when Alias = 'Epinephrine gtt' then 1 else 0 end) Epinephrine_gtt,
    sum(case when Alias = 'Epoprostenol' then 1 else 0 end) Epoprostenol,
    sum(case when Alias = 'Esmolol gtt' then 1 else 0 end) Esmolol_gtt,
    sum(case when Alias = 'Hydrocortisone' then 1 else 0 end) Hydrocortisone,
    sum(case when Alias = 'Hydroxychloroquine' then 1 else 0 end) Hydroxychloroquine,
    sum(case when Alias = 'Inhaled Nitric Oxide' then 1 else 0 end) Inhaled_Nitric_Oxide,
    sum(case when Alias = 'intravenous immunoglobulin' then 1 else 0 end) intravenous_immunoglobulin,
    sum(case when Alias = 'Isoproterenol gtt' then 1 else 0 end) Isoproterenol_gtt,
    sum(case when Alias = 'Levosimendan gtt' then 1 else 0 end) Levosimendan_gtt,
    sum(case when Alias = 'Lopinavir' then 1 else 0 end) Lopinavir,
    sum(case when Alias = 'Lopinavir/Ritonavir combination' then 1 else 0 end) Lopinavir_Ritonavir_combination,
    sum(case when Alias = 'Methylprednisolone' then 1 else 0 end) Methylprednisolone,
    sum(case when Alias = 'Milrinone gtt' then 1 else 0 end) Milrinone_gtt,
    sum(case when Alias = 'Norepinephrine gtt' then 1 else 0 end) Norepinephrine_gtt,
    sum(case when Alias = 'Phenylephrine gtt' then 1 else 0 end) Phenylephrine_gtt,
    sum(case when Alias = 'Plasma' then 1 else 0 end) Plasma,
    sum(case when Alias = 'Prednisolone' then 1 else 0 end) Prednisolone,
    sum(case when Alias = 'Prednisone' then 1 else 0 end) Prednisone,
    sum(case when Alias = 'Remdesivir' then 1 else 0 end) Remdesivir,
    sum(case when Alias = 'Ritonavir' then 1 else 0 end) Ritonavir,
    sum(case when Alias = 'Tocilizumab' then 1 else 0 end) Tocilizumab,
    sum(case when Alias = 'Vasopressin gtt' then 1 else 0 end) Vasopressin_gtt,
    sum(case when Alias = 'ACE' then 1 else 0 end) ACE,
    sum(case when Alias = 'ARB' then 1 else 0 end) ARB,
    sum(case when Alias = 'Fluvoxamine' then 1 else 0 end) Fluvoxamine,
    sum(case when Alias = 'Clarithromycin' then 1 else 0 end) Clarithromycin,
    sum(case when Alias = 'Erythromycin' then 1 else 0 end) Erythromycin,
    sum(case when Alias = 'Josamycin' then 1 else 0 end) Josamycin,
    sum(case when Alias = 'Macrolides' then 1 else 0 end) Macrolides,
    sum(case when Alias = 'Midecamycin' then 1 else 0 end) Midecamycin,
    sum(case when Alias = 'Miocamycin' then 1 else 0 end) Miocamycin,
    sum(case when Alias = 'Prednisone/Methylprednisolone' then 1 else 0 end) Prednisone_Methylprednisolone,
    sum(case when Alias = 'Rokitamycin' then 1 else 0 end) Rokitamycin,
    sum(case when Alias = 'Roxithromycin' then 1 else 0 end) Roxithromycin,
    sum(case when Alias = 'Spiramycin' then 1 else 0 end) Spiramycin,
    sum(case when Alias = 'troleandomycin' then 1 else 0 end) troleandomycin
FROM cov_med_365_both
group by person_id

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.cef5e6f6-2e13-4ab2-adb0-38439afa5b32"),
    neg_earliest=Input(rid="ri.foundry.main.dataset.19034815-c3a7-43e6-85f3-95725a63b850")
)
SELECT min(date_of_earliest_covid_diagnosis) as first_diag
FROM neg_earliest

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.d3b1dfd4-9a7c-4fa1-b29a-de58ed191d1b"),
    Covid_positive_persons_LDS=Input(rid="ri.foundry.main.dataset.a917b9e7-b891-4ac5-93dc-de1a3b93c7ad")
)
SELECT min(date_of_earliest_covid_diagnosis) as first_diag
FROM Covid_positive_persons_LDS

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.19034815-c3a7-43e6-85f3-95725a63b850"),
    Covid_positive_persons_LDS=Input(rid="ri.foundry.main.dataset.a917b9e7-b891-4ac5-93dc-de1a3b93c7ad"),
    covid_19_negative_patients=Input(rid="ri.foundry.main.dataset.6dbbbc52-960c-4e0d-954b-d343eca1ecff")
)
SELECT neg.person_id, Min(neg.measurement_date) as date_of_earliest_covid_diagnosis
FROM 

covid_19_negative_patients neg 
left join Covid_positive_persons_LDS pos 
on neg.person_id = pos.person_id

where pos.person_id is null
group by neg.person_id 

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.69fda5b2-1ddc-4e20-9d18-ad1c44019fe6"),
    co=Input(rid="ri.foundry.main.dataset.900fa2ad-87ea-4285-be30-c6b5bab60e86"),
    cs=Input(rid="ri.foundry.main.dataset.e670c5ad-42ca-46a2-ae55-e917e3e161b6")
)
SELECT
distinct
x.person_id , cast( x.hypertension as INTEGER), cast( x.upper_gi_bleed as INTEGER), cast( x.MI as INTEGER) , cast( x.CHF as INTEGER) , cast( x.PVD as INTEGER) , cast( x.stroke as INTEGER) , cast( x.dementia as INTEGER) , cast( x.pulmonary as INTEGER) ,
cast( x.rheumatic as INTEGER) , cast( x.PUD as INTEGER) , cast( x.liver_mild as INTEGER) , cast( x.liversevere as INTEGER) , cast( x.diabetes as INTEGER) ,  cast( x.dmcx as INTEGER) ,
 cast( x.paralysis as INTEGER) , cast( x.renal as INTEGER) , cast( x.cancer as INTEGER) ,  cast( x.mets as INTEGER) ,  cast( x.hiv as INTEGER) , cast( x.multiple as INTEGER) 
 --- , c.CCI_INDEX
--(x.MI*1 + x.CHF*1 + x.PVD*1 + x.stroke*1 + x.dementia*1 + x.pulmonary*1 + x.rheumatic*1 + x.PUD*1 + x.liver_mild*1 + x.diabetes*1 + x.dmcx*2 + x.paralysis*2 + x.renal*2 + x.cancer*2 + x.liversevere*3 + x.mets*6 + x.hiv*6) CCI_INDEX

FROM
(
SELECT
distinct
    person_id, 
    sum(case when comorbidity = 'MI' then 1 else 0 end) MI ,
    sum(case when comorbidity = 'CHF' then 1 else 0 end) CHF ,
    sum(case when comorbidity = 'PVD' then 1 else 0 end) PVD ,
    sum(case when comorbidity = 'Stroke' then 1 else 0 end) stroke ,
    sum(case when comorbidity = 'Dementia' then 1 else 0 end) dementia ,
    sum(case when comorbidity = 'Pulmonary' then 1 else 0 end) pulmonary ,
    sum(case when comorbidity = 'Rheumatic' then 1 else 0 end) rheumatic ,
    sum(case when comorbidity = 'PUD' then 1 else 0 end) PUD ,
    sum(case when comorbidity = 'LiverMild' then 1 else 0 end) liver_mild ,
    sum(case when comorbidity = 'DM' then 1 else 0 end) diabetes ,
    sum(case when comorbidity = 'DMcx' then 1 else 0 end) dmcx ,
    sum(case when comorbidity = 'Paralysis' then 1 else 0 end) paralysis ,
    sum(case when comorbidity = 'Renal' then 1 else 0 end) renal ,
    sum(case when comorbidity = 'Cancer' then 1 else 0 end) cancer ,
    sum(case when comorbidity = 'LiverSevere' then 1 else 0 end) liversevere ,
    sum(case when comorbidity = 'Mets' then 1 else 0 end) mets ,   
    sum(case when comorbidity = 'hypertension' then 1 else 0 end) hypertension,
    sum(case when comorbidity = 'Upper GI Bleed' then 1 else 0 end) upper_gi_bleed,
  --  sum(case when comorbidity = 'HIV' then 1 else 0 end) hiv, 
    sum(case when comorbidity = 'hiv infection' then 1 else 0 end) hiv,    
    case when count(*) > 1 then 1 else 0 end multiple
FROM (
SELECT 
distinct
co.person_id ,
--replace(cs.concept_set_name, 'Charlson - ','') comorbidity 
case when cs.codeset_id = 382527336 then 'hiv infection' when cs.codeset_id = 779214702 then 'hypertension' when cs.codeset_id = 438154696 then 'Upper GI Bleed' else replace(cs.concept_set_name, 'Charlson - ','') end as comorbidity
FROM 
co left outer join  cs on ( cs.concept_id = co.condition_source_concept_id or cs.concept_id = co.condition_concept_id )
and cs.is_most_recent_version = true
and cs.codeset_id in ( 535274723, 359043664, 78746470, 719585646, 403438288, 
    494981955, 248333963, 378462283, 259495957, 489555336, 510748896, 514953976, 376881697, 
    220495690, 765004404, 652711186, 779214702, 382527336, 438154696 
    )
) t
group by t.person_id
) x
-- INNER JOIN COVID_Positive_CCI_INDEX c on c.person_id = x.person_id

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.7487019f-58ba-4147-a1b3-952fcca48c10"),
    complete_patient_table_with_derived_scores=Input(rid="ri.foundry.main.dataset.d467585c-53b3-4f10-b09e-8e86a4796a1a"),
    cov_med_365_neg=Input(rid="ri.foundry.main.dataset.cd6a4e5d-4cdc-4438-be6c-33aa604a600b"),
    preexisting_comorbidities_both=Input(rid="ri.foundry.main.dataset.69fda5b2-1ddc-4e20-9d18-ad1c44019fe6"),
    selected_critical_visits=Input(rid="ri.foundry.main.dataset.a8a5404a-760c-4f8e-bee2-76ccc7b5b781")
)
SELECT C.*, 
       M.Amiodarone_gtt, 
       M.Anakinra, 
       M.Azithromycin, 
       M.Chloroquine, 
       M.Dexamethasone, 
       M.dialysis_CRRT_HD, 
       M.Dobutamine_gtt, 
       M.Dopamine_gtt, 
       M.Epinephrine_gtt, 
       M.Epoprostenol, 
       M.Esmolol_gtt, 
       M.Hydrocortisone, 
       M.Hydroxychloroquine, 
       M.Inhaled_Nitric_Oxide, 
       M.intravenous_immunoglobulin, 
       M.Isoproterenol_gtt, 
       M.Levosimendan_gtt, 
       M.Lopinavir, 
       M.Lopinavir_Ritonavir_combination, 
       M.Methylprednisolone, 
       M.Milrinone_gtt, 
       M.Norepinephrine_gtt, 
       M.Phenylephrine_gtt, 
       M.Plasma, 
       M.Prednisolone, 
       M.Prednisone, 
       M.Remdesivir, 
       M.Ritonavir, 
       M.Tocilizumab, 
       M.Vasopressin_gtt, 
       M.ACE, 
       M.ARB, 
       M.Fluvoxamine, 
       M.Clarithromycin, 
       M.Erythromycin, 
       M.Josamycin, 
       M.Macrolides, 
       M.Midecamycin, 
       M.Miocamycin, 
       M.Prednisone_Methylprednisolone, 
       M.Rokitamycin, 
       M.Roxithromycin, 
       M.Spiramycin, 
       M.troleandomycin, 
       CP.visit_concept_id, 
       CP.visit_start_date, 
       CP.visit_concept_name, 
       CP.AKI_in_hospital, 
       CP.ECMO, 
       CP.Invasive_Ventilation, 
       CP.in_death_table, 
       CP.age_at_visit_start_in_years_int, 
       CP.length_of_stay, 
       CP.Race, 
       CP.Ethnicity, 
       CP.gender_concept_name, 
       CP.smoking_status, 
       CP.blood_type, 
       CP.Severity_Type, 
       CP.InpatientOrED, 
       CP.Q_Score, 
       CP.BMI, 
       S.readmission --, 
    --    P.RUCA1, 
    --    P.RUCA2, 
    --    P.RUCA_Zip_Type, 
    --    P.RUCC_2013 
FROM   preexisting_comorbidities_both C 
       LEFT JOIN cov_med_365_neg M 
              ON C.person_id = M.person_id 
       LEFT JOIN complete_patient_table_with_derived_scores CP 
              ON C.person_id = CP.person_id 
       LEFT JOIN selected_critical_visits S 
              ON C.person_id = S.person_id 
    --    LEFT JOIN Person_RUCA_RUCC P 
    --           ON C.person_id = P.person_id 

@transform_pandas(
    Output(rid="ri.vector.main.execute.9abfc8f5-82e0-442b-8163-9d9d0dc0eb04"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7"),
    Covid_positive_persons_LDS=Input(rid="ri.foundry.main.dataset.a917b9e7-b891-4ac5-93dc-de1a3b93c7ad"),
    compiled_gfr=Input(rid="ri.foundry.main.dataset.76d20b33-17ef-4566-9c80-2d19872d35cc"),
    summary_Table=Input(rid="ri.foundry.main.dataset.5dd50df8-79ce-461c-9397-0a429352f52a")
)
(SELECT count(*) as GFR_pos_earliest
FROM compiled_gfr as cg1
inner join  Covid_positive_persons_LDS pose on cg1.person_id = pose.person_id ) --   --- 

-- (SELECT count(*) as GFR_neg_summary
-- FROM compiled_gfr as cg2
-- inner join  summary_Table summ on cg2.person_id = summ.person_id )

-- (SELECT count(*) as GFR_neg_ACE
-- FROM compiled_gfr as cg3
-- inner join  ACE_ARB_Pull pull_pos on cg3.person_id = pull_pos.person_id) 

@transform_pandas(
    Output(rid="ri.vector.main.execute.4edca9ef-3c10-4c36-8cb0-4d262ff5eee3"),
    ACE_ARB_Pull_negative=Input(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4"),
    compiled_gfr_neg=Input(rid="ri.foundry.main.dataset.de0b8ba5-9bca-47f1-9e40-8e233ed3c6ad"),
    neg_earliest=Input(rid="ri.foundry.main.dataset.19034815-c3a7-43e6-85f3-95725a63b850"),
    summary_Table_neg=Input(rid="ri.foundry.main.dataset.7487019f-58ba-4147-a1b3-952fcca48c10")
)
(SELECT count(*) as GFR_neg_earliest
FROM compiled_gfr_neg as cg1
inner join  neg_earliest nege on cg1.person_id = nege.person_id) 

-- (SELECT count(*) as GFR_neg_summary
-- FROM compiled_gfr as cg2
-- inner join  summary_Table_neg summ on cg2.person_id = summ.person_id )

-- (SELECT count(*) as GFR_neg_ACE
-- FROM compiled_gfr as cg3
-- inner join  ACE_ARB_Pull_negative pull_neg on cg3.person_id = pull_neg.person_id )

@transform_pandas(
    Output(rid="ri.vector.main.execute.c85924a2-9e59-412a-a174-e64d00be6f80"),
    covid_19_negative_patients=Input(rid="ri.foundry.main.dataset.6dbbbc52-960c-4e0d-954b-d343eca1ecff")
)
SELECT count(distinct person_id)
FROM covid_19_negative_patients

