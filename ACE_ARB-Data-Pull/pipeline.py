

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.6dbbbc52-960c-4e0d-954b-d343eca1ecff"),
    measurement=Input(rid="ri.foundry.main.dataset.d6054221-ee0c-4858-97de-22292458fa19")
)
def covid_19_negative_patients(measurement):

    df_measurement = measurement

    covid_measurement_concept_ids = [

        '757680', '757679', '757678', '757677', '723459', '715262', '715261', '715260', '706181',
        '706180', '706179', '706178', '706177', '706176', '706175', '706174', '706173', '706172',
        '706171', '706170', '706169', '706168', '706167', '706166', '706165', '706163', '706161',
        '706160', '706159', '706158', '706157', '706156', '706155', '706154', '586526', '586523',
        '586522', '586521', '586520', '586519', '586518', '586517', '586516', '586515'

    ]

    covid_measurement_value_as_concept_ids = ['45878583', '45880296', '9189', '9190']

    persons_with_no_covid_measurement = df_measurement.where(

        (df_measurement.measurement_concept_id.isin(covid_measurement_concept_ids))

        & (df_measurement.value_as_concept_id.isin(covid_measurement_value_as_concept_ids))

    ).selectExpr("person_id", "data_partner_id", "measurement_date" , "measurement_concept_id", "measurement_concept_name as concept" , "measurement_source_value" , "value_as_concept_id" , "value_source_value" , "value_as_concept_name" ).distinct().withColumn("covid_diagnosis", F.lit(1))

    return persons_with_no_covid_measurement

######################
## Template imports ##
######################
from pyspark.sql import functions as F
import pandas as pd
import numpy as np
import math
import statistics

