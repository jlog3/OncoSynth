import json
import os 
import csv


# Extract patient resources
# patient_resources = [entry['resource'] for entry in fhir_json.get('entry', []) if entry.get('resource', {}).get('resourceType') == 'Patient']

# PATIENT
# Function to format patient data
def format_patient_data(patient):
    # Extracting and formatting basic patient information
    name_info = patient.get('name', [{}])[0]
    name = f"{name_info.get('family', '')}, {' '.join(name_info.get('given', []))}"
    gender = patient.get('gender', 'unknown')
    birth_date = patient.get('birthDate', 'unknown')

    # Address information
    address_info = patient.get('address', [{}])[0]
    address = address_info.get('line', [''])[0]
    city_state = f"{address_info.get('city', '')}, {address_info.get('state', '')}"
    postal_code = address_info.get('postalCode', 'unknown')

    # Extensions (Race, Ethnicity, etc.)
    race = 'unknown'
    ethnicity = 'unknown'
    for ext in patient.get('extension', []):
        if ext.get('url') == 'http://hl7.org/fhir/us/core/StructureDefinition/us-core-race':
            race = ext.get('extension', [{}])[0].get('valueCoding', {}).get('display', 'unknown')
        elif ext.get('url') == 'http://hl7.org/fhir/us/core/StructureDefinition/us-core-ethnicity':
            ethnicity = ext.get('extension', [{}])[0].get('valueCoding', {}).get('display', 'unknown')

    # Language
    language = patient.get('communication', [{}])[0].get('language', {}).get('text', 'unknown')

    # Blood type is not typically included in the standard FHIR patient resource, so we set it as 'unknown'
    blood_type = 'unknown'

    # Format output
    formatted_data = f"PATIENT\nName\n{name}\nGender\n{gender}\nDate of Birth\n{birth_date}\nAddress\n{address}\nCity, State\n{city_state}\nPostal Code\n{postal_code}\nRace\n{race}\nEthnicity\n{ethnicity}\nLanguage\n{language}\nBlood Type\n{blood_type}"
    
    # return formatted_data
    return name, gender, birth_date, address, city_state, postal_code, race, ethnicity, language, blood_type

import csv
import os

def write_patient_to_csv(name, gender, birth_date, address, city_state, postal_code, race, ethnicity, language, blood_type, dead_t, write_dir):
    # Construct the data dictionary
    patient_data = {
        'Name': name,
        'Gender': gender,
        'Date of Birth': birth_date,
        'Address': address,
        'City, State': city_state,
        'Postal Code': postal_code,
        'Race': race,
        'Ethnicity': ethnicity,
        'Language': language,
        'Blood Type': blood_type,
        # "Death": dead_t  # or some other placeholder
        "Death": str(dead_t) if dead_t is not None else "Alive"  # or some other placeholder
    }

    # File path
    file_path = os.path.join(write_dir, 'patient.csv')

    # Check if file exists to append or write headers
    file_exists = os.path.isfile(file_path)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    with open(file_path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=patient_data.keys())

        # If file does not exist, write the header
        if not file_exists:
            writer.writeheader()

        # Write the patient data
        writer.writerow(patient_data)

# Example usage
# write_patient_to_csv('John Doe', 'Male', '1990-01-01', '123 Main St', 
#                      'Anytown, State', '12345', 'Caucasian', 'Non-Hispanic', 
#                      'English', 'A+', 'path/to/your/directory')



####
# CONDITIONS
# Extract condition resources
# condition_resources = [entry['resource'] for entry in fhir_json.get('entry', []) if entry.get('resource', {}).get('resourceType') == 'Condition']

# Function to format condition data
def format_condition_data(condition):
    # Extracting and formatting condition information
    code_info = condition.get('code', {}).get('coding', [{}])[0]
    snomed_code = code_info.get('code', 'N/A')
    condition_display = code_info.get('display', 'N/A')
    onset_date = condition.get('onsetDateTime', 'N/A').split('T')[0] if 'onsetDateTime' in condition else 'N/A'
    # Assuming 'resolved' date is not provided in the sample, set as 'N/A'
    resolved_date = 'N/A'

    # Format output
    formatted_data = f"{snomed_code}\t{condition_display}\t{onset_date}\t{resolved_date}"
    
    # return formatted_data
    return snomed_code, condition_display, onset_date, resolved_date

import csv
import os

def write_condition_to_csv(name, snomed_code, condition_display, onset_date, resolved_date, write_dir):
    # Construct the data dictionary with 'name' as the first key
    condition_data = {
        'Name': name,
        'SNOMED Code': snomed_code,
        'Condition': condition_display,
        'Onset Date': onset_date,
        'Resolved Date': resolved_date
    }

    # File path
    file_path = os.path.join(write_dir, 'condition.csv')

    # Check if file exists to append or write headers
    file_exists = os.path.isfile(file_path)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    with open(file_path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=condition_data.keys())

        # If file does not exist, write the header
        if not file_exists:
            writer.writeheader()

        # Write the condition data
        writer.writerow(condition_data)

# Example usage
# write_condition_to_csv('John Doe', '123456', 'Sample Condition', '2023-01-01', 
                    #    '2023-06-01', '/home/john/Desktop/cancer_vis/tcga_synthea_test')


# MEDICATIONS
# Function to format medication request data
def format_medication_request_data(medication_request):
    # Extracting medication information from the MedicationRequest resource
    medication_info = medication_request.get('medicationCodeableConcept', {}).get('coding', [{}])[0]
    rxnorm_code = medication_info.get('code', 'N/A')
    medication_display = medication_info.get('display', 'N/A')
    date_prescribed = medication_request.get('authoredOn', 'N/A').split('T')[0]
    status = medication_request.get('status', 'unknown')

    # Format output
    formatted_data = f"{rxnorm_code}\t{medication_display}\t{date_prescribed}\t{status}"
    # return formatted_data
    return rxnorm_code, medication_display, date_prescribed, status

def write_medication_to_csv(name, rxnorm_code, medication_display, date_prescribed, status, write_dir):
    # Construct the data dictionary
    medication_data = {
        'Name': name,
        'RxNorm Code': rxnorm_code,
        'Medication': medication_display,
        'Date Prescribed': date_prescribed,
        'Status': status
    }

    # File path
    file_path = os.path.join(write_dir, 'medications.csv')

    # Check if file exists to append or write headers
    file_exists = os.path.isfile(file_path)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    with open(file_path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=medication_data.keys())

        # If file does not exist, write the header
        if not file_exists:
            writer.writeheader()

        # Write the medication data
        writer.writerow(medication_data)

# OBSERVATIONS
# Function to format observation data
def format_observation_data(observation):
    # Extracting observation information
    code_info = observation.get('code', {}).get('coding', [{}])[0]
    loinc_code = code_info.get('code', 'N/A')
    observation_display = code_info.get('display', 'N/A')
    
    # Extracting value and unit
    value_quantity = observation.get('valueQuantity', {})
    value = value_quantity.get('value', 'N/A')
    unit = value_quantity.get('unit', '')
    observation_value = f"{value} {unit}".strip()

    # Extracting the date recorded
    date_recorded = observation.get('effectiveDateTime', 'N/A').split('T')[0]

    # Format output
    formatted_data = f"{loinc_code}\t{observation_display}\t{observation_value}\t{date_recorded}"
    # return formatted_data
    return loinc_code, observation_display, observation_value, date_recorded

def write_observation_to_csv(name, loinc_code, observation_display, observation_value, date_recorded, write_dir):
    # Construct the data dictionary with 'name' as the first key
    observation_data = {
        'Name': name,
        'LOINC Code': loinc_code,
        'Observation': observation_display,
        'Value': observation_value,
        'Date Recorded': date_recorded
    }

    # File path
    file_path = os.path.join(write_dir, 'observations.csv')

    # Check if file exists to append or write headers
    file_exists = os.path.isfile(file_path)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    with open(file_path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=observation_data.keys())

        # If file does not exist, write the header
        if not file_exists:
            writer.writeheader()

        # Write the observation data
        writer.writerow(observation_data)


# REPORTS
import base64

def format_report_data(report):
    # Extracting the LOINC code
    code_info = report.get('code', {}).get('coding', [{}])[0]
    loinc_code = code_info.get('code', '')

    # Extracting the report text (assuming it's base64 encoded in the 'data' field of 'presentedForm')
    presented_form = report.get('presentedForm', [{}])[0]
    data_encoded = presented_form.get('data', '')
    report_text = ''
    if data_encoded:
        report_text = base64.b64decode(data_encoded).decode('utf-8')

    # Extracting the effective date/time
    effective_date = report.get('effectiveDateTime', '')
    if effective_date:
        effective_date = effective_date.replace('T', ' ').replace('-04:00', '')

    return loinc_code, report_text, effective_date


def write_report_to_csv(name, loinc_code, report_text, effective_date, write_dir):
    # Construct the data dictionary with 'name' as the first key
    report_data = {
        'Name': name,
        'LOINC': loinc_code,
        'Report/Observation': report_text,
        'Date': effective_date
    }

    # File path
    file_path = os.path.join(write_dir, 'reports.csv')

    # Check if file exists to append or write headers
    file_exists = os.path.isfile(file_path)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    with open(file_path, 'a', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=report_data.keys())

        # If file does not exist, write the header
        if not file_exists:
            writer.writeheader()

        # Write the report data
        writer.writerow(report_data)


# CAREPLANS
def format_careplan_data(careplan):
    # Extracting the category code (SNOMED if available)
    category = careplan.get('category', [{}])
    snomed_code = ''
    for cat in category:
        coding = cat.get('coding', [{}])
        for code in coding:
            if code.get('system', '') == 'http://snomed.info/sct':
                snomed_code = code.get('code', '')
                break

    # Extracting care plan activities
    # Extracting care plan activities
    activities = careplan.get('activity', [])
    activity_details = []
    for act in activities:
        detail = act.get('detail', {})
        code_info = detail.get('code', {}).get('coding', [{}])[0]
        activity_display = code_info.get('display', '')
        if activity_display:
            activity_details.append(activity_display)

    activities_display = "Activity: " + "\nActivity: ".join(activity_details) if activity_details else "No activities"

    print(f'activities_display;{activities_display}')

    # Extracting the start date
    period = careplan.get('period', {})
    start_date = period.get('start', '')
    if start_date:
        # Formatting the date
        start_date = start_date.replace('T', ' ').replace('-04:00', '')

    return snomed_code, activities_display, start_date

# # Example usage
# careplan_json = {
#     "resource": {
#         # ... (your JSON data)
#     }
# }

# snomed_code, activities_display, start_date = format_careplan_data(careplan_json['resource'])
# formatted_data = f"{snomed_code}\t{activities_display}\t{start_date}"


def write_careplan_to_csv(name, snomed_code, activities_display, start_date, write_dir):
    # Construct the data dictionary with 'name' as the first key
    careplan_data = {
        'Name': name,
        'SNOMED Code': snomed_code,
        'Activities': activities_display,
        'Date': start_date
    }

    # File path
    file_path = os.path.join(write_dir, 'careplans.csv')

    # Check if file exists to append or write headers
    file_exists = os.path.isfile(file_path)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    with open(file_path, 'a', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=careplan_data.keys())

        # If file does not exist, write the header
        if not file_exists:
            writer.writeheader()

        # Write the careplan data
        writer.writerow(careplan_data)


# PROCEDURES
def format_procedure_data(procedure):
    # Extracting the SNOMED code and the display name
    code_info = procedure.get('code', {}).get('coding', [{}])[0]
    snomed_code = code_info.get('code', '')
    procedure_display = code_info.get('display', '')

    # Extracting the date performed
    performed_period = procedure.get('performedPeriod', {})
    start_date = performed_period.get('start', '')
    
    # Formatting the date if it exists
    if start_date:
        # Assuming date is in ISO format, converting to desired format
        start_date = start_date.replace('T', ' ').replace('-04:00', ' pm')
    
    return snomed_code, procedure_display, start_date


def write_procedure_to_csv(name, snomed_code, procedure_display, start_date, write_dir):
    # Construct the data dictionary with 'name' as the first key
    procedure_data = {
        'Name': name,
        'SNOMED Code': snomed_code,
        'Procedure': procedure_display,
        'Date Performed': start_date
    }

    # File path
    file_path = os.path.join(write_dir, 'procedures.csv')

    # Check if file exists to append or write headers
    file_exists = os.path.isfile(file_path)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    with open(file_path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=procedure_data.keys())

        # If file does not exist, write the header
        if not file_exists:
            writer.writeheader()

        # Write the procedure data
        writer.writerow(procedure_data)

# Example usage
write_procedure_to_csv('John Doe', '265764009', 'Renal dialysis (procedure)', '2023-12-04 1:52:33 pm', '/home/john/Desktop/cancer_vis/tcga_synthea_test')




# ENCOUNTERS
from datetime import datetime

def format_encounter_data(encounter):
    # Extracting the SNOMED code and encounter description
    type_info = encounter.get('type', [{}])[0].get('coding', [{}])[0]
    snomed_code = type_info.get('code', '')
    encounter_display = type_info.get('display', '')

    # Extracting start and end times
    period = encounter.get('period', {})
    start_time = period.get('start', '')
    end_time = period.get('end', '')

    # Calculating duration
    duration = ''
    if start_time and end_time:
        start_dt = datetime.fromisoformat(start_time)
        end_dt = datetime.fromisoformat(end_time)
        duration_seconds = (end_dt - start_dt).total_seconds()
        duration_minutes = duration_seconds / 60
        duration = f"{int(duration_minutes)} minutes"

    # Formatting the start time
    if start_time:
        start_time = start_time.replace('T', ' ').replace('-05:00', '')

    return snomed_code, encounter_display, start_time, duration

# snomed_code, encounter_display, start_time, duration = format_encounter_data(encounter_json['resource'])
# formatted_data = f"{snomed_code}\t{encounter_display}\t{start_time}\t{duration}"

def write_encounter_to_csv(name, snomed_code, encounter_display, start_time, duration, write_dir):
    # Construct the data dictionary with 'name' as the first key
    encounter_data = {
        'Name': name,
        'SNOMED': snomed_code,
        'Encounter': encounter_display,
        'Start Time': start_time,
        'Duration': duration
    }

    # File path
    file_path = os.path.join(write_dir, 'encounters.csv')

    # Check if file exists to append or write headers
    file_exists = os.path.isfile(file_path)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    with open(file_path, 'a', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=encounter_data.keys())

        # If file does not exist, write the header
        if not file_exists:
            writer.writeheader()

        # Write the encounter data
        writer.writerow(encounter_data)



# ALLERGY
def format_allergy_data(allergy):
    # Extracting the SNOMED code and allergy substance
    allergy_info = allergy.get('code', {}).get('coding', [{}])[0]
    snomed_code = allergy_info.get('code', '')
    allergy_substance = allergy_info.get('display', '')

    return snomed_code, allergy_substance

# Example usage with your provided JSON structure
# allergy_json = {
#     "resource": {
#         # ... (your JSON data)
#     }
# }

# snomed_code, allergy_substance = format_allergy_data(allergy_json['resource'])
# formatted_data = f"{snomed_code}\t{allergy_substance}"


def write_allergy_to_csv(name, snomed_code, allergy_substance, write_dir):
    # Construct the data dictionary with 'name' as the first key
    allergy_data = {
        'Name': name,
        'SNOMED': snomed_code,
        'Allergy': allergy_substance
    }

    # File path
    file_path = os.path.join(write_dir, 'allergies.csv')

    # Check if file exists to append or write headers
    file_exists = os.path.isfile(file_path)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    with open(file_path, 'a', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=allergy_data.keys())

        # If file does not exist, write the header
        if not file_exists:
            writer.writeheader()

        # Write the allergy data
        writer.writerow(allergy_data)

# Example usage
# write_allergy_to_csv('John Doe', snomed_code, allergy_substance, '/path/to/directory')


# VACCINATIONS
def format_vaccination_data(immunization):
    # Extracting the CVX code and vaccine name
    vaccine_info = immunization.get('vaccineCode', {}).get('coding', [{}])[0]
    cvx_code = vaccine_info.get('code', '')
    vaccine_name = vaccine_info.get('display', '')

    # Extracting the date given
    occurrence_date = immunization.get('occurrenceDateTime', '')
    if occurrence_date:
        occurrence_date = occurrence_date.replace('T', ' ').replace('-04:00', '')

    return cvx_code, vaccine_name, occurrence_date

# Example usage with your provided JSON structure
# immunization_json = {
#     "resource": {
#         # ... (your JSON data)
#     }
# }

# cvx_code, vaccine_name, occurrence_date = format_vaccination_data(immunization_json['resource'])
# formatted_data = f"{cvx_code}\t{vaccine_name}\t{occurrence_date}"
def write_vaccination_to_csv(name, cvx_code, vaccine_name, occurrence_date, write_dir):
    # Construct the data dictionary with 'name' as the first key
    vaccination_data = {
        'Name': name,
        'CVX': cvx_code,
        'Vaccine': vaccine_name,
        'Date Given': occurrence_date
    }

    # File path
    file_path = os.path.join(write_dir, 'vaccinations.csv')

    # Check if file exists to append or write headers
    file_exists = os.path.isfile(file_path)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    with open(file_path, 'a', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=vaccination_data.keys())

        # If file does not exist, write the header
        if not file_exists:
            writer.writeheader()

        # Write the vaccination data
        writer.writerow(vaccination_data)

# Example usage
# write_vaccination_to_csv('John Doe', cvx_code, vaccine_name, occurrence_date, '/path/to/directory')




# DOCUMENTS
def format_document_data(document):
    # Extracting the date
    date = document.get('date', '')
    if date:
        date = date.replace('T', ' ').replace('-04:00', '')

    # Extracting the document content (assuming it's base64 encoded)
    content_info = document.get('content', [{}])[0].get('attachment', {})
    data_encoded = content_info.get('data', '')
    document_content = ''
    if data_encoded:
        document_content = base64.b64decode(data_encoded).decode('utf-8')

    return date, document_content

def write_document_to_csv(name, date, document_content, write_dir):
    # Construct the data dictionary with 'name' as the first key
    document_data = {
        'Name': name,
        'Date': date,
        'Content': document_content
    }

    # File path
    file_path = os.path.join(write_dir, 'documents.csv')

    # Check if file exists to append or write headers
    file_exists = os.path.isfile(file_path)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    with open(file_path, 'a', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=document_data.keys())

        # If file does not exist, write the header
        if not file_exists:
            writer.writeheader()

        # Write the document data
        writer.writerow(document_data)


# For each category of data in the Synthea FHIR output (conditions, observations, etc.), create separate CSV files.

# Requirements
# Have already generated patient FHIR data using Synthea

# Use
# change fhir_dir to the location of your FHIR files
def main(): 
    # make keep file to get BCRA

    # synth_dir = '~/Desktop/cancer_vis/synthea/synthea/'
    fhir_dir = '/home/john/Desktop/cancer_vis/synthea/synthea/output/fhir'
    write_dir = '/home/john/Desktop/cancer_vis/tcga_synthea_test'

    # pt_filename = 'Bernardina456_Klocko335_db52edeb-f696-c792-0c06-10a9eb74b5bb.json'  #bedubg
    # pt_filename = 'Hulda44_Eusebia552_Runolfsdottir785_0d3cdf04-e9d4-2a77-ccd4-ae01b6f3d0ad.json'  #bedubg
    # pt_filename = 'Trinh357_Kilback373_c3abc213-3b98-97cc-2479-38f3069827f1.json'  #bedubg
    # if True: 
    # Iterate over files in the directory
    for pt_filename in os.listdir(fhir_dir):
        # Check if the file does not start with the specified strings
        if not pt_filename.startswith("practitionerInformation") and not pt_filename.startswith("hospitalInformation"):
            # Process the file here
            print(f"Processing {pt_filename}")

        pt_fhir = os.path.join(fhir_dir, pt_filename)

        # Load the FHIR JSON file
        with open(pt_fhir, 'r') as file:
            fhir_data = json.load(file)


        # Initialize a dictionary to hold grouped resources
        grouped_resources = {}

        # Iterate through each entry in the 'entry' array
        for entry in fhir_data.get('entry', []):
            # Extract the resource
            resource = entry.get('resource', {})

            # Get the resource type
            resource_type = resource.get('resourceType', 'Unknown')

            # Add the resource to the corresponding list in the dictionary
            if resource_type not in grouped_resources:
                grouped_resources[resource_type] = []
            grouped_resources[resource_type].append(resource)

        # Now you have a dictionary where each key is a resource type
        # and each value is a list of resources of that type

        # put Encounter first to check for death certificate and get T/F and Datetime and put that in the write_pt csv param to help us with matching to TCGA
        # dead = 'alive'
        dead_t = None
        # dead_t = 'alive'
        if 'Encounter' in grouped_resources:
            print(f'ENCOTUNER')
            for encounter in grouped_resources['Encounter']:
                snomed_code, encounter_display, start_time, duration = format_encounter_data(encounter)
                if encounter_display == "Death Certification": 
                    # dead = True
                    dead_t = start_time
                    cdt = type(dead_t)
                    print(f'\nthey died {dead_t}:{cdt}\n')
                    break
        print('<>' * 50)  


        # Example: Print all Patient resources
        if 'Patient' in grouped_resources:
            for patient in grouped_resources['Patient']: # ONLY ever one patient since iterating over FHIR files which are for each patient, important because using the returned "name" later to identify patients, replaced later by submitter_id
                name, gender, birth_date, address, city_state, postal_code, race, ethnicity, language, blood_type = format_patient_data(patient)
                print(f"PATIENT\nName\n{name}\nGender\n{gender}\nDate of Birth\n{birth_date}\nAddress\n{address}\nCity, State\n{city_state}\nPostal Code\n{postal_code}\nRace\n{race}\nEthnicity\n{ethnicity}\nLanguage\n{language}\nBlood Type\n{blood_type}\nDead Time\n{dead_t}")
                write_patient_to_csv(name, gender, birth_date, address, city_state, postal_code, race, ethnicity, language, blood_type, dead_t, write_dir)
                print('<>' * 50)  # Separator for readability


        if 'Encounter' in grouped_resources:
            print(f'ENCOTUNER')
            for encounter in grouped_resources['Encounter']:
                snomed_code, encounter_display, start_time, duration = format_encounter_data(encounter)

                print(f"{snomed_code}\t{encounter_display}\t{start_time}\t{duration}")
                # Now you can write formatted_data to a CSV
                write_encounter_to_csv(name, snomed_code, encounter_display, start_time, duration, write_dir)
        print('<>' * 50)  



        if 'Condition' in grouped_resources:
            # Print header for conditions
            print("CONDITIONS\nSNOMED\tCondition\tDate of Onset\tDate Resolved")

            # Process each condition resource
            for condition in grouped_resources['Condition']:
                snomed_code, condition_display, onset_date, resolved_date = format_condition_data(condition)
                # print(f"{snomed_code}\t{condition_display}\t{onset_date}\t{resolved_date}")
                # print('-' * 50)  # Separator for readability
                write_condition_to_csv(name, snomed_code, condition_display, onset_date, resolved_date, write_dir)
        print('<>' * 50)  # Separator for readability


        if 'MedicationRequest' in grouped_resources:
            # Print header for medications
            print("MEDICATIONS\nRxNorm\tMedication\tDate Prescribed\tStatus")

            # Process each medication request resource
            for medication_request in grouped_resources['MedicationRequest']:
                # print(format_medication_request_data(medication_request))
                rxnorm_code, medication_display, date_prescribed, status = format_medication_request_data(medication_request)
                # print(f"{rxnorm_code}\t{medication_display}\t{date_prescribed}\t{status}")
                # print('-' * 50)  # Separator for readability
                write_medication_to_csv(name, rxnorm_code, medication_display, date_prescribed, status, write_dir)
        print('<>' * 50)  # Separator for readability


        if 'Observation' in grouped_resources:
            # Print header for observations
            print("OBSERVATIONS\nLOINC\tObservation\tValue\tDate Recorded")

            # Process each observation resource
            for observation in grouped_resources['Observation']:
                # print(format_observation_data(observation))
                loinc_code, observation_display, observation_value, date_recorded = format_observation_data(observation)
                # print(f"{loinc_code}\t{observation_display}\t{observation_value}\t{date_recorded}")
                # print('-' * 50)  # Separator for readability
                write_observation_to_csv(name, loinc_code, observation_display, observation_value, date_recorded, write_dir)
        print('<>' * 50)  # Separator for readability


        if 'CarePlan' in grouped_resources:
            for careplan in grouped_resources['CarePlan']:
                snomed_code, activities_display, start_date = format_careplan_data(careplan)
                # print(f"{snomed_code}\t{activities_display}\t{start_date}")
                # Now you can write formatted_data to a CSV or perform other actions
                write_careplan_to_csv(name, snomed_code, activities_display, start_date, write_dir)
        print('<>' * 50)  


        if 'Procedure' in grouped_resources:
            for procedure in grouped_resources['Procedure']:
                snomed_code, procedure_display, start_date = format_procedure_data(procedure)
                # print(f"{snomed_code}\t{procedure_display}\t{start_date}")
                # Now you can write formatted_data to a CSV or perform other actions
                write_procedure_to_csv(name, snomed_code, procedure_display, start_date, write_dir)
        print('<>' * 50)  


        if 'DiagnosticReport' in grouped_resources:
            for report in grouped_resources['DiagnosticReport']:
                loinc_code, report_text, effective_date = format_report_data(report)
                # print(f"{loinc_code}\t{report_text}\t{effective_date}")
                # Now you can write formatted_data to a CSV
                write_report_to_csv(name, loinc_code, report_text, effective_date, write_dir)
        print('<>' * 50)  


    
        if 'DocumentReference' in grouped_resources:
            for document in grouped_resources['DocumentReference']:
                date, document_content = format_document_data(document)
                # print(f"{date}\t{document_content}")
                # Now you can write formatted_data to a CSV
                write_document_to_csv(name, date, document_content, write_dir)
        print('<>' * 50)  


        if 'AllergyIntolerance' in grouped_resources:
            # print("ALLREGY\nLOINC\tObservation\tValue\tDate Recorded")

            for allergy in grouped_resources['AllergyIntolerance']:
                snomed_code, allergy_substance = format_allergy_data(allergy)
                # print(f"{snomed_code}\t{allergy_substance}")
                # Now you can write formatted_data to a CSV
                write_allergy_to_csv(name, snomed_code, allergy_substance, write_dir)
        print('<>' * 50)  


        if 'Immunization' in grouped_resources:
            for immunization in grouped_resources['Immunization']:
                cvx_code, vaccine_name, occurrence_date = format_vaccination_data(immunization)
                # print(f"{cvx_code}\t{vaccine_name}\t{occurrence_date}")
                # Now you can write formatted_data to a CSV
                write_vaccination_to_csv(name, cvx_code, vaccine_name, occurrence_date, write_dir)
        print('<>' * 50)  



if __name__ == "__main__":
    main()
