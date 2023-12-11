import pandas as pd
import os



"""

Data Matching and Integration:

BASH

# create TCGA data in R 

cd ~/Desktop/cancer_vis
./run_synthea -p 10 -g F -a  -k ../../keep.json
python3 fhir_wrangle.py  # create csvs using fhir
python3 tcga_synthea_matching.py   # use the tcga files from the same dir to get submitter_ids and assign to a synthea patient

(manually added the normal rawcounts file to the zip )

zip -r tcga_synthea_test.zip tcga_synthea_test/
# this file i suploaded to shiny

"""

def get_matching_synthea_patient(tcga_row, synthea_df, used_patient_names):
    is_alive_tcga = tcga_row['vital_status'] == 'Alive'
    gender_tcga = tcga_row['gender'].lower()

    for _, patient_row in synthea_df.iterrows():
        # Check if patient is alive in Synthea data
        is_alive_synthea = patient_row['Death'] == 'Alive'
        gender_synthea = patient_row['Gender'].lower()

        if (is_alive_synthea == is_alive_tcga and 
            gender_synthea == gender_tcga and 
            patient_row['Name'] not in used_patient_names):
            return patient_row['Name']
    return None


def main(): 
    filedir = '/home/john/Desktop/cancer_vis/tcga_synthea_test'

    # Load the TCGA clinical data CSV file
    # tcga_clinical_path = os.path.join(filedir, 'tcga_clinical.csv')  # Replace with your actual file path
    tcga_clinical_df = pd.read_csv(os.path.join(filedir, 'tcga_clinical.csv'), index_col=0)  # Assuming the first column is the row name
    synthea_patients_df = pd.read_csv(os.path.join(filedir, 'patient.csv'))
    print(synthea_patients_df.columns)

    # Iterate through each row (sample)
    # for index, row in tcga_clinical_df.iterrows():
    #     # Extract needed information
    #     sample_id = index  # Row name as sample ID
    used_patient_names = set()  # Set to keep track of used Synthea patient names
    updated_data = {}  # Dictionary to store updated data for each file

    # Process each TCGA sample
    for index, tcga_row in tcga_clinical_df.iterrows():
        sample_id = index
        matching_patient_name = get_matching_synthea_patient(tcga_row, synthea_patients_df, used_patient_names)
        if matching_patient_name:
            used_patient_names.add(matching_patient_name)  # Add the used patient name to the set

            print('got amtching')
            print(matching_patient_name)
            # Process each Synthea file
            for filename in os.listdir(filedir):
                if filename.startswith("tcga") or filename.startswith(".") or filename.startswith("new_"):
                    continue  # Skip TCGA files and hidden files

                # Read the Synthea file
                file_path = os.path.join(filedir, filename)
                synthea_df = pd.read_csv(file_path)
                print(f'csv:{file_path}')

                                # Filter the Synthea data to only include rows for the matching patient
                filtered_df = synthea_df[synthea_df['Name'] == matching_patient_name].copy()
                filtered_df['TCGA_Sample_ID'] = sample_id  # Add the TCGA sample_id

                # Accumulate the updated data
                if filename not in updated_data:
                    updated_data[filename] = filtered_df
                else:
                    updated_data[filename] = pd.concat([updated_data[filename], filtered_df])

        else: 
            print('got no matching')

    # Save the accumulated updated data to new files
    for filename, data_df in updated_data.items():
        new_file_path = os.path.join(filedir, f'new_{filename}')
        data_df.to_csv(new_file_path, index=False)

    # Optionally, replace old Synthea files with new ones
    for filename in os.listdir(filedir):
        if filename.startswith('new_'):
            old_file_path = os.path.join(filedir, filename[4:])
            new_file_path = os.path.join(filedir, filename)
            os.remove(old_file_path)
            os.rename(new_file_path, old_file_path)



if __name__ == "__main__":
    main()