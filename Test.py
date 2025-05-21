import pandas as pd

def correct_cds_codes_csv(input_file, column_name, output_file):
    # Read the CSV file, ensuring the column is read as a string
    df = pd.read_csv(input_file, dtype={column_name: str})

    # Ensure all CDS codes are strings and add leading zero if length is 13
    df[column_name] = df[column_name].apply(lambda x: x.zfill(14) if len(x) == 13 else x)

    # Save the corrected data to a new CSV file
    df.to_csv(output_file, index=False)
    print(f"Corrected file saved as {output_file}")

input_file = r"C:\Users\huang\Documents\Deepsolar (Broken)\Results\first half\school_results.csv"  # Replace with your input file path
column = 'cds_code'  # Replace with the actual column name in your file
output_file = r'C:\Users\huang\Documents\Deepsolar (Broken)\Results\first half\corrected_cds_codes.csv'


correct_cds_codes_csv(input_file, column, output_file)
