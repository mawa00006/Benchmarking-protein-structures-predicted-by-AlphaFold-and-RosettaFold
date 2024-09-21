import os
import zipfile

def clean_name(name):
    """
    Helper function to remove '_longest_chain' from the folder or file name.
    """
    name = name.replace("_1", "")
    return name.replace('_longest_chain', '')

def extract_zip_files(zip_folder, destination_folder):
    # Get a list of all zip files in the source folder
    for file_name in os.listdir(zip_folder):
        if file_name.endswith(".zip"):
            # Path for the zip file
            zip_path = os.path.join(zip_folder, file_name)

            # Clean the folder name by removing '_longest_chain' and 'fold_ if present
            folder_name = clean_name(os.path.splitext(file_name)[0]).replace("fold_", "").upper()
            new_folder_path = os.path.join(destination_folder, folder_name)

            # Create the new directory if it doesn't exist
            os.makedirs(new_folder_path, exist_ok=True)

            # Extract the zip file into a temporary folder
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(new_folder_path)

            # Now rename the files inside the extracted folder, removing '_longest_chain'
            for root, dirs, files in os.walk(new_folder_path):
                for original_file_name in files:
                    # Clean the file name by removing '_longest_chain'
                    cleaned_file_name = clean_name(original_file_name).upper()

                    # If the file name changes, rename the file
                    if cleaned_file_name != original_file_name:
                        original_file_path = os.path.join(root, original_file_name)
                        cleaned_file_path = os.path.join(root, cleaned_file_name)
                        os.rename(original_file_path, cleaned_file_path)

            print(f"Extracted {file_name} to {new_folder_path} (cleaned folder and file names)")

# Example usage:
zip_folder = "alphafold"           # Folder containing zip files
destination_folder = "alphafold_extracted"  # Folder where you want to extract the files

extract_zip_files(zip_folder, destination_folder)
