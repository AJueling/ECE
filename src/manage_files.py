# this file collects functions to handle file, e.g., moving or renaming them
import os
import shutil

def move_files(source_dir, destination_dir, file_extensions=None):
    """
    Move files from source directory to destination directory.

    Parameters:
        source_dir (str): Path to the source directory.
        destination_dir (str): Path to the destination directory.
        file_extensions (list): List of file extensions to move. If None, move all files.

    Returns:
        None
    """
    # Ensure both directories exist
    if not os.path.exists(destination_dir):
        print(f"Destination directory {destination_dir} does not exist.")
        return
    if not os.path.exists(source_dir):
        print("Source directory does not exist.")
        return

    # Create destination directory if it doesn't exist
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    # Get list of files in the source directory
    files = os.listdir(source_dir)

    # Filter files by extensions if provided
    if file_extensions:
        files = [file for file in files if os.path.isfile(os.path.join(source_dir, file)) and
                 file.endswith(tuple(file_extensions))]

    # Move each file to the destination directory
    for file in files:
        source_file_path = os.path.join(source_dir, file)
        destination_file_path = os.path.join(destination_dir, file)
        try:
            shutil.move(source_file_path, destination_file_path)
            print(f"Moved {source_file_path} to {destination_file_path}")
        except Exception as e:
            print(f"Failed to move {source_file_path}: {str(e)}")


def rename_files(directory, old_part, new_part):
    """
    Rename files in the directory by replacing old_part with new_part.

    Parameters:
        directory (str): Path to the directory containing the files.
        old_part (str): The part of the filename to be replaced.
        new_part (str): The string to replace old_part with.

    Returns:
        None
    """
    # Ensure the directory exists
    if not os.path.exists(directory):
        print("Directory does not exist.")
        return

    # Iterate over files in the directory
    for filename in os.listdir(directory):
        # Construct the new filename
        new_filename = filename.replace(old_part, new_part)

        # Check if the filename actually changed
        if new_filename != filename:
            try:
                # Rename the file
                os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))
                print(f"Renamed {filename} to {new_filename}")
            except Exception as e:
                print(f"Failed to rename {filename}: {str(e)}")


def create_and_write_file(file_path, lines):
    """
    Create a file and write lines to it.

    Parameters:
        file_path (str): Path to the file to be created.
        lines (list): List of strings representing lines to be written to the file.

    Returns:
        None
    """
    try:
        # Open the file in write mode
        with open(file_path, 'w') as file:
            # Write each line to the file
            for line in lines:
                file.write(line + '\n')
        print(f"File '{file_path}' created and lines added successfully.")
    except Exception as e:
        print(f"Failed to create file or add lines: {str(e)}")


def remove_empty_folders(path_abs):
    walk = list(os.walk(path_abs))
    for path, _, _ in walk[::-1]:
        if len(os.listdir(path)) == 0:
            os.rmdir(path)