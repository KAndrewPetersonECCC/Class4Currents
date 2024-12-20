import os
def execute_python_file(file_path):
    if ( os.path.isfile(file_path) ):
        with open(file_path, 'r') as file:
            python_code = file.read()
            exec(python_code)
    else:   
        print(f"Error: The file '{file_path}' does not exist.")
