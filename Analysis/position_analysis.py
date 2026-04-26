import csv
import pandas as pd
import os

 
def main():
    xvg_file = 'RUN_position.xvg'  # Replace with your actual .xvg file
    start_atom_index = 0  # Replace with the start of the desired atom index range
    end_atom_index = 14753 # Replace with the end of the desired atom index range 125991
    timeframes = [i/1.0 for i in range(0, 200001, 10)]
    calculate_averages_position(xvg_file, timeframes, start_atom_index, end_atom_index)

def extract_coordinates_from_xvg(xvg_file, time, start_atom_index, end_atom_index):
    """
    Extracts X, Y, and Z coordinates for all atoms from a specified start index to an end index at a specific time from an XVG file.

    Args:
        xvg_file (str): Path to the XVG file.
        time (float): The specific time to extract the coordinates.
        start_atom_index (int): The start index of the atoms.
        end_atom_index (int): The end index of the atoms.

    Returns:
        List of tuples (time, atom, X, Y, Z) representing the coordinates for all atoms at the specified time.
    """
    coordinates = []
    with open(xvg_file, 'r') as xvg:
        for line in xvg:
            if line.startswith(("#", "@")):
                continue  # Skip comments and metadata lines
            tokens = line.split()
            if len(tokens) >= 4 and float(tokens[0]) == time:
                for atom_index in range(start_atom_index, end_atom_index + 1):
                    x, y, z = float(tokens[atom_index * 3 + 1]), float(tokens[atom_index * 3 + 2]), float(tokens[atom_index * 3 + 3])
                    coordinates.append((time, atom_index, x, y, z))
                break
    DF = pd.DataFrame(coordinates, columns=["time", "index", "Px", "Py", "Pz"])        
    return DF

def calculate_averages_position(xvg_file, timeframes, start_atom_index, end_atom_index):
    # Initialize an empty dataframe to store the results
    results = pd.DataFrame(columns=['time', 'avg_Px', 'avg_Py', 'avg_Pz'])

    for timeframe in timeframes:
        print(timeframe)
        # Construct the file name
        df = extract_coordinates_from_xvg(xvg_file, timeframe, start_atom_index, end_atom_index)

        # Calculate the averages
        avg_Px = df['Px'].mean()
        avg_Py = df['Py'].mean()
        avg_Pz = df['Pz'].mean()

        # Append the results to the dataframe
        new_row = pd.DataFrame({'time': [timeframe], 'avg_Px': [avg_Px], 'avg_Py': [avg_Py], 'avg_Pz': [avg_Pz]})
        results = pd.concat([results, new_row], ignore_index=True)

    # Save the results to a csv file
    results.to_csv(f"Position_averages.csv", index=False)

    print(f"Averages saved to Position_averages.csv")

if __name__ == "__main__":
    main()


