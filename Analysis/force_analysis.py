
import pandas as pd

def main():
    xvg_file = 'RUN_force.xvg'  # Replace with your actual .xvg file
    start_atom_index = 0  # Replace with the start of the desired atom index range
    end_atom_index = 14753 # Replace with the end of the desired atom index range 125991
    timeframes = [i/1.0 for i in range(0, 200001, 10)]
    calculate_averages(xvg_file, timeframes, start_atom_index, end_atom_index)


def extract_forces_at_time(xvg_file, start_atom_index, end_atom_index, time_input):
    """
    Extracts forces (X, Y, Z components) for a range of atoms at time=0 from a GROMACS .xvg file and writes them to a CSV file.

    Args:
        xvg_file (str): Path to the GROMACS .xvg file.
        start_atom_index (int): Start index of the atom range (e.g., 1 for Atom 1).
        end_atom_index (int): End index of the atom range (inclusive).
        output_file (str): Path to the output CSV file.

    Returns:
        None
    """
    forces_at_time_zero = []
    with open(xvg_file, 'r') as file:  
        for line in file:
            if not line.startswith(('#', '@')):  # Skip comment and parameter lines
                tokens = line.split()
                time = float(tokens[0])
                if time == time_input:  # Check if time is 0
                    if len(tokens) >= (end_atom_index * 3 + 4):  # Check if tokens has enough elements
                        for atom_index in range(start_atom_index, end_atom_index + 1):
                            fx, fy, fz = float(tokens[atom_index * 3 + 1]), float(tokens[atom_index * 3 + 2]), float(tokens[atom_index * 3 + 3])
                            forces_at_time_zero.append([time, atom_index, fx, fy, fz])
                        break  # Exit the loop once forces at time=0 are found

    DF = pd.DataFrame(forces_at_time_zero, columns=["time", "index", "Fx", "Fy", "Fz"])
    #df.to_csv(output_file, index=False)
    return DF

def calculate_averages(xvg_file, timeframes, start_atom_index, end_atom_index):
    # Initialize an empty dataframe to store the results
    results = pd.DataFrame(columns=['time', 'avg_Fx', 'avg_Fy', 'avg_Fz'])

    for time in timeframes:
        print(time)
        # Check if the file exists
        # Read the csv file
        df = extract_forces_at_time(xvg_file, start_atom_index, end_atom_index, time)

        # Calculate the averages
        avg_Fx = df['Fx'].mean()
        avg_Fy = df['Fy'].mean()
        avg_Fz = df['Fz'].mean()

        # Append the results to the dataframe
        new_row = pd.DataFrame({'time': [time], 'avg_Fx': [avg_Fx], 'avg_Fy': [avg_Fy], 'avg_Fz': [avg_Fz]})
        results = pd.concat([results, new_row], ignore_index=True)

    # Save the results to a csv file
    results.to_csv(f"Force_averages.csv", index=False)

    print(f"Averages saved to Force_averages.csv")

if __name__ == "__main__":
    main()

