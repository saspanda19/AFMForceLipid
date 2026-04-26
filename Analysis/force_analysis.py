import pandas as pd

def main():
    xvg_file = 'RUN_force.xvg'  # Replace with your actual .xvg file
    start_atom_index = 0  # Replace with the start of the desired atom index range
    end_atom_index = 14753  # Replace with the end of the desired atom index range
    timeframes = [i / 1.0 for i in range(0, 200001, 10)]  # 0 to 100000 with step 10
    calculate_sums(xvg_file, timeframes, start_atom_index, end_atom_index)

def extract_forces_at_time(xvg_file, start_atom_index, end_atom_index, time_input):
    """
    Extracts forces (X, Y, Z components) for a range of atoms at a given time from a GROMACS .xvg file.

    Args:
        xvg_file (str): Path to the GROMACS .xvg file.
        start_atom_index (int): Start index of the atom range (0-based).
        end_atom_index (int): End index of the atom range (inclusive).
        time_input (float): Simulation time to extract forces for.

    Returns:
        DataFrame: DataFrame containing time, atom index, and force components (Fx, Fy, Fz).
    """
    forces_at_time = []
    with open(xvg_file, 'r') as file:
        for line in file:
            if not line.startswith(('#', '@')):
                tokens = line.split()
                time = float(tokens[0])
                if time == time_input:
                    if len(tokens) >= (end_atom_index * 3 + 4):
                        for atom_index in range(start_atom_index, end_atom_index + 1):
                            fx = float(tokens[atom_index * 3 + 1])
                            fy = float(tokens[atom_index * 3 + 2])
                            fz = float(tokens[atom_index * 3 + 3])
                            forces_at_time.append([time, atom_index, fx, fy, fz])
                    break

    return pd.DataFrame(forces_at_time, columns=["time", "index", "Fx", "Fy", "Fz"])

def calculate_sums(xvg_file, timeframes, start_atom_index, end_atom_index):
    results = pd.DataFrame(columns=['time', 'sum_Fx', 'sum_Fy', 'sum_Fz'])

    for time in timeframes:
        print(f"Processing time = {time} ps")
        df = extract_forces_at_time(xvg_file, start_atom_index, end_atom_index, time)

        if not df.empty:
            sum_Fx = df['Fx'].sum()
            sum_Fy = df['Fy'].sum()
            sum_Fz = df['Fz'].sum()

            new_row = pd.DataFrame({
                'time': [time],
                'sum_Fx': [sum_Fx],
                'sum_Fy': [sum_Fy],
                'sum_Fz': [sum_Fz]
            })
            results = pd.concat([results, new_row], ignore_index=True)

    results.to_csv("Force_sums.csv", index=False)
    print("Sums saved to Force_sums.csv")

if __name__ == "__main__":
    main()
