import tomlkit
import sys


def parse_old_clusterfile(file_path):
    # Read the content of the file
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # A dictionary to store parsed data
    data = {}
    
    current_section = None

    # Loop through each line
    for line in lines:
        line = line.strip()  # Remove leading/trailing whitespace

        # If the line starts with "[", it's a section header
        if line.startswith('['):
            current_section = line[1:-1]
            data[current_section] = {}
        else:
            # Check for monomer keyword
            if line == "monomer":
                data[current_section]["isMonomer"] = True
            elif "composition" in line:
                comps = line.split()[1:]
                data[current_section]["composition"] = [int(comp) for comp in comps]
            elif "sigma" in line:
                data[current_section]["sigma"] = int(line.split()[1])
            elif "energy" in line:
                data[current_section]["energy"] = float(line.split()[1])
            elif "volume" in line:
                data[current_section]["volume"] = float(line.split()[1])
            elif "coordinates" in line:
                data[current_section]["coordinates"] = line.split()[1]
            elif "frequencies" in line:
                data[current_section]["frequencies"] = line.split()[1]

    return data

def write_to_toml(data, output_path):
    with open(output_path, 'w') as f:
        f.write(tomlkit.dumps(data))

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <clusterset-file> <toml_clusterset-file>")
        sys.exit(1)

    old_file = sys.argv[1]
    new_file = sys.argv[2]


    parsed_cluster_data = parse_old_clusterfile(old_file)
    write_to_toml(parsed_cluster_data, new_file)

if __name__ == "__main__":
    main()