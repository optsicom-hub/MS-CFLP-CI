import argparse
import os
import random
import re


def build_output_file_path(input_file_path, suffix):
    """Build an output file name using the requested variant suffix."""
    base_name, extension = os.path.splitext(input_file_path)
    return f"{base_name}_{suffix}{extension}"


def get_instance_stats(input_file_path):
    """Read the instance and return the number of stores and original incompatible pairs."""
    with open(input_file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    stores_match = re.search(r'Stores\s*=\s*(\d+)\s*;', content)
    if not stores_match:
        raise ValueError("Could not find the 'Stores' definition in the file.")

    stores = int(stores_match.group(1))

    pairs_match = re.search(r'IncompatiblePairs\s*=\s*\[\|\s*(.*?)\s*\|\]\s*;', content, re.DOTALL)
    if not pairs_match:
        raise ValueError("Could not find the 'IncompatiblePairs = [| ... |];' matrix in the file.")

    raw_pairs_block = pairs_match.group(1)
    original_pairs = []
    pairs_found = re.findall(r'(\d+)\s*,\s*(\d+)', raw_pairs_block)
    for u, v in pairs_found:
        pair = (min(int(u), int(v)), max(int(u), int(v)))
        original_pairs.append(pair)

    return stores, len(set(original_pairs))


def parse_and_modify_instance(input_file_path, output_file_path, target_count):
    """
    Reads an instance text file, extracts the current incompatible pairs,
    and generates a new file with a target number of incompatibilities.
    """
    with open(input_file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    stores_match = re.search(r'Stores\s*=\s*(\d+)\s*;', content)
    if not stores_match:
        raise ValueError("Could not find the 'Stores' definition in the file.")

    stores = int(stores_match.group(1))
    max_possible = int((stores * (stores - 1)) / 2)

    pairs_match = re.search(r'IncompatiblePairs\s*=\s*\[\|\s*(.*?)\s*\|\]\s*;', content, re.DOTALL)
    if not pairs_match:
        raise ValueError("Could not find the 'IncompatiblePairs = [| ... |];' matrix in the file.")

    raw_pairs_block = pairs_match.group(1)
    original_pairs = []
    pairs_found = re.findall(r'(\d+)\s*,\s*(\d+)', raw_pairs_block)
    for u, v in pairs_found:
        pair = (min(int(u), int(v)), max(int(u), int(v)))
        original_pairs.append(pair)

    original_set = set(original_pairs)
    all_possible_pairs = set((i, j) for i in range(1, stores + 1) for j in range(i + 1, stores + 1))
    non_existing_pairs = list(all_possible_pairs - original_set)

    current_count = len(original_set)
    target_count = max(0, min(target_count, max_possible))

    if target_count < current_count:
        new_pairs = random.sample(list(original_set), target_count)
    else:
        needed = target_count - current_count
        needed = min(needed, len(non_existing_pairs))
        if needed > 0:
            added_pairs = random.sample(non_existing_pairs, needed)
            new_pairs = list(original_set) + added_pairs
        else:
            new_pairs = list(original_set)

    new_pairs.sort()
    new_count = len(new_pairs)

    # Format all pairs in a single line
    formatted_pairs = "IncompatiblePairs = [| "
    line_elements = []
    for p in new_pairs:
        line_elements.append(f"{p[0]}, {p[1]}")
    formatted_pairs += " | ".join(line_elements) + " |];"

    modified_content = re.sub(
        r'(Incompatibilidades|Incompatibilities)\s*=\s*\d+\s*;',
        lambda m: f'{m.group(1)} = {new_count};',
        content,
        flags=re.IGNORECASE
    )
    modified_content = re.sub(
        r'IncompatiblePairs\s*=\s*\[\|.*?\|\]\s*;',
        formatted_pairs,
        modified_content,
        flags=re.DOTALL
    )

    with open(output_file_path, 'w', encoding='utf-8') as f:
        f.write(modified_content)

    print(f"Instance successfully generated!")
    print(f"-> Stores: {stores}")
    print(f"-> Target incompatibilities: {target_count}")
    print(f"-> Actual generated incompatibilities: {new_count}")
    print(f"-> File saved to: '{output_file_path}'")

    return new_count


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate half and double incompatibility variants from an instance.")
    parser.add_argument("instance_name", help="Input instance file, for example: wlp01.txt")
    args = parser.parse_args()

    original_instance = args.instance_name
    stores, original_count = get_instance_stats(original_instance)
    max_possible = int((stores * (stores - 1)) / 2)

    half_target = max(0, original_count // 2)
    double_target = min(max_possible, original_count * 2)

    half_output = build_output_file_path(original_instance, "half")
    double_output = build_output_file_path(original_instance, "double")

    print(f"Generating variants for: {original_instance}")
    print(f"-> Original incompatibilities: {original_count}")
    print(f"-> Half target: {half_target}")
    print(f"-> Double target: {double_target}")

    parse_and_modify_instance(original_instance, half_output, half_target)
    parse_and_modify_instance(original_instance, double_output, double_target)
