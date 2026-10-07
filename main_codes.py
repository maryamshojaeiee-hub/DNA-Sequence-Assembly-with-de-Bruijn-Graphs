"""
author : [ Maryam Shojaei Shahrokhabadi ]
studentnumber : [ 2469734 ]
"""

from pathlib import Path
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from re import findall


def read_csv(name):
    "This function reads the file"
    name = f"{name}.csv"
    # Define the folder path
    my_path = Path(__file__).parent / "csv-files"
    file_path = my_path / name  # indicate the name of file in the folder
    df = pd.read_csv(file_path, header=None)
    df.columns = ['SegmentNr', 'Position', 'A', 'C', 'G', 'T']
    return df


def deal_Duplicated_position(data):
    for segment in data['SegmentNr'].unique():
        filter = (data['SegmentNr'] == segment)
        data[filter].sort_values(by=['Position'])
        numpy_ar = data[filter].to_numpy()
        ss = numpy_ar.shape
        flag = 0
        if ss[0] > max(data[filter]['Position']):
            rows_to_drop = []
            for i in range(1, ss[0] - 1):
                if np.array_equal(numpy_ar[i][1], numpy_ar[i - 1][1]):
                    if np.array_equal(numpy_ar[i][2:], numpy_ar[i - 1][2:]):
                        rows_to_drop.append(data[filter].index[i])
                    else:
                        flag = 1
            if flag == 1:
                a = data[filter].index
                for value in a:
                    rows_to_drop.append(value)
            data = data.drop(rows_to_drop)
    return data


def deal_Missing_position(data):
    for segment in data['SegmentNr'].unique():
        b = 0
        filter = (data['SegmentNr'] == segment)
        for i in range(1, max(data[filter]['Position']) + 1):
            if i in data[filter]['Position'].values:
                b += 1
        if b != max(data[filter]['Position']):
            indexes_drop = data[filter].index
            data = data.drop(indexes_drop)
    return data


def deal_Wrong_position(data):
    for segment in data['SegmentNr'].unique():
        b = 0
        filter = (data['SegmentNr'] == segment)
        data[filter].sort_values(by=['Position'])
        numpy_array = data[filter].to_numpy()
        ss = numpy_array.shape
        for i in range(ss[0]):
            number = sum(numpy_array[i][2:])
            if number != 1:
                b += 1
        if b != 0:
            indexes_drop = data[filter].index
            data = data.drop(indexes_drop)
    return data


def deal_Duplicated_segments(data):
    seen_sequences = {}
    extra_segments = []
    for segment in data['SegmentNr'].unique():
        filter = (data['SegmentNr'] == segment)
        data[filter].sort_values(by=['Position'])
        numpy_array = data[filter].to_numpy()
        new_array = numpy_array[:, 1:]
        sequence_key = []
        for i in range(new_array.shape[0]):
            row_tuple = tuple(row.item() for row in new_array[i])
            sequence_key.append(row_tuple)
        if tuple(sequence_key) not in seen_sequences:
            seen_sequences[tuple(sequence_key)] = segment.item()
        else:
            extra_segments.append(segment.item())
    data = data[~data['SegmentNr'].isin(extra_segments)]
    return data


def clean_data(data):
    step1 = deal_Duplicated_position(data)
    step2 = deal_Duplicated_segments(step1)
    step3 = deal_Wrong_position(step2)
    final_data = deal_Missing_position(step3)
    return final_data


def decode(segment):
    segment.sort_values(by=['Position'])
    bases = ['A', 'C', 'G', 'T']
    sequence = ''
    for i in range(segment.shape[0]):
        d = segment.iloc[i]
        for item in bases:
            if d[item] == 1:
                sequence += item
    return sequence


def generate_sequences(data):
    JSON_sequences = {}
    for segment in data['SegmentNr'].unique():
        filter = (data['SegmentNr'] == segment)
        nn = decode(data[filter])
        JSON_sequences[segment.item()] = nn
    return JSON_sequences


def kmers(JSON_sequences, k):
    kmer_dict = {}
    for segment, sequence in JSON_sequences.items():
        kmers_list = []
        for i in range(len(sequence) - k + 1):
            kmers = sequence[i:i+k]
            kmers_list.append(kmers)
        kmer_dict[segment] = kmers_list
    return kmer_dict


def L_R_segments(kmer_dict, k):
    all_L_R = []
    kk = k - 1
    for kmers in kmer_dict.values():
        for value in kmers:
            L = value[:kk]
            R = value[1:]
            all_L_R.append(L)
            all_L_R.append(R)
    return all_L_R


def construct_graph(json_data, k):
    graph = nx.MultiDiGraph()
    kmer_dict = kmers(json_data, k)
    L_R = L_R_segments(kmer_dict, k)
    for i in range(0, len(L_R), 2):
        L = L_R[i]
        R = L_R[i+1]
        graph.add_edge(L, R)
    return graph


def count_parallel_edges(graph):
    edge_counts = {}
    for u, v, _ in graph.edges(keys=True):
        if (u, v) in edge_counts:
            edge_counts[(u, v)] += 1
        else:
            edge_counts[(u, v)] = 1
    return edge_counts


def plot_graph(graph, filename):
    nodes = list(graph.nodes)
    pos = {}
    nodes_per_row = 5
    for i, node in enumerate(nodes):
        row = i // nodes_per_row
        col = i % nodes_per_row
        pos[node] = (col, -row)
    plt.figure(figsize=(8, 8))
    nx.draw(graph, pos, with_labels=True, arrows=True,
            node_size=1000, font_size=10)
    edge_counts = count_parallel_edges(graph)
    edge_labels = {edge: str(count) for edge, count in edge_counts.items()}
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels,
                                 font_color='red')
    plt.savefig(f"{filename}.png", format='png')


def check_connected(graph):
    undirected_graph = nx.Graph()
    for u, v in graph.edges():
        undirected_graph.add_edge(u, v)
    visited = []
    stack = []
    nodes = list(undirected_graph)
    stack.append(nodes[0])
    while len(stack) != 0:
        i = stack.pop()
        if i not in visited:
            visited.append(i)
            neighbours = list(undirected_graph[i])
            for w in neighbours:
                stack.append(w)
    return len(visited) == len(nodes)


def startnode_eulerian_path(graph):
    for node in graph.nodes():
        if graph.out_degree(node) - graph.in_degree(node) == 1:
            return node
    for node in sorted(graph.nodes()):
        return node


def is_valid_graph(graph):
    diff_count = 0
    nodes = []
    if check_connected(graph):
        for node in graph.nodes():
            in_deg = graph.in_degree(node)
            out_deg = graph.out_degree(node)
            if in_deg != out_deg:
                diff_count += 1
                nodes.append((node, in_deg - out_deg))
        if diff_count == 0:
            return True
        elif diff_count == 2:
            diff_vector = [diff for _, diff in nodes]
            if sorted(diff_vector) == [-1, 1]:
                return True
        else:
            return False
    else:
        return False


def out_edges_multidigraph(graph):
    out_edges = {}
    for u, v in graph.edges():
        if u not in out_edges:
            out_edges[u] = []
        out_edges[u].append(v)
    return out_edges


def path_to_stuck(start_node, out_edges):
    path = []
    current = start_node
    while out_edges[current]:
        path.append(current)
        next_node = out_edges[current].pop(0)
        current = next_node
    path.append(current)
    return path, out_edges


def extend_euler_path(graph, euler_path, out_edges):
    for i, node in enumerate(euler_path):
        if out_edges[node]:
            new_path, out_edges = path_to_stuck(node, out_edges)
            euler_path = euler_path[:i] + new_path + euler_path[i+1:]
            return extend_euler_path(graph, euler_path, out_edges)
    return euler_path


def find_eulerian_path(graph, start_node):
    out_edges = out_edges_multidigraph(graph)
    initial_path, out_edges = path_to_stuck(start_node, out_edges)
    full_path = extend_euler_path(graph, initial_path, out_edges)
    return full_path


def construct_dna_sequence(graph):
    out_str = ' '
    if is_valid_graph(graph):
        start_node = startnode_eulerian_path(graph)
        Euler_path = find_eulerian_path(graph, start_node)
        for_print = ' - '.join(Euler_path)
        print(for_print)
        for i in range(len(Euler_path) - 1):
            out_str += Euler_path[i][0]
        out_str += Euler_path[-1]
    else:
        out_str = "DNA sequence can not be constructed."
    return out_str


def save_output(s, filename):
    with open(f"{filename}.txt", "w") as file:
        file.write(s)


file_name = 'DNA_2_5'
pattern = r"DNA_(\d+)_(\d+)"
matches = findall(pattern, file_name)
if matches:
    out_file = f"DNA_{matches[0][0]}"
    k = int(matches[0][1])
else:
    print("Filename format is invalid.")
a = read_csv(file_name)
data = clean_data(a)
JSON_sequences = generate_sequences(data)
graph = construct_graph(JSON_sequences, k)
plot_graph(graph, out_file)
string = construct_dna_sequence(graph)
save_output(string, out_file)
