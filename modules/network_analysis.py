# modules/network_analysis.py

import pandas as pd
import networkx as nx
import os

def main():
    validated_path = "results/validated_triplets.csv"
    network_path = "results/cerna_network.graphml"
    cytoscape_path = "results/cerna_network.sif"
    centrality_path = "results/centrality_scores.csv"
    nodes_path = "results/cerna_network_nodes.csv"  # New: nodes export
    edges_path = "results/cerna_network_edges.csv"  # New: edges export

    # Load validated triplets
    validated = pd.read_csv(validated_path)

    # Build network
    G = nx.Graph()
    for _, row in validated.iterrows():
        G.add_edge(row['lncRNA'], row['mRNA'], miRNA=row['miRNA'], score=row['score'])

    # Add gene names as node attributes
    all_genes = set(validated['lncRNA']) | set(validated['mRNA'])
    for gene in all_genes:
        G.nodes[gene]['name'] = gene

    # Save as GraphML
    nx.write_graphml(G, network_path)
    print(f"Network saved to {network_path}")

    # Save as SIF for Cytoscape
    with open(cytoscape_path, "w") as f:
        for u, v, data in G.edges(data=True):
            f.write(f"{u} interacts {v}\n")
    print(f"Cytoscape SIF file saved to {cytoscape_path}")

    # Compute and save centrality scores
    centrality = nx.degree_centrality(G)
    centrality_df = pd.DataFrame(list(centrality.items()), columns=['gene', 'degree_centrality'])
    centrality_df = centrality_df.sort_values('degree_centrality', ascending=False)
    centrality_df.to_csv(centrality_path, index=False)
    print(f"Centrality scores saved to {centrality_path}")

    # Export nodes to CSV
    nodes = pd.DataFrame(G.nodes(data=True), columns=['gene', 'attributes'])
    nodes.to_csv(nodes_path, index=False)
    print(f"Nodes exported to {nodes_path}")

    # Export edges to CSV
    edges = pd.DataFrame(G.edges(data=True), columns=['source', 'target', 'attributes'])
    edges.to_csv(edges_path, index=False)
    print(f"Edges exported to {edges_path}")

if __name__ == "__main__":
    main()