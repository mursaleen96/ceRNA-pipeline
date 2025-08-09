# modules/generate_report.py

import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from plotly.offline import plot
import os

def main():
    validated_path = "results/validated_triplets.csv"
    network_path = "results/cerna_network.graphml"
    report_path = "results/cerna_analysis_report.html"

    # Load validated triplets
    validated = pd.read_csv(validated_path)

    # Load network
    G = nx.read_graphml(network_path)

    # Generate interactive network visualization with gene names as labels
    pos = nx.spring_layout(G)
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(x=edge_x, y=edge_y, line=dict(width=0.5, color='#888'), hoverinfo='none', mode='lines')

    node_x = []
    node_y = []
    node_text = []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(node)  # Gene name as label

    node_trace = go.Scatter(x=node_x, y=node_y, mode='markers+text', text=node_text, textposition='top center',
                            hoverinfo='text', marker=dict(size=10, color='skyblue'))  # Simplified marker without colorbar

    fig = go.Figure(data=[edge_trace, node_trace], layout=go.Layout(showlegend=False, hovermode='closest', margin=dict(b=20, l=5, r=5, t=40), xaxis=dict(showgrid=False, zeroline=False, showticklabels=False), yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
    
    # Embed full Plotly JS
    network_html = plot(fig, output_type='div', include_plotlyjs=True)

    # Generate HTML report
    html_content = """
    <html>
    <head><title>ceRNA Analysis Report</title></head>
    <body>
    <h1>ceRNA Analysis Report</h1>
    <h2>Validated Triplets</h2>
    {table}
    <h2>ceRNA Network</h2>
    {network}
    </body>
    </html>
    """
    table_html = validated.to_html(index=False)
    full_html = html_content.format(table=table_html, network=network_html)

    with open(report_path, "w") as f:
        f.write(full_html)
    print(f"Report generated at {report_path}")

if __name__ == "__main__":
    main()