"""
Updated by Sofiane
Generate toy graphs with a groundtruth entity graph.
"""
import random
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import collections as mc
import os

def generate_trajectory_graphs(trajectories_path, n_subgraphs, bbox=None):
	"""

	:param trajectories_path:
	:param bbox:
	:param n_subgraphs:
	:return: list of subgraphs corresponding to traces, bbox of the traces considered.
	"""
	max_lon, min_lon, max_lat, min_lat = float("-inf"), float("inf"), float("-inf"), float("inf")
	SGs = []
	for sg_idx, fname in enumerate(os.listdir(trajectories_path)[:n_subgraphs]):
		with open(os.path.join(trajectories_path, fname)) as f:
			edge_lst = []
			lines = f.readlines()
			for i in range(len(lines) - 1):
				s = map(float, lines[i].split(','))
				t = map(float, lines[i+1].split(','))
				# extract node as (lon, lat)
				edge_lst.append(((s[2], s[1]), (t[2], t[1])))

				# update bbox:
				if max_lon < s[2]: max_lon = s[2]
				if min_lon > s[2]: min_lon = s[2]
				if max_lat < s[1]: max_lat = s[1]
				if min_lat > s[1]: min_lat = s[1]
			if len(edge_lst) > 0:
				rn = nx.DiGraph(edge_lst)
				for node in rn.nodes():
					rn.node[node]['subgraph'] = sg_idx
				SGs.append(rn)
	return SGs, (max_lon, min_lon, max_lat, min_lat)


def draw_subgraphs(SGs, bbox, showNodes=False):
	fig, ax = plt.subplots(nrows=len(SGs)/3+1, ncols=3, figsize=(30, 40))

	union_lines = []
	for i, G in enumerate(SGs):
		lines = [[s, t] for s, t in G.edges()]
		union_lines += lines
		lc = mc.LineCollection(lines, colors='black', linewidths=1)
		ax1 = plt.subplot(len(SGs)/3+1, 3, i+1)
		ax1.add_collection(lc)
		ax1.autoscale()
		ax1.set_xlim(bbox[1], bbox[0])
		ax1.set_ylim(bbox[3], bbox[2])
		if showNodes == True:
			plt.scatter([node[0] for node in G.nodes()], [node[1] for node in G.nodes()], s=10)

	# Print Union of all subgraphs:
	lc = mc.LineCollection(union_lines, colors='red', linewidths=2)
	ax1 = plt.subplot(len(SGs) / 3 + 1, 3, len(SGs) + 1)
	ax1.add_collection(lc)
	ax1.autoscale()
	ax1.set_xlim(bbox[1], bbox[0])
	ax1.set_ylim(bbox[3], bbox[2])
	if showNodes == True:
		plt.scatter([node[0] for G in SGs for node in G.nodes()], [node[1] for G in SGs for node in G.nodes()], s=10)
	plt.show()

if __name__ == '__main__':
	# sGs_, G_ = generate_graphs(20, 2, 10, 1, 3)
	# draw_graphs([G_]+sGs_)
	traj_path = '/home/sofiane/projects/data/GPS_Dataprivatebus_trips_20s_sample/'
	SGs, actual_bbox = generate_trajectory_graphs(trajectories_path=traj_path, n_subgraphs=20)
	for node in SGs[0].nodes():
		print node, SGs[0].node[node]['subgraph']
	print actual_bbox
	draw_subgraphs(SGs, bbox=actual_bbox, showNodes=True)